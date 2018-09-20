#!/usr/bin/env python3
import os
import sys
import time
import random
import subprocess

from ete3 import Tree

from .utils import data_sources as ds

##
CHUNK_FOLDER_SUFFIX="fasta_chunks"
CHUNK_INPUT_PREFIX="fasta_samples"
CHUNK_OUTPUT_PREFIX=""

#### FASTTREE BINARIES

FASTTREE_WEB = os.path.join(ds.BINARIES,'FastTree')
FASTTREE_FULL_OPTIONS = ' -nopr -pseudo -mlacc 3 -slownni'
# Options for slow nni version from ete build workflow "fasttree_full":
# -nopr [noprogress]
# -slownni to turn off heuristics to avoid constant subtrees (affects both ML and ME NNIs)
# -mlacc 3 to always optimize all 5 branches at each NNI, and to optimize all 5 branches in 3 rounds
# -pseudo [weight] -- Use pseudocounts to estimate distances between sequences with little or no overlap.
#   (Off by default.) Recommended if analyzing the alignment has sequences with little or no overlap.
#   If the weight is not specified, it is 1.0

#### MAFFT BINARIES

MAFFT_WEB = os.path.join(ds.BINARIES,'mafft-linux64/mafft.bat')
MAFFT_AUTO_OPTION = " --auto -"
MAFFT_LEAVEGAPS_OPTION = " --auto --leavegappyregion -"

QSUB_TEMPLATE="""
#!/bin/bash

### SGE VARIABLES

#$ -S /bin/bash
#$ -cwd
#$ -sync y
#$ -N tree_task
#$ -l 'h_rt=3:00:00'
#$ -q short.virtual.q 
#$ -e {base_dir}/{higher_level}.fasta_chunks/
#$ -o {base_dir}/{higher_level}.fasta_chunks/
#$ -t 1:{chunk_no}

# SGE_TASK_ID=1 # testing only

DATA_DIR="{base_dir}/{higher_level}.fasta_chunks"

MEMSAVE_SAMPLES="{memsave_sample_ids}"
"""

TREE_TEMPLATE="""

### TREE BUILDING PIPELINE

MAFFT='/mnt/eris/davide/programs/mafft/bin/mafft --auto --amino'
FASTTREE='/mnt/eris/davide/programs/FastTree -nopr -pseudo -mlacc 3 -slownni'

TASK_SUFFIX=$(printf "%04d" ${SGE_TASK_ID})

INFILE="${DATA_DIR}/fasta_samples.${TASK_SUFFIX}.tsv"
TMPFILE=$(mktemp)
OUTFILE="${DATA_DIR}/fasta_samples.${TASK_SUFFIX}.nw"

# Go through tab file [NOG_ID]tab[SAMPLE_ID]tab[FASTA_SEQ]
while IFS=$'\\t' read -r -a myArray
do
    NOG_ID=${myArray[0]}
    SAMPLE_ID=${myArray[1]}
    FASTA_STR=${myArray[2]}
    
    FASTA=$(echo "$FASTA_STR"| tr -d "'")
    
    # M1 check if any line of the fasta file surpasses 10k chars
    # echo -e "$FASTA" | awk '{ if (length($0)>10000) exit 1 }'
    # M2 check list membership:
    echo $MEMSAVE_SAMPLES | grep -q -w $SAMPLE_ID
    if [ $? -eq 0 ]; then
        echo "$NOG_ID $SAMPLE_ID memsave" >&2
        NW=$($MAFFT --memsave <(echo -e "$FASTA") | $FASTTREE);
    else
        echo "$NOG_ID $SAMPLE_ID" >&2
        NW=$($MAFFT <(echo -e "$FASTA") | $FASTTREE);
    fi
    echo -e "${myArray[0]}\\t${myArray[1]}\\t$NW" >> $TMPFILE
done < $INFILE

mv $TMPFILE $OUTFILE

"""

def get_tree_chunks(tree_jobs,n):
    random.shuffle(tree_jobs)
    return [tree_jobs[i:i+n] for i in range(0,len(tree_jobs),n)]
    # does this cover all?

def write_tree_tasks(higher_level,tree_jobs,chunk_no = 1,min_size=1,max_size=50):
    if chunk_no > 1:
        
        memsave_sample_ids=[]
        
        tot_job_no = len(tree_jobs)
        chunk_avg = tot_job_no / chunk_no
        chunk_size = max(min_size,min(chunk_avg,max_size))
        
        chunk_dir = '%d.%s'%(higher_level,CHUNK_FOLDER_SUFFIX)
        os.mkdir(chunk_dir)
        
        # divide jobs into chunks
        tree_chunks = get_tree_chunks(tree_jobs, chunk_size)
        
        for i,tree_chunk in enumerate(tree_chunks):
            chunk_file = os.path.join(chunk_dir,'%s.%04d.tsv'%(CHUNK_INPUT_PREFIX,i+1)) 
            with open(chunk_file,'w') as f:             
                for nog_id, sample_counter, fasta_sequences in tree_chunk:
                    f.write("%s\t%d\t%s\n"%(
                        nog_id,sample_counter,repr(fasta_sequences)))
                    
                    # check if any sequence is beyond 5k aa
                    if max(map(len,fasta_sequences.split('>'))) > 4000:
                        memsave_sample_ids.append(sample_counter)
                    
        return len(tree_chunks),memsave_sample_ids
    else:
        with open('%d.fasta_samples.tsv'%higher_level,'w') as f:             
            for nog_id, sample_counter, fasta_sequences in tree_jobs:
                f.write("%s\t%d\t%s\n"%(nog_id,sample_counter,repr(fasta_sequences)))
        return 1

def read_tree_tasks(higher_level):
    results_dir = "%d.%s"%(higher_level,CHUNK_FOLDER_SUFFIX)
    assert os.path.exists(results_dir),'No directory found at %s'%results_dir
    input_tsv_wildcard = '%s/*.tsv'%results_dir

    print('concatenating %s'%input_tsv_wildcard)

    p = subprocess.Popen(['cat %s'%input_tsv_wildcard],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    fasta_input, err = p.communicate()
    
    tree_lines = [x.split('\t') for x in fasta_input.split('\n')]
    fasta_jobs = [(int(x[0]),int(x[1]),eval(x[2])) for x in tree_lines if len(x) == 3]
    
    return sorted(fasta_jobs,key=lambda x: x[1])

def submit_cluster_arrayTask(higher_level,tree_jobs):
    
    granularity = 1000
    chunk_no,memsave_sample_ids = write_tree_tasks(higher_level,tree_jobs,granularity)
    #chunk_no = 1456
    
    qsub_cmd = QSUB_TEMPLATE.format(
        chunk_no=chunk_no,
        higher_level=higher_level,
        base_dir=os.getcwd(),
        memsave_sample_ids=memsave_sample_ids)
    
    qsub_cmd += TREE_TEMPLATE
    
    qsub_script='tree_task.%d.qsub'%higher_level
    with open(qsub_script,'wb') as f:
        f.write(qsub_cmd)
        
    subprocess.call(['qsub',qsub_script])
    
def collect_cluster_results(higher_level):
    results_dir = "%d.%s"%(higher_level,CHUNK_FOLDER_SUFFIX)
    assert os.path.exists(results_dir),'No directory found at %s'%results_dir
    results_nw_wildcard = '%s/*.nw'%results_dir

    print('concatenating %s'%results_nw_wildcard)

    p = subprocess.Popen(['cat %s'%results_nw_wildcard],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    results, err = p.communicate()
    
    tree_lines = [x.split('\t') for x in results.split('\n')]
    tree_computations = [(int(x[0]),int(x[1]),x[2]) for x in tree_lines if len(x) == 3]
    
    return sorted(tree_computations,key=lambda x: x[1])

def ete_build_workflow(tree_job):
    
    nog_id, i, input_fasta = tree_job
    
    # check if sample folder exists, if not create
    input_dir = '%d.samples'%nog_id
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)
    
    # write fasta file to disk to use ete build workflows
    input_file = os.path.join(input_dir,"%d.fa"%i)
    with open(input_file,'w') as f:
        f.write(input_fasta)
    
    # Save results in directory with .ete extension    
    tree_dir = os.path.join(input_dir,'%d.ete'%i)
    
    time_start = time.time()
    # execute ete build workflow
    ete_workflow = "mafft_default-none-none-fasttree_full"
    subprocess.call(
        "ete3 build -w %s -a %s -o %s --noimg -v 0 --compress --cpu 1"%(
            ete_workflow,input_file,tree_dir),shell=True)
    time_end = time.time()
    
    tree_file = os.path.join(tree_dir,
                             "%s/%d.fa.final_tree.nw"%(
                                ete_workflow,i))
    with open(tree_file) as f:
        tree_nw = f.readline()
        
    return (nog_id, time_end - time_start, tree_nw)

def direct(mafft_cmd, fasttree_cmd, tree_job):
    
    nog_id, sample_no, input_fasta = tree_job
    
    p = subprocess.Popen("%s | %s"%(mafft_cmd,fasttree_cmd),
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    
    tree_nw, err = p.communicate(input_fasta.encode())
    
    return (nog_id, sample_no, tree_nw.decode())

def direct_fast(tree_job):
    mafft_cmd = MAFFT_WEB + MAFFT_AUTO_OPTION
    fasttree_cmd = FASTTREE_WEB
    return direct(mafft_cmd, fasttree_cmd, tree_job)
  
def direct_full(tree_job):
    mafft_cmd = MAFFT_WEB + MAFFT_AUTO_OPTION
    fasttree_cmd = FASTTREE_WEB + FASTTREE_FULL_OPTIONS
    return direct(mafft_cmd, fasttree_cmd, tree_job)

def backup_alg(mafft_cmd, fasttree_cmd, tree_job):
    
    nog_id, i, input_fasta = tree_job
    
    # cool stuff: tee `mktemp %s.XXXX`
    
    t_start = time.time()
    p = subprocess.Popen("%s | tee %d.%d.fa | %s"%(mafft_cmd,nog_id,i,fasttree_cmd),
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    
    tree_nw, err = p.communicate(input_fasta)
    t_end = time.time()
    
    return (nog_id, t_end - t_start, tree_nw)
    
def reroot_trees(tree_computations, species_tree_polytomies):
    rerooted = []
    for nog_id,tree_time,tree_nw in tree_computations:
        assert tree_nw, "Tree newick is non existant for %d %d"%(nog_id,tree_time)
        t = Tree(tree_nw)
        if species_tree_polytomies:
            # reconciliation algorithm can only have one input with multifurcations/polytomies
            t.resolve_polytomy(recursive=True)
        
        node = t.get_midpoint_outgroup()
        if node:
            t.set_outgroup(node)
            rerooted_job = (nog_id,tree_time,t.write())
            rerooted.append(rerooted_job)
        else:
            sys.stderr.write('Problems in rerooting %s %s %s'%(nog_id,tree_time,tree_nw))
            rerooted.append((nog_id,tree_time,tree_nw))
            
    return rerooted
    
