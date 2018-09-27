#!/usr/bin/env python3
# Module containing functions for gene Tree - species Tree reconciliation
# author:   davide.heller@imls.uzh.ch
# date:     06-03-2017

import os
import sys
import time
import random
import tempfile
import argparse
import subprocess

from ete3 import Tree
from ete3.parser.newick import NewickError
from multiprocessing import Pool

from .utils import data_sources as ds

EGGNOGv4_SPECIES_TREE = os.path.join(ds.EGGNOG_OUTPUT,"eggNOG_species_tree.nw")

QSUB_TEMPLATE="""
#!/bin/bash

### SGE VARIABLES

#$ -S /bin/bash
#$ -cwd
#$ -l 'h_vmem=2G'
#$ -sync y
#$ -q short.virtual.q 
#$ -l 'h_rt=3:00:00'
#$ -N reconciliation_task
#$ -e {base_dir}/{higher_level}.tree_chunks/
#$ -o {base_dir}/{higher_level}.tree_chunks/
#$ -t 1:{chunk_no}

#SGE_TASK_ID=1 # testing only

DATA_DIR="{base_dir}/{higher_level}.tree_chunks"

NOTUNG='java -Xmx1G -XX:+UseSerialGC -jar /mnt/eris/davide/programs/Notung-2.9.jar'
NOTUNG_OPT='{notung_options} --speciestag prefix --parsable --edgeweights name'
NOTUNG_EXT='{notung_ext}'

"""
REARRANGE_OPTIONS="--rearrange --threshold 0.9"
REARRANGE_EXT="rearrange"

ROOTING_OPTION="--root"
ROOTING_EXT="rooting"

TRANSFER_OPTION="--reconcile --infertransfers true"
TRANSFER_EXT="reconciled"

NOTUNG_TEMPLATE="""

### RECONCILIATION PIPELINE

TASK_SUFFIX=$(printf "%04d" ${SGE_TASK_ID})

INFILE="${DATA_DIR}/tree_samples.${TASK_SUFFIX}.tsv"
TMPFILE=$(mktemp)
OUTFILE="${DATA_DIR}/tree_samples.${TASK_SUFFIX}.ntg"

# Go through tab file [NOG_ID]tab[SAMPLE_ID]tab[FASTA_SEQ]
while IFS=$'\\t' read -r -a myArray
do
    # prepare input
    TREENAME="${myArray[0]}.${myArray[1]}"
    GENETREE="${TREENAME}.gTree"
    echo ${myArray[2]} > $GENETREE
    SPECIESTREE="${TREENAME}.sTree"
    echo ${myArray[3]} > $SPECIESTREE

    # run reconciliation
    $NOTUNG -g $GENETREE -s $SPECIESTREE $NOTUNG_OPT

    # collect output
    RECONCILED="${GENETREE}.${NOTUNG_EXT}.0"
    PARSABLE="${RECONCILED}.parsable.txt"

    echo "#TREESAMPLE $TREENAME" >> $TMPFILE
    cat $RECONCILED >> $TMPFILE
    cat $PARSABLE >> $TMPFILE

    # remove output
    rm $GENETREE
    rm $SPECIESTREE
    rm $RECONCILED
    rm $PARSABLE

done < $INFILE

mv $TMPFILE $OUTFILE

"""


def get_tree_chunks(tree_jobs,n):
    random.shuffle(tree_jobs)
    return [tree_jobs[i:i+n] for i in range(0,len(tree_jobs),n)]

def write_reconciliation_tasks(higher_level,tree_computations,chunk_no = 1,min_size=1,max_size=50):

    if chunk_no > 1:
    
        # compute granularity and divide jobs into chunks
        tot_job_no = len(tree_computations)
        chunk_avg = tot_job_no / chunk_no
        chunk_size = max(min_size,min(chunk_avg,max_size))
        
        chunk_dir = '%d.tree_chunks'%higher_level
        os.mkdir(chunk_dir)
        
        tree_chunks = get_tree_chunks(tree_computations, chunk_size)
        
        # write chunks and compute species tree for each gene tree
        for i,tree_chunk in enumerate(tree_chunks):
            chunk_file = os.path.join(chunk_dir,'tree_samples.%04d.tsv'%(i+1)) 
            with open(chunk_file,'w') as f:             
                for computation in tree_chunk:
                    nog_id,sample_id,gene_nw, species_nw = computation
                    f.write("%d\t%d\t%s\t%s\n"%(nog_id,sample_id,gene_nw,species_nw))
                    
        return len(tree_chunks)
    else:
        # alternatively save all jobs into one file
        with open('%d.tree_samples.tsv'%higher_level,'w') as f:
            for computation in tree_computations[1:10]:
                nog_id,sample_id,gene_nw, species_nw = computation
                f.write("%d\t%d\t%s\t%s\n"%(nog_id,sample_id,gene_nw,species_nw))
        return 1

def submit_taskArray(higher_level,tree_computations,keep_polytomies,root_notung,infer_transfers):
    granularity = 1000
    chunk_no = write_reconciliation_tasks(higher_level,tree_computations,granularity)

    if root_notung:
        notung_options=ROOTING_OPTION
        notung_ext=ROOTING_EXT
    elif infer_transfers:
        notung_options=TRANSFER_OPTION
        notung_ext=TRANSFER_EXT
    else:
        notung_options=REARRANGE_OPTIONS
        notung_ext=REARRANGE_EXT

    qsub_cmd = QSUB_TEMPLATE.format(
        chunk_no=chunk_no,
        higher_level=higher_level,
        base_dir=os.getcwd(),
        notung_options=notung_options,
        notung_ext=notung_ext)
    
    qsub_cmd += NOTUNG_TEMPLATE
    
    qsub_script='reconciliation_task.%d.qsub'%higher_level
    with open(qsub_script,'wb') as f:
        f.write(qsub_cmd)
        
    subprocess.call(['qsub',qsub_script])

def extract_reconciled_tree(nwk):
    nwk_splitted = nwk.split('[')
    nwk_cleaned = nwk_splitted[0]
    for nwk_seq in nwk_splitted[1:]:
        nwk_cleaned += nwk_seq.split(']')[1]
    return Tree(nwk_cleaned, format=1)

def annotate_with_parsable(t, parsable):
    events = {}
    
    for line in parsable:
        if line.startswith('#'):
            #if letter is not dict insert header first
            event_type = line[1]
            
            if event_type not in ['D','T','C']:
                continue
            
            if event_type in events:
                event_node,event_details = node_mapping(t,event_type,line)
                events[event_type][event_node] = event_details[2:]
            else:
                #create category & skip header line
                events[event_type] = {}
    
    # annotate all remaining nodes with speciation events                
    for node in t.traverse():
        if node.is_leaf():
            continue
        elif not hasattr(node,'event'):
            node.add_features(event="S")
    
    return t

def parse_taskArray(reconciliation,legacy_return=True):
    lines = reconciliation.split('\n')
    sample_id = lines[0].strip()
    sample_nwk = lines[1]
    
    try:
        sample_t = extract_reconciled_tree(sample_nwk)
    except NewickError:
        print('error in reconciliation:\n%s'%reconciliation)
        return (sample_id,None)

    parsable = lines[2:]
    annotated_tree = annotate_with_parsable(sample_t,parsable)
    
    if legacy_return:
        root_event = annotated_tree.event
        split_index = compute_split_index(annotated_tree)
        sample_solution = (root_event,split_index)
    else:
        sample_solution = annotated_tree
        
    return (sample_id,sample_solution)

def collect_taskArray(higher_level,legacy_return=True,
                      sampling=None,alt_dir=None, cpu_cores=1):
    
    if sampling is None:
        results_dir = "%d.tree_chunks"%higher_level
    else:
        results_dir = "%d.%s.tree_chunks"%(higher_level,sampling)
        
    if alt_dir is not None:
        results_dir = os.path.join(alt_dir,results_dir)
        
    assert os.path.exists(results_dir),'No directory found at %s'%results_dir
    results_nw_wildcard = '%s/*.ntg'%results_dir
    
    p = subprocess.Popen(['cat %s'%results_nw_wildcard],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    results, err = p.communicate()
    
    reconciliations = results.split('#TREESAMPLE')[1:]
    
    print('Starting to collect %d tree samples'%(len(reconciliations))) 
    
    if cpu_cores > 1:
        p = Pool(cpu_cores)
        reconciled_trees = dict(
            p.map(parse_taskArray,reconciliations))
        p.close()
        p.join()
    else:
        reconciled_trees = {}
        for reconciliation in reconciliations:
            sample_id,sample_solution = parse_taskArray(reconciliation)
            reconciled_trees[sample_id] = sample_solution
    
    max_id = 0  
    for sample_id in reconciled_trees:
        nog_id,sample_no = sample_id.split('.')
        if max_id < int(sample_no):
            max_id = int(sample_no)
        
    # generate legacy output (touple list with preserved order)
    if legacy_return:
        legacy_output = [None] * (max_id  + 1)
        for sample_id in reconciled_trees:
            nog_id, sample_no = [int(x) for x in sample_id.split('.')]
            # '10096041102931.714' -> 10096041102931, 714
            sample_solution = reconciled_trees[sample_id]
            legacy_output[sample_no] = (nog_id,sample_solution)
        return legacy_output
    else:
        return reconciled_trees

def load_reconciliations(higher_level, results_dir):
    # load precomputed reconciliation
    level_reconciliations = os.path.join(results_dir, "%d.tree_reconciliations.tsv"%higher_level)
    
    reconciliations = []
    with open(level_reconciliations) as f:
        for line in f:
            l = line.rstrip().split()
            assert len(l) == 3, "Non standart line %s"%line
            nog_id = int(l[0])
            root_event = l[1]
            split_index = float(l[2])
            reconciliations.append((nog_id,(root_event,split_index)))
        
    # aggregate reconciliations by nog_id
    reconciliation_dict = {}
    for nog_id, result in reconciliations:
        assert result, 'No reconciliation found for %d'%nog_id
        root_event,split_index = result
        if nog_id not in reconciliation_dict:
            reconciliation_dict[nog_id] = {'D':[],'S':[]}
        reconciliation_dict[nog_id][root_event].append(split_index)
    
    return reconciliation_dict

def prune_species_tree(gene_tree, cached_species_tree=None, keep_polytomies=False):
    
    gTree = Tree(gene_tree)
    
    #species reading
    
    #leaf names should be of the type [speciesID_ProteinName]
    leaf_names = gTree.get_leaf_names()
    
    species_list = {x.split('_')[0] for x in leaf_names}
    species_list = list(species_list)

    species_ids = {''.join(filter(str.isdigit,x)):x for x in species_list}

    #big species tree
    if cached_species_tree:
        s = cached_species_tree
    else:
        s = Tree(EGGNOGv4_SPECIES_TREE)
    
    #get lca for core
    common_ancestor = s.get_common_ancestor(list(species_ids.keys())).copy()
    
    #prune to subset
    # common_ancestor.prune(species_ids) # slower method
    leaves = {x.name:x for x in common_ancestor.get_leaves()}
    to_remove = leaves.keys() - species_ids.keys()
    for species_id in to_remove:
        if species_id in leaves:
            leaves[species_id].delete()
    assert(len(common_ancestor.get_leaf_names()) == len(species_ids))
    
    #binarize
    if not keep_polytomies:
        common_ancestor.resolve_polytomy(recursive=True)
    
    #change names
    for leaf in common_ancestor.get_leaves():
        leaf.name = species_ids[leaf.name]
    
    #write out reconciliation_job
    species_nw = common_ancestor.write(format=5)

    return species_nw

def prepare_reconciliation_job(tree_job):
    
    nog_id,sample_id,gene_nw,cached_species_tree,keep_polytomies,root_notung = tree_job
    
    species_nw = prune_species_tree(gene_nw,cached_species_tree,keep_polytomies)
    
    reconciliation_job = (nog_id,sample_id,gene_nw,species_nw,root_notung)
    return reconciliation_job

def compute_split_index(t, annotate_tree=False):
    
    splits = {}
    
    for leaf_name in t.get_leaf_names():
        if "*LOST" in leaf_name:
            continue
        else:
            split_id = leaf_name.split("_")[1][0:3]
            if split_id != "xxx":
                # don't compute split index for singletons
                if split_id not in splits:
                    splits[split_id] = [leaf_name]
                else:
                    splits[split_id].append(leaf_name)
    
    clade_errors = []
    
    for split_id,split_proteins in splits.items():
        
        if len(split_proteins) > 1:
            
            lca = t.get_common_ancestor(split_proteins)
            
            lca_leaves = lca.get_leaf_names()
            lost_leaves = [x for x in lca_leaves if "*LOST" in x]
            singleton_leaves = [x for x in lca_leaves if "_xxx" in x]
            errors = [x for x in lca_leaves if x not in split_proteins+lost_leaves+singleton_leaves]
            
            if errors:           
                error_rate =  len(errors) / float(len(split_proteins) + len(errors)) 
            else:
                error_rate = 0
                
            clade_errors.append( error_rate )
            
            if annotate_tree:
                split_label = "%s[%d|%d|%d|%d]"%(split_id,
                                                 len(lca_leaves),len(split_proteins),
                                                 len(errors),len(lost_leaves))
                
                if hasattr(lca,'split'):
                    lca.split += split_label
                else:
                    lca.add_features(split=split_label)
                    
    if clade_errors:
        return sum(clade_errors)/len(clade_errors)
    else:
        sys.stderr.write('No clade errors found in\n%s\n'%t.get_ascii(attributes=["name","event","split"]))
        return -7.0

def node_mapping(gTree, event_type, line):
    
    event_details = line.split()
    
    # identify node to annotate
    gNode_name = event_details[1]
    nodes = gTree.search_nodes(name=gNode_name)
    assert(nodes),"No node found with name: %s\nExtracted from %s"%(gNode_name,line)
    assert(len(nodes)==1),"More than 1 node found with the same name: %s"%gNode_name
    gNode = nodes[0]
    
    # annotate node with event description
    gNode.add_features(event=event_type)
    
    # special annotation for (T)ransfer events
    if event_type == 'T':
        acceptor = event_details[2]
        acceptor_species = event_details[4]
        donor_species = event_details[3]
        
        donor = None
        for child in gNode.get_children():
            if child.name == acceptor:
                child.add_features(is_acceptor=True)
            else:
                donor = child.name
                child.add_features(is_acceptor=False)
            
        assert(donor),"No donor found!"
        gNode.add_features(donor=donor_species,acceptor=acceptor_species)
    
    return gNode,event_details

def annotate_notung(tree_file, parable_file):
    
    # remove squared brackets from notung nwk format
    with open(tree_file) as f:
        nwk = f.readline().rstrip()
    
    nwk_splitted = nwk.split('[')
    nwk_cleaned = nwk_splitted[0]
    for nwk_seq in nwk_splitted[1:]:
        nwk_cleaned += nwk_seq.split(']')[1]
        
    t = Tree(nwk_cleaned, format=1)
    
    # parsable file
    parsable_file = tree_file + ".parsable.txt"
    events = {}
    
    with open(parsable_file) as f:
        for line in f:
            if line.startswith('#'):
                #if letter is not dict insert header first
                event_type = line[1]
                
                if event_type not in ['D','T']:
                    continue
                
                if event_type in events:
                    event_node,event_details = node_mapping(t,event_type,line)
                    events[event_type][event_node] = event_details[2:]
                else:
                    #create category & skip header line
                    events[event_type] = {}
    
    # annotate all remaining nodes with speciation events                
    for node in t.traverse():
        if node.is_leaf():
            continue
        elif not hasattr(node,'event'):
            node.add_features(event="S")
            
    return t

# @profile    
def reconcile(job, return_tree=False):
    
    nog_id, sample_id, geneTree, speciesTree, root_notung = job
    
    annotated_tree = None
    
    time_start = time.time()
    
    # [TODO] evaluate use of RAM disk for temporary files
    f1, f1_name = tempfile.mkstemp(dir='./')
    f2, f2_name = tempfile.mkstemp(dir='./')
    try:
        # add info 
        os.write(f1, geneTree.encode())
        os.close(f1)
        os.write(f2, speciesTree.encode())
        os.close(f2)
        
        ### [TODO] separate NOTUNG specific code into separate function
        
        # run notung on tmp files
        subprocess.check_output(
            'java -jar %s'%os.path.join(ds.BINARIES,'Notung-2.9.jar')+
            ' -g %s -s %s'%(f1_name, f2_name)+
            ' --rearrange --speciestag prefix --threshold 0.9 --edgeweights name --parsable',
            shell=True)
        
        # find the reconciled files
        notung_tree = f1_name + ".rearrange.0"
        notung_parsable = notung_tree + ".parsable.txt"
        
        if os.path.exists(notung_tree):
            # extract the relevant information
            annotated_tree = annotate_notung(notung_tree,notung_parsable)
            # 2. event dictionary [event][node.name]=[event.details]
            
            # clean up
            os.remove(notung_tree)
            os.remove(notung_parsable)
        else:
            # reconciliation failed
            pass
        
        ### NOTUNG END
    finally:
        os.remove(f1_name)
        os.remove(f2_name)
        
    time_end = time.time()
    
    if annotated_tree:
        if return_tree:
            split_index = compute_split_index(annotated_tree, True)
            return nog_id,(annotated_tree,split_index)
        else:
            root_event = annotated_tree.event
            split_index = compute_split_index(annotated_tree)
            return nog_id,(root_event,split_index) # [TODO] add time_end - time_start
    else:
        return nog_id,None
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Reconcile batch file")
    parser.add_argument("batch_file")
    parser.add_argument("-c","--cpu",type=int, default=1)
    parser.add_argument("-s","--test_species_trees", action="store_true", default=False)
    parser.add_argument("-l","--limit_computation", type=int, default=-1)
    parser.add_argument("-t","--test_annotation", action="store_true", default=False)
    parser.add_argument("-o","--output_trees", action="store_true", default=False)
    args = parser.parse_args()
    
    assert os.path.exists(args.batch_file), "No batch file found at: %s"%args.batch_file
    lim = args.limit_computation
    reconciliation_jobs = []
    
    eggNOG_speciesTree = Tree(EGGNOGv4_SPECIES_TREE)
    
    with open(args.batch_file) as f:
        c = 1
        for line in f:
            l = line.rstrip().split('\t')
            assert len(l) == 3, 'Unexpected split length: %s'%l
            nog_id = l[0]
            timing = l[1]
            geneTree = l[2]
            
            # prune eggNOG species tree to fit sample gene tree
            speciesTree = prune_species_tree(geneTree,eggNOG_speciesTree)
            
            if args.output_trees:
                with open('%d.geneTree.nw'%c,'w') as f:
                    f.write(geneTree)
                with open('%d.speciesTree.nw'%c,'w') as f:
                    f.write(speciesTree)
            c+=1
            
            # define reconciliation job
            if args.test_species_trees:
                job = (nog_id,timing,speciesTree)
            else:
                job = (nog_id,geneTree,speciesTree)
            reconciliation_jobs.append(job)
            
            # limit jobs to compute
            if len(reconciliation_jobs) == lim:
                break
    
    if args.test_species_trees:
        for job in reconciliation_jobs:
            nog_id, timing, speciesTree = job
            print("%s\t%s\t%s"%(nog_id,timing,speciesTree))
    else:
        # schedule workers
        if args.cpu > 1:
            p = Pool(args.cpu)
            reconciliations = p.map(reconcile,reconciliation_jobs)
            
            # wait for all workers to finish
            p.close()
            p.join()
        else:
            reconciliations = [ reconcile(
                job,args.test_annotation) for job in reconciliation_jobs]
        
        for nog_id, results in reconciliations:
            if args.test_annotation:
                annotated_tree,split_index = results
                print('%s\t%s\t%.3f'%(nog_id,annotated_tree.event,split_index))
                print(annotated_tree.get_ascii(
                    attributes=["name","event","split"]))
            else:
                root_event,split_index = results
                print("%s\t%s\t%.5f"%(nog_id, root_event, split_index))
