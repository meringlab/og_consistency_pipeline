#!/usr/bin/env python3

import sys
import time
from multiprocessing import Pool

from tqdm import tqdm

from methods import tree_building


def read_samples(tsv_path):
    tree_jobs = []
    with open(tsv_path) as f:
        for line in f:
            l = line.rstrip().split()
            assert len(l) == 3, "Found unexpected number of elements in %s"%line
            inconsistent_nog = int(l[0])
            nog_sample_no = int(l[1])
            nog_fasta_samples = eval(l[2])
            tree_job = (inconsistent_nog,nog_sample_no,nog_fasta_samples)
            tree_jobs.append(tree_job)
    return tree_jobs

def write_trees(tsv_path,tree_computations):
    with open(tsv_path, 'w') as f:
        for nog_id,sample_no,tree_nw in tree_computations:
            f.write('%d\t%d\t%s\n'%(nog_id,sample_no,tree_nw))

def build_trees(higher_level,input_samples, tree_method,
                alignment_software,tree_software, keep_polytomies,
                root_notung, cpu_cores, output_unrooted, output_rooted):
    
    tree_building.MAFFT_WEB = alignment_software
    tree_building.FASTTREE_WEB = tree_software
    
    tree_jobs = read_samples(input_samples)
    
    if not tree_jobs:
        # no tree jobs found => create two empty files for snakemake
        open(output_unrooted, 'a').close()
        open(output_rooted, 'a').close()
        return
    
    if tree_method == 'ete':
        tree_method = tree_building.ete_build_workflow
    elif tree_method == 'website':
        tree_method = tree_building.direct_fast
    elif tree_method == 'website_full':
        tree_method = tree_building.direct_full
    elif tree_method == 'cluster':
        pass
    else:
        sys.stderr.write("Unknown tree building method: %s\n"%tree_method)
        sys.exit()
        
    sys.stderr.write('Starting %d tree computations..\n'%len(tree_jobs))
            
    start = time.time()
    
    if tree_method == 'cluster':
        # 1. submit <- tree_jobs
        tree_building.submit_cluster_arrayTask(higher_level,tree_jobs)
        
        # 2. collect -> tree_computations
        tree_computations = tree_building.collect_cluster_results(higher_level)
        
    elif cpu_cores > 1:
        with Pool(cpu_cores) as p:
            tree_computations = list(tqdm(p.imap(tree_method, tree_jobs, 10),total=len(tree_jobs)))
    else:
        tree_computations = [tree_method(x) for x in tqdm(tree_jobs)]
    
    stop = time.time()
    total = stop - start
    avg_tree_time = total * cpu_cores / len(tree_jobs)
    
    sys.stderr.write('Generated %d trees in %.2f seconds with %d cores\n'%(
        len(tree_jobs), total, cpu_cores))
    sys.stderr.write('Avg.comp.time: %.5f\n'%avg_tree_time)
    
    write_trees(output_unrooted,tree_computations)
    
    ## 3b. rerooting
    sys.stderr.write('Applying mid-point rooting to all %d tree computations..\n'%len(tree_jobs))
    if not root_notung:
        tree_computations = [tree_building.reroot_trees(x,keep_polytomies) for x in tqdm(tree_computations)]

    write_trees(output_rooted,tree_computations)
    
if __name__ == '__main__':
    build_trees(higher_level=int(snakemake.wildcards.level_id),
                input_samples=snakemake.input.samples,
                tree_method=snakemake.params.tree_method,
                alignment_software=snakemake.input.alignment_software,
                tree_software=snakemake.input.tree_software,
                keep_polytomies=snakemake.params.keep_polytomies,
                root_notung=snakemake.params.root_notung,
                cpu_cores=snakemake.threads,
                output_unrooted=snakemake.output.trees_unrooted,
                output_rooted=snakemake.output.trees_rooted)

