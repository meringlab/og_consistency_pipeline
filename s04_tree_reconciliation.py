#!/usr/bin/env python3

import sys
import time
from multiprocessing import Pool

from tqdm import tqdm
from ete3 import Tree

from methods import reconciliation
from methods.utils import eggNOG_utils as eu

def read_trees(tree_tsv):
    tree_computations = []
    with open(tree_tsv) as f:
        c = 0
        for line in f:
            c+=1
            l = line.rstrip().split()
            if len(l) == 3:
                nog_id = int(l[0])
                tree_time = float(l[1])
                tree_nw = l[2]
                tree_computations.append((nog_id,tree_time,tree_nw))
            else:
                sys.stderr.write("Non standard line @ %d: %s\n"%(c,line))
                sys.stderr.write('Completed loading of %d tree compuations from %s\n'%(
                    len(tree_computations),args.input_trees))
    return tree_computations

def write_reconciliations(reconciliations_tsv,reconciliations):
    with open(reconciliations_tsv,'w') as f:
        for nog_id, result in reconciliations:
            assert result, 'No reconciliation found for %d'%nog_id
            root_event,split_index = result
            f.write('%d\t%s\t%.3f\n'%(nog_id,root_event,split_index))

def reconcile_trees(higher_level,input_trees,
                    computation_method,cpu_cores,
                    keep_polytomies,root_notung,infer_transfers,
                    output_reconciliations):
    
    tree_computations = read_trees(input_trees)
    eggNOG_level_species = eu.read_level_species()
    
    # species tree generation
    eggNOG_speciesTree = Tree(reconciliation.EGGNOGv4_SPECIES_TREE)
    # do intial pruning to exclude all non-euNOG species
    higher_node = eggNOG_speciesTree.get_common_ancestor([str(x) for x in eggNOG_level_species[higher_level]])
    eggNOG_speciesTree = higher_node.detach()
    
    sys.stderr.write('Generating species trees reconciliation...\n')
    if cpu_cores > 1:
        cached_jobs = [(x,y,z,eggNOG_speciesTree,keep_polytomies,root_notung) for x,y,z in tree_computations]
        with Pool(cpu_cores) as p:
            reconciliation_jobs = list(tqdm(
                p.imap(reconciliation.prepare_reconciliation_job,cached_jobs),total=len(cached_jobs)))
    else:
        reconciliation_jobs = []
        for nog_id,sample_no,tree_nw in tqdm(tree_computations):
            species_nw = reconciliation.prune_species_tree(tree_nw,eggNOG_speciesTree,keep_polytomies)
            job = (nog_id,sample_no,tree_nw,species_nw,root_notung)
            reconciliation_jobs.append(job)

    sys.stderr.write('Starting reconciliation for %d jobs...\n'%len(reconciliation_jobs))
    reconciliation_method = reconciliation.reconcile
    
    # reconciliation
    if computation_method == 'cluster':
        reconciliation.submit_taskArray(higher_level,reconciliation_jobs,
                                        keep_polytomies,root_notung,infer_transfers)
        reconciliations = reconciliation.collect_taskArray(higher_level,cpu_cores=cpu_cores)
    elif cpu_cores > 1:
        with Pool(cpu_cores) as p:
            reconciliations = list(tqdm(
                p.imap(reconciliation_method,reconciliation_jobs),total=len(reconciliation_jobs)))
    else:
        reconciliations = [reconciliation_method(x) for x in tqdm(reconciliation_jobs)]
        
    # flag incomplete reconciliations
    to_delete = []
    for i in range(len(reconciliations)):
        if reconciliations[i] is None:
            to_delete.append(i)
        else:
            nog_id, result = reconciliations[i]
            if not result:
                # reconciliation failed
                reconciliations[i] = (nog_id,('S',-5.0))

    # delete incomplete 
    if to_delete:
        for i in sorted(to_delete, reverse=True):
            sys.stderr.write('Deleting reconciliation entry %d.%d because empty\n'%(higher_level,i))
            del reconciliations[i]
    
    write_reconciliations(output_reconciliations,reconciliations)   

if __name__ == '__main__':
    
    reconcile_trees(higher_level=int(snakemake.wildcards.level_id),
                    input_trees=snakemake.input.trees,
                    computation_method=snakemake.params.computation_method,
                    cpu_cores=snakemake.threads,
                    keep_polytomies=snakemake.params.keep_polytomies,
                    root_notung=snakemake.params.root_notung,
                    infer_transfers=snakemake.params.infer_transfers,
                    output_reconciliations=snakemake.output.reconciliations)
    