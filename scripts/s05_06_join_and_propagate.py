#!/usr/bin/env python3

import os
import sys
import time
from collections import Counter,defaultdict,OrderedDict
from multiprocessing import Pool

from tqdm import tqdm
from methods import expansion, join, test, hgt_utils
from methods.utils import eggNOG_utils as eu

def read_reconciliations(path_tsv):
    reconciliations = []
    with open(path_tsv) as f:
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
            
        if root_event != 'D':
            # hack for C (co-Divergence) and T (Transfer)
            root_event = 'S'
            
        reconciliation_dict[nog_id][root_event].append(split_index)
    
    return reconciliation_dict

def read_inconsistencies(inconsistent_tsv, default_tsv):
    inconsistencies = {}
    with open(inconsistent_tsv) as f:
        for line in f:
            nog_id, inconsistent_id = [int(x) for x in line.rstrip().split()]
            if nog_id in inconsistencies:
                inconsistencies[nog_id][inconsistent_id] = {}
            else:
                inconsistencies[nog_id] = {inconsistent_id:{}}
    
    # populate default decisions
    with open(default_tsv) as f:
        for line in f:
            nog_id, inconsistent_id, root_event, value = line.rstrip().split()
            nog_id = int(nog_id)
            inconsistent_id = int(inconsistent_id)
            value = float(value)
            assert nog_id in inconsistencies
            assert inconsistent_id in inconsistencies[nog_id]
            
            # form solution vector
            solution = {'decision':'','S':[],'D':[]}
            solution[root_event].append(value)
            if root_event == 'S':
                solution['decision'] = 'merge'
            else:
                assert root_event == 'D'
                solution['decision'] = 'split'
            
            inconsistencies[nog_id][inconsistent_id] = solution
    
    return inconsistencies

def write_protein_nogs(output_dir,protein_nogs):
    for level_id in protein_nogs:
        level_tsv = os.path.join(output_dir,'%d.tsv'%level_id)
        sys.stderr.write('Writing %s\n'%level_tsv)
        with open(level_tsv,'w') as f:
            for protein_id, nog_id in sorted(protein_nogs[level_id].items(),
                                             key=lambda x: '%d %d'%(x[1],x[0])):
                f.write('%d\t%d\n'%(nog_id,protein_id))
        
        # pickle version
        hgt_utils.save_eggNOG_nog_mapping(level_id,protein_nogs[level_id],output_dir)
                
def write_singletons(path_tsv, new_singletons):
    with open(path_tsv,'w') as f:
        for level_id in new_singletons: # TODO sort
            for nog_id in new_singletons[level_id]: # TODO sort
                f.write('%d\t%d\n'%(level_id,nog_id))

def join_solutions(higher_level, output_dir, input_definition,
                   input_reconciliations, input_inconsistencies, input_default,
                   cpu_cores, majority_vote_threshold,
                   output_consistent, output_singletons):
    
    eggNOG_children = eu.read_eggNOG_treeRev()
    eggNOG_names = eu.read_eggNOG_names()
    eggNOG_children_ids = {y:x for x,y in eggNOG_names.items() if x in eggNOG_children}
    eggNOG_level_species = eu.read_level_species()
    
    # pre-loading data
    protein_names = hgt_utils.load_eggNOG_protein_names_pickle()

    # load og definition
    protein_nogs, nog_proteins = hgt_utils.load_join_data(higher_level, output_dir, input_definition)
    
    reconciliation_dict = read_reconciliations(input_reconciliations)
    inconsistencies = read_inconsistencies(input_inconsistencies,input_default)
    
    consistency_jobs = []
    for nog_id in inconsistencies:
                
        # rebuild complete exNOG (all-sublevels)
        exNOG = expansion.ExpandedNOG(nog_id)
        exNOG.set_cache(protein_nogs, nog_proteins, protein_names)
        exNOG.build(complete_expansion=True)
        exNOG.remove_cache() # remove big objects for jobs
        
        cog_intersection = set()
        # if higher_level in [2,2759,2157]:
        #     # do not merge COG/KOGs
        #     cog_intersection = cogs & exNOG.higher_nogs.keys()
                
        # extract relevant reconciliations
        for inconsistent_nog in inconsistencies[nog_id]:
            if 'decision' not in inconsistencies[nog_id][inconsistent_nog]:
                # no solution (i.e. no simple merge)
                # assert inconsistent_nog in reconciliation_dict, 'No reconciliation for %d'%inconsistent_nog
                if inconsistent_nog not in reconciliation_dict:
                    print('No reconciliation for %d'%inconsistent_nog)
                    continue
                
                nog_solution = reconciliation_dict[inconsistent_nog]

                tot_solutions = len(nog_solution['D']) + len(nog_solution['S'])
                duplication_threshold = majority_vote_threshold * tot_solutions

                if len(nog_solution['D']) > duplication_threshold:
                    nog_solution['decision'] = 'split'
                else:
                    nog_solution['decision'] = 'merge'
                
                inconsistencies[nog_id][inconsistent_nog] = nog_solution
            
        # define jobs as the exNOG w/o cache and the relative solution dict
        job = (exNOG, inconsistencies[nog_id], cog_intersection)
        consistency_jobs.append(job)
            
    # reconcilations -> join
    sys.stderr.write('Starting joining %d reconciliations...\n'%len(consistency_jobs))
    start = time.time()
        
    if cpu_cores > 1:
        with Pool(cpu_cores) as p:
            consistent_nogs = list(tqdm(
                p.imap(join.apply_solutions,consistency_jobs),total=len(consistency_jobs)))
    else:
        consistent_nogs = [join.apply_solutions(x) for x in tqdm(consistency_jobs)]
    
    stop = time.time()
    total = stop - start
    sys.stderr.write('Completed %d joins in %.2fs \n'%(len(consistency_jobs),total))
        
    new_singletons = defaultdict(set) # dict[level]:set(singletons)
    singleton_reference = defaultdict(dict)
    
    sys.stderr.write('Updating old definition.. \n')
    # 1. remove new singletons for every level in old definiton
    for consistent_id, consistent_definiton, consistent_singletons, nog_history in tqdm(consistent_nogs):
        for level_id in consistent_singletons:
            old_definition = protein_nogs[level_id]
            for protein_id in consistent_singletons[level_id]:
                
                # safety check for double removal of singletons
                if level_id in singleton_reference:
                    if protein_id in singleton_reference[level_id]:
                        sys.stderr.write("Singleton %d was already removed previously by: %d\n"%(
                            protein_id,singleton_reference[level_id][protein_id]))
                        sys.stderr.write("Current removal attempt by: %d\n"%consistent_id)
                        sys.exit()
                        
                assert protein_id in old_definition, "Singleton %d not found @ %d"%(protein_id, level_id)
                del old_definition[protein_id]
                
                # testing double singleton removal
                singleton_reference[level_id][protein_id] = consistent_id
            
            new_singletons[level_id].update(consistent_singletons[level_id])
            
        # 2. update the old protein mapping with consistent nog definitions
        #    and identify the NOGs that need a new id:
        #     - higher level: merges (3) -> (3)[just coordinate new merges]
        #     - lower level: splits (6) -> (5) [previous 5 must be considered!]
        
        # 2a. initialize nog counter to not collide with previous nog ids
        nog_counter = Counter()
        for level_id in nog_proteins:
            
            for nog_id in nog_proteins[level_id]:
                if str(nog_id).startswith('5'):
                    
                    if level_id not in nog_counter:
                        nog_counter[level_id] = 0
                        
                    _, level_counter = eu.decompose_nog_id(nog_id)
                    if nog_counter[level_id] < level_counter:
                        nog_counter[level_id] = level_counter
    
        assert higher_level not in nog_counter, "Found higher level %d in counter: %s"%(
            higher_level,nog_counter[higher_level])
        
        # 3a. assign new nog_ids and substitute proteins
        level_history = defaultdict(OrderedDict)
        
        for consistent_id, consistent_definiton, consistent_singletons, nog_history in consistent_nogs:
            nog_mapper = {}
            for level_id in consistent_definiton:
                
                old_definition = protein_nogs[level_id]
                
                if level_id == higher_level:
                    prefix = {'3':3,'6':5}
                else:
                    prefix = {'6':5}
                    
                for protein_id in sorted(consistent_definiton[level_id]):
                    nog_id = consistent_definiton[level_id][protein_id]
                    nog_str = str(nog_id)
                    if nog_str[0] in prefix:
                        if nog_id not in nog_mapper:
                            nog_counter[level_id] += 1
                            new_nog_id = eu.get_copy_id(nog_id,nog_counter[level_id],prefix[nog_str[0]])
                            nog_mapper[nog_id] = int(new_nog_id)
                                
                        old_definition[protein_id] = nog_mapper[nog_id]
            
            # update the ids in the nog_history
            for old_id, new_ids in nog_history.items():
                level_history[consistent_id][old_id] = [nog_mapper[x] if x in nog_mapper else x for x in new_ids]
    
    # write level id history to file
    # history_path = '%d.nog_history.tsv'%higher_level
    # with open(history_path,'w') as f:
    #     for consistent_id in level_history:
    #         for old_id, new_ids in level_history[consistent_id].items():
    #             f.write('%d\t%d\t%s\n'%(consistent_id,old_id,','.join([str(x) for x in new_ids])))
    
    sys.stderr.write('Writing new definition.. \n')
    write_protein_nogs(output_consistent, protein_nogs)
    write_singletons(output_singletons, new_singletons)
    
    # test consistency
    sys.stderr.write('Testing new definition.. \n')
    lower_levels = list(eggNOG_children[higher_level])
    test.test_consistency([higher_level],lower_levels,protein_nogs)
    while lower_levels:
        lower_level = lower_levels.pop()
        if lower_level in eggNOG_children:
            sub_levels = list(eggNOG_children[lower_level])
            test.test_consistency([lower_level],sub_levels,protein_nogs)
            lower_levels.extend(sub_levels)

if __name__ == '__main__':
    join_solutions(
        higher_level=int(snakemake.wildcards.level_id),
        output_dir=os.path.join(snakemake.config['output_dir'],'new_definition'),
        input_definition=os.path.dirname(snakemake.input.parent),
        input_reconciliations=snakemake.input.reconciliations,
        input_inconsistencies=snakemake.input.inconsistencies,
        input_default=snakemake.input.default_solutions,
        cpu_cores=snakemake.threads,
        majority_vote_threshold=snakemake.params.majority_vote_threshold,
        output_consistent=os.path.dirname(snakemake.output.consistent_ogs),
        output_singletons=snakemake.output.new_singletons)
