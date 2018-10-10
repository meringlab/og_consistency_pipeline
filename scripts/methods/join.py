#!/usr/bin/env python3
# module to join solutions from the reconciliation procedure to make the
# expanded NOGs consistent

import time
from .utils import eggNOG_utils as eu
from collections import defaultdict, Counter

NOG_PREFIXES = {
    'old':1,
    'old_singleton':2,
    'new_merge':3,
    'new_singleton':4,
    'old_split':5,
    'new_split':6 # [TODO] description
}

def get_new_id(nog_type,level_id,new_nog_counter):
    assert nog_type in NOG_PREFIXES
    
    new_nog_counter[level_id] += 1
    level_counter = new_nog_counter[level_id]
    new_nog_id = eu.get_nog_id(level_id, level_counter, NOG_PREFIXES[nog_type])

    return new_nog_id

def apply_solutions(job, debug=False):
    
    exNOG, solutions, cog_intersection = job
    eggNOG_tree = exNOG.eggNOG_dict
    higher_level = exNOG.higher_level

    new_singletons = defaultdict(set)
    new_nog_counter = Counter() 
    complex_decisions = {}
    
    # 0a. reduce exNOG complexity by splitting weak connections in advance (<transfer_thr)
    split_results = exNOG.pre_split(cogs=cog_intersection)
    for nog_id, nog_proteins in split_results.items():
        if len(nog_proteins) > 1:
            new_nog_counter[higher_level] += 1
        else:
            singleton_protein_id = int(list(nog_proteins)[0])
            new_singletons[higher_level].add(singleton_protein_id)
    
    # 0b. reduce singleton and smaller nog scattering by pre-merging
    merge_results = exNOG.pre_merge()
    new_nog_counter[higher_level] += len(merge_results)
    
    if len(solutions) == 1 and exNOG.input_id in solutions:
        if solutions[exNOG.input_id]['S'] == -9.0:
            solutions = {}
    
    if debug:
        if split_results:
            print("%d: applied pre-splits %s"%(exNOG.input_id,exNOG.get_nog_id_changes()))
        
        if merge_results:
            print("%d: applied pre-merges %s"%(exNOG.input_id,exNOG.get_nog_id_changes()))
        
        print("Applying simple solutions:")
        t_start = time.clock()
    
    # 1. Apply singleton merges (simple decisions) and collect complex decisions
    for inconsistent_nog,solution in solutions.items():
        merge_values = solution['S']
        if any([ x < 0 for x in merge_values ]):
            new_nog_id = get_new_id("new_merge",higher_level, new_nog_counter)
            exNOG.merge(inconsistent_nog,new_nog_id)
        else:
            # record reconciliation decision, i.e. only normal split or merge
            complex_decisions[inconsistent_nog] = solution
    
    if debug:
        print(time.clock()-t_start)
        print("Applying complex solutions:")
        t_start = time.clock()
    
    # 2. Define order of application for complex decisions         
    if len(complex_decisions) > 1 and debug:
        # [TODO] one coherent action should be taken (merge or split)
        # i.e. all remaining reconciliations must agree
        print(len(complex_decisions), exNOG.input_id)
        
    #else:
    decision_order = sorted(complex_decisions.keys())
    
    # 3. Apply final split/merges if there are
    for inconsistent_nog in decision_order:
        
        solution = complex_decisions[inconsistent_nog]
        lower_level = eu.decompose_nog_id(inconsistent_nog)[0]
        # [TODO] check if still needed: higher_level = eggNOG_tree[lower_level]
        
        if solution['decision'] == 'split':
            new_nog_id = get_new_id("new_split",lower_level, new_nog_counter)
            splitted_nog = exNOG.split(inconsistent_nog, new_nog_id)
            
            # filter out new singletons and update counter (i.e. only one id was created)
            for nog_id, nog_proteins in splitted_nog.items():
                if len(nog_proteins) == 1:
                    singleton_protein_id = int(list(nog_proteins)[0])
                    new_singletons[lower_level].add(singleton_protein_id)
                else:
                    if nog_id != new_nog_id:
                        new_nog_counter[lower_level] += 1
        else:
            new_nog_id = get_new_id("new_merge",higher_level, new_nog_counter)
            exNOG.merge(inconsistent_nog, new_nog_id)
    
    if debug:
        print(time.clock()-t_start)
        print("Propagating solutions:")
        
    # 4. Check that the reconciliation has introduced new inconsistencies
    #    in lower levels. If yes, propagate split until consistency is met
    
    level_iter = 0
    new_inconsistencies = exNOG.find_inconsistencies()
    while new_inconsistencies:
        
        if debug:
            t_start = time.clock()
            print("\nit.%d (%d): "%(level_iter,len(new_inconsistencies)), end=' ')
            level_iter += 1
        
        for inconsistent_nog in new_inconsistencies:
            
            # split inconsistent_nog
            level_id, nog_number = eu.decompose_nog_id(inconsistent_nog)
            new_nog_id = get_new_id("new_split", level_id, new_nog_counter)
            splitted_nog = exNOG.split(inconsistent_nog, new_nog_id)
            
            # filter out new singletons and update counter (i.e. only one id was created)
            for nog_id, nog_proteins in splitted_nog.items():
                if len(nog_proteins) == 1:
                    singleton_protein_id = int(list(nog_proteins)[0])
                    new_singletons[level_id].add(singleton_protein_id)
                else:
                    if nog_id != new_nog_id:
                        new_nog_counter[level_id] += 1
        
        if debug:
            print(time.clock()-t_start, end=' ')
            t_start = time.clock()
        
        new_inconsistencies = exNOG.find_inconsistencies()
        
        if debug:
            print(' (r.%.3f)'%(time.clock()-t_start), end=' ')
    
    # 5. Obtain the new consistent definition of the expanded NOG without singletons     
    nog_definitions = exNOG.get_protein_mapping_wo_singletons(int) #[TODO] int ref
    
    return (exNOG.input_id, nog_definitions, new_singletons, exNOG.get_nog_id_changes())