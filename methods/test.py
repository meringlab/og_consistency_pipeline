import os
import sys
from collections import defaultdict

def read_nog_files(nog_file_list):
    
    nog_definitions = defaultdict(set)
    nog_proteins = dict()
    
    for nog_file in nog_file_list:
        with open('new_definition/%d.tsv'%nog_file) as f:
            for line in f:
                nog_id, protein_id = line.rstrip().split()
                nog_id = int(nog_id)
                protein_id = int(protein_id)
                nog_definitions[nog_id].add(protein_id)
                nog_proteins[protein_id] = nog_id
                
    return nog_definitions, nog_proteins

def read_cache(levels,protein_nogs):
    
    nog_definitions = defaultdict(set)
    nog_proteins = dict()
    
    for level_id in levels:
        for protein_id,nog_id in protein_nogs[level_id].items():
            nog_definitions[nog_id].add(protein_id)
            nog_proteins[protein_id] = nog_id
    
    return nog_definitions, nog_proteins

def fast_join(levels, protein_nogs, nog_proteins):
    
    # join multiple levels in one definition
    nog_definitions = {} # nog_id => set(protein_ids)
    protein_definitions = {} # protein_id => nog_id
    
    for level_id in levels:
        nog_definitions.update(nog_proteins[level_id])
        protein_definitions.update(protein_nogs[level_id])
        
    return nog_definitions, protein_definitions


def test_consistency(higher_levels,lower_levels,protein_nogs=None):
    
    # read definitions
    if protein_nogs:
        lower_nogs,lower_proteins = read_cache(lower_levels,protein_nogs)
        higher_nogs,higher_proteins = read_cache(higher_levels,protein_nogs)
    else:
        lower_nogs,lower_proteins = read_nog_files(lower_levels)
        higher_nogs,higher_proteins = read_nog_files(higher_levels)

    for lower_nog,lower_nog_proteins in lower_nogs.items():
        
        higher_nog = -1        
        
        for protein_id in lower_nog_proteins:
            
            assert protein_id in higher_proteins, "No higher level match found for %d [NOG %d]"%(
                protein_id,lower_nog)
            
            if higher_nog == -1:
                # Initialize NOG (i.e. first protein of the lower_nog)
                higher_nog = higher_proteins[protein_id]
            else:
                if higher_nog != higher_proteins[protein_id]:
                    
                    sys.stderr.write(
                        "Inconsistency detected for %d, NOG should be [%d] but is [%d]\n"%(
                        protein_id, higher_nog, higher_proteins[protein_id]))
                    sys.stderr.write("Lower nog: %d\n"%lower_nog)
                    
                    # Write out were individual proteins went
                    for nog_protein in lower_nog_proteins:
                        if nog_protein in higher_proteins: # for singleton mode
                            sys.stderr.write("%d\t%d\n"%(nog_protein,higher_proteins[nog_protein]))
                        else:
                            sys.stderr.write("%d\tsingleton\n"%nog_protein)
                            
                    sys.exit(1)
                    
    sys.stderr.write('Successfully tested consistency between %s <> %s\n'%(
        higher_levels,lower_levels))
    
def find_inconsistencies(higher_levels,lower_levels,protein_nogs,nog_proteins):
    
    inconsistencies = []
    
    # read definitions
    lower_nogs,lower_proteins = fast_join(lower_levels,protein_nogs, nog_proteins)
    higher_nogs,higher_proteins = fast_join(higher_levels,protein_nogs, nog_proteins)
    
    for lower_nog,lower_nog_proteins in lower_nogs.items():
        
        higher_nog = -1        
        
        for protein_id in lower_nog_proteins:
            
            if protein_id not in higher_proteins:
                inconsistencies.append(lower_nog)
                break
            
            if higher_nog == -1:
                # Initialize NOG (i.e. first protein of the lower_nog)
                higher_nog = higher_proteins[protein_id]
            else:
                if higher_nog != higher_proteins[protein_id]:
                    inconsistencies.append(lower_nog)
                    break
                    
    return inconsistencies

def find_singletons(higher_levels,lower_levels,protein_nogs,nog_proteins):
    
    singletons_inconsistencies = []
    
    # read definitions
    lower_nogs,lower_proteins = fast_join(lower_levels,protein_nogs, nog_proteins)
    higher_nogs,higher_proteins = fast_join(higher_levels,protein_nogs, nog_proteins)
    
    for lower_nog,lower_nog_proteins in lower_nogs.items():
        
        higher_nog = -1        
        
        singleton_inconsistency = False
        for protein_id in lower_nog_proteins:
            
            if protein_id not in higher_proteins:
                singleton_inconsistency = True
                continue
            
            if higher_nog == -1:
                # Initialize NOG (i.e. first protein of the lower_nog)
                higher_nog = higher_proteins[protein_id]
            else:
                if higher_nog != higher_proteins[protein_id]:
                    singleton_inconsistency = False
                    break
                
        if singleton_inconsistency:
            singletons_inconsistencies.append(lower_nog)
                    
    return singletons_inconsistencies


def test_all_levels(protein_nogs,nog_proteins,eggNOG_children,test_function=find_inconsistencies):
    higher_level = 2759
    
    # set test for consistency before writing new definition
    lower_levels = list(eggNOG_children[higher_level])
    global_inconsistencies = test_function([higher_level],lower_levels,protein_nogs,nog_proteins)
    print("Found %d inconsistencies @ %d"%(len(global_inconsistencies),higher_level))
    
    while lower_levels:
        lower_level = lower_levels.pop()
        if lower_level in eggNOG_children:
            sub_levels = list(eggNOG_children[lower_level])
            inconsistencies = test_function([lower_level],
                                            sub_levels,
                                            protein_nogs,nog_proteins)
            print("Found %d inconsistencies @ %d"%(len(inconsistencies),lower_level))
            lower_levels.extend(sub_levels)
            global_inconsistencies.extend(inconsistencies)
            
    return global_inconsistencies