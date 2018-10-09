#!/usr/bin/env python3
import os
import sys
import time
import random
from tqdm import tqdm

from methods import hgt_utils, expansion, sampling
from methods.utils import eggNOG_utils as eu


def write_tree_tasks(path_tsv,tree_jobs):
    with open(path_tsv,'w') as f:             
        for nog_id, sample_counter, fasta_sequences in tree_jobs:
            f.write("%s\t%d\t%s\n"%(nog_id,sample_counter,repr(fasta_sequences)))

def write_default_solutions(path_tsv, inconsistencies):
    with open(path_tsv,'w') as f:
        for nog_id in sorted(inconsistencies):
            for inconsistent_nog in sorted(inconsistencies[nog_id]):
                solution = inconsistencies[nog_id][inconsistent_nog]
                if solution:
                    if solution['decision'] == 'split':
                        root = 'D'
                    else:
                        root = 'S'
                    
                    f.write('%d\t%d\t%s\t%.1f\n'%(
                        nog_id, inconsistent_nog, root, solution[root][0]))

def write_inconsistencies(path_tsv, inconsistencies):
    with open(path_tsv,'w') as f:
        for nog_id in sorted(inconsistencies):
            for inconsistent_nog in sorted(inconsistencies[nog_id]):
                f.write('%d\t%d\n'%(nog_id,inconsistent_nog))

def expand_and_sample(higher_level,output_dir,
                      input_definition,input_fasta,input_clades,
                      sample_no,sample_size,sample_method,random_seed,
                      default_action, tree_limit, verbose,
                      tree_tsv, default_tsv, inconsistent_tsv):
    
    eggNOG_children = eu.read_eggNOG_treeRev()
    eggNOG_names = eu.read_eggNOG_names()
    eggNOG_children_ids = {y:x for x,y in eggNOG_names.items() if x in eggNOG_children}
    eggNOG_level_species = eu.read_level_species()
    
    # pre-loading data
    protein_names = hgt_utils.load_eggNOG_protein_names_pickle()
    protein_fasta = hgt_utils.get_level_fasta(higher_level, input_fasta)

    # define sampling parameters
    random.seed(random_seed)
    sampler = sampling.InconsistencySampler(sample_no,sample_size,sample_method,
                                            protein_fasta)
    
    # load og definition
    protein_nogs, nog_proteins = hgt_utils.load_join_data(higher_level, output_dir, input_definition)

    # save nogs you already encountered
    discovered_nogs = set()
    
    # dictionary of inconsistencies to process
    inconsistencies = {}
    
    tree_jobs = []
    sample_counter = 0
    total_pre_splits = 0
    total_pre_merges = 0
    
    for nog_level in eggNOG_children[higher_level]:
        
        sys.stderr.write("Started expanding %s [%s] @ %s\n"%(
            nog_level,eggNOG_names[nog_level],time.strftime("%Y%m%d_%H%M%S",time.localtime())))
        
        nog_names = sorted(nog_proteins[nog_level])
        
        sys.stderr.write('Analyzing %d nogs for consistency\n'%len(nog_names))
        sys.stderr.write('Covering a total of %d proteins\n'%len(protein_nogs[nog_level]))  

        for nog_name in tqdm(nog_names):
            if nog_name not in discovered_nogs:
                
                # un-comment for testing
                if tree_limit > 0 and len(tree_jobs) > tree_limit:
                    break
                
                ############## 1. NOG EXPANSION ##############
                
                # compute a new exNOG in master script and pass
                exNOG = expansion.ExpandedNOG(str(nog_level),str(nog_name),verbose=verbose)
                exNOG.set_cache(protein_nogs, nog_proteins, protein_names)
                exNOG.build()#complete_expansion = True)
        
                discovered_nogs.update(exNOG.get_expansion_keys())
                
                nog_id = exNOG.input_id
                if exNOG.has_inconsistencies():
                    inconsistencies[nog_id] = {}
                    # exNOG.write_dependency('%d.graph'%nog_id)
                    
                    has_cogs = False
                    cog_intersection = set()
                    # if higher_level in [2,2759,2157]:
                    #     # do not merge COG/KOGs
                    #     cog_intersection = cogs & exNOG.higher_nogs.keys()
                    #     has_cogs = len(cog_intersection) > 1
                            
                    # reduce problem complexity by pre-splitting based on threshold
                    # cogs are never split by design
                    pre_splits = exNOG.pre_split(cogs=cog_intersection)
                    assert exNOG.has_inconsistencies(), 'Presplitting of %d failed'%nog_id
                    assert len(pre_splits) % 2 == 0
                    total_pre_splits += len(pre_splits) / 2
                    
                    # pre-aggregate singletons and smaller NOGs
                    pre_merges = exNOG.pre_merge()
                    total_pre_merges += len(pre_merges)
                    if not exNOG.has_inconsistencies():
                        inconsistencies[nog_id][nog_id] = {'decision':'merge','S':[-9.0], 'D':[]}
                        
                    if has_cogs:
                        # update the cog_ids in case they were changed by the pre_merge step    
                        current_cogs = {exNOG.get_cog_id_correspondence(x) for x in cog_intersection} 
                
                ############## 2. INCONSISTENT NOG SAMPLING ##############
                
                for inconsistent_nog in exNOG.find_inconsistencies():
                    
                    if default_action:
                        if default_action == 'split':
                            nog_solution = {'decision':default_action, 'S':[], 'D':[-1.0]}
                        elif default_action == 'merge':
                            nog_solution = {'decision':default_action, 'S':[-1.0], 'D':[]}
                        elif default_action == 'random':
                            merge = bool(random.getrandbits(1))
                            if merge:
                                nog_solution = {'decision':'merge', 'S':[-1.0], 'D':[]}
                            else:
                                nog_solution = {'decision':'split', 'S':[], 'D':[-1.0]}
                        else:
                            sys.stderr.write('Unknown default action "%s", please review'%default_action)
                            sys.exit()
                        inconsistencies[nog_id][inconsistent_nog] = nog_solution
                        continue
                    
                    # check if problem contains at least 2 non-sigleton higher nogs
                    # otherwise perform simple merge operation
                    higher_nogs = exNOG.lower_nogs[inconsistent_nog]
                    higher_nogs_size = { x:len(exNOG.expanded_nog[x]) for x in higher_nogs }
                    plus2_sizes = [ size for size in higher_nogs_size.values() if size > 1 ]
                    
                    if len(plus2_sizes) < 2:
                        nog_solution = {'decision':'merge', 'S':[-2.0 + len(plus2_sizes)], 'D':[]} # -2.0 for 0; -1.0 for 1
                        inconsistencies[nog_id][inconsistent_nog] = nog_solution
                        continue
                    
                    if has_cogs:
                        # check how many of the higher nogs are COGs
                        if len(higher_nogs & current_cogs) > 1:
                            # do not allow COGs to be merged
                            nog_solution = {'decision':'split','S':[],'D':[2.0]}
                            inconsistencies[nog_id][inconsistent_nog] = nog_solution
                            continue
                    
                    # Sample inconsistency
                    proteins, levels, splits, paralogs = exNOG.get_composition(inconsistent_nog,input_clades)     
                    
                    # levels check, lca must be direct child of higher_level
                    # empty levels dictionary implies a disconnected child_level
                    if len(levels) < 2:
                        nog_solution = {'decision':'merge','S':[-6.0],'D':[]}
                        inconsistencies[nog_id][inconsistent_nog] = nog_solution
                        continue
                    
                    # species check (sample needs at least 2 species)
                    species = set()
                    for level_species in levels.values():
                        species.update(level_species)
                    if len(species) < 2:
                        nog_solution = {'decision':'split','S':[],'D':[1.0]}
                        inconsistencies[nog_id][inconsistent_nog] = nog_solution
                        continue
                    
                    if verbose: sys.stderr.write('%d: sampling %d: '%(nog_id,inconsistent_nog))
                    samples = sampler.sample_inconsistency(proteins,
                                                           levels,
                                                           splits,
                                                           paralogs)
                    if verbose: sys.stderr.write('%d samples\n'%len(samples))
                    
                    for fasta_sequences in samples:
                        tree_job = (inconsistent_nog, sample_counter, fasta_sequences)
                        tree_jobs.append(tree_job)
                        sample_counter += 1

                    inconsistencies[nog_id][inconsistent_nog] = {}
    
    sys.stderr.write('Completed expansion with a total of %d pre-splits and %d pre-merges\n'%(
        total_pre_splits,total_pre_merges))
    
    # save tree fasta samples and default solutions
    write_tree_tasks(tree_tsv,tree_jobs)
    write_default_solutions(default_tsv,inconsistencies)
    write_inconsistencies(inconsistent_tsv,inconsistencies)

if __name__ == '__main__':
    expand_and_sample(higher_level=int(snakemake.wildcards.level_id),
                      output_dir=snakemake.config['output_dir'],
                      input_definition=os.path.dirname(snakemake.input.parent),
                      input_fasta=snakemake.config['fasta_dir'],
                      input_clades=snakemake.config['clades_dir'],
                      random_seed=snakemake.params.random_seed,
                      sample_no=snakemake.params.sample_no,
                      sample_size=snakemake.params.sample_size,
                      sample_method=snakemake.params.sample_method,
                      default_action=snakemake.params.default_action,
                      tree_limit=snakemake.params.tree_limit,
                      verbose=snakemake.params.verbose,
                      tree_tsv=snakemake.output.samples,
                      default_tsv=snakemake.output.default_solutions,
                      inconsistent_tsv=snakemake.output.inconsistencies)
