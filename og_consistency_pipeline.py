# Launch script to make a hierarchy of Orthologous Groups consistent
# author:   Davide Heller
# email:    davide.heller@imls.uzh.ch
# version:  0.1 [2018-09-13]

from methods import hgt_utils, expansion, sampling
from methods import tree_building, reconciliation, join
from methods import test

from methods.utils import eggNOG_utils as eu

from ete3 import Tree

import time
import sys
import os
import random
import argparse
import subprocess

#from ipdb import set_trace as bp
from collections import defaultdict, OrderedDict
from multiprocessing import Queue, Pool
from collections import Counter

def get_level_order(level_id):
    
    t = Tree(eu.get_eggNOG_nhx(),format=8)
    
    # define order in which to visit the nodes (lower levels first)
    level_node = t.search_nodes(name=str(level_id))
    assert len(level_node) == 1, "Ambiguous level_id %d: %s"%(level_id, level_node)
    level_node = level_node[0]
    
    # traverse subtree by level order and leave out leaves
    level_order = [ int(x.name) for x in level_node.traverse(
        strategy='levelorder') if not x.is_leaf() ]
    
    # return reversed order to build up hierarchy bottom-up
    return list(reversed(level_order))

def write_default_solutions(higher_level, inconsistencies):
    with open('%d.default_solutions.tsv'%higher_level,'w') as f:
        for nog_id in inconsistencies:
            for inconsistent_nog, solution in inconsistencies[nog_id].items():
                if solution:
                    if solution['decision'] == 'split':
                        root = 'D'
                    else:
                        root = 'S'
                    
                    f.write('%d\t%s\t%.1f\n'%(
                        inconsistent_nog, root, solution[root][0]))

def filter_fusions(fusions,fusion_level,children,protein_mapping,nog_mapping):

    to_keep = {}
    to_remove = {}

    for fusion_id, fusion_nogs in fusions.items():
    
        assert fusion_id in protein_mapping[fusion_level],'fusion %d not found @ %d'%(fusion_id,fusion_level)
    
        higher_proteins = {x:nog_mapping[fusion_level][x] for x in fusion_nogs}
    
        for level_id in children[fusion_level]:
            if fusion_id in protein_mapping[level_id]:
                lower_nog = protein_mapping[level_id][fusion_id]
                lower_proteins = nog_mapping[level_id][lower_nog]
                
                overlap = Counter({x:len(lower_proteins & y) for x,y in higher_proteins.items()})
                sorted_by_overlap = overlap.most_common()
                
                to_keep[fusion_id] = sorted_by_overlap[0][0]
                to_remove[fusion_id] = {x[0] for x in sorted_by_overlap[1:]}
                break
        else:
            sorted_by_size = Counter({x:len(y) for x,y in higher_proteins.items()}).most_common()
            to_keep[fusion_id] = sorted_by_size[0][0]
            to_remove[fusion_id] = {x[0] for x in sorted_by_size[1:]}
    
    return to_keep, to_remove        
                
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Queue Reconciliation')
    parser.add_argument("final_level",)
    parser.add_argument("-a1","--input_fasta",help="line format: [nog_id]\t[sample_id]\t[fasta_samples]",default=None)
    parser.add_argument("-a2","--input_trees",help="line format: [nog_id]\t[sample_id]\t[unrooted_nw]",default=None)
    parser.add_argument("-a3","--input_reconciliation",help="line format: [nog_id]\t[root_event]\t[split_index]",default=None)
    parser.add_argument("-o","--output_dir",default=None)
    parser.add_argument("-c","--cpu",type=int, default=1)
    parser.add_argument("-t","--tree_method",type=str, default='website')
    parser.add_argument("-s","--seed",type=int, default=1)
    parser.add_argument("-l","--limit",type=int, default=-1)
    parser.add_argument("-n","--sample_number", help="number of samples per inconsistency", type=int, default=20)
    parser.add_argument("-m","--sample_size", help="minimum number of sequences per sample", type=int, default=10)
    parser.add_argument("--sample_method",help="sampling strategy to use",type=str,default='combined')
    parser.add_argument("-d","--dry_run", help="don't write out files", action="store_true", default=False)
    parser.add_argument("-i","--input_definition", help="supply input definition to use instead of 4.0", default=None)
    parser.add_argument("-k","--keep_polytomies",help="don't remove polytomies from the species tree", action="store_true", default=False)
    parser.add_argument("--default_action",help="apply a always default action (merge or split)",default=None)
    parser.add_argument("-r","--root_notung",help="perform the rooting with NOTUNG instead of ete",action="store_true",default=False)
    parser.add_argument("-v","--majority_vote_threshold",help="bias D/S vote",type=float,default=0.5)
    parser.add_argument('--infer_transfers', help="use default DTL model for reconciliation", action="store_true", default=False)
    parser.add_argument("--verbose", help="expansion verbosity", action="store_true", default=False)
    args = parser.parse_args()
    
    if args.infer_transfers:
        # DTL model requires no polytomies in geneTree, polytomies in speciesTree can therefore be kept
        args.keep_polytomies = True
    
    sys.stderr.write('Parameters for run started at %s\n'%time.strftime("%Y%m%d_%H%M%S",time.localtime()))
    args_dict = vars(args)
    for parameter_name in sorted(args_dict):
        if args_dict[parameter_name]:
            sys.stderr.write('%s:\t%s\n'%(parameter_name,args_dict[parameter_name]))
    
    assert args.tree_method in ['website','website_full','ete','cluster']
    assert args.majority_vote_threshold < 1.0 and args.majority_vote_threshold > 0.0
    assert args.sample_method in ['combined','random','splits','levels','paralogs'], 'Unknown sampling method %s'%args.sample_method
    
    if args.tree_method == 'cluster':
        assert 'mnt' in os.getcwd(), "Reposition to mnt to launch cluster command"
    
    if args.default_action:
        assert args.default_action in ['merge','split','random'], 'default action %s not available!'%args.default_action
    default_action = args.default_action
    
    if args.output_dir:
        assert os.path.exists(args.output_dir),"Output dir doesn't exist %s"%args.output_dir
        output_dir = args.output_dir
    else:
        output_dir = 'new_definition'
        if not args.dry_run:
            os.mkdir(output_dir)
    
    eggNOG_children = eu.read_eggNOG_treeRev()
    eggNOG_names = eu.read_eggNOG_names()
    eggNOG_children_ids = {y:x for x,y in eggNOG_names.items() if x in eggNOG_children}
    eggNOG_level_species = eu.read_level_species()
    
    if args.final_level.isdigit():
        final_level = int(args.final_level)
        assert final_level in eggNOG_children, '%d is not higher eggnog level!'%final_level
        final_name = eggNOG_names[final_level]
    else:
        final_name = args.final_level
        assert final_name in eggNOG_children_ids, '%s is not higher eggnog level!'%final_name
        final_level = eggNOG_children_ids[final_name]
    
    # check if user supplied specific definition to work on
    if args.input_definition:
        input_definition = args.input_definition
    else:
        input_definition = hgt_utils.V4_PICKLE_DIR
    assert os.path.exists(input_definition), "No input definition found to work on @ %s"%input_definition
    
    sys.stderr.write("Task to make %d[%s] consistent started at %s\n"%(
        final_level,
        final_name,
        time.strftime("%Y%m%d_%H%M%S",time.localtime())))
    
    # pre-loading data
    protein_names = hgt_utils.load_eggNOG_protein_names_pickle()
    protein_fasta = hgt_utils.get_level_fasta(final_level)

    # define sampling parameters
    sample_no = args.sample_number
    sample_size = args.sample_size
    sample_method = args.sample_method
    sampler = sampling.InconsistencySampler(sample_no,sample_size,sample_method,
                                            protein_fasta)
    


    # set random seed for all sessions, source: http://stackoverflow.com/a/11527011
    if args.seed == -1:
        random.seed(None)
    else:
        random.seed(args.seed)
    
    # set cpu cores
    cpu_cores = args.cpu

    level_order = get_level_order(final_level)
    sys.stderr.write("Level order: %s\n"%(
        " -> ".join(["%d[%s]"%(x,eggNOG_names[x]) for x in level_order])))
    
    with open('timings.log','a') as f:
        f.write('#level\t0load\t1exp\t2tree\t3root\t4species\t5recon\t6join\t7merge\n')
    
    for higher_level in level_order:
        
        timings = [time.time()]
        
        higher_path = os.path.join(output_dir,'%d.tsv'%higher_level)
        if os.path.exists(higher_path):
            sys.stderr.write("Found NOG definition for %d [%s], skipping!\n"%(
                higher_level,eggNOG_names[higher_level]))
            continue
        
        sys.stderr.write("Started consistency pipeline for %d [%s] @ %s\n"%(
                    higher_level,eggNOG_names[higher_level],
                    time.strftime("%Y%m%d_%H%M%S",time.localtime())))
        
        # [TODO] need to reload this to the latest definition? every time?
        protein_nogs, nog_proteins = hgt_utils.load_join_data(higher_level, output_dir, input_definition)
        
        if higher_level in [2,2759,2157]:
            sys.stderr.write('Starting to remove fusions from %d\n'%higher_level)
            
            # load fusions, i.e. dict contains the fusion_id:list(fusion_nogs)
            fusions = eu.read_eggNOG_fusions(higher_level) 
            fusions_to_keep, fusions_to_remove = filter_fusions(fusions, higher_level, eggNOG_children, protein_nogs, nog_proteins)
            
            # apply filtering        
            for fusion_id, fusion_nog in fusions_to_keep.items():
                    
                # ensure correct mapping of protein
                protein_nogs[higher_level][fusion_id] = fusion_nog
            
                # remove other associations
                for nog_id in fusions_to_remove[fusion_id]:
                    nog_proteins[higher_level][nog_id].remove(fusion_id)
                    
            # save fusions to be re-added later
            if not args.dry_run:
                with open('%d.removed_fusions.tsv'%higher_level,'w') as f:
                    for fusion_id, fusion_nogs in fusions_to_remove.items():
                        f.write('%d\t%s\n'%(fusion_id,','.join([str(x) for x in fusion_nogs])))
            
            sys.stderr.write('Completed removal of %d fusion proteins\n'%len(fusions_to_remove))
            
            cogs = hgt_utils.load_COG_mapping(higher_level,eggNOG_names[higher_level])
            cogs = set(cogs.values())
            
        timings.append(time.time())

        # save nogs you already encountered
        discovered_nogs = set()
    
        # dictionary of inconsistencies to process
        inconsistencies = {}
        
        tree_jobs = []
        if args.input_fasta:
            assert os.path.exists(args.input_fasta), 'No fasta file found @ %s'%args.input_fasta
            with open(args.input_fasta) as f:
                for line in f:
                    # tree_job = (inconsistent_nog, fasta_sequences)
                    l = line.rstrip().split()
                    assert len(l) == 3, "Found unexpected number of elements in %s"%line
                    inconsistent_nog = int(l[0])
                    nog_sample_no = int(l[1])
                    nog_fasta_samples = eval(l[2])
                    tree_job = (inconsistent_nog,nog_sample_no,nog_fasta_samples)
                    tree_jobs.append(tree_job)
                    
                    if args.limit > 0 and len(tree_jobs) == args.limit:
                            break
                    
            sys.stderr.write('Successfully imported %d samples from %s\n'%(
                len(tree_jobs),args.input_fasta))
            
            # remove flag for successive runs
            args.input_fasta = None
            
        else:
            # [TODO] define sample_fasta(args) or sample_fasta(higher_level...) function
            sample_counter = 0
            total_pre_splits = 0
            total_pre_merges = 0
            for nog_level in eggNOG_children[higher_level]:
                
                sys.stderr.write("Started expanding %s [%s] @ %s\n"%(
                    nog_level,eggNOG_names[nog_level],time.strftime("%Y%m%d_%H%M%S",time.localtime())))
                
                nog_names = sorted(nog_proteins[nog_level])
                
                sys.stderr.write('Analyzing %d nogs for consistency\n'%len(nog_names))
                sys.stderr.write('Covering a total of %d proteins\n'%len(protein_nogs[nog_level]))  
        
                for nog_name in nog_names:
                    if nog_name not in discovered_nogs:
                        
                        # un-comment for testing
                        if args.limit > 0 and len(tree_jobs) > args.limit:
                            break
                        
                        ############## 1. NOG EXPANSION ##############
                        
                        # compute a new exNOG in master script and pass
                        exNOG = expansion.ExpandedNOG(str(nog_level),str(nog_name),verbose=args.verbose)
                        exNOG.set_cache(protein_nogs, nog_proteins, protein_names)
                        exNOG.build()#complete_expansion = True)
                
                        discovered_nogs.update(exNOG.get_expansion_keys())
                        
                        nog_id = exNOG.input_id
                        if exNOG.has_inconsistencies():
                            inconsistencies[nog_id] = {}
                            # exNOG.write_dependency('%d.graph'%nog_id)
                            
                            has_cogs = False
                            cog_intersection = set()
                            if higher_level in [2,2759,2157]:
                                # do not merge COG/KOGs
                                cog_intersection = cogs & exNOG.higher_nogs.keys()
                                has_cogs = len(cog_intersection) > 1
                                    
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
                            proteins, levels, splits, paralogs = exNOG.get_composition(inconsistent_nog)     
                            
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
                            
                            if args.verbose: sys.stderr.write('%d: sampling %d: '%(nog_id,inconsistent_nog))
                            samples = sampler.sample_inconsistency(proteins,
                                                                   levels,
                                                                   splits,
                                                                   paralogs)
                            if args.verbose: sys.stderr.write('%d samples\n'%len(samples))
                            
                            if args.dry_run:
                                tree_jobs.append(inconsistent_nog)
                            else:
                                for fasta_sequences in samples:
                                    tree_job = (inconsistent_nog, sample_counter, fasta_sequences)
                                    tree_jobs.append(tree_job)
                                    sample_counter += 1
    
                            inconsistencies[nog_id][inconsistent_nog] = {}
            
            sys.stderr.write('Completed expansion with a total of %d pre-splits and %d pre-merges\n'%(
                total_pre_splits,total_pre_merges))
            
            # save fasta samples
            if args.input_trees or args.input_reconciliation:
                pass 
            elif args.dry_run:
                #write_default_solutions(higher_level,inconsistencies)
                sys.stderr.write('Completed assembly of %d fasta samples\n'%len(tree_jobs))
                sys.exit()
            else:
                write_default_solutions(higher_level,inconsistencies)
                #tree_building.write_tree_tasks(higher_level,tree_jobs)
        
        timings.append(time.time())        
        ############## 3. TREE BUILDING ##############
        
        if args.input_trees:
            # load precomputed trees
            tree_computations = []
            with open(args.input_trees) as f:
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
            
            # remove flag for successive runs
            args.input_trees = None
            
        elif args.input_reconciliation:
            pass
        elif not tree_jobs:
            sys.stderr.write('WARNING: no tree jobs found to process!\n')
            tree_computations = []
            pass
        else:
            # map the workers to tree tasks via pool
            # if not tree_jobs:
            #     sys.exit()
            
            if args.tree_method == 'ete':
                tree_method = tree_building.ete_build_workflow
            elif args.tree_method == 'website':
                tree_method = tree_building.direct_fast
            elif args.tree_method == 'website_full':
                tree_method = tree_building.direct_full
            elif args.tree_method == 'cluster':
                pass
            else:
                sys.stderr.write("Unknown tree building method: %s\n"%args.tree_method)
                sys.exit()
            
            sys.stderr.write('Starting %d tree computations..\n'%len(tree_jobs))
            
            start = time.time()
            
            if args.tree_method == 'cluster':
                # 1. submit <- tree_jobs
                tree_building.submit_cluster_arrayTask(higher_level,tree_jobs)
                
                # 2. collect -> tree_computations
                tree_computations = tree_building.collect_cluster_results(higher_level)
                
            elif cpu_cores > 1:
                p = Pool(cpu_cores)
                tree_computations = p.map(tree_method, tree_jobs, 10)
    
                # wait for all workers to finish
                p.close()
                p.join()
            else:
                tree_computations = list(map(tree_method, tree_jobs))
            
            stop = time.time()
            total = stop - start
            avg_tree_time = total * cpu_cores / len(tree_jobs)
            
            sys.stderr.write('Generated %d trees in %.2f seconds with %d cores\n'%(
                len(tree_jobs), total, cpu_cores))
            sys.stderr.write('Avg.comp.time: %.5f\n'%avg_tree_time)
            
            if args.dry_run:    
                sys.stderr.write('Completed building %d trees\n'%len(tree_computations))
            else:
                with open('%d.tree_computations_unrooted.tsv'%higher_level, 'w') as f:
                    for nog_id,sample_no,tree_nw in tree_computations:
                        f.write('%d\t%d\t%s'%(nog_id,sample_no,tree_nw))
            
            timings.append(time.time())
            
            ## 3b. rerooting
            sys.stderr.write('Applying mid-point rooting to all %d tree computations..\n'%len(tree_jobs))
            if not args.root_notung:
                tree_computations = tree_building.reroot_trees(tree_computations,args.keep_polytomies)
            
            # save tree building results
            if args.dry_run:    
                sys.stderr.write('Completed building %d trees\n'%len(tree_computations))
            else:
                with open('%d.tree_computations.tsv'%higher_level, 'w') as f:
                    for nog_id,sample_no,tree_nw in tree_computations:
                        f.write('%d\t%d\t%s\n'%(nog_id,sample_no,tree_nw))
        
            timings.append(time.time())
        ############## 4. TREE RECONCILIATION ##############
        
        if args.input_reconciliation:
            # load precomputed reconciliation
            reconciliations = []
            with open(args.input_reconciliation) as f:
                for line in f:
                    l = line.rstrip().split()
                    assert len(l) == 3, "Non standart line %s"%line
                    nog_id = int(l[0])
                    root_event = l[1]
                    split_index = float(l[2])
                    reconciliations.append((nog_id,(root_event,split_index)))
            
            # remove flag for successive runs
            args.input_reconciliation = None
            
        elif not tree_computations:
            sys.stderr.write('WARNING: no tree computations found to process!\n')
            reconciliations = []
            pass
        else:
            sys.stderr.write('Generating species trees reconciliation...\n')
            
            # species tree generation
            eggNOG_speciesTree = Tree(reconciliation.EGGNOGv4_SPECIES_TREE)
            # do intial pruning to exclude all non-euNOG species
            higher_node = eggNOG_speciesTree.get_common_ancestor([str(x) for x in eggNOG_level_species[higher_level]])
            eggNOG_speciesTree = higher_node.detach()
            
            if cpu_cores > 1:
                p = Pool(cpu_cores)
                cached_jobs = [(x,y,z,eggNOG_speciesTree,args.keep_polytomies) for x,y,z in tree_computations]
                reconciliation_jobs = p.map(reconciliation.prepare_reconciliation_job,cached_jobs)
                p.close()
                p.join()
            else:
                reconciliation_jobs = []
                for nog_id,sample_no,tree_nw in tree_computations:
                    species_nw = reconciliation.prune_species_tree(tree_nw,eggNOG_speciesTree,args.keep_polytomies)
                    job = (nog_id,sample_no,tree_nw,species_nw,args.root_notung)
                    reconciliation_jobs.append(job)
        
            sys.stderr.write('Starting reconciliation for %d jobs...\n'%len(reconciliation_jobs))
            reconciliation_method = reconciliation.reconcile
            
            timings.append(time.time())
            
            # reconciliation
            if args.tree_method == 'cluster':
                reconciliation.submit_taskArray(higher_level,reconciliation_jobs,
                                                args.keep_polytomies,args.root_notung,args.infer_transfers)
                reconciliations = reconciliation.collect_taskArray(higher_level,cpu_cores=cpu_cores)
            elif cpu_cores > 1:
                p = Pool(cpu_cores)
                reconciliations = p.map(reconciliation_method,reconciliation_jobs)
                p.close()
                p.join()
            else:
                reconciliations = list(map(reconciliation_method,reconciliation_jobs))
                
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
            
            # save reconciliation results
            if args.dry_run:
                sys.stderr.write('Completed %d reconciliations...\n'%len(reconciliations))
            else:
                with open('%d.tree_reconciliations.tsv'%higher_level,'w') as f:
                    for nog_id, result in reconciliations:
                        assert result, 'No reconciliation found for %d'%nog_id
                        root_event,split_index = result
                        f.write('%d\t%s\t%.3f\n'%(nog_id,root_event,split_index))
        
        timings.append(time.time())
        ############## 5. SOLUTION JOINING ##############
                
        if not inconsistencies:
            sys.stderr.write('No inconsistencies found to resolve!')
            sys.exit()
        else:
            
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
            
            consistency_jobs = []
            
            #bp()
            
            # assign reconciliation to exNOGs
            for nog_id in inconsistencies:
                
                # rebuild complete exNOG (all-sublevels)
                exNOG = expansion.ExpandedNOG(nog_id)
                exNOG.set_cache(protein_nogs, nog_proteins, protein_names)
                exNOG.build(complete_expansion=True)
                exNOG.remove_cache() # remove big objects for jobs
                
                cog_intersection = set()
                if higher_level in [2,2759,2157]:
                    # do not merge COG/KOGs
                    cog_intersection = cogs & exNOG.higher_nogs.keys()
                
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
                        duplication_threshold = args.majority_vote_threshold * tot_solutions

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
                p = Pool(cpu_cores)
                consistent_nogs = p.map(join.apply_solutions,consistency_jobs)
                p.close()
                p.join()
            else:
                consistent_nogs = list(map(join.apply_solutions,consistency_jobs))
            
            stop = time.time()
            total = stop - start
            sys.stderr.write('Completed %d joins in %.2fs \n'%(len(consistency_jobs),total))
            
            # merge all the consistent nogs with the old definition
            if args.dry_run:
                pass
            else:
                pass
                # [TODO] maybe output the order of decision application?
            
            timings.append(time.time())
            ############## 6. SOLUTION MERGING ##############
            
            new_singletons = defaultdict(set) # dict[level]:set(singletons)
            singleton_reference = defaultdict(dict)
            
            # 1. remove new singletons for every level in old definiton
            for consistent_id, consistent_definiton, consistent_singletons, nog_history in consistent_nogs:
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
                
                #for old_id,new_id in nog_mapper.items():
                #    print("%d\t%d\t%d"%(consistent_id,old_id,new_id))
            
            # write level id history to file
            history_path = '%d.nog_history.tsv'%higher_level
            with open(history_path,'w') as f:
                for consistent_id in level_history:
                    for old_id, new_ids in level_history[consistent_id].items():
                        f.write('%d\t%d\t%s\n'%(consistent_id,old_id,','.join([str(x) for x in new_ids])))
                
            # write a copy of the definition to disk before testing
            copy_dir = '%d.new_definition'%higher_level
            os.makedirs(copy_dir)
            for level_id in protein_nogs:
                with open(os.path.join(copy_dir,'%d.tsv'%level_id),'w') as f:
                    for protein_id, nog_id in protein_nogs[level_id].items():
                        f.write('%d\t%d\n'%(nog_id,protein_id))
                        
                new_definition_pickle = hgt_utils.save_eggNOG_nog_mapping(
                    level_id,protein_nogs[level_id],copy_dir)
            
            # set test for consistency before writing new definition
            lower_levels = list(eggNOG_children[higher_level])
            test.test_consistency([higher_level],lower_levels,protein_nogs)
            while lower_levels:
                lower_level = lower_levels.pop()
                if lower_level in eggNOG_children:
                    sub_levels = list(eggNOG_children[lower_level])
                    test.test_consistency([lower_level],sub_levels,protein_nogs)
                    lower_levels.extend(sub_levels)
            
            # write to disk new definitions
            for level_id in protein_nogs:
                with open(os.path.join(output_dir,'%d.tsv'%level_id),'w') as f:
                    for protein_id, nog_id in protein_nogs[level_id].items():
                        f.write('%d\t%d\n'%(nog_id,protein_id))
                        
                new_definition_pickle = hgt_utils.save_eggNOG_nog_mapping(
                    level_id,protein_nogs[level_id],output_dir)
                
            sys.stderr.write("Wrote consistent NOG definitions to %s\n"%output_dir)
            
            with open('%d.new_singletons.txt'%higher_level,'w') as f:
                for level_id in new_singletons:
                    for nog_id in new_singletons[level_id]:
                        f.write('%d\t%d\n'%(level_id,nog_id))
            
            timings.append(time.time())
    
            with open('timings.log','a') as f:
                step_times = timings
                #assert len(step_times) == 8, "Unexpected number of step times: %s"%step_times
                step_durations = ["%.3f"%(step_times[i] - step_times[i-1]) for i in range(1,len(step_times))]
                f.write('%d\t'%higher_level + '\t'.join(step_durations) + '\n')
            
        # sys.stderr.write('Inconsistencies to process: %s\n'%inconsistencies)
        # sys.stderr.write('Simple merges to process: %s\n'%simple_merges)
        # sys.stderr.write('Total # of inconsistencies: %s\n'%sum(inconsistencies.values()))
        # sys.stderr.write('Total # of simple merges: %s\n'%sum(simple_merges.values()))
        # # [TODO] add stats about split and join 
        # 
