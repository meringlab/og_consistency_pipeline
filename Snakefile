from os.path import join
from methods.utils import eggNOG_utils as eu

configfile: 'config.yaml'

level_hierarchy = eu.read_eggNOG_treeRev()

def get_children_paths(wildcards):
    level_id = int(wildcards.level_id)
    assert level_id in level_hierarchy
    children = level_hierarchy[level_id]
    children_paths = []
    for child_id in children:
        if child_id in level_hierarchy:
            children_paths.append(join(config['output_dir'],'consistent_ogs/%d.tsv'%child_id))
        else:
            # leaf
            children_paths.append(join(config['input_dir'],'%d.tsv'%child_id))
    return children_paths
    
rule all:
    input:
        join(config['output_dir'],'consistent_ogs/{level_id}.tsv')

rule join:
    input:
        parent=join(config['input_dir'],'{level_id}.tsv'),
        children=get_children_paths,    
        reconciliations=join(config['output_dir'],'reconciliations/{level_id}.tsv'),
        #default_solutions=join(config['output_dir'],'default_solutions/{level_id}.tsv'),
    output:
        join(config['output_dir'],'consistent_ogs/{level_id}.tsv')
    shell:
        'touch {output}'

rule tree_reconciliation:
    input:
        join(config['output_dir'],'trees/{level_id}.tsv')
    output:
        join(config['output_dir'],'reconciliations/{level_id}.tsv')
    shell:
        'touch {output}'

rule tree_building:
    input:
        join(config['output_dir'],'samples/{level_id}.tsv')
    output:
        join(config['output_dir'],'trees/{level_id}.tsv')
    shell:
        "touch {output}"

rule expansion:
    input:
        join(config['input_dir'],'{level_id}.tsv'),
        children=get_children_paths
    output:
        samples=join(config['output_dir'],'samples/{level_id}.tsv')
    params:
        random_seed = 1,
        sample_no = 20,
        sample_size = 10,
        sample_method = 'combined',
        verbose = False
    run:
        from tqdm import tqdm
        from methods import expansion, sampling, hgt_utils
        from methods.utils import eggNOG_utils as eu
        
        import random
        random.seed(params.random_seed) # local only? maybe level_id multiplication
        
        level_id = int(wildcards.level_id)
        
        protein_names = hgt_utils.load_eggNOG_protein_names_pickle()
        protein_fasta = hgt_utils.get_level_fasta(level_id)
        
        sampler = sampling.InconsistencySampler(params.sample_no,
                                                params.sample_size,
                                                params.sample_method,
                                                protein_fasta)
            
        protein_nogs, nog_proteins = hgt_utils.load_join_data(level_id, None,'test_data')
        
        # dictionary of inconsistencies to process
        discovered_nogs = set()
        sample_counter = 0
        total_pre_splits = 0
        total_pre_merges = 0
        buffer = ""
        for children_id in level_hierarchy[level_id]:
            for nog_id in tqdm(nog_proteins[children_id]):
                
                if nog_id in discovered_nogs:
                    continue
                
                exNOG = expansion.ExpandedNOG(nog_id,verbose=params.verbose)
                exNOG.set_cache(protein_nogs, nog_proteins, protein_names)
                exNOG.build()
                
                discovered_nogs.update(exNOG.get_expansion_keys())
        
                if exNOG.has_inconsistencies():
                    
                    pre_splits = exNOG.pre_split()
                    total_pre_splits += len(pre_splits) / 2
                    pre_merges = exNOG.pre_merge()
                    total_pre_merges += len(pre_merges)
                    
                    for inconsistent_nog in exNOG.find_inconsistencies():
                        
                        higher_nogs = exNOG.lower_nogs[inconsistent_nog]
                        higher_nogs_size = { x:len(exNOG.expanded_nog[x]) for x in higher_nogs }
                        plus2_sizes = [ size for size in higher_nogs_size.values() if size > 1 ]
                        
                        if len(plus2_sizes) < 2:
                            continue
                     
                        proteins, levels, splits, paralogs = exNOG.get_composition(inconsistent_nog)     
                        samples = sampler.sample_inconsistency(proteins, levels, splits, paralogs)
                        
                        for fasta_sequences in samples:
                            buffer += "%s\t%d\t%s\n"%(nog_id,sample_counter,repr(fasta_sequences))
                            sample_counter += 1
                            
        with open(output.samples,'w') as f:
            f.write(buffer)
        
        print('Completed expansion with a total of %d pre-splits and %d pre-merges'%(
                total_pre_splits,total_pre_merges))