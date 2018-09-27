from os.path import join
from scripts.methods.utils import eggNOG_utils as eu

configfile: 'config.yaml'

level_hierarchy = eu.read_eggNOG_treeRev() # TODO replace with config

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
        default_solutions=join(config['output_dir'],'default_solutions/{level_id}.tsv'),
        inconsistencies=join(config['output_dir'],'inconsistencies/{level_id}.tsv')
    output:
        consistent_ogs=join(config['consistent_ogs'],'{level_id}.tsv'),
        new_singletons=join(config['output_dir'],'new_singletons/{level_id}.tsv')
    params:
        majority_vote_threshold=0.5
    threads:
        20 # max=20, i.e. threads = min(threads, cores)
    script:
        'scripts/s05_06_join_and_propagate.py'

rule tree_reconciliation:
    input:
        trees = join(config['output_dir'],'trees/{level_id}.tsv'),
        reconciliation_software = 'bin/Notung-2.9.jar'
    output:
        reconciliations = join(config['output_dir'],'reconciliations/{level_id}.tsv')
    threads:
        20 # max=20, i.e. threads = min(threads, cores)
    params:
        computation_method = 'multicore',
        root_notung=False,
        keep_polytomies=False,
        infer_transfers=False
    script:
        'scripts/s04_tree_reconciliation.py'

rule tree_building:
    input:
        samples=join(config['output_dir'],'samples/{level_id}.tsv'),
        alignment_software = 'bin/mafft-linux64/mafft.bat',
        tree_software = 'bin/FastTree'
    output:
        trees_rooted=join(config['output_dir'],'trees/{level_id}.tsv'),
        trees_unrooted=join(config['output_dir'],'unrooted_trees/{level_id}.tsv')
    threads:
        20 # max=20, i.e. threads = min(threads, cores)
    params:
        tree_method='website',
        root_notung=False,
        keep_polytomies=False,
    script:
        "scripts/s03_tree_building.py"

rule expansion:
    input:
        join(config['input_dir'],'{level_id}.tsv'),
        children=get_children_paths
    output:
        samples=join(config['output_dir'],'samples/{level_id}.tsv'),
        default_solutions=join(config['output_dir'],'default_solutions/{level_id}.tsv'),
        inconsistencies=join(config['output_dir'],'inconsistencies/{level_id}.tsv')
    params:
        random_seed = 1,
        sample_no = 20,
        sample_size = 10,
        sample_method = 'combined',
        default_action = None,
        tree_limit = -1, # no limit
        verbose = False
    script:
        'scripts/s01_02_expand_and_sample.py'
        
rule test_data:
    input:
        'data/pickles/9443.nogINT.setProteinINT.pkl2'    
    output:
        'test_data/9443.tsv',
        'test_data/9604.tsv',
        'test_data/314294.tsv'
    script:
        'scripts/setup.py'
        
rule expand_data:
    input:
        data='data.tar.gz'
    output:
        'data/pickles/9443.nogINT.setProteinINT.pkl2'
    shell:
        "tar -xzf data.tar.gz"
        
rule expand_binaries:
    input:
        binaries='bin.tar.gz'
    output:
        alignment_software = 'bin/mafft-linux64/mafft.bat',
        tree_software = 'bin/FastTree',
        reconciliation_software = 'bin/Notung-2.9.jar'
    shell:
        "tar -xzf bin.tar.gz"