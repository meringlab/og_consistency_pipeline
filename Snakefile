# Snakemake workflow to make a hierarchy of Orthologous Groups (OGs) consistent
# Copyright (C) 2018  Davide Heller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# author:   Davide Heller
# email:    davide.heller@imls.uzh.ch
# version:  0.3 [2018-09-27]

from os.path import join
from collections import defaultdict
# to download tool binaries
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: 'config.yaml'

level_hierarchy = None

def read_level_hierarchy():
    # returns level hierarchy as dictionary {parent -> children}
    species = set()
    with open(config['species_names']) as f:
        for line in f:
            l = line.rstrip().split('\t')
            species.add(l[0])
    
    hierarchy = defaultdict(list)
    with open(config['level_hierarchy']) as f:
        for line in f:
            node, parent = line.rstrip().split()
            if node in species:
                continue
            hierarchy[parent].append(node)
    return hierarchy

def get_children_paths(wildcards):
    # returns location of children levels, using input_dir if leaf (no sub-levels)
    children_paths = []

    global level_hierarchy
    if level_hierarchy is None:
        level_hierarchy = read_level_hierarchy()
    assert wildcards.level_id in level_hierarchy, 'level_id %s not found in level hiearchy!'%wildcards.level_id

    # return children location
    for child_id in level_hierarchy[wildcards.level_id]:
        if child_id in level_hierarchy:
            # inner level
            children_paths.append(join(config['consistent_ogs'],'%s.tsv'%child_id))
        else:
            # leaf level
            children_paths.append(join(config['input_dir'],'%s.tsv'%child_id))

    return children_paths

rule all:
    input:
        join(config['consistent_ogs'],'%s.tsv'%config['target'])

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
        reconciliation_software = 'bin/Notung-2.9/Notung-2.9.jar'
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

include: 'rules/tools.smk'
rule download_tools:
    # using included rules/tools.smk
    input:
        'bin/FastTree',
        'bin/mafft-linux64/mafft.bat',
        'bin/Notung-2.9/Notung-2.9.jar'

include: 'rules/preprocessing.smk'
rule preprocess_data:
    # using included rules/preprocessing.smk
    input:
        tree_tsv="preprocessed_data/eggNOG_tree.tsv",
        levels_only_tsv='preprocessed_data/eggNOG_tree.levels_only.tsv',
        members_tsv="preprocessed_data/eggNOG_level_members.tsv",
        species_txt='preprocessed_data/eggNOG_species.txt',
        tree_nhx="preprocessed_data/eggNOG_tree.levels_only.nhx",
        species_tree = 'preprocessed_data/eggNOG_species_tree.nw',
        protein_names_pickle='preprocessed_data/proteinINT.tupleSpeciesINT_ShortnameSTR.pkl'

rule expand_data:
    input:
        data='data.tar.gz'
    output:
        'data/orthologous_groups/9443.tsv',
        'data/orthologous_groups/9604.tsv',
        'data/orthologous_groups/314294.tsv'
    shell:
        "tar -xzf data.tar.gz"
        
rule generate_test_data:
    input:
        rules.preprocess_data.input,
        rules.expand_data.output
    params:
        random_samples=100
    output:
        'test_data/9443.tsv',
        'test_data/9604.tsv',
        'test_data/314294.tsv'
    script:
        'scripts/setup.py'
