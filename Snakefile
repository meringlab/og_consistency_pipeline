#!/usr/bin/env python3
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
from ete3 import Tree

configfile: 'config.yaml'

def get_children_paths(wildcards):
    children_paths = []

    # read level hierarchy
    t = Tree(config['level_hierarchy'],format=8)
    node = t.search_nodes(name=wildcards.level_id)
    assert len(node), 'level_id %s not found in level hiearchy!'%wildcards.level_id
    node = node[0]

    # return children location
    for child in node.get_children():
        child_id = int(child.name)
        if child.is_leaf():
            children_paths.append(join(config['input_dir'],'%d.tsv'%child_id))
        else:
            children_paths.append(join(config['consistent_ogs'],'%d.tsv'%child_id))

    return children_paths

rule all:
    input:
        join(config['consistent_ogs'],config['target'])

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
