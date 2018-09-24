from os.path import join
from methods.utils import eggNOG_utils as eu

ODIR = 'test_output'
INCONSISTENT_OGS = '/mnt/gaia/davide/eggnog/eggnog5/subluca_consistent/'

level_hierarchy = eu.read_eggNOG_treeRev()

def get_children_paths(level_id):
    level_id = int(level_id)
    assert level_id in level_hierarchy
    children = level_hierarchy[level_id]
    children_paths = []
    for child_id in children:
        if child_id in level_hierarchy:
            children_paths.append(join(ODIR,'consistent_ogs/%d.tsv'%child_id))
        else:
            # leaf
            children_paths.append(join(INCONSISTENT_OGS,'%d.tsv'%child_id))
    return children_paths
    
rule all:
    input:
        join(ODIR,'consistent_ogs/{level_id}.tsv')

rule join:
    input:
        parent=join(INCONSISTENT_OGS,'{level_id}.tsv'),
        children=lambda wildcards: get_children_paths(wildcards.level_id),    
        solutions=join(ODIR,'reconciliations/{level_id}.tsv')
    output:
        join(ODIR,'consistent_ogs/{level_id}.tsv')
    shell:
        'touch {output}'

rule tree_reconciliation:
    input:
        join(ODIR,'trees/{level_id}.tsv')
    output:
        join(ODIR,'reconciliations/{level_id}.tsv')
    shell:
        'touch {output}'

rule tree_building:
    input:
        join(ODIR,'samples/{level_id}.tsv')
    output:
        join(ODIR,'trees/{level_id}.tsv')
    shell:
        "touch {output}"

rule expansion:
    input:
        join(INCONSISTENT_OGS,'{level_id}.tsv'),
        children=lambda wildcards: get_children_paths(wildcards.level_id)
    output:
        join(ODIR,'samples/{level_id}.tsv')
    shell:
        "touch {output}"

