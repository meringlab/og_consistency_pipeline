from methods.utils import eggNOG_utils as eu

level_hierarchy = eu.read_eggNOG_treeRev()

def get_children_paths(level_id):
    level_id = int(level_id)
    assert level_id in level_hierarchy
    children = level_hierarchy[level_id]
    children_paths = []
    for child_id in children:
        if child_id in level_hierarchy:
            children_paths.append('consistent_ogs/%d.tsv'%child_id)
        else:
            # leaf
            children_paths.append('/mnt/gaia/davide/eggnog/eggnog5/subluca_consistent/%d.tsv'%child_id)
    return children_paths
    
rule make_consistent:
    input:
        lambda wildcards: get_children_paths(wildcards.level_id)
    output:
        consistent_level='consistent_ogs/{level_id}.tsv'
    shell:
        "touch {output.consistent_level}"