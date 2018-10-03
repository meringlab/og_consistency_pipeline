from collections import Counter, defaultdict
from ete3 import Tree
from os.path import join

rule build_eggNOG_species_tree:
    input:
        guiding_tree = config['guiding_tree'],
        species_txt = 'preprocessed_data/eggNOG_species.txt'
    output:
        species_tree = 'preprocessed_data/eggNOG_species_tree.nw'
    run:
        t = Tree(input.guiding_tree)
        leaf_names = set(t.get_leaf_names())
        
        with open(input.species_txt) as f:
            species = {line.rstrip() for line in f}
            
        assert species.issubset(leaf_names), 'Not all eggNOG species found in guide tree: '%(species - leaf_names)
        
        t.write(outfile=output.species_tree,format=5)
    
rule build_eggNOG_nhx:
    """build tree in NHX format of the eggNOG level hierarchy 
    NHX (ete export format=8 + root) includes:
        node.name = tax_id of level (e.g. 33213)
        node.nog_prefix = level_name (e.g. 'maNOG') 
    """
    input:
        tree_tsv="preprocessed_data/eggNOG_tree.levels_only.tsv",
        name_tsv="preprocessed_data/eggNOG_names.tsv"
    output:
        tree_nhx="preprocessed_data/eggNOG_tree.levels_only.nhx"
    run:
        levels = dict([ tuple( line.rstrip().split('\t') ) for line in open(input.tree_tsv) ])
        names = dict([tuple(line.rstrip().split('\t')) for line in open(input.name_tsv)])
        
        # revert levels dictionary => parent = set(children)
        children = defaultdict(set)
        for child,parent in levels.items():
            children[parent].add(child)
        
        # initialize tree
        t = Tree()

        # insert root
        t.name = '1'
        for child in children[t.name]:
            t.add_child(name=child)
        del children[t.name]

        # insert children at leaves
        while children:
            leaves = t.get_leaves()
            for n in leaves:
                if n.name in children:
                    for child in children[n.name]:
                        n.add_child(name=child)
                    del children[n.name]
                    # print('Removed %s (%d left)'%(n.name,len(children)))
        
        for n in t.traverse():
            n.add_feature("nog_prefix",names[n.name])
        
        # write out tree
        print(t.get_ascii(attributes=['name','nog_prefix']))
        t.write(format=8,outfile=output.tree_nhx,format_root_node=True,features=['nog_prefix'])

rule build_eggNOG_names:
    input:
        level_names=config['level_names'],
        species_names=config['species_names'],
        tree_tsv="preprocessed_data/eggNOG_tree.tsv"
    output:
        name_tsv="preprocessed_data/eggNOG_names.tsv"
    run:
        names = {}
        
        # load level names
        with open(input.level_names) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                l = line.rstrip().split('\t')
                level_id = l[0]
                level_name = l[1]
                level_size = int(l[3])
                if level_size == 1:
                    assert level_name == 'cryNOG', 'another level with a single species found:%s'%l
                    continue
                
                if 'newSublevelOf' in level_name:
                    level_name = level_name[len('newSublevelOf'):]
                names[level_id] = level_name
        
        # load compact species names
        with open(input.species_names) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                l = line.rstrip().split('\t')
                species_id = l[0]
                names[species_id] = l[3]

        # check id match
        with open(input.tree_tsv) as f:
            tree = dict((line.rstrip().split('\t')) for line in f)
            for tax_id in names:
                if tax_id != '1':
                    assert tax_id in tree, 'Missing name in tree: %s [%s]'%(tax_id,names[tax_id])
        
        # write out
        with open(output.name_tsv,'w') as f:
            for tax_id, tax_name in names.items():
                f.write('%s\t%s\n'%(tax_id,tax_name))

rule build_eggNOG_tree:
    input:
        level_definition=config['level_definition']
    output:
        tree_tsv="preprocessed_data/eggNOG_tree.tsv",
        levels_only_tsv='preprocessed_data/eggNOG_tree.levels_only.tsv',
        members_tsv="preprocessed_data/eggNOG_level_members.tsv",
        species_txt='preprocessed_data/eggNOG_species.txt'
    run:
        # read in all the level species
        levels = {}
        with open(input.level_definition) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                l = line.rstrip().split('\t')
                tax_id = l[0]
                tax_name = l[1]
                member_species = {x for x in l[-1].split(' ')}
                
                if len(member_species) == 1:
                    assert tax_name == 'cryNOG', 'another level with a single species found:%s'%l
                    continue
                
                levels[tax_id] = member_species
        
        # find minimum superset to construct tree
        level_tree = {}
        for tax_id, member_species in levels.items():
            overlap = {}
            for other_id in levels:
                if tax_id != other_id:
                    other_species = levels[other_id]
                    if member_species.issubset(other_species):
                        difference = other_species - member_species
                        overlap[other_id] = len(difference)
            
            #overlap = Counter({x:len(member_species & y) for x,y in levels.items() if tax_id != y and member_species.issubset(y)})
            if overlap:
                overlap = Counter(overlap)
                minimum_superset = overlap.most_common()[-1]
                level_tree[tax_id] = minimum_superset[0] # link to higher level, e.g. biNOG -> meNOG
        
        level_tree_w_species = dict(level_tree)
        
        # attach the species_ids to tree
        all_species = set.union(*levels.values())
        for species_id in all_species:
            species_levels = Counter({x:len(y) for x,y in levels.items() if species_id in y})
            leaf_level = species_levels.most_common()[-1]
            level_tree_w_species[species_id] = leaf_level[0]
        
        # write out
        with open(output.tree_tsv,'w') as f:
            for low_id, high_id in level_tree_w_species.items():
                f.write('%s\t%s\n'%(low_id,high_id))
                
        with open(output.levels_only_tsv,'w') as f:
            for low_id, high_id in level_tree.items():
                f.write('%s\t%s\n'%(low_id,high_id))
                
        with open(output.members_tsv,'w') as f:
            for tax_id, tax_members in levels.items():
                f.write('%s\t%d\t%s\n'%(tax_id,len(tax_members),','.join(tax_members)))
                
        with open(output.species_txt,'w') as f:
            for tax_id in all_species:
                f.write('%s\n'%tax_id)
