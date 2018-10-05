from collections import Counter, defaultdict
from ete3 import Tree
from os.path import join

rule build_eggNOG_species_tree:
    input:
        guiding_tree = config['species_tree'],
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

rule build_eggNOG_files:
    input:
        level_tree=config['level_hierarchy'],
        species_names=config['species_names'],
        level_names=config['level_names']
    output:
        tree_tsv="preprocessed_data/eggNOG_tree.tsv",
        levels_only_tsv='preprocessed_data/eggNOG_tree.levels_only.tsv',
        members_tsv="preprocessed_data/eggNOG_level_members.tsv",
        species_txt='preprocessed_data/eggNOG_species.txt',
        names='preprocessed_data/eggNOG_names.tsv'
    run:
        def read_tsv_dict(tsv_file):
            with open(tsv_file) as f:
                return dict(line.rstrip().split('\t') for line in f)
        
        #load data
        eggNOG_species = read_tsv_dict(input.species_names)
        eggNOG_tree = read_tsv_dict(input.level_tree)
        eggNOG_levels= read_tsv_dict(input.level_names)
        
        #build eggNOG level sets of memeber species
        eggNOG_level_sets = defaultdict(set)
        eggNOG_level_sets['1'].update(eggNOG_species)

        for species_id in eggNOG_species:
            assert species_id in eggNOG_tree, 'Species {id}[{name}] not found in hierarchy!'.format(
                id=species_id,
                name=eggNOG_species[species_id]
            )
            
            next_level = eggNOG_tree[species_id]
            assert next_level in eggNOG_levels, 'Unknown level {id} not found in level names'.format(
                id=next_level
            )
            
            while next_level != '1':
                eggNOG_level_sets[next_level].add(species_id)
                next_level = eggNOG_tree[next_level]
        
        # write out
        with open(output.tree_tsv,'w') as f:
            for low_id, high_id in eggNOG_tree.items():
                f.write('%s\t%s\n'%(low_id,high_id))
                
        with open(output.levels_only_tsv,'w') as f:
            for low_id, high_id in eggNOG_tree.items():
                if low_id in eggNOG_levels:
                    f.write('%s\t%s\n'%(low_id,high_id))
                
        with open(output.members_tsv,'w') as f:
            for tax_id, tax_members in eggNOG_level_sets.items():
                f.write('%s\t%d\t%s\n'%(tax_id,len(tax_members),','.join(tax_members)))
                
        with open(output.species_txt,'w') as f:
            for tax_id in eggNOG_species:
                f.write('%s\n'%tax_id)
        
        with open(output.names,'w') as f:
            for tax_id,tax_name in eggNOG_species.items():
                f.write('%s\t%s\n'%(tax_id,tax_name))
            for tax_id,tax_name in eggNOG_levels.items():
                f.write('%s\t%s\n'%(tax_id,tax_name))
            
rule pickle_protein_names:
    input:
        protein_names=config['protein_names']
    output:
        protein_names_pickle='preprocessed_data/proteinINT.tupleSpeciesINT_ShortnameSTR.pkl'
    run:
        import pickle
        protein_dict = {}
        
        with open(input.protein_names, 'r') as f:
            for line in f:
                # e.g. 394.NGR_c00010	1
                protein_name, protein_id = line.strip().split('\t')
                
                # get protein species and shortname
                separator_idx = protein_name.find('.')
                assert(separator_idx != -1)
                
                # filter speciesID, e.g. 394
                protein_species = protein_name[:separator_idx]
                assert(protein_species.isdigit())
                protein_species = int(protein_species)
                
                # if args.filter_species is not None:
                #     if protein_species not in filter_species:
                #         continue
                
                # filter shortname, e.g. NGR_c00010
                protein_short = protein_name[separator_idx+1:]
                
                # filter proteinID, e.g. 1
                assert(protein_id.isdigit())
                protein_id = int(protein_id)
                
                # fill dictionary,  e.g. [1] = (394,NGR_c00010)
                protein_dict[protein_id] = (protein_species,protein_short)
        
        # 2. Pickle protein dictionary
        with open(output.protein_names_pickle,'wb') as pfile:
            pickle.dump(protein_dict,pfile)