#!/usr/bin/env python3
# generate simple test data set with 100 randomly expanded NOGs from hominidae level

import os
import random
from collections import defaultdict
from tqdm import tqdm
from methods import hgt_utils, expansion

if __name__ == '__main__':
    
    output_dir = 'test_data'
    assert not os.path.exists(output_dir), '[test_data] already present, remove to generate new'
    os.makedirs(output_dir)
    
    print('# Loading definition')
    protein_nogs, nog_proteins = hgt_utils.load_join_data(9443, output_dir , 'data/pickles')
    protein_names = hgt_utils.load_eggNOG_protein_names_pickle()

    print('# Sampling nogs')
    nog_level = 9604 # start with homNOG
    random.seed(1)
    random_nogs = random.sample(list(nog_proteins[9604]),100)

    print('# Expanding nogs')
    visited_nogs = set()
    for nog_id in tqdm(sorted(random_nogs)):
        exNOG = expansion.ExpandedNOG(nog_id)
        exNOG.set_cache(protein_nogs, nog_proteins, protein_names)
        exNOG.build()
        visited_nogs.update(exNOG.get_expansion_keys())
        # print(exNOG)

    print('# Writing out test definition')
    for level_id in tqdm(protein_nogs):
        test_definition = {}
        with open(os.path.join(output_dir,'%d.tsv'%level_id),'w') as f:
            for protein_id, nog_id in protein_nogs[level_id].items():
                if nog_id in visited_nogs:
                    test_definition[protein_id] = nog_id
                    f.write('%d\t%d\n'%(nog_id,protein_id))
                    
        hgt_utils.save_eggNOG_nog_mapping(
            level_id,test_definition,output_dir)
