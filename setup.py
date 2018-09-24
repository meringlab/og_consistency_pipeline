#!/usr/bin/env python3
# generate simple test data set with 100 randomly expanded NOGs from hominidae level

import os
import random
from collections import defaultdict
from tqdm import tqdm
from methods import hgt_utils
from methods.utils import eggNOG_utils as eu

if __name__ == '__main__':
    
    print('# Expanding example data and program binaries..')
    assert os.path.exists('data.tar.gz'), 'Example data not found!'
    assert os.path.exists('bin.tar.gz'), 'Program binaries not found!'
    
    os.system('tar -xzf data.tar.gz')
    os.system('tar -xzf bin.tar.gz')

    print('# Testing extraction and module loading')    
    assert os.path.exists('data/eggNOG_names.tsv'), "no data dir found, please check if data.tar.gz was extracted correctly!"
    
    assert eu.read_eggNOG_names()[9606] == 'Homo sapiens', "Expected 'Homo sapiens' but obtained: %s"%eu.read_eggNOG_names()[9606]
    print('Successfully read eggNOG species names!')
    
    assert 'CL4' in hgt_utils.load_v4clades('homNOG'), "Expected clade CL4 not found in homNOG clades: %s"%hu.load_v4clades('homNOG')
    print('Successfully read hominidae [homNOG] clades!')
    
    print("# Starting test data generation")     
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
    from methods import expansion
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
