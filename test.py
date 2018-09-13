# basic testing of script modules

import os
assert os.path.exists('data/eggNOG_names.tsv'), "no data dir found, please check if data.tar.gz was extracted correctly!"

from methods.utils import eggNOG_utils as eu
assert eu.read_eggNOG_names()[9606] == 'Homo sapiens', "Expected 'Homo sapiens' but obtained: %s"%eu.read_eggNOG_names()[9606]
print('Successfully read eggNOG species names!')

from methods import hgt_utils as hu
assert 'CL4' in hu.load_v4clades('homNOG'), "Expected clade CL4 not found in homNOG clades: %s"%hu.load_v4clades('homNOG')
print('Successfully read hominidae [homNOG] clades!')

