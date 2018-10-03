#!/usr/bin/env python3

#Script defining where raw data files can be found

#author:    Davide Heller
#email:     davide.heller@imls.uzh.ch
#date:      2015-08-25

#TODO: possibly replace constant string with constant class field:
#      see: http://stackoverflow.com/a/2688086

from os.path import join

#computed eggNOG reference files
EGGNOG_OUTPUT = "preprocessed_data"

#eggNOG clades
EGGNOG_CLADES = join(EGGNOG_OUTPUT,"clades")

#tool binaries
BINARIES = "bin"