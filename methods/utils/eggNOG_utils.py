#!/usr/bin/env python

"""
Collection of helper function to ease interaction with the eggNOG data and general file handling.

author: Davide Heller
email:  davide.heller@imls.uzh.ch
date:   2014-03-04
"""

import os
import sys
import ast
import gzip

from . import file_utils
from . import data_sources as ds

def loadBigProteinIndexDict():
    
    input_file = ds.EGGNOGv4_PROT_MAP
    
    assert(os.path.isfile(input_file)),'No file at %s'%input_file
    
    protein_dict = {}
    
    pc = file_utils.ProgressCounter(9637435,'Loading protein integers..')
    
    with open(input_file, 'r') as f:
        for line in f:
            l = line.strip().split('\t')
                
            protein_name = l[0]
                
            protein_integer = l[1]
            
            protein_dict[protein_name] = protein_integer
            
            pc.count()
            
    print('Completed')        
    return protein_dict

def getSpecies2CladeDict(og_label):
    """Read a single eggNOG label and all contained species with their clade
    
    output:
    dict[Species_ID] = Clade_ID     
    """
    
    og_file = "%s/%s/species_and_clades.txt"%(
        ds.EGGNOGv4_COMPUTATIONS,og_label)
    
    if(not os.path.exists(og_file)):
        print("OG label %s not found in %s please verify correctness!"%(og_label,og_file))
        sys.exit()
        
    #dictionary to contain the species_and_clades.txt file
    id2clade = {}
    
    #read in the particular id2clade of veNOG - vertebrate NOGs
    for line in open(og_file):
        l = line.rstrip().split()
        if(len(l) > 1):
            for taxid in l[1:]:
                id2clade[taxid] = l[0]
    
    return id2clade

def read_NOG_species():
    """Return all available species at the highest eggNOG level, i.e. 'NOG_COG'"""
    specie2clade = getSpecies2CladeDict('NOG_COG')
    return list(specie2clade.keys())


def read_ids(og_label):
    """Return all available species at a given eggNOG level, e.g. 'biNOG'"""
    id2clade = getSpecies2CladeDict(og_label)
    return list(id2clade.keys())

def read_tsv_dict(
    file_name,field_no = 2,type1 = str, type2 = str,
    REVERSE_ORDER = False,second_idx = -1,
    skip_header = False,skip_comment = False,skip_assert = False):
    """
    see [file_utils.read_tsv_dict] Reference only maintained for
    backwards compatability.
    """
    return file_utils.read_tsv_dict(
        file_name,field_no,type1, type2,REVERSE_ORDER,second_idx,
        skip_header,skip_comment,skip_assert)
    

def save_tsv_dict(tsv_dict, output_file_name = "tmp.tsv", REVERSE_ORDER = False, header_str="none"):
    """
    see [file_utils.save_tsv_dict] Reference only maintained for
    backwards compatability.
    """
    return file_utils.save_tsv_dict(
        tsv_dict, output_file_name, REVERSE_ORDER, header_str)

def get_eggNOG_nhx():
    """Location of eggNOG tree in nhx format, with ete use format 8 to load"""
    return os.path.join(ds.EGGNOG_OUTPUT,"eggNOG_tree.levels_only.nhx")
    
def get_og_labels_map():
    """Returns dictionary with all eggNOG labels mapped to their lca species id"""
    og_label_file = ds.EGGNOGv4_DERIVED_LABELS
    og_label2int = {}
    
    for line in open(og_label_file):
        l = line.rstrip().split('\t')
        assert(len(l) == 3)
        og_label = l[0]
        og_tax_id = l[2]

        og_label2int[og_label] = int(og_tax_id)
        
    return og_label2int

def get_full_level_names():
    """Returns dictionary with all eggNOG labels mapped to their lca species id"""
    og_label_file = ds.EGGNOGv4_DERIVED_LABELS
    og_lable_full_name = {}
    
    for line in open(og_label_file):
        if line.startswith('#'):
            continue
        l = line.rstrip().split('\t')

        og_tax_id = int(l[0])
        og_full_name = l[2]

        og_lable_full_name[og_tax_id] = og_full_name
        
    return og_lable_full_name

def get_level_hierarchy(level_id,eggNOG_tree=None, eggNOG_names=None):
    if eggNOG_tree is None:
        eggNOG_tree = read_eggNOG_tree()
    assert level_id in eggNOG_tree
    level_hierarchy = []
    while level_id != 1:
        level_id = eggNOG_tree[level_id]
        level_hierarchy.append(level_id)
    
    if eggNOG_names is None:
        return level_hierarchy
    else:
        return [eggNOG_names[x] for x in level_hierarchy]

def get_nog_id(level_id,nog_number,nog_type = 1):
    """ Given level_id and nog_number returns nog_id,
    e.g. 33213 [biNOG], 1 returns 1000033213000000001 (default type == 1)"""
    
    nog_str = '%d%09d%09d'%(nog_type,level_id,nog_number)
    nog_id = int(nog_str)
    
    return nog_id
    
def decompose_nog_id(nog_id):
    """ Given 1000033213000602804 returns 33213 and 602804"""
    nog_id = str(nog_id)

    assert nog_id.isdigit(), 'nog id is not a number! %s'%nog_id
    assert len(nog_id) == 19, 'nog id is not 19 characters long! %s'%nog_id
    
    level_id = int(nog_id[1:10])
    nog_number = int(nog_id[10:])

    return level_id, nog_number

def get_copy_id(nog_id,nog_number,nog_type):
    """ Given an existing nog_id a new nog_id is computed by extracting the level from the old"""
    
    level_id,_ = decompose_nog_id(nog_id)
    
    return get_nog_id(level_id,nog_number,nog_type)
    
def read_core_species():
    """Returns list with all core species of eggNOG"""
    input_file = ds.EGGNOG_OUTPUT + "/core_periphery_species.tsv"
    field_no = 2
    core_periphery_dict = read_tsv_dict(input_file, field_no, int, str)
    
    core_species = []
    for tax_id in core_periphery_dict:
        if(core_periphery_dict[tax_id] == "core"):
            core_species.append(tax_id)
            
    return core_species

def read_core_species_set():
    """Duplicate of read_core_species, returns set instead of list"""
    input_file = ds.EGGNOG_OUTPUT + "/core_periphery_species.tsv"
    field_no = 2
    core_periphery_dict = read_tsv_dict(input_file, field_no, int, str)
    
    core_species = set()
    for tax_id in core_periphery_dict:
        if(core_periphery_dict[tax_id] == "core"):
            core_species.add(tax_id)
            
    return core_species

def read_all_species():
    input_file = ds.EGGNOG_OUTPUT + "/eggNOG_species.txt"
    with open(input_file) as f:
        species = [int(line.rstrip()) for line in f]
    
    return species

def read_level_species():
    """
    Return a dictionary with the set of species for each level
    
    return dict[level_id] = set([species_id])
    """
    #load eggNOG data
    eggNOG_species = read_all_species()
    eggNOG_tree = read_eggNOG_tree()
    
    #build eggNOG species level sets
    eggNOG_level_sets = {}
    eggNOG_level_sets[1] = set(eggNOG_species)

    for species_id in eggNOG_species:
        next_level = eggNOG_tree[species_id]
        while next_level != 1:
            if next_level not in eggNOG_level_sets:
                eggNOG_level_sets[next_level] = set()
            eggNOG_level_sets[next_level].add(species_id)
            next_level = eggNOG_tree[next_level]
    
    return eggNOG_level_sets

def read_eggNOG_tree():
    """Returns the tree structure of eggNOG (species + levels) as dictionary"""
    input_file = ds.EGGNOG_OUTPUT + "/eggNOG_tree.tsv"
    field_no = 2
    eggNOG_tree = read_tsv_dict(input_file, field_no, int, int)
    
    return eggNOG_tree

def read_eggNOG_treeRev():
    """
    Reverses the tree dictionary with sets for lower levels
    """
    reverse_dict = {}
    eggNOG_tree = read_eggNOG_tree()
    eggNOG_species = read_all_species()
    
    for lower_level,higher_level in eggNOG_tree.items():
        
        if lower_level in eggNOG_species:
            continue
        
        if higher_level in reverse_dict:
            reverse_dict[higher_level].add(lower_level)
        else:
            reverse_dict[higher_level] = {lower_level}
            
    return reverse_dict

def read_eggNOG_names():
    """Returns a dictionary that maps al the eggNOG integer ids to their string name"""
    input_file = ds.EGGNOG_OUTPUT + "/eggNOG_names.tsv"
    eggNOG_names = read_tsv_dict(input_file, 2, int, str)

    return eggNOG_names

def read_eggNOG_namesRev():
    """Returns a dictionary that maps al the eggNOG integer ids to their string name"""
    input_file = ds.EGGNOG_OUTPUT + "/eggNOG_names.tsv"
    eggNOG_names = read_tsv_dict(input_file, 2, int, str,REVERSE_ORDER = True)

    return eggNOG_names

def read_eggNOG_tree_extended():
    """Same functionality as read_eggNOG_tree but includes all
    intermediate levels in between eggNOG levels where species
    are added to the tree
    e.g. between homNOG and prNOG, Catarrhini is added for the entry of Macaca Mulatta
    """
    input_file = ds.EGGNOG_OUTPUT + "/eggNOG_tree.extended.tsv"
    field_no = 2
    eggNOG_tree = read_tsv_dict(input_file, field_no, int, int)
    
    return eggNOG_tree

def read_eggNOG_names_extended():
    """Same functionality as read_eggNOG_name but for exteded version
    see read_eggNOG_names() and read_eggNOG_tree_extended()
    """
    input_file = ds.EGGNOG_OUTPUT + "/eggNOG_names.extended.tsv"
    eggNOG_names = read_tsv_dict(input_file, 2, int, str)

    return eggNOG_names

def read_eggNOG_fusions(level_id):
    """ read in fusion proteins, i.e. proteins that appear in at least two nogs in the same level """
    fusion_file = os.path.join(ds.EGGNOG_OUTPUT,'eggNOG_fusions.tsv')
    
    with open(fusion_file) as f:
        for line in f:
            if line.startswith("%d\t"%level_id):
                fusions = ast.literal_eval(line.rstrip().split('\t')[1])
                return fusions
    
    return {}

def compute_first_names():
    """
    Return a dictionary with only the first name of eggNOG species
    Names starting with special characters have been custom renamed
    
    dict[9606] = "Homo" instead of "Homo sapiens"
    """
    
    eggNOG_names = read_eggNOG_names()
    
    first_names = {str(x):(y.split(' ')[0]) for x,y in eggNOG_names.items()}

    bad_names = {"'Nostoc":'Nostoc',
                 '[Cellvibrio]':'Cellvibrio',
                 'butyrate-producing':'butyrate',
                 '[Bacteroides]':'pectinophilus',
                 'SAR116':'SAR'}
    
    for species_id,first_name in first_names.items():
        for bad_name in bad_names:
            if first_name.startswith(bad_name):
                first_names[species_id] = bad_names[bad_name]
                
    return first_names

def compute_short_names(short_length=5):
    """
    Return a dictionary with short version of eggNOG species names
    by extracting the first [short_length] alphabetic characters. E.g.
    
    dict[9606] = "HomoS" instead of "Homo sapiens"
    """
    eggNOG_names = read_eggNOG_names()
    
    short_names = {}
    
    for species_id,species_name in eggNOG_names.items():
        if 'NOG' not in species_name:
            alpha_name = ''.join(filter(str.isalpha, species_name))
            assert len(alpha_name) > short_length, 'Chosen length not compatible with %s'%species_name
            short_names[str(species_id)] = alpha_name[:short_length]
            
    return short_names

def read_og_names(og_label_id):
    """Returns dict[NOG_integer_id] =  original_NOG_name_str"""
    input_file = ds.EGGNOG_OUTPUT + "/og_integer_mapping.updated/"+\
                 str(og_label_id)+".og_prot.stringent.tsv"
    
    return read_tsv_dict(input_file,2, str, int, REVERSE_ORDER = True)

def read_og_ids(og_label_id):
    """Returns NOG_integer_id for original_NOG_name_str"""
    input_file = ds.EGGNOG_OUTPUT + "/og_integer_mapping.updated/"+\
                 str(og_label_id)+".og_prot.stringent.tsv"
    
    # [nog_str_name]:[nog_integer_id]
    nog_ids = {}
    
    with open(input_file) as nog_file:
        for line in nog_file:
            nog_name, nog_id = line.rstrip().split('\t')
            nog_ids[nog_name] = int(nog_id)
    
    return nog_ids           

def read_og_prot(og_label_id,eggNOG_version = 40):
    """Retruns dict[protein_integer_id] =  NOG_integer_id
    
    Stringency level: stringent
    
    DISCLAIMER: Please note that the following methods do not take into account the
    existance of fusion proteins, i.e. the same protein id assigned to several NOGs
    but only attempt to check the existance of a protein at a given level.
    """
    if(eggNOG_version == 41):
        input_file = ds.EGGNOG_OUTPUT + "/og_protein_mapping.41/"+\
            str(og_label_id)+".og_id_prot.tsv"            
     
        return read_tsv_dict(input_file, 2, str, int, REVERSE_ORDER = True)
    elif(eggNOG_version == 40):
        input_file = ds.EGGNOG_OUTPUT + "/og_protein_mapping/"+\
            str(og_label_id)+".og_id_prot.stringent.tsv";
     
        return read_tsv_dict(input_file, 2, int, int, REVERSE_ORDER = True)
    else:
        print('Requested unknown eggnog version: %d'%eggNOG_version)
        sys.exit()

def semi_stringent_exist(og_label_id):
    """Retruns True if the semi-stringent level is available for the input level"""
    file_name = ds.EGGNOG_OUTPUT + "/og_protein_mapping.semi_stringent/"+\
                 str(og_label_id)+".og_id_prot.stringent.tsv";
    
    if(os.path.exists(file_name)):
        return True
    else:
        return False

def read_og_prot_semi_stringent(og_label_id):
    """Retruns dict[protein_integer_id] =  NOG_integer_id
    
    Stringency level: semi-stringent
    
    DISCLAIMER about fusion proteins see read_og_prot() 
    """
    input_file = ds.EGGNOG_OUTPUT + "/og_protein_mapping.semi_stringent/"+\
                 str(og_label_id)+".og_id_prot.stringent.tsv";
    return read_tsv_dict(input_file, 2, int, int, REVERSE_ORDER = True)

def read_og_prot_non_stringent(og_label_id):
    """Retruns dict[protein_integer_id] =  NOG_integer_id
    
    Stringency level: non-stringent
    
    DISCLAIMER about fusion proteins see read_og_prot() 
    """
    input_file = ds.EGGNOG_OUTPUT + "/og_protein_mapping.non_stringent/"+\
                 str(og_label_id)+".og_id_prot.stringent.tsv";
    return read_tsv_dict(input_file, 2, int, int, REVERSE_ORDER = True)

def read_nog_proteins(nog_id):
    """Returns the set of proteins from a nog given the nog INTEGER id"""
    nog_protein_dir = '%s/%s'%(ds.EGGNOG_OUTPUT,'og_protein_mapping')
    nog_str = str(nog_id)
    level_id = int(nog_str[1:7])
    nog_protein_file = "%s/%d.og_id_prot.stringent.tsv"%(nog_protein_dir,level_id)
    
    f = open(nog_protein_file)
    nog_proteins = { l.rstrip().split('\t')[1]
                    for l in f if l.startswith(nog_str) }
    f.close()
    
    return nog_proteins

# Read species proteins

def read_species_prot(specie_id):
    """duplicate of read_species_prot_ids(), please use latter"""
    input_file = ds.EGGNOG_OUTPUT + "/divided_protein.shorthands/"+\
                 str(specie_id)+".tsv"
    
    species_proteins = read_tsv_dict(input_file, 2, str, int, REVERSE_ORDER = True)
    
    return list(species_proteins.keys())

def read_species_prot_full(specie_id):
    """read full name as key e.g. dict["9606.ENSP00000001146"] = 1842113"""
    input_file = ds.EGGNOG_OUTPUT + "/divided_protein.shorthands/"+\
                 str(specie_id)+".tsv"
    
    species_proteins = read_tsv_dict(input_file, 2, str, int)
    
    return species_proteins

def read_species_prot_ids(specie_id):
    """read integer id as key e.g. dict[1842113] = "9606.ENSP00000001146"""
    input_file = ds.EGGNOG_OUTPUT + "/divided_protein.shorthands/"+\
                 str(specie_id)+".tsv"
    
    species_proteins = read_tsv_dict(input_file, 2, str, int, REVERSE_ORDER = True)
    
    return list(species_proteins.keys())

def read_species_prot_int_2_str(specie_id):
    """read integer id as key e.g. dict[1842113] = "9606.ENSP00000001146"""
    input_file = ds.EGGNOG_OUTPUT + "/divided_protein.shorthands/"+\
                 str(specie_id)+".tsv"
    
    species_proteins = read_tsv_dict(input_file, 2, str, int, REVERSE_ORDER = True)
    
    return species_proteins

def read_species_prot_names(specie_id):
    """read alternative names of proteins in a given species(ALERT: coverage not complete!)
    key is integer id, e.g. dict[1842113] = "CYP26B1"
    """
    input_file = ds.EGGNOG_OUTPUT + "/divided_protein.names/"+\
                 str(specie_id)+".tsv"
    
    species_protein_names = read_tsv_dict(input_file, 2, int, str)
    
    return species_protein_names

def read_rev_species_prot_names(specie_id):
    """reverse of read_species_prot_names method, i.e.alternative name is
    the dictionary key, e.g. dict[1842113] = "CYP26B1"
    """
    input_file = ds.EGGNOG_OUTPUT + "/divided_protein.names/"+\
                 str(specie_id)+".tsv"
    
    species_protein_names = read_tsv_dict(input_file, 2, int, str, REVERSE_ORDER = True)
    
    return species_protein_names

def printListTree(node, tree, depth = 0):
    """Recursive function to print a tree with increasing indentation,
    starting from the input node level
    
    inputs:
        node = single id
        tree = dictionary representing a tree, i.e. dict[child] = parent
    """
    
    print("\t" * depth, node)

    if node in tree:
        node_list = tree[node]
        for val in node_list:
            assert(val != node),"Double occurrence of %s in tree"%str(val)
            printListTree(val,tree, depth+1)
            
def printListTree_withNames(node, tree, name_dict, depth = 0):
    """Recursive function to print a tree with increasing indentation
    using the name supplied by name_dict (if available) instead of the
    integer id.
    """
    if(node in name_dict):
        print("\t" * depth, name_dict[node], "[", node, "]")
    else:
        print("\t" * depth, node)

    if node in tree:
        node_ids = tree[node]
        for val in tree[node]:
            printListTree_withNames(val,tree,name_dict, depth+1)
            
def save_tree(tree, output_file_name = "tmp.tsv"):
    """Save a dictionary tree in tsv format
    i.e. dict[child] = parent -> str(child[\t]parent[\n])
    """
    output_file = open(output_file_name, 'w')
    
    for child in tree:
        parent = tree[child]
        line = "%d\t%s\n"%(child,str(parent))
        output_file.write(line)
    
    output_file.close()

def find_lca(nodes_to_visit, ncbi_parent_dict, visited_nodes={}, rev_nodes={}):
    """Returns the last common ancestor of the input_nodes in a given tree dictionary
    
    Inputs:
        nodes_to_visit = list of node ids for which to find the lca
        ncbi_partent_dict = dictionary representing the tree, i.e. dict[child] = parent
        visited_nodes = dictionary of nodes that have been visited so far, value is a list of source nodes
        rev_nodes = reverse of visited nodes, i.e. key is the source node, value is a list of destination nodes
    
    Requirements:
        - Tree root should be 1
        - id type should be str or int
        - nodes should have a single parent
    
    NOTE: Recursive function which finds the LCA of the inputed id's
    Logic behind this function is that by adding parent ids to
    the visited_nodes dictionary the search level of the different input
    members is independent. Every parent id is checked for presence
    before being inserted. In the case that the id is already present
    it means that another input ID already visited that node, thus
    the common ancestor is that node. In case of more than two inputs
    the search continues till the search list is reduced to two.
    """    
    
    #if input list has only one element either input
    #was a single species or a multifurcation has occurred
    #and root <> 1 was reached
    if(len(nodes_to_visit) == 1):
        lca = nodes_to_visit[0]
        if(str(lca) == '1'):
            previous_nodes = visited_nodes[lca]
        
            #look for first fork with more than 1 element
            while(len(previous_nodes) == 1):
                lca = previous_nodes[0]
                if lca in visited_nodes:
                    previous_nodes = visited_nodes[lca]
                else:
                    #this should not be possible! LOST multifurcation point
                    break    
            return lca 
    
    new_nodes_to_visit = []
    for k in nodes_to_visit:
        parent_id = ncbi_parent_dict[k]
        if(parent_id == k):
            continue
        rev_nodes[k] = [parent_id]
        if(parent_id not in visited_nodes):
            new_nodes_to_visit.append(parent_id)
            visited_nodes[parent_id] = []
        visited_nodes[parent_id].append(k)
                   
    if(not new_nodes_to_visit):
        if(not nodes_to_visit):
            print("No elements to visit, abort!")
            sys.exit()
            
        lca = nodes_to_visit[0]
        while(lca in rev_nodes):
            #Go up the tree till the highest explored node
            lca = rev_nodes[lca][0]
            
        new_nodes_to_visit.append(lca)
       
    return find_lca(new_nodes_to_visit, ncbi_parent_dict, visited_nodes, rev_nodes)

class eggNOG_lca_finder:
    
    def __init__(self):
        
        #load data
        eggNOG_species = read_all_species()
        self.eggNOG_tree = read_eggNOG_tree()
        eggNOG_names = read_eggNOG_names()
        
        #build custom eggNOG sets
        self.eggNOG_level_sets = {}
        self.eggNOG_level_sets[1] = set(eggNOG_species)

        for species_id in eggNOG_species:
            next_level = self.eggNOG_tree[species_id]
            while(next_level != 1):
                if(next_level not in self.eggNOG_level_sets):
                    self.eggNOG_level_sets[next_level] = set()
                self.eggNOG_level_sets[next_level].add(species_id)
                next_level = self.eggNOG_tree[next_level]
    
    def find_level(self, input_species):
        
        species_in_eggNOG = set(input_species)
        assert(species_in_eggNOG.issubset(self.eggNOG_level_sets[1])), "input species are not fully contained in eggNOG!"
        
        lca_level_id = 1
        
        for species_id in species_in_eggNOG:
            next_level = self.eggNOG_tree[species_id]
            next_level_set = self.eggNOG_level_sets[next_level]
            
            while not species_in_eggNOG.issubset(next_level_set):
                next_level = self.eggNOG_tree[next_level]
                next_level_set = self.eggNOG_level_sets[next_level]
            
            lca_level_id = next_level
            break
        
        return lca_level_id
    
    def get_intermediateLevels(self, parent_level, child_level):
        if(parent_level != 1):
            assert(parent_level in self.eggNOG_tree),'Parent level %d is not valid'%parent_level
            
        if(child_level == 1):
            return []
        assert(child_level in self.eggNOG_tree),'Child Level %d is not valid'%child_level
        
        intermediate_levels = []
        
        next_level = self.eggNOG_tree[child_level]
        while(next_level != parent_level):
            intermediate_levels.append(next_level)
            if(next_level == 1):
                return []
            else:
                next_level = self.eggNOG_tree[next_level]
        
        return intermediate_levels
    
    def get_eggNOG_species_set(self,level_id):
        return self.eggNOG_level_sets[level_id]
    