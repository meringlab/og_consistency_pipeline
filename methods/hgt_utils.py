#!/usr/bin/env python

import os
import sys
import time
from glob import glob
import shelve
import cPickle as pickle

from ete3 import Tree
from collections import defaultdict, Counter

from utils import eggNOG_utils as eu
from utils import data_sources as ds
from utils import file_utils as fu

V4_PICKLE_DIR = os.path.join(ds.EGGNOG_OUTPUT,'pickles')
V4_FASTA_DIR = ds.EGGNOGv4_FASTA_DIR

PICKLE_PROTOCOL=2
PICKLE_NOG_EXT = "nogINT.setProteinINT.pkl2"

def load_v4clades(level_name,species_type=str):
    
    clade_file = glob("%s/[1-9]*.%s/species_and_clades.txt"%(ds.EGGNOG_CLADES,level_name))
    
    assert len(clade_file) == 1,"Clade file doesn't exist: %s"%level_name
    
    clade_file = clade_file[0]
    
    clade_species = defaultdict(set)
    with open(clade_file) as f:
        for line in f:
            l = line.rstrip().split()
            clade_id = l[0]
            species = l[1:]
            
            if species_type != str:
                species = [species_type(x) for x in species]
                
            clade_species[clade_id].update(species)
    
    return clade_species

def get_level_order(level_id,with_leaves=False):
    
    level_id = str(level_id)
    
    t = Tree(eu.get_eggNOG_nhx(),format=8)
    
    # define order in which to visit the nodes (lower levels first)
    level_node = t.search_nodes(name=level_id)
    assert len(level_node) == 1, "Ambiguous level_id %s: %s"%(level_id, level_node)
    level_node = level_node[0]
    
    # leaf out leaves
    if with_leaves:
        level_order = list(reversed(
            [ int(x.name) for x in level_node.traverse(strategy='levelorder') ]))
    else:
        level_order = list(reversed(
            [ int(x.name) for x in level_node.traverse(strategy='levelorder') if not x.is_leaf() ]))
    
    return level_order

def read_pickle(file_name):
    assert os.path.exists(file_name), 'Pickle file not found: %s'%file_name
    with open(file_name, 'rb') as pfile:
        return pickle.load(pfile)

def save_pickle(data_to_pickle, file_name):
    with open(file_name, 'wb') as pfile:
        pickle.dump(data_to_pickle, pfile, protocol=PICKLE_PROTOCOL)
        
def load_eggNOG_protein_names_pickle():
    # Returns a complete dictionary for proteins in eggNOG
    # i.e. dict[protein_id]:(species_id,protein_shortname)
    
    sys.stderr.write("Starting to load protein names\n")
    t_start = time.clock()
    names_pickle_file = os.path.join(
        V4_PICKLE_DIR,"proteinINT.tupleSpeciesINT_ShortnameSTR.pkl2")
    protein_names = read_pickle(names_pickle_file)
    t_end = time.clock()
    sys.stderr.write("Finished loading protein names in %.2f s\n"%(t_end - t_start))
    
    return protein_names

def load_only_nog_mapping(level_id, pickle_dir):
    pickle_file = "%d.%s"%(level_id,PICKLE_NOG_EXT)
    pickle_path = os.path.join(pickle_dir,pickle_file)
    
    # [nog_id]:set([protein_id])
    nog_mapping = read_pickle(pickle_path)
    
    return nog_mapping

def load_eggNOG_nog_mapping(level_id, pickle_dir=V4_PICKLE_DIR, pickle_ext=PICKLE_NOG_EXT):
    
    pickle_file = "%d.%s"%(level_id,pickle_ext)
    pickle_path = os.path.join(pickle_dir,pickle_file)
    
    # [nog_id]:set([protein_id])
    nog_mapping = read_pickle(pickle_path)
    
    protein_mapping = defaultdict(dict)
    for nog_id, nog_proteins in nog_mapping.items():
        # protein_mapping.update({x:y for x,y in zip(nog_proteins,[nog_id]*len(nog_proteins))
        for protein_id in nog_proteins:
            protein_mapping[protein_id] = nog_id
    
    # close for further modifications
    protein_mapping.default_factory = None
    if type(nog_mapping) == defaultdict:
        nog_mapping.default_factory = None
    
    return protein_mapping, nog_mapping

def load_COG_mapping(level_id,level_name):
    """ load all COGs/KOGs/arCOGs ids based on the NOG prefix"""
    
    cog_prefix = {2:'COG', 2759:'KOG', 2157:'arCOG'}
    
    assert level_id in cog_prefix, 'no cogs available for requested level %d'%level_id
    
    cog_dir = os.path.join(V4_PICKLE_DIR,'converted_groups')
    cog_path = os.path.join(cog_dir,'%d.%s.conversion.tsv'%(level_id,level_name))
    assert os.path.exists(cog_path), 'No conversion file found at %s'%cog_path
    
    cogs = {}
    with open(cog_path) as f:
        for line in f:
            if line.startswith(cog_prefix[level_id]):
                cog_id, integer_id = line.rstrip().split()
                cogs[cog_id] = int(integer_id)
    
    return cogs

def save_eggNOG_nog_mapping(level_id, protein_mapping, output_dir):
    new_definition_pickle = "%d.%s"%(level_id,PICKLE_NOG_EXT)
    new_definition_pickle_path = os.path.join(output_dir,new_definition_pickle)
    
    nog_mapping = defaultdict(set)
    for protein_id,nog_id in protein_mapping.items():
        nog_mapping[nog_id].add(protein_id)
    
    save_pickle(nog_mapping, new_definition_pickle_path)
    
    return new_definition_pickle_path
    
def load_join_data(higher_level, new_definition, old_definition=V4_PICKLE_DIR):
    """
    Loads the [OLD] nog definitions for the level itself
    and the immidiate leaf children. Loads the [NEW] definition 
    for all the remaining sublevels
    
    e.g. when loading higher_level = 314146 [spriNOG]
    
                    /9443[prNOG][NEW]---9604[homNOG][NEW]
    -314146[spriNOG][OLD]
                    \-9989[roNOG][OLD]
                    
    returns:
    protein_mapping [level_id][proteinID]:[nogID]
    nog_mapping     [level_id][nog_id]:set([protein_id])
    
    """
    eggNOG_tree = Tree(eu.get_eggNOG_nhx() ,format=8)
    higher_node = eggNOG_tree.search_nodes(name=str(higher_level))
    assert len(higher_node) == 1
    higher_node = higher_node[0]
    
    nog_mapping = dict()
    protein_mapping = dict()
    
    # load yourself from old definition
    proteins, nogs = load_eggNOG_nog_mapping(higher_level, old_definition)
    protein_mapping[higher_level] = proteins
    nog_mapping[higher_level] = nogs
    
    # load leaf levels from old and intermediate from new
    for node in higher_node.get_children():
        
        level_id = int(node.name)
        
        if node.is_leaf():
            proteins, nogs = load_eggNOG_nog_mapping(level_id, old_definition)
            protein_mapping[level_id] = proteins
            nog_mapping[level_id] = nogs
        else:
            for child in node.traverse():
                child_id = int(child.name)
                proteins, nogs = load_eggNOG_nog_mapping(child_id, new_definition)
                protein_mapping[child_id] = proteins
                nog_mapping[child_id] = nogs    
            
    return protein_mapping, nog_mapping

def load_subtree_data(higher_level,new_definition):
    eggNOG_tree = Tree(eu.get_eggNOG_nhx() ,format=8)
    higher_node = eggNOG_tree.search_nodes(name=str(higher_level))
    assert len(higher_node) == 1
    higher_node = higher_node[0]
    
    nog_mapping = dict()
    protein_mapping = dict()
    
    # load new definitions
    for node in higher_node.traverse():
        level_id = int(node.name)
        proteins, nogs = load_eggNOG_nog_mapping(level_id, new_definition)
        protein_mapping[level_id] = proteins
        nog_mapping[level_id] = nogs    
            
    return protein_mapping, nog_mapping

def load_reconciliation_data(higher_level,new_definition):
    """
    Loads the [OLD] nog definitions for the level itself
    and the immidiate leaf children. A [NEW] definition is
    instead used in case of children that are inner nodes,
    assuming that these levels were computed in a previous
    reconciliation step.
    
    e.g. when loading higher_level = 314146 [spriNOG]
    
                    /9443[prNOG][NEW]---9604[homNOG][NOT_LOADED]
    -314146[spriNOG][OLD]
                    \-9989[roNOG][OLD]
                    
    returns:
    protein_mapping [level_id][proteinID]:[nogID]
    nog_mapping     [level_id][nog_id]:set([protein_id])
    
    """
    eggNOG_tree = Tree(eu.get_eggNOG_nhx() ,format=8)
    higher_node = eggNOG_tree.search_nodes(name=str(higher_level))
    assert len(higher_node) == 1
    higher_node = higher_node[0]
    
    nog_mapping = dict()
    protein_mapping = dict()
    
    # load yourself from old definition
    proteins, nogs = load_eggNOG_nog_mapping(higher_level)
    protein_mapping[higher_level] = proteins
    nog_mapping[higher_level] = nogs
    
    # load leaf levels from old and intermediate from new
    for node in higher_node.get_children():
        
        level_id = int(node.name)
        
        if node.is_leaf():
            proteins, nogs = load_eggNOG_nog_mapping(level_id)
        else:
            proteins, nogs = load_eggNOG_nog_mapping(level_id, new_definition)
            
        protein_mapping[level_id] = proteins
        nog_mapping[level_id] = nogs
            
    return protein_mapping, nog_mapping  
    
def test_mapping_conversion(protein_mapping):
    # testing conversion vs saved pickle
    
    # Read reverse eggNOG tree, parent -> set(children)
    eggNOG_tree_rev = eu.read_eggNOG_treeRev()
    
    t_start = time.clock()
    protein_mapping2 = dict()
    protein_mapping2[max_level] = read_pickle(
        '%s/%d.proteinINT.nogINT.pkl2'%(V4_PICKLE_DIR,max_level))
    for child_level in eggNOG_tree_rev[max_level]:
        protein_mapping2[child_level] = read_pickle(
            '%s/%d.proteinINT.nogINT.pkl2'%(V4_PICKLE_DIR,child_level))
    t_end = time.clock()
    sys.stderr.write("Pickle conversion in %.2f s\n"%(t_end - t_start))
    
    #test equality
    for level_id in protein_mapping:
        mapping1 = set(protein_mapping[level_id])
        mapping2 = set(protein_mapping2[level_id])
        if not mapping1.issubset(mapping2) or not mapping2.issubset(mapping1):
            sys.stderr.write("Difference found @ level %d: %d\n"%(
                level_id, len(mapping1 - mapping2)))
    
    return protein_mapping, protein_mapping1
    
def load_eggNOG_data(max_level):
    """
    Precache eggnog data in memory:
    protein_nog_mapping[level_id][protein_id] = nog_id
    nog_protein_mapping[level_id][nog_id] = set(protein_id)
    protein_name_mapping[protein_id] = [species_id]_[protein_str]
    
    """
    
    protein_nog_mapping = {}
    nog_protein_mapping = defaultdict(dict)
    protein_name_mapping = {}
    
    eggNOG_tree = Tree(eu.get_eggNOG_nhx(),format=8)
    
    #load protein_id -> nog_id mapping
    max_node = eggNOG_tree.search_nodes(name=str(max_level))[0]
    
    levels_to_load = 0
    for n in max_node.traverse():
        levels_to_load+=1
        
    sys.stderr.write("Started loading %d levels..\n"%levels_to_load)
    t_start = time.clock()
    pc = fu.ProgressCounter(levels_to_load)
    
    for n in max_node.traverse():
        
        # load protein -> nog_id mapping
        
        level_id = int(n.name)
        shelve_file = '%s/%s.40/%s.og_id_prot.db'%(
            ds.EGGNOG_OUTPUT,'og_protein_mapping',level_id)
        
        try:
            s = shelve.open(shelve_file)
            protein_nog_mapping[level_id] = dict(s)
        except Exception as e:
            sys.stderr.write("Exception occurred!\n %s"%e)
        finally:
            s.close()
        
        #load nog_id -> protein_id list mapping
        
        nog_protein_dir = '%s/%s'%(
            ds.EGGNOG_OUTPUT,'og_protein_mapping')        
        nog_protein_file = "%s/%s.og_id_prot.stringent.tsv"%(
            nog_protein_dir,level_id)

        current_id = None
        with open(nog_protein_file) as f:
            for line in f:
                l = line.rstrip().split()
                nog_id = int(l[0])
                if nog_id != current_id:
                    nog_protein_mapping[level_id][nog_id] = set()
                    current_id = nog_id
                nog_protein_mapping[level_id][nog_id].add(l[1])
        
        pc.count()
        pc.show_progress
    
    t_end = time.clock()
    sys.stderr.write("Finished loading NOG data in %.2f s\n"%(t_end-t_start))
    
    #load protein names shelve
    protein_id_shelve = "%s/%s"%(ds.EGGNOG_OUTPUT,'eggNOG.pInt.pSpecies_pShort.shelve')
    sys.stderr.write("Starting to load protein names\n")
    t_start = time.clock()
    try:
        s = shelve.open(protein_id_shelve)
        protein_name_mapping = dict(s)
    except Exception as e:
        sys.stderr.write("Exception occurred!\n %s\n"%e)
    finally:
        s.close()
    t_end = time.clock()
    sys.stderr.write("Finished loading protein names in %.2f s\n"%(t_end - t_start))
    
    return  (protein_nog_mapping , nog_protein_mapping , protein_name_mapping)

def get_protein_fasta(protein_query_dict,header_format=None):
    """
    method to retrieve the fasta sequences of all proteins in the
    query dictionary.
    
    Input should be of format:
    dict[species_id] = set(protein_short_names)
    
    Returns dictionary of format:
    dict[protein_name] = fasta_entries
    
    """

    fasta_entries = {}
    
    for species_id in protein_query_dict:
        
        fasta_file = "%s/%s.fa"%(V4_FASTA_DIR,species_id)
        
        protein_set = protein_query_dict[species_id]
        
        current_protein = None
        
        with open(fasta_file) as f:
            for line in f:
                
                if(current_protein):
                    if(line.startswith(">")):
                        protein_set.remove(current_protein)
                        current_protein = None
                    else:
                        fasta_entries[current_protein] += line
                        
                if(protein_set):
                    for protein_name in protein_set:    
                        if(line.startswith(">%s"%protein_name)):
                            
                            #safe fasta with new id: [speciesID]_[proteinName]
                            new_start = len(protein_name) + 1 #plus one for starting ">"
                            
                            if header_format == '.':
                                new_header = ">%s.%s\n"%(species_id,protein_name)
                            else:
                                new_header = ">%s_%s%s"%(species_id,protein_name,line[new_start:])
                            
                            current_protein = protein_name
                            
                            # safety check to not override something
                            assert current_protein not in fasta_entries, 'WARNING! Attempting to override %s by %s.%s'%(
                                fasta_entries[current_protein].split()[0],species_id,current_protein)
                            
                            fasta_entries[current_protein] = new_header
                            break
                else:
                    break
                    
    return fasta_entries

def get_species_fasta(species_id):
    """
    get all fasta sequences from the input species
    
    return dict[protein_id]=[fasta_seq]
    
    """
    fasta_entries = {}
    fasta_file = "%s/%s.fa"%(V4_FASTA_DIR,species_id)
    
    current_protein = None
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                current_protein = line.split()[0][1:]
                fasta_entries[current_protein] = ">%s\n"%current_protein
            else:
                assert current_protein
                # fix for Selenocysteine
                fix_line = line.replace('U','X')
                # fix for Pyrrolysine encoded by the 'amber' stop codon UAG in bacteria
                fix_line = fix_line.replace('O','X')
                fasta_entries[current_protein] += fix_line
                
    return fasta_entries

def get_level_fasta(level_id):
    """
    get all fasta sequences from all species in at the input level
    
    return dict[species_id][protein_id]=[fasta_seq]
    """
    
    level_species = eu.read_level_species()
    
    assert level_id in level_species, 'level_id %s not found in eggNOG'%level_id
    
    sys.stderr.write("Started loading FASTA data for %d species..\n"%len(level_species[level_id]))
    t_start = time.clock()
    
    level_fasta = {}
    for species_id in level_species[level_id]:
        level_fasta[species_id] = get_species_fasta(species_id)

    t_end = time.clock()
    sys.stderr.write("Finished loading FASTA data in %.2f s\n"%(t_end-t_start))
        
    return level_fasta

def filter_fusions(fusions,fusion_level,protein_mapping,nog_mapping):
    """ filter fusions based on a size criteria, the protein fusion_id is kept
    for the biggest NOG keeps while it will be removed from the other, smaller NOGs"""
    to_keep = {}
    to_remove = {}

    for fusion_id, fusion_nogs in fusions.items():
    
        assert fusion_id in protein_mapping[fusion_level],'fusion %d not found @ %d'%(fusion_id,fusion_level)
    
        higher_proteins = {x:nog_mapping[fusion_level][x] for x in fusion_nogs}
    
        sorted_by_size = Counter({x:len(y) for x,y in higher_proteins.items()}).most_common()
        to_keep[fusion_id] = sorted_by_size[0][0]
        to_remove[fusion_id] = {x[0] for x in sorted_by_size[1:]}
    
    return to_keep, to_remove 

def remove_fusions(higher_level,protein_nogs,nog_proteins):
    # load fusions, i.e. dict contains the fusion_id:list(fusion_nogs)
    fusions = eu.read_eggNOG_fusions(higher_level) 
    fusions_to_keep, fusions_to_remove = filter_fusions(fusions, higher_level, protein_nogs, nog_proteins)
    
    # apply filtering        
    for fusion_id, fusion_nog in fusions_to_keep.items():
            
        # ensure correct mapping of protein
        protein_nogs[higher_level][fusion_id] = fusion_nog
    
        # remove other associations
        for nog_id in fusions_to_remove[fusion_id]:
            nog_proteins[higher_level][nog_id].remove(fusion_id)
    
    return protein_nogs, nog_proteins

class HgtTreeClassifier:

    def __init__(self,
                 eggNOG_lca_finder,
                 eggNOG_names,
                 verbose=False):

        #fields
        self.lca_finder = eggNOG_lca_finder
        self.names = eggNOG_names
        self.level_ids = {str_id: int_id for (int_id,str_id) in eggNOG_names.items()}
        self.visit_order = 1;
        self.locked_nodes = set();
        self.geneTree = None

        #classifications
        self._classification = {}
        self._hgtFusions = {}
        
        #flags
        self.verbose = verbose
        self.debug = False
    
    #read-only property shortcut: https://docs.python.org/2/library/functions.html#property
    
    @property
    def classification(self):
        return self._classification
    
    @property
    def hgtFusions(self):
        return self._hgtFusions

    def addLevelMapping(self,eggNOG_level,tree_node):
        """
        Creates a new correspondence to an eggNOG level
        and adds the level to the classification if
        it is detected for the first time
        
        If the list stays empty tree_node will serve
        as NOG
        """
        
        if eggNOG_level not in self._classification:
            self._classification[eggNOG_level] = {}
        
        #List will be empty if no duplication downstream
        if tree_node not in self._classification[eggNOG_level]:
            self._classification[eggNOG_level][tree_node] = []
        
    def addLevelNOGs(self,lca_level,tree_node,nog_nodes):
        """
        Adds a new NOG to the specified lca_level+node
        """
        self._classification[lca_level][tree_node].extend(nog_nodes)
        
    def updateLevelNOGs(self,lca_level,tree_node,old_nog_node,new_nog_nodes):
        """
        Check if old_nog_nodes were inserted (e.g. as branch of a speciation)
        in lca_level, if yes remove it and insert the new new_nog_nodes
        """
        if(old_nog_node in self._classification[lca_level][tree_node]):
            self._classification[lca_level][tree_node].remove(old_nog_node)
        self.addLevelNOGs(lca_level,tree_node,new_nog_nodes)
        
    def hasSpeciationLock(self,lca_level,tree_node,ancestors_of_nog_candidate):
        """
        The locked nodes are checked to avoid the insertion of a
        one sided duplication (i.e. on only one branch after the a speciation)
        because of insufficient coverage with eggNOG levels
        """
        #check if there is a more recent locked ancestor within
        for sub_node in self._classification[lca_level][tree_node]:
            if(sub_node in ancestors_of_nog_candidate):
                if sub_node.event == 'S':
                    if(self.verbose):
                        print("\trefOG insertion in %s blocked by speciation %s"%(
                            self.names[lca_level],sub_node.name))
                    return True
        
        #No locked speciation node found
        return False
    
    def mapChildrenLevelThatDiffer(self,level_lca, node, children):
        """
        If the children's lca level is different then the parental
        one they are inserted as TreeNode correspondences for their level
        
        Children with equal eggNOG lca are not mapped because
        they do not bring additional information that can be captured
        by eggNOG. An eggNOG level with higher resolution would be needed.
        """
        def map_levels(child_lca):
            self.addLevelMapping(child_lca, child)
            
            level_lca_id = self.level_ids[level_lca]
            child_lca_id = self.level_ids[child_lca]
            
            intermediate_levels = self.lca_finder.get_intermediateLevels(level_lca_id,child_lca_id)        
            
            for level_id in intermediate_levels:
                level_str = self.names[level_id]
                self.addLevelMapping(level_str,child)
            
        for child in children:
            if level_lca != child.lca:
                map_levels(child.lca)

    
    def analyze_complete_geneTree(self, gTree):
        """
        complete pipeline to transform an annotated gene Tree
        (i.e. every node has .event=D|T|S)
        """
        
        #1. annotate the gene Tree with lca levels
        self.annotate_lca_levels(gTree)
        
        #2. form new nogs
        new_nogs, unassigned_acceptors = self.form_new_nogs(gTree)
        
        #3. resolve unassigned acceptors
        self.assign_acceptors(unassigned_acceptors,new_nogs,gTree)
        
        return new_nogs
        
    def form_new_nogs(self,tree_node):
           
        # Transition elements to be excluded in initial nog assignment    
        acceptors = set()
        
        # Leaf function that stops on acceptor nodes
        def is_acceptor(node):
            if hasattr(node,'is_acceptor'):
                return node.is_acceptor
            else:
                return False
        
        
        
        # Assign nogs by traversing from top to bottom
        # (likewise to eggnas.utils.classification_utils.py)
        # NOTE: unless root, a node always assigns the children and not itself!
        for node in tree_node.traverse('preorder', is_leaf_fn=is_acceptor):
            
            if(self.verbose):
                print("%d Visiting: %s "%(self.visit_order,node.name))
            self.visit_order += 1
            
            if node.is_leaf():
                continue
            
            if is_acceptor(node):
                acceptors.add(node)
                continue
            
            lca_level = node.lca
            
            if node.is_root():
                self.addLevelMapping(lca_level,node)
        
            #Classify children based on event type of node
            children = node.get_children()
            
            if node.event == 'D':
                #Legacy algorithm:
                
                #To insert the refOG at the right position the
                #ancestor match has to be found
                matched_ancestor = None
                
                if node in self._classification[lca_level]:
                    self.addLevelNOGs(lca_level,node,children)
                    matched_ancestor = node
                else:
                    ancestors = node.get_ancestors()
                    for mapping_node in self._classification[lca_level]:
                        if mapping_node in ancestors:
                            matched_ancestor = mapping_node
                            if matched_ancestor.event == 'D':
                                #find if there was a speciation in between
                                #matched_ancestor and current dup.node that
                                #would block the insertion
                                if not self.hasSpeciationLock(lca_level,matched_ancestor,ancestors):
                                    self.updateLevelNOGs(lca_level,matched_ancestor,node,children)
                            break
                
                assert(matched_ancestor is not None),"Missing ancestor in %s (%d)"%(
                    self.names[lca_level],
                    len(self._classification[lca_level]))
                
                #assign children levels that differ
                self.mapChildrenLevelThatDiffer(lca_level,node,children)
                
            elif node.event == 'T':
                
                #only assign transfer donor 
                donor = [child for child in children if not child.is_acceptor]
                
                #assert(len(donor) == 1), "%d donor where found: %s"%(len(donor),donor)
                #
                ##detect whether node was inserted a NOG and in case
                ##insert donor
                #matched_ancestor = None
                #
                ##find ancestor inserted for the level
                #if node in self._classification[lca_level]:
                #    self.addLevelNOGs(lca_level,node,donor)
                #    matched_ancestor = node
                #else:
                #    ancestors = node.get_ancestors()
                #    for mapping_node in self._classification[lca_level]:
                #        if mapping_node in ancestors:
                #            matched_ancestor = mapping_node
                #            self.updateLevelNOGs(lca_level,matched_ancestor,node,donor)
                #            break
                #
                #assert(matched_ancestor is not None),"Missing ancestor in %s (%d)"%(
                #    self.names[lca_level],
                #    len(self._classification[lca_level]))
                #
                ##assign donor level that differs and intermediate ones
                #self.mapChildrenLevelThatDiffer(lca_level,node,donor)
                    
            elif node.event == 'S':
                self.mapChildrenLevelThatDiffer(lca_level,node,children)
            else:
                print("UNKNOWN event: %s!"%node.event)
                sys.exit()      
            
            #could add a flag to the node that no further exploration
            #should be performed (with is_leaf_fn option)
    
        return acceptors

    def assign_acceptors(self, acceptors, gTree):
        #perform duplications checks for acceptors with lca levels already
        #present in new_nogs
        #also consider intermdeiate levels wrt up [TODO]
        for acceptor in acceptors:
            
            if(self.verbose):
                print("Assigning acceptor: %s"%acceptor.name)
            
            if acceptor.lca in self._classification:
                
                #Logic: try fusing the acceptor to the closest compatible node
                nearest_node = None
                nearest_node_dist = 10000
                
                for nog_node in self._classification[acceptor.lca]:
                    
                    ancestor = gTree.get_common_ancestor([nog_node,acceptor])
                    
                    node_distance = gTree.get_distance(ancestor, nog_node, topology_only=True)
                    
                    if node_distance < nearest_node_dist:
                        nearest_node_dist = node_distance
                        nearest_node = ancestor
                    
                if nearest_node:
                    
                    ancestor = nearest_node
                    if(self.verbose):
                        print("Common ancestor %s has level %s"%(ancestor.name,ancestor.lca))
                    
                    if ancestor.lca == acceptor.lca:
                        if(self.verbose):
                            print("Will fuse to nog_node: %s"%nog_node.name)
                        
                        if self._classification[acceptor.lca][nog_node]:
                            #i.e. there is already a NOG sub-definition, e.g. duplications
                            self.addLevelNOGs(acceptor.lca, nog_node, [acceptor])
                        else:
                            if acceptor.lca not in self._hgtFusions:
                                self._hgtFusions[acceptor.lca] = {}
                            
                            #i.e. no pre-existing NOG sub-definition. Add nog_node itself to create fusion
                            if nog_node in self._hgtFusions[acceptor.lca]:
                                self._hgtFusions[acceptor.lca][nog_node].append(acceptor)
                            else:
                                self._hgtFusions[acceptor.lca][nog_node] = [acceptor]
                            #could this be moved to ancestor node rather than nog_node?
                        #todo fuse operation
                    else:
                        self.addLevelMapping(acceptor.lca,acceptor)
                        if(self.verbose):
                            print("Will be a new NOG")
                        
            else:
                #Lca level wasn't detected yet (e.g. transfer to distant kingdom)
                self.addLevelMapping(acceptor.lca,acceptor)
                if(self.verbose):
                    print("Will be a new NOG")
            
            #can the node be further classified? (i.e. acceptor not terminal leaf)
            #unblock head node (i.e. is_acceptor)
            acceptor.is_acceptor = False
            new_acceptors = self.form_new_nogs(acceptor)
            
            if(self.verbose):
                self.report_classification()
                print "acceptors to resolve:\t%s"%new_acceptors
            
            if new_acceptors:
                self.assign_acceptors(new_acceptors,gTree)

    def duplication_on_path(self,node, other, gTree):
        """
        verifies whether between node and other a duplication exists
        duplications will stop Transfer events from fusing
        """
        #TODO: identify the location of the lca of the fusion partners, should
        #      lie in the same level to comply to the orthologous group definition
        #TODO: here we need the original geneTree or at least recover it!
        
        ancestor = gTree.get_common_ancestor([node,other])
        
        for parent in [node,other]:
            #test if path from parent to lca contains duplications
            while ancestor != parent:
                if parent.event == 'D':
                    return True
                parent = parent.up
                
        return False
    
    def annotate_proteins_wrt_lca(self,gTree):
        """
        Annotates the nodes with the respective proteins
        they would contain should they be taken as a NOG
        """
        for descendant in gTree.traverse("postorder"):
            
            if descendant.is_leaf():
                protein_name = descendant.name #for real name: .replace('_','.',1)
                descendant.add_features(proteins=[protein_name])
            else:
                descendant.add_features(proteins=[])
                for child in descendant.get_children():
                    
                    #do not consider proteins coming from an acceptor branch
                    if descendant.event == 'T':
                        if child.is_acceptor:
                            continue
                        
                    descendant.proteins.extend(child.proteins)
                    
    def report_classification(self):
        for level in self._classification:
            print(level)
            for node in self._classification[level]:
                print("\t%s %s"%(node.name,[x.name for x in self._classification[level][node]]))
            if level in self._hgtFusions:
                for node in self._hgtFusions[level]:
                    print("\tfusion: %s %s"%(node.name,[x.name for x in self._hgtFusions[level][node]]))
                
    
    def annotate_lca_levels(self,gTree,include_transfers=False):
        """
        Annotate the nodes of the tree with the minimal level of eggnog
        """
    
        eggNOG_names = self.names
        lca_level_finder = self.lca_finder
        
        detected_lca_levels = set()
        
        for node in gTree.traverse("postorder"):
           
            if node.is_leaf():
                #initialize leaves
                leaf_species_id = node.name.split('_')[0]
                node.add_features(species=[int(leaf_species_id)],
                                  mapping=leaf_species_id)
            else:
                
                assert(hasattr(node,'event')), "Something is wrong, node is missing event"   
                node.add_features(species=[])
                
                for child in node.get_children():
                    #exclude species from acceptor branch for transitions
                    if node.event == 'T':
                        if hasattr(child,'is_acceptor'):
                            if child.is_acceptor:
                                continue
                        else:
                            assert(hasattr(child,'mapping')), "Something is wrong, child %s has no mapping"%child.name
                            if child.mapping == node.acceptor:
                                child.add_features(is_acceptor=True)
                                continue
                            else:
                                child.add_features(is_acceptor=False)
                            
                        #else: 
                        #    assert(child.mapping == node.donor),'Mismatch in donor/child.mapping @ %s ! [child]%s != %s[donor]'%(node.name,child.mapping, node.donor)
                        #currently off because of: AssertionError: Mismatch in donor/child.mapping @ m90 ! [child]n60 != n45[donor]
                    node.species.extend(child.species)
            
            
            #Compute the to last common ancestor eggNOG levels
            #one including the species in the transfer acceptor (wT)
            #one without the species in the transfer acceptor (woT)
            
            if node.is_leaf():
                all_species = node.species
            else:
                all_species = [x.species[0] for x in node.get_leaves()]
            
            if include_transfers:
                lca_level = lca_level_finder.find_level(all_species)
            else:
                lca_level = lca_level_finder.find_level(node.species)
                
            lca_level_name = eggNOG_names[lca_level]
            
            detected_lca_levels.add(lca_level_name)
            
            node.add_features(lca=lca_level_name)
        
        return detected_lca_levels
    
