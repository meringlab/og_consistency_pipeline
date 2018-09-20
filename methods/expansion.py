import sys
import os
import time
import argparse
import shelve

from collections import defaultdict, OrderedDict, Counter
from ete3 import Tree

from .utils import data_sources as ds
from .utils import eggNOG_utils as eu
from .utils import file_utils as fu

from .hgt_utils import get_protein_fasta, load_v4clades

class ExpandedNOG:
    # 'Class representing an expanded NOG'
    
    eggNOG_dict = eu.read_eggNOG_tree()
    eggNOG_dict_rev = eu.read_eggNOG_treeRev()
    eggNOG_names = eu.read_eggNOG_names()
    eggNOG_tree = Tree(eu.get_eggNOG_nhx(),format=8)
    eggNOG_short_names = eu.compute_first_names()
    protein_id_shelve = "%s/%s"%(ds.EGGNOG_OUTPUT,'eggNOG.pInt.pSpecies_pShort.shelve')
    
    use_permanence = True
    protein_nog_mapping = {}
    nog_protein_mapping = defaultdict(dict)
    protein_name_mapping = {}
    
    # [TODO] switch position of nog_level and nog_name
    def __init__(self, nog_level, nog_name=None, max_level=-1, verbose=False):
        
        self.verbose=verbose
        
        nog_level = str(nog_level)
        
        if nog_name is None: 
            #nog_level is a nog_int identifier. Extract level id and reassign
            nog_name = nog_level
            level_id, nog_number = eu.decompose_nog_id(nog_name)
            nog_level = level_id
        
        #initialize expanded NOG
        if self.verbose:
            print('Initializing new ExpandedNOG based on %s from %s'%(
                nog_name, nog_level))
        
        self.input_name = nog_name
        self.input_level = int(nog_level)
        if max_level == -1:
            self.higher_level = self.eggNOG_dict[self.input_level]
        else:
            self.higher_level = int(max_level)
            
            if self.higher_level != 1:
                assert self.higher_level in self.eggNOG_dict, "Unknown max_level chosen: %d"%self.higher_level
            
            #Check that self.higher_level is parent of input_level
            parent_level = self.eggNOG_dict[self.input_level]
            while parent_level in self.eggNOG_dict:
                if parent_level == self.higher_level:
                    break
                else:
                    parent_level = self.eggNOG_dict[parent_level]
            if parent_level != self.higher_level:
                print('Chosen max_level is not parent of input nog_level!')
                
        
        #initialize containers
        self.expanded_nog = defaultdict(set)
        self.protein_species = {}
        self.higher_nogs = defaultdict(set)
        self.lower_nogs = defaultdict(set)
        self.singletons = defaultdict(set)
        self.directly_attached = defaultdict(set)
        self.already_reversed = set()
        self.already_expanded = set()
        self.nogs_history = OrderedDict() # e.g. split 10066560306208 -> set{100665600000001,100665600000002}
        
        #0. resolve nog name to nog integer (6656_NOG05039 -> 10066560306208)
        if str(nog_name).isdigit():
            self.input_id = int(nog_name)
        else:
            nog_index_shelve = '%s/%s'%(ds.EGGNOG_OUTPUT,'eggNOGv40.levelIDnogSTR.nogINT.shelve')
            shelve_id = "%s%s"%(nog_level,nog_name)
            self.input_id = fu.shelve_lookup(shelve_id,nog_index_shelve)
    
    def __str__(self):
        return (
            "ExNOG:\t%s_%s\n"%(self.input_level,self.input_name)
            + "Start size:\t%d\n"%len(self.expanded_nog[self.input_id])
            + "Total size:\t%d\n"%self.get_size()
            + "H-nogs:\t%d\n"%len(self.higher_nogs)
            + "L-nogs:\t%d\n"%len(self.lower_nogs))
    
    def pretty_print(self):
        def ex_len(x):
            return "%d|%d"%(x,len(self.expanded_nog[x]))

        print('\n'.join(
            ["%s => %s"%(ex_len(x),
                         [ex_len(z) for z in y])
             for x,y in self.lower_nogs.items()]))
        
    def set_cache(self, protein_nog_mapping, nog_protein_mapping, protein_name_mapping):
        """
        Optional method to speed up computations by preloading the data needed
        to expand the NOG
        """
        self.protein_nog_mapping = protein_nog_mapping
        self.nog_protein_mapping = nog_protein_mapping
        self.protein_name_mapping = protein_name_mapping
        self.use_permanence = False
        
    def remove_cache(self):
        self.protein_nog_mapping = None
        self.nog_protein_mapping = None
        self.protein_name_mapping = None
        self.use_permanence = True
    
    def build(self, complete_expansion = False):
        
        #1. get proteins
        nog_proteins = self.get_nog_proteins(self.input_id)
        if self.verbose:
            print("Initial size: %d"%len(self.get_full_protein_set()))
        
        #2. initial expansion
        to_expand = set([self.input_id]) # self.expand(self.input_id)
        to_reverse = set()
        if complete_expansion:
            to_reverse.add(self.input_id)

        # if self.verbose:
        #    print("Size after first expansion: %d"%len(self.get_full_protein_set()))
        
        #3. Reverse and expand util all nogs are fully explored
        while(to_reverse or to_expand):
            
            if self.verbose:
                print("New nogs to reverse expand: %s"%to_reverse)
                print("New nogs to reverse expand: %d"%len(to_reverse))
            
            if complete_expansion:
                # Check for new nogs to expand among nogs to reverse, i.e.
                # add expansion to H level [in order to obtain H+1]
                for nog_id in to_reverse:
                    if nog_id not in self.already_expanded:
                        level_id = self.get_level_id(nog_id)
                        if level_id != self.higher_level:
                            to_expand.add(nog_id)
            
            # Do reverse expansion    
            while to_reverse:
                nog_id = to_reverse.pop()
                new_nogs = self.reverse(nog_id)
                to_expand.update(new_nogs)
        
            if self.verbose:
                print("Size after reverse expansion: %d"%len(self.get_full_protein_set()))
        
            if self.verbose:
                print("New nogs to expand: %s"%to_expand)
                print("New nogs to expand: %d"%len(to_expand))
            
            if complete_expansion:
                # Check for new nogs to reverse among nogs to expand, i.e.
                # add reverse expansion to L level [in order to obtain L-1]
                for nog_id in to_expand:
                    if nog_id not in self.already_reversed:
                        level_id = self.get_level_id(nog_id)
                        if level_id in self.eggNOG_dict_rev: # [TODO] change safety condition to minimum lower level (like higher)
                            to_reverse.add(nog_id)

            while to_expand:
                nog_id = to_expand.pop()
                new_nogs = self.expand(nog_id)
                to_reverse.update(new_nogs)
            
            if self.verbose:    
                print("Size after expansion: %d"%len(self.get_full_protein_set()))
    
    def pre_split(self,transfer_thr=0.1, cogs=set()):
        """ pre split expanded nog network based on transfer threshold (percentage)"""
        split_counter = 1
        split_results = {}
        for nog_id in sorted(self.find_inconsistencies()):
            nog_proteins = self.expanded_nog[nog_id]
            transfer_min = max(len(nog_proteins) * transfer_thr,2) # split always singletons
            for higher_id in sorted(self.lower_nogs[nog_id]):
                
                if higher_id in cogs:
                    # do not split COGs
                    continue
                
                higher_proteins = self.expanded_nog[higher_id]
                if nog_proteins.issuperset(higher_proteins):
                    # no further action required since end-point-nog
                    continue
                
                if len(nog_proteins & higher_proteins) < transfer_min:
                    split_id = eu.get_copy_id(higher_id,split_counter,6)
                    if self.verbose:
                        sys.stderr.write('%d:%d PRE-SPLIT: %d|%d -> %d|%d\n'%(self.input_id,nog_id,
                                                                        higher_id,
                                                                        len(self.expanded_nog[higher_id]),
                                                                        split_id,
                                                                        len(nog_proteins & higher_proteins)))
                    
                    split_nogs = self.partial_split(higher_id,nog_id,split_id)
                    split_counter += len([x for x,y in split_nogs.items() if len(y) > 1])
                    split_results.update(split_nogs)
        
        return split_results
    
    def pre_merge(self,transfer_thr=0.1):
        """ pre merge all singletons or contributions below 10% to the NOG with the largest overlap"""
        merge_counter = 1
        merge_results = {}
        for nog_id in sorted(self.find_inconsistencies()):
            nog_proteins = self.expanded_nog[nog_id]
            transfer_min = max(len(nog_proteins) * transfer_thr,2) # 2 be the minimum to merge singletons
            higher_ids = sorted(self.lower_nogs[nog_id])
            transfers = Counter({x:len(nog_proteins & self.expanded_nog[x]) for x in higher_ids})
            top_id, top_size = transfers.most_common(1)[0]
            org_id = top_id
            #if top_size > transfer_min:
            # given a top transfer nog with sufficiently big overlap, merge too small nogs (e.g. singletons) into it
            for higher_id in sorted(transfers.keys()):
                
                if higher_id == org_id:
                    # don't add yourself
                    continue
                
                if transfers[higher_id] >= transfer_min:
                    # don't merge larger groups that transfer the minimum or more
                    continue
                
                if len(self.higher_nogs[higher_id]) > 1:
                    # don't merge groups that have more than 1 input
                    continue
                                
                merged_id = eu.get_copy_id(higher_id,merge_counter,3)
                if self.verbose:
                    sys.stderr.write('%d:%d PRE-MERGE: %d|%d + %d|%d -> %d\n'%(self.input_id,nog_id,
                                                                            higher_id,len(self.expanded_nog[higher_id]),
                                                                            top_id,len(self.expanded_nog[top_id]),
                                                                            merged_id))
                merged_nog = self.merge_all([higher_id,top_id],merged_id)
                merge_counter += 1
                merge_results.update(merged_nog)
                top_id = merged_id
    
        return merge_results
    
    def contains(self, nog_id):
        return nog_id in self.expanded_nog
    
    def get_expansion_keys(self):
        return set(self.expanded_nog.keys())
    
    def get_name(self):
        return "%s_%s"%(self.input_level,self.input_name)
    
    def get_size(self):
        return len(self.get_full_protein_set())
    
    def get_level_id(self,nog_id):
        level_id, nog_number = eu.decompose_nog_id(nog_id)
        return level_id
    
    def get_full_protein_set(self):
        full_protein_set = set()
        for nog_proteins in self.expanded_nog.values():
            full_protein_set.update(nog_proteins)
        return full_protein_set
    
    def get_protein_mapping_wo_singletons(self,output_type=str):
        # Output proteins mapping to nogs for all levels (w/o singletons)
        protein_mapping = defaultdict(dict)
        for nog_id, nog_proteins in self.expanded_nog.items():
            if len(nog_proteins) > 1:
                level_id = self.get_level_id(nog_id)
                #output conversion
                nog_proteins = {output_type(x) for x in nog_proteins}
                for protein_id in nog_proteins:
                    assert protein_id not in protein_mapping[level_id], 'Double insertion for %s @ %d: %d vs %d'%(
                        protein_id,level_id,nog_id,protein_mapping[level_id][protein_id] )
                    protein_mapping[level_id][protein_id] = nog_id
                    
        return protein_mapping

    def get_nog_proteins(self,nog_id):
        nog_proteins = set()
        level_id = self.get_level_id(nog_id)
    
        if self.use_permanence:
            nog_protein_dir = '%s/%s'%(
                ds.EGGNOG_OUTPUT,'og_protein_mapping')        
            nog_protein_file = "%s/%s.og_id_prot.stringent.tsv"%(
                nog_protein_dir,level_id)
        
            with open(nog_protein_file) as f:
                for line in f:
                    if line.startswith(str(nog_id)):
                        l = line.rstrip().split()
                        nog_proteins.add(l[1])
        
        else:
            nog_proteins = self.nog_protein_mapping[level_id][nog_id]
        
        self.expanded_nog[nog_id] = nog_proteins
        return nog_proteins
    
    def get_protein_nog(self,protein_id, level_id):
        
        if self.use_permanence:
            og_db_file = '%s/%s.40/%s.og_id_prot.db'%(
                ds.EGGNOG_OUTPUT,'og_protein_mapping',level_id)
                
            nog_id = fu.shelve_lookup(protein_id,og_db_file) # Returns None if not present
        else:
            if protein_id in self.protein_nog_mapping[level_id]:
                nog_id = self.protein_nog_mapping[level_id][protein_id]
            else:
                nog_id = None
            
        return nog_id
    
    def get_nog_id_changes(self):
        return self.nogs_history
    
    def get_cog_id_correspondence(self,cog_id):
        
        # check if the id is currently used
        if cog_id in self.expanded_nog:
            return cog_id
        else:
            assert cog_id in self.nogs_history, "COG %d has no reference in this exNOG.%d"%(cog_id,self.input_id)
        
        # Find current denomination by searching in the nog_history
        new_id = cog_id
        while new_id in self.nogs_history:
            new_id = self.nogs_history[new_id]
            assert len(new_id) == 1, "WARNING: COG %d was split into %s"%(cog_id,new_id)
            new_id = list(new_id)[0]
        
        return new_id
    
    def expand(self, nog_id):
        
        #expansion step
        level_id = self.get_level_id(nog_id)
        to_reverse = set()
        self.already_expanded.add(nog_id)
        
        higher_level = self.eggNOG_dict[level_id]
        for protein_id in self.expanded_nog[nog_id]:
            
            higher_nog = self.get_protein_nog(protein_id,higher_level)
            
            if higher_nog:
                self.higher_nogs[higher_nog].add(nog_id)
                self.lower_nogs[nog_id].add(higher_nog)
                
                if higher_nog not in self.expanded_nog:
                    self.get_nog_proteins(higher_nog)
                    to_reverse.add(higher_nog)
            else:
                higher_nog = eu.get_nog_id(higher_level,protein_id,nog_type=2) # [TODO] rename to singleton_nog (starts with 2)
                self.higher_nogs[higher_nog].add(nog_id)
                self.lower_nogs[nog_id].add(higher_nog)
                
                if protein_id not in self.expanded_nog[higher_nog]:
                    self.expanded_nog[higher_nog].add(protein_id)
                    self.singletons[higher_level].add(protein_id)
                    to_reverse.add(higher_nog)
            
            #debug if 'singletons' in str(nog_id):
            #    print("%s -> %s"%(nog_id,higher_nog))
        
        return to_reverse

    def reverse(self, nog_id):
        #reverse expansion
        level_id = self.get_level_id(nog_id)
        to_expand = set()
        self.already_reversed.add(nog_id)
        
        for protein_id in self.expanded_nog[nog_id]:
            
            # 0. Find out the species of the protein to reduce search space (to 1 out of n-children levels)
            if protein_id not in self.protein_species:
                if self.use_permanence:
                    protein_species,protein_short = fu.shelve_lookup(protein_id, self.protein_id_shelve)
                else:
                    protein_species,protein_short = self.protein_name_mapping[protein_id]
                    
                self.protein_species[protein_id] = protein_species
                
            species_id = self.protein_species[protein_id]
            
            # Check if species is directly attached to level, i.e. no L level
            if self.eggNOG_dict[species_id] == level_id:
                self.directly_attached[level_id].add(protein_id)
            else:
                # 1. find closest level (between species leaf and current level)
                species_level = species_id
                while self.eggNOG_dict[species_level] != level_id:
                    species_level = self.eggNOG_dict[species_level]
                
                # 2. find nog or assign singleton
                lower_level = species_level
                lower_nog = self.get_protein_nog(protein_id,lower_level)
                
                if lower_nog:
                    self.higher_nogs[nog_id].add(lower_nog)
                    self.lower_nogs[lower_nog].add(nog_id)
                    
                    if lower_nog not in self.expanded_nog:
                        self.get_nog_proteins(lower_nog)
                        to_expand.add(lower_nog)
                else:
                    lower_nog = eu.get_nog_id(lower_level,protein_id, nog_type=2) # [TODO] change to singleton_nog (starts with 2)
                    self.higher_nogs[nog_id].add(lower_nog)
                    self.lower_nogs[lower_nog].add(nog_id)
                    
                    if protein_id not in self.expanded_nog[lower_nog]:
                        self.expanded_nog[lower_nog].add(protein_id)
                        self.singletons[lower_level].add(protein_id)
                        to_expand.add(lower_nog)
        
        # 3. double check by rebuilding all assumed origins
        new_proteins = self.expanded_nog[nog_id] & self.directly_attached[level_id] # [TODO] review utility of global directly_attached
        old_proteins = set()
        for lower_nog in self.higher_nogs[nog_id]:
            old_proteins.update(self.expanded_nog[nog_id] & self.expanded_nog[lower_nog])
        old_singletons = self.expanded_nog[nog_id]
        joint_set = new_proteins | old_proteins | old_singletons #union
        
        #test symmetric difference
        assert len(self.expanded_nog[nog_id] ^ joint_set) == 0, "Origins not matched: %d vs %d[%d,%d,%d]"%(
            len(self.expanded_nog[nog_id]),len(joint_set),len(new_proteins),len(old_proteins),len(old_singletons))
        
        return to_expand
        
    def merge(self, nog_id, new_nog_id, apply_changes = True):
        
        # Join higher nogs (n) into a single (1) nog
        assert nog_id in self.lower_nogs, "%d does not exist in exNOG.%d's lower nogs"%(nog_id,self.input_id)
            
        return self.merge_all(self.lower_nogs[nog_id], new_nog_id)

    def merge_all(self, higher_nogs, new_nog_id, apply_changes = True):
        
        merged_proteins = set()
        
        to_eliminate = []
        for higher_nog in higher_nogs:
            higher_proteins = self.expanded_nog[higher_nog]
            merged_proteins.update(higher_proteins)
            to_eliminate.append(higher_nog)
        
        if apply_changes:
            
            self.expanded_nog[new_nog_id] = merged_proteins
            self.higher_nogs[new_nog_id] = set()
            
            # remove old references
            for higher_nog in to_eliminate:
                
                self.nogs_history[higher_nog] = {new_nog_id}
                
                # [TODO] either introduce singleton removal step or delete self.singleton use
                del self.expanded_nog[higher_nog]
                
                # top to bottom [higher_nogs]
                self.higher_nogs[new_nog_id].update(self.higher_nogs[higher_nog])
                del self.higher_nogs[higher_nog]
                
                # bottom to top [lower_nogs]        
                for lower_nog in self.higher_nogs[new_nog_id]:
                    if higher_nog in self.lower_nogs[lower_nog]:
                        self.lower_nogs[lower_nog].remove(higher_nog)
                        self.lower_nogs[lower_nog].add(new_nog_id) # set guarantees single occurrence
                
        return { new_nog_id:merged_proteins }
    
    def partial_split(self, higher_id, lower_id, new_nog_id, apply_changes = True ):
        """ partially split higher nog based on lower nog to simplify the problem """
        
        assert higher_id in self.higher_nogs
        assert lower_id in self.lower_nogs
        
        overlap = self.expanded_nog[higher_id] & self.expanded_nog[lower_id]
        assert overlap
        
        remaining = self.expanded_nog[higher_id] - overlap
        assert remaining
        
        splitted_proteins = defaultdict(set)
        minus_one = self.higher_nogs[higher_id]
        
        if apply_changes:
            self.nogs_history[higher_id] = set()
        
        for split_proteins in [overlap,remaining]:
            
            split_id = new_nog_id
            if len(split_proteins) == 1:
                # use nog ids starting with 4 for singletons w/o increasing split_id counter
                split_id = eu.get_copy_id(higher_id,list(split_proteins)[0],nog_type=4)
                splitted_proteins[split_id] = split_proteins
            else:
                splitted_proteins[split_id] = split_proteins
                new_nog_id += 1
            
            if apply_changes:
                    
                self.nogs_history[higher_id].add(split_id)
                
                # update references only if user specified a new nog_id
                self.higher_nogs[split_id] = set()
                self.expanded_nog[split_id] = split_proteins
                
                for lower_nog in minus_one:
                    if self.expanded_nog[lower_nog] & split_proteins:
                        self.lower_nogs[lower_nog].add(split_id)
                        self.higher_nogs[split_id].add(lower_nog)
                
        if apply_changes:  
            # delete references
            del self.expanded_nog[higher_id]
            
            # delete L -> L-1
            del self.higher_nogs[higher_id]
            
            # delete nog_id from L-1 -> L
            for lower_nog in minus_one:
                self.lower_nogs[lower_nog].remove(higher_id)
        
        return splitted_proteins
        
    def split(self, nog_id, new_nog_id, apply_changes = True ):
        # Divide lower nog (1) according to higher nogs (n)
        assert nog_id in self.lower_nogs
        
        plus_one = self.lower_nogs[nog_id]
        splitted_proteins = defaultdict(set) # list of sets
        
        lower_proteins = self.expanded_nog[nog_id]
    
        # L - 1 (if complete expansion was applied)        
        minus_one = self.higher_nogs[nog_id]
        
        if apply_changes:
            self.nogs_history[nog_id] = set()
        
        for higher_nog in plus_one:
            higher_proteins = self.expanded_nog[higher_nog]
            split_proteins = higher_proteins & lower_proteins
            
            split_id = new_nog_id
            if len(split_proteins) == 1:
                # use nog ids starting with 4 for singletons w/o increasing split_id counter
                split_id = eu.get_nog_id(self.get_level_id(nog_id),list(split_proteins)[0], nog_type=4)
                splitted_proteins[split_id] = split_proteins
            else:
                splitted_proteins[split_id] = split_proteins
                new_nog_id += 1 # [TODO] rewrite split condition by decreasing new_nog_id by one befor the for-loop?
                
            if apply_changes:
                
                self.nogs_history[nog_id].add(split_id)
                
                # update references only if user specified a new nog_id
                self.lower_nogs[split_id] = set([higher_nog])
                self.higher_nogs[higher_nog].add(split_id)
                self.expanded_nog[split_id] = split_proteins
                
                # test for L-1, if overlap exists, insert split_id (old nog_id will be removed at the end)
                for lower_nog in minus_one:
                    if self.expanded_nog[lower_nog] & split_proteins:
                        self.lower_nogs[lower_nog].add(split_id)
                        self.higher_nogs[split_id].add(lower_nog)
                
        if apply_changes:  
            # delete references
            del self.expanded_nog[nog_id]
            
            # delete L -> L+1
            del self.lower_nogs[nog_id]
            
            # delete L -> L-1
            del self.higher_nogs[nog_id]
            
            # delete nog_id from L+1 -> L
            for higher_nog in plus_one:
                self.higher_nogs[higher_nog].remove(nog_id)
            
            # delete nog_id from L-1 -> L
            for lower_nog in minus_one:
                self.lower_nogs[lower_nog].remove(nog_id)
        
        return splitted_proteins
    
    def write_ids(self, output_file): # [TODO] consider to kick it out, redundant and not complete (see split_char)
        
        with open(output_file,'w') as f:
            for protein_id in self.get_full_protein_set():
                
                if protein_id not in self.protein_species: # [TODO] shouldn't this be an assert?
                    protein_species,protein_short = fu.shelve_lookup(protein_id, self.protein_id_shelve)
                    self.protein_species[protein_id] = protein_species
                
                species_id = self.protein_species[protein_id]
                
                split_char = 'A'
                species_short = 'test'
                
                simplified_id = "%s%s_%s%s"%(species_short,species_id,split_char,protein_id)
                
                f.write('%s\t%s\n'%(simplified_id,'fasta_name'))
                
        #print("Wrote ids to: %s"%output_file)
    
    def write_dependency(self, output_file, split_minimum=0, transfer_minimum=0, write_mode='w'): # [TODO] add default output %d.graph
        
        with open(output_file, write_mode) as f:
            
            #header
            if write_mode == 'w': # or os.path.getsize(output_file) == 0
                f.write('#source\ttarget\ttransfer\n')
            
            for l_nog in self.lower_nogs:
                
                transfers = {}
                
                for h_nog in self.lower_nogs[l_nog]:
                    transfer_size = len(self.expanded_nog[l_nog] & self.expanded_nog[h_nog])
                    if transfer_size >= len(self.expanded_nog[l_nog])*transfer_minimum:
                        transfers[h_nog] = transfer_size
            
                if len(transfers) >= split_minimum:
                    for h_nog, transfer_size in transfers.items():      
                        f.write('%s\t%s\t%d\n'%(l_nog,h_nog,transfer_size))
            
        #print("Wrote dependency graph to: %s"%output_file)
        
    def write_nogs(self, output_file, write_mode='w', filtered_nogs = None, filtered_species=None): # [TODO] add default output %d.nogs
        
        if filtered_species is not None:
            assert filtered_species in self.eggNOG_names, "Species to filter for is not known: %s"%filtered_species
        
        with open(output_file,write_mode) as f:
            if write_mode == 'w': # or os.path.getsize(output_file) == 0:
                f.write('#id\tlevel\tsize\n')
            for nog_id in self.expanded_nog:
                level_id = self.get_level_id(nog_id)
                if filtered_nogs:
                    if nog_id in filtered_nogs[level_id]:
                        nog_size = len(filtered_nogs[level_id][nog_id])
                    else:
                        assert len(self.expanded_nog[nog_id]) == 1
                        nog_size = 1
                else:
                    if filtered_species is None:
                        nog_size = len(self.expanded_nog[nog_id])
                    else:
                        nog_size = len([x for x in self.expanded_nog[nog_id] if self.protein_name_mapping[x][0] == filtered_species])
                        
                f.write('%s\t%s\t%d\n'%(nog_id,level_id,nog_size))
                
        #print("Wrote nog's level and size to: %s"%output_file)
    
    #Simple solution, protein tracking would be better    
    def get_level_species(self, nog_species, test_lca=False):
        lca_level_finder = eu.eggNOG_lca_finder()
        lca_level = lca_level_finder.find_level(list(nog_species.keys()))
        
        if test_lca:
            # testing the LCA implies ensuring that the species
            # composition's root should lie in the higher NOG level,
            # exception: root is a leaf level (=> use clade definition)
            if lca_level != self.higher_level and lca_level in self.eggNOG_dict_rev:
                # assumption: the lca level is a grandchild of the higher level
                parent_level = self.eggNOG_dict[lca_level]
                while parent_level != self.higher_level:
                    parent_level = self.eggNOG_dict[parent_level]
                    if parent_level == 1:
                        break
                assert parent_level == self.higher_level, "lca %d is not grandchild of %d"%(
                    parent_level,higher_level)
                return {} 
        
        #level_species[level_id] = {species_ids}
        level_species = defaultdict(set)
        
        nog_species_set = set(nog_species)
        
        # Check if lca_level is a leaf level. If yes use clade definitions
        if lca_level not in self.eggNOG_dict_rev:
            level_clades = load_v4clades(self.eggNOG_names[lca_level],species_type=int)
            for clade_id, clade_species in level_clades.items():
                if clade_species & nog_species_set:
                    level_species[clade_id] = clade_species & nog_species_set
        else:
            for species_id in nog_species:
                intermediates = lca_level_finder.get_intermediateLevels(
                    lca_level,species_id)
                
                if intermediates:
                    for intermediate in intermediates:
                        #Filter only the level immediately under the LCA_level
                        if self.eggNOG_dict[intermediate] == lca_level:
                            level_species[intermediate].add(species_id)
                            #level_species[lca_level].add(intermediate)
                else:
                    #meaning no intermediates are available, species is directly attached to H_level
                    level_species[lca_level].add(species_id)
            
        return level_species
        
    def find_inconsistencies(self):
        
        inconsistent_nogs = set()
        
        for L_nog in sorted(self.lower_nogs.keys()):
            H_nogs = self.lower_nogs[L_nog]
            if len(H_nogs) > 1:
                inconsistent_nogs.add(L_nog)
        
        return inconsistent_nogs
    
    def has_inconsistencies(self):
        
        return len(self.find_inconsistencies()) > 0
    
    def get_composition(self,nog_id,test_lca=True):
        
        # assert that input is compatible with object
        assert nog_id in self.lower_nogs, "nog_id %d is not present in exNOG object"%nog_id
        assert nog_id in self.expanded_nog, "nog_id %d is present in exNOG but lacks a definition"%nog_id
        
        # prepare proteins and create alphabetical enumeration for splits
        nog_proteins = set()
        nog_proteins.update(self.expanded_nog[nog_id]) # exluded for top_nogs selection
        nog_char = {}
        nog_labels = {}
        
        # extract proteins from upper nogs and assign labels for id
        for h_nog in sorted(self.lower_nogs[nog_id]):
            nog_proteins.update(self.expanded_nog[h_nog])
            
            if str(h_nog)[0] in ('2','4'): # [TODO] Define better starting labels for singletons
                nog_label = 'xxx'
            else:
                # triplet label system, up to 26**3 splits, reversed for better overview
                nog_label = ''
                split_size = len(nog_labels)
                for i in reversed(range(3)):
                    i_size = split_size//(26 ** i)
                    nog_label += chr(ord('A') + i_size)
                    split_size-= i_size * (26 ** i)
                
                assert split_size == 0, 'Split labelling failed for %d: %d splits with outcome %d'%(nog_id,len(nog_labels),split_size)
                nog_label = nog_label[::-1] 
                nog_labels[nog_label] = h_nog
                
            level_id = int(str(h_nog)[1:7])
                
            for protein_id in self.expanded_nog[h_nog]:
                nog_char[protein_id] = nog_label
        
        # obtain species_ids and protein_shortname (fasta retrieval)
        protein_shortnames = {}                 # dict[protein_id]=[protein_shortname]
        protein_species = {}                    # dict[protein_id]=[species_id] 
        final_names = {}                        # dict[protein_id]=[final_protein_name]
        proteins = {}
        
        paralogs = defaultdict(set)         # dict[species_id]=set([final_protein_name])
        split_species = defaultdict(set)
        
        for protein_id in nog_proteins:

            if self.use_permanence:
                species_id,protein_short = fu.shelve_lookup(protein_id, self.protein_id_shelve)
            else:
                species_id,protein_short = self.protein_name_mapping[protein_id]
                
            protein_shortnames[protein_id] = protein_short
            protein_species[protein_id] = species_id

            species_short = self.eggNOG_short_names[str(species_id)]
            split_char = nog_char[protein_id]
            split_species[split_char].add(species_id)
            
            simplified_id = "%s%s_%s%s"%(species_short,species_id,
                                         split_char,protein_id)
            paralogs[species_id].add(simplified_id)
            final_names[protein_id] = simplified_id
            
            #final output
            proteins[simplified_id] =  (species_id, protein_short)
        
        # 1. prepare level_species dict[level_id]=set([species_id])
        # [TODO] if level_cache:
        level_species = self.get_level_species(paralogs, test_lca=test_lca)
        
        return proteins, level_species, split_species, paralogs

