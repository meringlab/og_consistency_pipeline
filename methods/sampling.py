import os
import sys
import random
import argparse
from collections import defaultdict
from collections import Counter

def random_sampling(headers,sample_size):
    return random.sample(headers,sample_size)

def splits_sampling(split_species, paralogs, sample_size, sample=set(), sample_species=set()):
    species_with_paralogs = {x for x,s in paralogs.items() if len(s) > 1}
    
    #shuffle splits to avoid biasing sampling by order
    shuffled_split_ids = list(split_species.keys())
    random.shuffle(shuffled_split_ids)

    for split_id in shuffled_split_ids:
        if not split_id == 'xxx':
            if not sample_species.intersection(split_species[split_id]):
                
                split_species_with_paralogs = split_species[split_id].intersection(species_with_paralogs)
                if split_species_with_paralogs:
                    species_id = random.sample(split_species_with_paralogs,1)[0]
                else:
                    species_id = random.sample(split_species[split_id],1)[0]
                
                sampled_proteins = paralogs_sampling(paralogs[species_id])
                sample.update(sampled_proteins)
            
                sample_species.add(species_id)
                
        #print("Updated sampling for split %s with %d proteins from %s"%(split_id,len(paralogs[species_id]),species_id))
            
    return sample

def paralogs_sampling(paralogs):

    paralog_sample = set()
    
    # divide paralogs (same species) into their splits
    paralog_splits = defaultdict(set)
    for protein_str in paralogs:
        assert len(protein_str.split("_")) == 2, "Problem with %s"%protein_str
        species_str, protein_id = protein_str.split("_")
        split_id = protein_id[0:3]
        paralog_splits[split_id].add(protein_str)
    
    # for each split, sample one protein
    for split_id,split_proteins in paralog_splits.items():
        paralog_sample.add(random.sample(split_proteins,1)[0])

    return paralog_sample

def level_sampling(level_species, paralogs, sample_size, sample=set(),sample_species=set()):

    species_with_paralogs = {x for x,s in paralogs.items() if len(s) > 1}
    
    for level_id in level_species:
        # don't sample from level if any of it's species is already present
        if not sample_species.intersection(level_species[level_id]):
        
            level_species_with_paralogs = level_species[level_id].intersection(species_with_paralogs)
            if level_species_with_paralogs:
                species_id = random.sample(level_species_with_paralogs,1)[0]
            else:
                species_id = random.sample(level_species[level_id],1)[0]
            
            sampled_proteins = paralogs_sampling(paralogs[species_id])
            sample.update(sampled_proteins)
            sample_species.add(species_id)
    
    attempts = sample_size        
    while len(sample) < sample_size and attempts:
        attempts -= 1
        for level_id in level_species:
            species_id = random.sample(level_species[level_id],1)[0]
            
            sampled_proteins = paralogs_sampling(paralogs[species_id])
            sample.update(sampled_proteins)
            sample_species.add(species_id)
        
    #print("Updated sampling for level %s with %d proteins from %s"%(level_id,len(paralogs[species_id]),species_id))
    return sample

def combined_sampling(split_species,paralogs,level_species,sample_size):
    sample = set()
    sample_species = set()
    
    # fix a number between 1/4 and 1/3 to be slightly above
    species_no = sample_size / 3
    
    # 1. choose randomly two levels
    for level_id in random.sample(list(level_species.keys()),2):
        # 2. choose up to 3 species from each level
        species = random.sample(level_species[level_id],min(len(level_species[level_id]),species_no))
        # 3. choose 1 protein from each species and level
        for species_id in species:
            protein_id = random.sample(paralogs[species_id],1)[0]
            sample.add(protein_id)
    
    # 1. choose randomly two splits
    for split_id in random.sample(list(split_species.keys()),2):
        # 2. choose up to 3 species from each split
        species = random.sample(split_species[split_id],min(len(split_species[split_id]),species_no))
        # 3. choose 1 protein from each species and split
        split_prefix = '_%s'%split_id
        for species_id in species:
            protein_id = random.sample([x for x in paralogs[species_id] if split_prefix in x],1)[0]
            sample.add(protein_id)
    
    attempts = sample_size
    while len(sample) < sample_size and attempts:
        species_id = random.sample(list(paralogs.keys()),1)[0]
        protein_id = random.sample(paralogs[species_id],1)[0]
        sample.add(protein_id)
        attempts -= 1
        
    return sample

def write_fasta_sample(i, sequences_to_extract, fasta_sequences, sample_path):
    with open(sample_path,'w') as sample_file:
        for protein_id in sequences_to_extract:
            fasta_seq = fasta_sequences[protein_id]
            assert fasta_seq.endswith('\n'), "Fasta sequence for protein %"
            sample_file.write(fasta_seq)

class InconsistencySampler:
    
    def __init__(self,sample_no,sample_size,sampling_method,fasta_sequences):
        
        self.sample_no = sample_no
        self.sample_size = sample_size # TODO add self... 
        self.sampling_method = sampling_method
        self.fasta_sequences = fasta_sequences
    
    def sample_inconsistency(self,
                             proteins,
                             level_species,
                             split_species,
                             paralogs,
                             output_dir=None):
        
        samples = set()
        headers = set(proteins.keys())
        fasta_cache = {}
        
        # get specific fasta sequences
        for protein_id in proteins:
            
            species_id,protein_short = proteins[protein_id]
            full_name = "%d.%s"%(species_id,protein_short)
            fasta_seq = self.fasta_sequences[species_id][full_name]
                
            lines = fasta_seq.split('\n')
            lines[0] = ">%s"%protein_id    
            fasta_seq = "\n".join(lines)
            fasta_cache[protein_id] = fasta_seq
        
        if self.sample_size > len(headers) - 2:
            # subtract 2 to make sampling possible,
            # i.e. sample 10 from at least 10 + 2
            samples.add(frozenset(headers))
            sample_no = 1
        else:
            sample_no = self.sample_no
        
        rejected_sampling_attempts = 0
        while len(samples) < sample_no:
            
            if rejected_sampling_attempts == 10:
                break
            
            #Get sequences by chosen sampling method
            if self.sampling_method == "random":
                sequences_to_extract = random_sampling(headers,self.sample_size)
            elif self.sampling_method == "splits":
                sequences_to_extract = splits_sampling(split_species,paralogs,self.sample_size)
            elif self.sampling_method == "paralogs":
                sequences_to_extract = paralogs_sampling(paralogs, self.sample_size)
            elif self.sampling_method == "levels":
                sequences_to_extract = level_sampling(level_species, paralogs, self.sample_size)
            elif self.sampling_method == "combined":
                sequences_to_extract = combined_sampling(split_species,paralogs,level_species,self.sample_size)
            else:
                print("Unknown method %s"%self.sampling_method)
                sys.exit()
            
            frozen_sample = frozenset(sequences_to_extract)
            if frozen_sample in samples:
                # test if sample was already found
                rejected_sampling_attempts += 1
                # sys.stderr.write('>>>>>>>> REJECTED [duplicate] %d: %s\n'%(rejected_sampling_attempts,sequences_to_extract))
            elif len({x.split('_')[0] for x in frozen_sample}) < 2:
                # test whether at least two species are present in the sample
                rejected_sampling_attempts += 1
                sys.stderr.write('>>>>>>>> REJECTED [1 species] %d: %s\n'%(rejected_sampling_attempts,sequences_to_extract))
            else:
                samples.add(frozen_sample)
        
        # retrieve actual fasta sequences
        if len(samples) == 1:
            samples = list(samples) * self.sample_no
        
        sample_sequences = [""] * len(samples)
        for i,sequences_to_extract in enumerate(samples):
            
            if output_dir:
                sample_path = os.path.join(output_dir,'%d.fa'%i)
                write_fasta_sample(i, sequences_to_extract, fasta_cache, sample_path)
                sample_sequences[i] = sample_path
            else:
                sequences_to_extract = list(sequences_to_extract)
                random.shuffle(sequences_to_extract)
                for protein_id in sequences_to_extract:
                    assert protein_id in fasta_cache, '%s not found in %s from %s'%(protein_id,list(fasta_cache.keys()),sequences_to_extract)
                    fasta_seq = fasta_cache[protein_id]
                    sample_sequences[i] += fasta_seq
        
        samples = {}    
        return sample_sequences