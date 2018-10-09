# Rules to generate pickle files for orthologous group definitions

rule convert_nogs:
    input:
        level_path = join(config['input_dir'],'{level_id}.tsv'),
        protein_names = config['protein_names']
    output:
        converted_path = "preprocessed_data/converted_groups/{level_id}.tsv",
        names_path = "preprocessed_data/converted_groups/{level_id}.conversion.tsv"
    run:
        with open(input.protein_names) as f:
            protein_names = dict( (line.rstrip().split('\t')) for line in f)
            
        level_str = wildcards.level_id
        assert len(level_str) < 10, "level string exeeded 9 chrs for %s"%level_path
        level_int = int(level_str)
        level_nog_counter = 1
        with open(input.level_path) as f, open(output.converted_path, 'w') as g, open(output.names_path,'w') as n:
            for line in f:
                l = line.rstrip().split('\t')
                nog_old = l[0]
                
                # create new integer based identifier
                nog_int = '1%09d%09d'%(level_int,level_nog_counter)
                protein_ints = [protein_names[x] for x in l[1].split(',')]
                
                # write out new mapping and new conversion
                g.write('%s\t%s\n'%(nog_int,','.join(protein_ints)))
                n.write('%s\t%s\n'%(nog_old,nog_int))
                
                level_nog_counter += 1
    
rule pickle_nogs:
    input:
        og_wide_tsv="preprocessed_data/converted_groups/{level_id}.tsv",
    output:
        og_pickle="preprocessed_data/orthologous_groups/{level_id}.nogINT.setProteinINT.pkl2",
        og_long_tsv="preprocessed_data/orthologous_groups/{level_id}.tsv"
    run:
        import pickle
        
        nogs = {}
        with open(input.og_wide_tsv) as f:
            for line in f:
                nog_id, protein_ids = line.rstrip().split('\t')
                nog_id = int(nog_id)
                protein_ids = {int(x) for x in protein_ids.split(',')}
                nogs[nog_id] = protein_ids
        
        
        with open(output.og_pickle,'wb') as pfile:
            pickle.dump(nogs,pfile)
        
        with open(output.og_long_tsv, 'w') as f:
            for nog_id in nogs:
                for protein_id in nogs[nog_id]:
                    f.write('%d\t%d\n'%(nog_id,protein_id))
    
rule reconvert_nogs:
    input:
        long_tsv=join(config['output_dir'],'new_definition/{level_id}.tsv'),
        protein_names_pickle='preprocessed_data/proteinINT.tupleSpeciesINT_ShortnameSTR.pkl'
    output:
        wide_tsv=join(config['consistent_ogs'],'{level_id}.tsv')
    run:
        import os
        import tqdm
        import glob
        import pickle
            
        with open(input.protein_names_pickle,'rb') as f:
            protein_names = pickle.load(f)
        
        input_dir = os.path.dirname(input.long_tsv)
        output_dir = os.path.dirname(output.wide_tsv)
        
        input_levels = glob.glob(input_dir + '/*.tsv')
        
        print('Converting final definition (including %d sublevels..)'%(len(input_levels)-1))
        for long_tsv in tqdm.tqdm(input_levels):
            
            # read in long format (nog_id -> protein_id list)
            nogs = defaultdict(set)
            with open(long_tsv) as f:
                for line in f:
                    nog_id,protein_id = line.rstrip().split()
                    nogs[int(nog_id)].add(int(protein_id))
            
            # write out wide format (nog_str -> protein_str list)
            wide_tsv = os.path.join(output_dir,os.path.basename(long_tsv))
            with open(wide_tsv,'w') as g:
                for i, nog_id in enumerate(sorted(nogs)):
                    nog_str = 'NOG%06d'%(i+1)
                    proteins_str = ('%d.%s'%protein_names[x] for x in nogs[nog_id])
                    g.write('%s\t%s\n'%(nog_str,','.join(proteins_str)))
