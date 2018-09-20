#!/usr/bin/env python3
#Utilities to work with files

#author:    Davide Heller
#email:     davide.heller@imls.uzh.ch
#date:      2015-08-19

import sys
import os
import subprocess
import shelve

def read_tsv_dict(file_name,
                  field_no = 2,
                  type1 = str, type2 = str,
                  REVERSE_ORDER = False,
                  second_idx = -1,
                  skip_header = False,
                  skip_comment = False,
                  skip_assert = False):
    
    """Versatile function to read a tab separated value file into a dictionary
    
    Keyword arguments:
    file_name -- path to the file to read
    field_no -- number tab separated fields in the file
    type1 -- variable type to be read for the first import field
    type2 -- variable type to be read for the second import field
    REVERSE_ORDER -- set True to reverse the default: dict[field1] = field2
    second_idx -- set other index if the second field is not the last, i.e. field_no - 1
    skip_header -- set True if first line should be skipped
    skip_commet -- set True if all the lines starting with # should be skipped
    """
    
    tsv_dict = {}
    
    with open(file_name) as f:
        for line in f:
            
            #skip header line if requested
            if(skip_header):
                skip_header = False
                continue
            
            if(skip_comment):
                if(line.startswith("#")):
                    continue
            
            l = line.rstrip().split('\t')
            
            if(skip_assert):
                second_idx = len(l) - 1
            else:
                assert len(l) == field_no, "%d text fields found: %s"%(len(l),str(l))
                
            if(second_idx == -1):
                second_idx = field_no-1
            else:
                assert(second_idx < field_no), "Field idx to extract is higher than field_no! (%d > %d)"%(second_idx,field_no)
            
            field1 = type1(l[0])
            field2 = type2(l[second_idx])
    
            if(REVERSE_ORDER):
                tsv_dict[field2] = field1
            else:
                tsv_dict[field1] = field2
        
    return tsv_dict

def save_tsv_dict(tsv_dict, output_file_name = "tmp.tsv", REVERSE_ORDER = False, header_str="none"):
    """Write a dictionary to a file.
    default: dict[key] = value <> str(key[\t]value[\n])
    
    Keyword arguments:
    tsv_dict -- dictionary to write
    output_file_name -- path where to to write the file
    REVERSE_ORDER -- Set True to invert the order between key and value
    """
    
    output_file = open(output_file_name, 'w')
    
    if(header_str != "none"):
        header_line = "#"+header_str+'\n'
        output_file.write(header_line)
    
    for key in tsv_dict:
        value = tsv_dict[key]
        
        if(REVERSE_ORDER):
            line = "%s\t%s\n"%(str(value),str(key))
        else:
            line = "%s\t%s\n"%(str(key),str(value))
            
        output_file.write(line)
    
    output_file.close()

def wccount(filename):
    """
    Function to count the lines in an input file by
    using the command line utility wc with option [-l]
    """
    assert(os.path.isfile(filename)),'Incorrect path:%s'%filename
    
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

class ProgressCounter:
    """
    Utility to visualize the progress on command line.
    The class is initialized with the total number of elements to
    be processed and with an optional starting message. The function
    count serves to increase the counter of processed elements and
    outputs a message to standard out whenever 10% of the elements
    were processed.
    """
    def __init__(self,tot_number,start_message='Counter initialized'):
        self.tot_number = float(tot_number)
        self.current_number = 0
        self.pct_number = max(1,(tot_number / 100) - 1);
        print("%s: %d entries to be analyzed"%(start_message,tot_number))
        
    def count(self):
        self.current_number += 1
        if(self.current_number % self.pct_number == 0):
            self.show_progress()
        
    def show_progress(self):
        sys.stdout.write("\r%.2f%%" %(self.current_number*100/self.tot_number))
        sys.stdout.flush()
        
    def get_count(self):
        return self.current_number
    
    def set_zero(self):
        self.current_number = 0
    
def shelve_lookup(key,shelve_file):
    """
    Returns entry for key in shelve_file
    """
    
    s = shelve.open(shelve_file)
    
    try:
        value = s[str(key)]
    except KeyError as e:
        value = None
    finally:
        s.close()
        
    return value