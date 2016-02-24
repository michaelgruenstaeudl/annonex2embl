#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'Submission Preparation Tool for Sequences of Phylogenetic '\
           'Datasets (SPTSPD)'
__version__ = '2016.02.24.1900'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

###########
# CLASSES #
###########

class Inp:
    ''' This class contains functions to conduct input and output operations.
    
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass


    def extract_fn(self, in_path):
        ''' This function splits a the path from path+filename. '''
        import os
        
        path, fn =  os.path.split(in_path)
        return fn

    
    def repl_fileend(self, fn, new_end):
        ''' This function replaces the file ending (e.g. ".csv") with a 
        different ending. '''
        
        return fn[:fn.rfind('.')] + '.' + new_end

    
    def parse_csv_file(self, path_to_csv):
        ''' This function parses a csv file. '''
        from csv import DictReader
        
        try:
            reader = DictReader(open(path_to_csv, 'rb'), delimiter=',', 
                quotechar='"', skipinitialspace=True)
            a_matrix = list(reader)
        except:
            raise ME.MyException('Parsing of .csv-file unsuccessful.')
        return a_matrix    


    def parse_nexus_file(self, path_to_nex):
        ''' This function parses a NEXUS file. '''
        from Bio.Nexus import Nexus
        
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_nex)
            charsets = aln.charsets
            matrix = aln.matrix
        except:
            raise ME.MyException('Parsing of .nex-file unsuccessful.')
        return (charsets, matrix)