#!/usr/bin/env python
'''
Custom operations for submission preparation tool
'''

#####################
# IMPORT OPERATIONS #
#####################

from copy import copy

###############
# AUTHOR INFO #
###############

__author__ = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
__copyright__ = "Copyright (C) 2016 Michael Gruenstaeudl"
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic Datasets (SPTSPD)"
__version__ = "2016.02.06.2000"

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

class DegapButMaintainAnno:
    ''' This class contains functions to degap DNA sequences while maintaining 
    annotations. 
    
    Specifically, the functions remove dashes from strings while 
    maintaining annotations on these strings. Only some of the implementations 
    work if the charsets are overlapping.
    
    Args:
        seq (str):      a string that represents an aligned DNA sequence;
                        example: "ATG-C"
        charsets (dict):a dictionary with gene names (str) as keys and lists 
                        of nucleotide positions (list) as values; example: 
                        {"gene_1":[0,1],"gene_2":[2,3,4]}
    
    Returns:
        tupl.   The return consists of the degapped sequence and the 
                corresponding degapped charsets; example: 
                (degapped_seq, degapped_charsets)
    
    Raises:
        currently nothing
    '''
    
    def __init__(self, seq, charsets):
        self.seq = seq
        self.charsets = charsets
        
    def intersection_exists(ranges_list):
        init_range = set(ranges_list[0])
        for r in ranges_list[1:]:
            if init_range.intersection(r):
                return True # Exits the entire function with 'True' (i.e., stops all of the loops)
        return False

    def degap_1(self):
        '''
        In its current implementation, this function only works if none of the 
        charsets are overlapping; hence, the initial check.
    
        Examples:
        
            Example 1: # Contains an internal gap
                >>> seq = "ATG-C"
                >>> annot = {"gene_1":[0,1],"gene_2":[2,3,4]}
                >>> DegapButMaintainAnno(seq, annot)
                Out: ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})
            
            Example 2: # Contains start and end gaps
                >>> seq = "AA----TT"
                >>> annot = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, annot)
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})
                    
            Example 3: # Entire genes missing
                >>> seq = "AA----TT"
                >>> annot = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
                >>> DegapButMaintainAnno(seq, annot)
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})
                    
            Example 4: # Overlapping genes with internal gaps
                >>> seq = "A--AT--T"
                >>> annot = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, annot)
                Out: ('AATTT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})
                    
            Example 5: # Overlapping genes with start and end gaps
                >>> seq = "AA----TT"
                >>> annot = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, annot)
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [1, 2]})
            
            Example 6: # Contains start and end gaps; incorrect charset order
                >>> seq = "AT----GC"
                >>> annot = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
                >>> DegapButMaintainAnno(seq, annot)
                Out: ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})
            
        TODO:
            (i)   Error in example 4: Out: AATT[sic!]T
            (ii)  Error in example 5: Out: {'gene1': [0, 1], 'gene2': [1[sic!], 2]})
        
        Notes:
            (i) In its current implementation, this function only works if 
            none of the charsets are overlapping; hence, the initial check.
            (ii) Order of charset in charsets seems to be irrelevant.
            
        '''
        seq = self.seq
        charsets = self.charsets
        
        if DegapButMaintainAnno.intersection_exists(charsets.values()):
            raise ValueError('MY ERROR: Character sets are overlapping.')
    
        degapped_seq = ''
        degapped_charsets = {}
        gaps_cumulative = 0
        for gene_name, index_list in charsets.items():
            gaps_within_gene = 0
            for pos, nucl in enumerate(seq):
                if pos in index_list and nucl == '-':
                    index_list.remove(pos)
                    gaps_within_gene += 1
                if pos in index_list and nucl != '-':
                    degapped_seq += nucl
                    index_list[index_list.index(pos)] = pos - gaps_within_gene
            index_list = [i-gaps_cumulative for i in index_list]
            degapped_charsets[gene_name] = index_list
            gaps_cumulative += gaps_within_gene
        return (degapped_seq, degapped_charsets)

    def degap_2(self):
        ''' This function works on overlapping charsets and is preferable over 
        "degap_1".
        Source: http://stackoverflow.com/questions/35233714/maintaining-overlapping-annotations-while-removing-dashes-from-string
        
            Examples:
                      
            Example 4: # Overlapping genes with internal gaps
                >>> seq = "A--AT--T"
                >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> rewriteGene(seq, charsets)
                Out: ('AATTT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})
                    
            Example 5: # Overlapping genes with start and end gaps
                >>> seq = "AA----TT"
                >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> rewriteGene(seq, charsets)
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [1, 2]})
           
        '''
        seq = self.seq
        charsets = self.charsets
        
        annotations = copy(charsets)
        index = seq.find('-')
        while index > -1:
            for gene_name, indices in annotations.items():
                if index in indices:
                    indices.remove(index)
                annotations[gene_name] = [e-1 if e > index else e for e in indices]
            seq = seq[:index] + seq[index+1:]
            index = seq.find('-')
        return seq, annotations

#############
# FUNCTIONS #
#############

