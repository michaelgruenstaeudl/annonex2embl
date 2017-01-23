#!/usr/bin/env python
'''
Unit Tests for the classes of the module `DegappingOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import unittest

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import DegappingOps as DgOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.01.21.2200'

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

class DegapButMaintainAnnoTestCases(unittest.TestCase):
    ''' Tests for class `DegapButMaintainAnno` '''

    def test_1_DegapButMaintainAnno(self):
        ''' This test evaluates the case where a gene contains an internal gap.
        '''
        seq = "ATG-C"
        rmchar = "-"
        charsets = {"gene_1":[0,1],"gene_2":[2,3,4]}
        out_ideal = ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_2_DegapButMaintainAnno(self):
        ''' This test evaluates the case where a gene contains gaps at the 
        start and at the end.
        '''
        seq = "AA----TT"
        rmchar = "-"
        charsets = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)
       
    def test_3_DegapButMaintainAnno(self):
        ''' This test evaluates the case where an entire gene is missing.
        '''
        seq = "AA----TT"
        rmchar = "-"
        charsets = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_4_DegapButMaintainAnno(self):
        ''' This test evaluates the case where genes with internal gaps are 
        overlapping.
        '''        
        seq = "A--AT--T"
        rmchar = "-"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_5_DegapButMaintainAnno(self):
        ''' This test evaluates the case where genes with start and end gaps 
        are overlapping.
        '''
        seq = "AA----TT"
        rmchar = "-"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_6_DegapButMaintainAnno(self):
        ''' This test evaluates the case where two genes are entirely 
        overlapping.
        '''
        seq = "ATG-C"
        rmchar = "-"
        charsets = {"gene1":[0,1,2], "gene2":[0,1,2], "gene3":[2,3,4]}
        out_ideal = ('ATGC', {'gene1': [0, 1, 2], 'gene2': [0, 1, 2], 
            'gene3': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_7_DegapButMaintainAnno(self):
        ''' This test evaluates the case where a gene contains start and end 
        gaps and the charset order is incorrect.
        '''
        seq = "AT----GC"
        rmchar = "-"
        charsets = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
        out_ideal = ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_8_DegapButMaintainAnno(self):
        ''' This test evaluates the case where three genes with internal gaps
        are overlapping.
        '''
        seq = "AA----TT"
        rmchar = "-"
        charsets = {"gene1":[0,1,2,3,4,5,6,7], "gene2":[0,1,2,3,4,5,6,7], 
            "gene3":[0,1,2,3,4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0,1,2,3], 'gene2': [0,1,2,3],
            'gene3': [0,1,2,3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, rmchar, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)


class RmAmbigsButMaintainAnnoTestCases(unittest.TestCase):
    ''' Tests for class `RmAmbigsButMaintainAnno` '''
    
    def test_1_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where only leading ambiguities 
        (n=2) exist and the annotations do NOT overlap. The annotations 
        enclose the first and the last nucleotides.
        '''
        seq = "NNATGCCCC"
        rmchar = "N"
        charsets = {"gene1":[0,1,2,3],"gene2":[4,5,6,7,8]}
        out_ideal_step1 = ('ATGCCCC', {'gene1': [0,1], 'gene2': [2,3,4,5,6]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step1)
    
    def test_2_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where only leading ambiguities
        (n=2) exist and the annotations DO overlap. The annotations 
        enclose the first and the last nucleotides.
        '''
        seq = "NNATGCCCC"
        rmchar = "N"
        charsets = {"gene1":[0,1,2,3,4],"gene2":[3,4,5,6,7,8]}
        out_ideal_step1 = ('ATGCCCC', {'gene1': [0,1,2], 'gene2': [1,2,3,4,5,6]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step1)
    
    def test_3_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where only leading ambiguities
        (n=2) exist and the annotations DO overlap. The annotations 
        do NOT enclose the first and the last nucleotides.
        '''
        seq = "NNATGCCCC"
        rmchar = "N"
        charsets = {"gene1":[2,3,4],"gene2":[3,4,5,6]}
        out_ideal_step1 = ('ATGCCCC', {'gene1': [0,1,2], 'gene2': [1,2,3,4]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step1)
    
    def test_4_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where only trailing ambiguities
        (n=2) exist and the annotations do NOT overlap. The annotations 
        enclose the first and the last nucleotides.
        '''
        seq = "AATGCCCNN"
        rmchar = "N"
        charsets = {"gene1":[0,1,2],"gene2":[3,4,5,6,7,8]}
        out_ideal_step1 = (seq, charsets)
        out_ideal_step2 = ('AATGCCC', {'gene1': [0,1,2], 'gene2': [3,4,5,6]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step2)
    
    def test_5_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where only trailing ambiguities
        (n=2) exist and the annotations DO overlap. The annotations 
        enclose the first and the last nucleotides.
        '''
        seq = "AATGCCCNN"
        rmchar = "N"
        charsets = {"gene1":[0,1,2,3,4],"gene2":[3,4,5,6,7,8]}
        out_ideal_step1 = (seq, charsets)
        out_ideal_step2 = ('AATGCCC', {'gene1': [0,1,2,3,4], 'gene2': [3,4,5,6]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step2)
    
    def test_6_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where only trailing ambiguities
        (n=2) exist and the annotations DO overlap. The annotations 
        do NOT enclose the first and the last nucleotides.
        '''
        seq = "AATGCCCNN"
        rmchar = "N"
        charsets = {"gene1":[3,4],"gene2":[3,4,5,6,7,8]}
        out_ideal_step1 = (seq, charsets)
        out_ideal_step2 = ('AATGCCC', {'gene1': [3,4], 'gene2': [3,4,5,6]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step2)
    
    def test_7_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where leading (n=2) and 
        trailing (n=3) ambiguities exist and the annotations do NOT 
        overlap. The annotations enclose the first and the last
        nucleotides.
        '''
        seq = "NNATGCNNN"
        rmchar = "N"
        charsets = {"gene1":[0,1,2,3],"gene2":[4,5,6,7,8]}
        out_ideal_step1 = ('ATGCNNN', {'gene1': [0,1], 'gene2': [2,3,4,5,6]})
        out_ideal_step2 = ('ATGC', {'gene1': [0,1], 'gene2': [2,3]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step2)
    
    def test_8_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where leading (n=1) and 
        trailing (n=3) ambiguities exist and the annotations DO overlap.
        The annotations enclose the first and the last nucleotides.
        '''
        seq = "NATGCNNN"
        rmchar = "N"
        charsets = {"gene1":[0,1,2],"gene2":[2,3,4,5,6,7]}
        out_ideal_step1 = ('ATGCNNN', {'gene1': [0,1], 'gene2': [1,2,3,4,5,6]})
        out_ideal_step2 = ('ATGC', {'gene1': [0,1], 'gene2': [1,2,3]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step2)
    
    def test_9_RmAmbigsButMaintainAnno(self):
        ''' This test evaluates the case where leading (n=1) and 
        trailing (n=3) ambiguities exist and the annotations DO overlap.
        The annotations do NOT enclose the first and the last 
        nucleotides.
        '''
        seq = "NATGCNNN"
        rmchar = "N"
        charsets = {"gene1":[1,2],"gene2":[2,3,4,5]}
        out_ideal_step1 = ('ATGCNNN', {'gene1': [0,1], 'gene2': [1,2,3,4]})
        out_ideal_step2 = ('ATGC', {'gene1': [0,1], 'gene2': [1,2,3]})
        
        out_actual_1 = DgOps.RmAmbigsButMaintainAnno().rm_leadambig(seq, rmchar, charsets)
        self.assertTupleEqual(out_actual_1, out_ideal_step1)
        out_actual_2 = DgOps.RmAmbigsButMaintainAnno().rm_trailambig(out_actual_1[0], rmchar, out_actual_1[1])
        self.assertTupleEqual(out_actual_2, out_ideal_step2)


#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
