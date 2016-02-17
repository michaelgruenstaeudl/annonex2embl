#!/usr/bin/env python
'''
Unit Tests for Custom operations module for EMBL submission preparation
tool
'''

#####################
# IMPORT OPERATIONS #
#####################

import unittest
import CustomOps as COps

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'Submission Preparation Tool for Sequences of Phylogenetic'\
           'Datasets (SPTSPD)'
__version__ = "2016.02.17.1900"

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

class AnnoChecksTestCases(unittest.TestCase):
    ''' Tests for class `AnnoChecks` '''
    
    def test_AnnoChecks_example_1(self):
        ''' Test to evaluate example 1 of AnnoChecks.check

        This test evaluates the default behaviour of the function. '''
        
        extract = Seq("ATGGCCTAA", generic_dna)
        location = FeatureLocation(0, 8)
        self.assertTrue(COps.AnnoChecks(extract,
            location).for_unittest())    
    
    def test_AnnoChecks_example_2(self):
        ''' Test to evaluate example 2 of AnnoChecks.check

        This test evaluates the situation where an internal stop codon is 
        present in the input sequence. '''
        
        extract = Seq("ATGTAATAA", generic_dna)
        location = FeatureLocation(0, 8)
        self.assertTrue(COps.AnnoChecks(extract,
            location).for_unittest())

    def test_AnnoChecks_example_3(self):
        ''' Test to evaluate example 3 of AnnoChecks.check

        This test evaluates the situation where the input sequence does not 
        start with a Methionine. '''
        
        extract = Seq("AAGTAA", generic_dna)
        location = FeatureLocation(0, 5)
        self.assertFalse(COps.AnnoChecks(extract,
            location).for_unittest())


class DegapButMaintainAnnoTestCases(unittest.TestCase):
    ''' Tests for class `DegapButMaintainAnno` '''

    def test_DegapButMaintainAnno_example_1(self):
        ''' Test to evaluate example 1 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains an internal gap.
        '''
        seq = "ATG-C"
        charsets = {"gene_1":[0,1],"gene_2":[2,3,4]}
        out_ideal = ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})

        out_actual = COps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_2(self):
        ''' Test to evaluate example 2 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains start and end gaps.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = COps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)
       
    def test_DegapButMaintainAnno_example_3(self):
        ''' Test to evaluate example 3 of DegapButMaintainAnno.degap

        This test evaluates the case where an entire gene is missing.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})

        out_actual = COps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_4(self):
        ''' Test to evaluate example 4 of DegapButMaintainAnno.degap
        
        This test evaluates the case where genes with internal gaps are 
        overlapping.
        '''
        
        seq = "A--AT--T"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})

        out_actual = COps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_5(self):
        ''' Test to evaluate example 5 of DegapButMaintainAnno.degap

        This test evaluates the case where genes with start and end gaps are 
        overlapping.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = COps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_6(self):
        ''' Test to evaluate example 6 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains start and end gaps 
        and the charset order is incorrect.
        '''
        seq = "AT----GC"
        charsets = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
        out_ideal = ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = COps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)    


#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
