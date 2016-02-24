#!/usr/bin/env python
'''
Unit Tests for Custom operations module for EMBL submission preparation
tool
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio # Do not remove; important for assertIsInstance
import unittest
import DegappingOps as DgOps

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio import SeqFeature

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'Submission Preparation Tool for Sequences of Phylogenetic '\
           'Datasets (SPTSPD)'
__version__ = '2016.02.18.1100'

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

    def test_DegapButMaintainAnno_example_1(self):
        ''' Test to evaluate example 1 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains an internal gap.
        '''
        seq = "ATG-C"
        charsets = {"gene_1":[0,1],"gene_2":[2,3,4]}
        out_ideal = ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_2(self):
        ''' Test to evaluate example 2 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains start and end gaps.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)
       
    def test_DegapButMaintainAnno_example_3(self):
        ''' Test to evaluate example 3 of DegapButMaintainAnno.degap

        This test evaluates the case where an entire gene is missing.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_4(self):
        ''' Test to evaluate example 4 of DegapButMaintainAnno.degap
        
        This test evaluates the case where genes with internal gaps are 
        overlapping.
        '''
        
        seq = "A--AT--T"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_5(self):
        ''' Test to evaluate example 5 of DegapButMaintainAnno.degap

        This test evaluates the case where genes with start and end gaps are 
        overlapping.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_6(self):
        ''' Test to evaluate example 6 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains start and end gaps 
        and the charset order is incorrect.
        '''
        seq = "AT----GC"
        charsets = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
        out_ideal = ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal) 

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
