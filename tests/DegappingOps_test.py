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
__version__ = '2016.02.24.2000'

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
        charsets = {"gene_1":[0,1],"gene_2":[2,3,4]}
        out_ideal = ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_2_DegapButMaintainAnno(self):
        ''' This test evaluates the case where a gene contains gaps at the 
        start and at the end.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)
       
    def test_3_DegapButMaintainAnno(self):
        ''' This test evaluates the case where an entire gene is missing.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_4_DegapButMaintainAnno(self):
        ''' This test evaluates the case where genes with internal gaps are 
        overlapping.
        '''        
        seq = "A--AT--T"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_5_DegapButMaintainAnno(self):
        ''' This test evaluates the case where genes with start and end gaps 
        are overlapping.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_6_DegapButMaintainAnno(self):
        ''' This test evaluates the case where two genes are entirely 
        overlapping.
        '''
        seq = "ATG-C"
        charsets = {"gene1":[0,1,2], "gene2":[0,1,2], "gene3":[2,3,4]}
        out_ideal = ('ATGC', {'gene1': [0, 1, 2], 'gene2': [0, 1, 2], 
            'gene3': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_7_DegapButMaintainAnno(self):
        ''' This test evaluates the case where a gene contains start and end 
        gaps and the charset order is incorrect.
        '''
        seq = "AT----GC"
        charsets = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
        out_ideal = ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = DgOps.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_8_DegapButMaintainAnno(self):
        ''' This test evaluates the case where three genes with internal gaps
        are overlapping.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3,4,5,6,7], "gene2":[0,1,2,3,4,5,6,7], 
            "gene3":[0,1,2,3,4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0,1,2,3], 'gene2': [0,1,2,3],
            'gene3': [0,1,2,3]})

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
