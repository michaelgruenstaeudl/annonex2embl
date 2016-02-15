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
__version__ = '2016.02.15.2000'

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

class AnnoChecksTestCase(unittest.TestCase):
    ''' Tests for class "AnnoChecks" '''
    
    def test_AnnoChecks_example_1(self):
        ''' Test to evaluate example 1 of AnnoChecks.check

        This test evaluates the default behaviour of the function.
        '''
        
        extract = Seq("ATGGCCTAA", generic_dna)
        location = FeatureLocation(0, 8)
        self.assertTrue(COps.AnnoChecks(extract,
            location).for_unittest())
    
    
    def test_AnnoChecks_example_2(self):
        ''' Test to evaluate example 2 of AnnoChecks.check

        This test evaluates the situation where an internal stop codon is 
        present in the input sequence.
        '''
        
        extract = Seq("ATGTAATAA", generic_dna)
        location = FeatureLocation(0, 8)
        self.assertTrue(COps.AnnoChecks(extract,
            location).for_unittest())

    def test_AnnoChecks_example_3(self):
        ''' Test to evaluate example 3 of AnnoChecks.check

        This test evaluates the situation where the input sequence does not 
        start with a Methionine.
        '''
        
        extract = Seq("AAGTAA", generic_dna)
        location = FeatureLocation(0, 5)
        self.assertFalse(COps.AnnoChecks(extract,
            location).for_unittest())


#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
