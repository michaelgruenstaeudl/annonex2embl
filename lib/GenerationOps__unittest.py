#!/usr/bin/env python
'''
Unit Tests for the classes of the module `GenerationOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio # Do not remove; important for assertIsInstance
import unittest
import GenerationOps as GnOps

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

class GenerateFeatLocTestCases(unittest.TestCase):
    ''' Tests for class `GenerateFeatLoc` '''

    def test_GenerateFeatLoc_example_1(self):
        ''' Test to evaluate example 1 of GenerateFeatLoc.exact

        This test evaluates an exact feature location.
        '''
        start_pos = 1
        stop_pos = 12

        out = GnOps.GenerateFeatLoc(start_pos, stop_pos).exact()
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation)

class GenerateSeqFeatureTestCases(unittest.TestCase):
    ''' Tests for class `GenerateSeqFeature` '''

    def test_GenerateSeqFeature_source_feat_example_1(self):
        ''' Test to evaluate example 1 of GenerateSeqFeature().source_feat()

        This test evaluates the correct generation of the SeqFeature `source`.
        '''
        feature_len = 500
        quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
        transl_table = 11

        out = GnOps.GenerateSeqFeature().source_feat(feature_len, quals, 
                                                  transl_table)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

    def test_GenerateSeqFeature_regular_feat_example_1(self):
        ''' Test to evaluate example 1 of GenerateSeqFeature().regular_feat()

        This test evaluates the correct generation of a regular, non-coding 
        SeqFeature.
        '''
        
        feature_name = 'psbI'
        feature_type = 'intron'
        feature_range = [2, 3, 4, 5]

        out = GnOps.GenerateSeqFeature().regular_feat(feature_name, feature_type,
            feature_range)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
