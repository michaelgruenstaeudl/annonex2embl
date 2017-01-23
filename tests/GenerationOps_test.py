#!/usr/bin/env python
'''
Unit Tests for the classes of the module `GenerationOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio # Do not remove; important for assertIsInstance
import unittest

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import GenerationOps as GnOps

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio import SeqFeature

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.06.13.1600'

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
        ''' Test to evaluate example 1 of GenerateFeatLoc.make_location

        This test evaluates the function 'make_location' using
        FeatureLocations.
        '''
        charset_range = range(1,8)

        out = GnOps.GenerateFeatLoc().make_location(charset_range)
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation) # FeatureLocation

    def test_GenerateFeatLoc_example_2(self):
        ''' Test to evaluate example 2 of GenerateFeatLoc.make_location

        This test evaluates the function 'make_location' using
        CompoundLocations.
        '''
        charset_range = [1,2,3,7,8]

        out = GnOps.GenerateFeatLoc().make_location(charset_range)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation

    def test_GenerateFeatLoc_example_2(self):
        ''' Test to evaluate if function `make_location` correctly
        generates a compound feature location if only a single base 
        separates them.
        '''
        charset_range = [1,2,3,5,6]

        out = GnOps.GenerateFeatLoc().make_location(charset_range)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation
        self.assertEqual(len(out.parts), 2)
    
    def test_GenerateFeatLoc_example_4(self):
        ''' Test to evaluate the function GenerateFeatLoc.make_start_fuzzy()

        This test evaluates if start FeatureLocations are made fuzzy.
        '''
        from Bio import SeqFeature
        start_pos = SeqFeature.ExactPosition(5)
        end_pos = SeqFeature.ExactPosition(9)
        location_object = SeqFeature.FeatureLocation(start_pos, end_pos)

        out = GnOps.GenerateFeatLoc().make_start_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation) # FeatureLocation
        self.assertIsInstance(out.start, Bio.SeqFeature.BeforePosition) # Fuzzy Start

    def test_GenerateFeatLoc_example_5(self):
        ''' Test to evaluate the function GenerateFeatLoc.make_start_fuzzy()

        This test evaluates if start FeatureLocations are made fuzzy.
        '''
        from Bio import SeqFeature
        charset_range = [1,2,3,7,8]
        location_object = GnOps.GenerateFeatLoc().make_location(charset_range)

        out = GnOps.GenerateFeatLoc().make_start_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation
        self.assertIsInstance(out.parts[0].start, Bio.SeqFeature.BeforePosition) # Fuzzy Start

    def test_GenerateFeatLoc_example_6(self):
        ''' Test to evaluate the function GenerateFeatLoc.make_end_fuzzy()

        This test evaluates if end FeatureLocations are made fuzzy.
        '''
        from Bio import SeqFeature
        start_pos = SeqFeature.ExactPosition(5)
        end_pos = SeqFeature.ExactPosition(9)
        location_object = SeqFeature.FeatureLocation(start_pos, end_pos)

        out = GnOps.GenerateFeatLoc().make_end_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation) # FeatureLocation
        self.assertIsInstance(out.end, Bio.SeqFeature.AfterPosition) # Fuzzy End

    def test_GenerateFeatLoc_example_7(self):
        ''' Test to evaluate the function GenerateFeatLoc.make_end_fuzzy()

        This test evaluates if end FeatureLocations are made fuzzy.
        '''
        from Bio import SeqFeature
        charset_range = [1,2,3,7,8]
        location_object = GnOps.GenerateFeatLoc().make_location(charset_range)

        out = GnOps.GenerateFeatLoc().make_start_fuzzy(location_object)
        self.assertIsInstance(out, Bio.SeqFeature.CompoundLocation) # CompoundLocation
        self.assertIsInstance(out.parts[0].start, Bio.SeqFeature.BeforePosition) # Fuzzy Start


class GenerateSeqFeatureTestCases(unittest.TestCase):
    ''' Tests for class `GenerateSeqFeature` '''

    def test_GenerateSeqFeature_source_feat_example_1(self):
        ''' Test to evaluate example 1 of GenerateSeqFeature().source_feat()

        This test evaluates the correct generation of the SeqFeature `source`.
        '''
        full_len = 509
        quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
        transl_table = 11

        out = GnOps.GenerateSeqFeature().source_feat(full_len,
            quals, transl_table)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

    def test_GenerateSeqFeature_regular_feat_example_1(self):
        ''' Test to evaluate example 1 of GenerateSeqFeature().regular_feat()

        This test evaluates the correct generation of a regular, non-coding 
        SeqFeature.
        '''
        feature_name = 'psbI'
        feature_type = 'intron'
        charset_range = [2, 3, 4, 5]
        feature_loc = GnOps.GenerateFeatLoc().make_location(charset_range)

        out = GnOps.GenerateSeqFeature().regular_feat(feature_name,
            feature_type, feature_loc)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
