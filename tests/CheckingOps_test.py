#!/usr/bin/env python
'''
Unit Tests for the classes of the module `CheckingOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import Bio # Do not remove; important for assertIsInstance

import unittest
import sys, os

# Add specific directory to sys.path in order to import its modules
# Note: This relative importing is amateurish; why can I not replace it with 'import annonex2embl'?
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'annonex2embl'))

import CheckingOps as CkOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2020 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2020.03.06.1800'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

###########
# CLASSES #
###########


class AnnoChecksTestCases(unittest.TestCase):
    ''' Tests for class `AnnoChecks` '''
    def setUp(self):
        warnings.simplefilter('ignore')


    def test_AnnoChecks_example_1(self):
        ''' Test to evaluate function `check` of class `AnnoChecks`.
        This test evaluates the default behaviour of the function. '''
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        from Bio.SeqFeature import FeatureLocation
        from Bio import SeqFeature
        extract = Seq("ATGATATAA", generic_dna)
        loc = FeatureLocation(0, 9) # Stop_pos must be +1
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='CDS')
        record_id = 'foobar'
        self.assertTrue(CkOps.AnnoCheck(extract, feature,
            record_id).for_unittest())


    def test_AnnoChecks_example_2(self):
        ''' Test to evaluate function `check` of class `AnnoChecks`.
        This test evaluates the situation where an internal stop 
        codon exists in the input sequence. '''
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        from Bio.SeqFeature import FeatureLocation
        from Bio import SeqFeature
        extract = Seq("ATGATATAATAA", generic_dna)
        loc = FeatureLocation(0, 12) # Stop_pos must be +1
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='CDS')
        record_id = 'foobar'
        self.assertTrue(CkOps.AnnoCheck(extract, feature,
            record_id).for_unittest())


    def test_AnnoChecks_example_3(self):
        ''' Test to evaluate function `check` of class `AnnoChecks`.
        This test evaluates the situation where the input sequence
        does not start with a methionin. '''
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        from Bio.SeqFeature import FeatureLocation
        from Bio import SeqFeature
        extract = Seq("AAGATATAA", generic_dna)
        loc = FeatureLocation(0, 9) # Stop_pos must be +1
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='CDS')
        record_id = 'foobar'
        self.assertTrue(CkOps.AnnoCheck(extract, feature,
            record_id).for_unittest())


class QualifierCheckTestCases(unittest.TestCase):
    ''' Tests for class `QualifierCheck` '''
    def setUp(self):
        warnings.simplefilter('ignore')


    def test_QualifierCheck__label_present__1(self):
        ''' Test to evaluate the static method `_label_present` of class 
            `QualifierCheck`.
            This test evaluates the situation where the label is present in ALL 
            key lists. '''
        lst_of_dcts = [
            {'foo': 'foobarqux', 'bar': 'foobarqux', 'qux': 'foobarqux'},
            {'foo': 'foobarbaz', 'bar': 'foobarbaz', 'baz': 'foobarbaz'}]
        label = 'foo'
        self.assertTrue(CkOps.QualifierCheck._label_present(lst_of_dcts, label))
    

    def test_QualifierCheck__label_present__2(self):
        ''' Test to evaluate the static method `_label_present` of class 
            `QualifierCheck`.
            This test evaluates the situation where the label is not present in 
            EACH key list. '''
        lst_of_dcts = [
            {'foo': 'foobarqux', 'bar': 'foobarqux', 'qux': 'foobarqux'},
            {'foo': 'foobarbaz', 'bar': 'foobarbaz', 'baz': 'foobarbaz'}]
        label = 'qux'
        with self.assertRaises(Exception):
            CkOps.QualifierCheck._label_present(lst_of_dcts, label)
    

    def test_QualifierCheck__label_present__3(self):
        ''' Test to evaluate the static method `_label_present` of class 
            `QualifierCheck`.
            This test evaluates the situation where the label is not present in 
            ANY key list. '''
        lst_of_dcts = [
            {'foo': 'foobarqux', 'bar': 'foobarqux', 'qux': 'foobarqux'},
            {'foo': 'foobarbaz', 'bar': 'foobarbaz', 'baz': 'foobarbaz'}]
        label = 'norf'
        with self.assertRaises(Exception):
            CkOps.QualifierCheck._label_present(lst_of_dcts, label)
    

    def test_QualifierCheck__rm_empty_qual__1(self):
        ''' Test to evaluate the static method `_rm_empty_modifier` of class 
            `QualifierCheck`.
            This test evaluates the situation where empty modifiers are in
            qualifiers lists. '''
        lst_of_dcts = [
            {'foo': 'foo', 'bar': 'bar', 'qux': ''},
            {'foo': '', 'bar': 'bar', 'baz': 'baz'}]
        out_ideal = [
            {'foo': 'foo', 'bar': 'bar'},
            {'bar': 'bar', 'baz': 'baz'}]
        out_actual = CkOps.QualifierCheck._rm_empty_qual(lst_of_dcts)
        self.assertEqual(out_actual, out_ideal)
    

    def test_QualifierCheck__valid_INSDC_quals__1(self):
        ''' Test to evaluate example 1 of function `_valid_INSDC_quals` of 
            class `QualifierCheck`.
            This test evaluates the situation where only valid qualifiers are
            present. '''
        lst_of_dcts = [
            {'allele': 'foobar', 'altitude': 'foobar', 'anticodon': 'foobar'},
            {'trans_splicing': 'foobar', 'type_material': 'foobar', 'variety': 'foobar'}]
        self.assertTrue(CkOps.QualifierCheck._valid_INSDC_quals(lst_of_dcts))


    def test_QualifierCheck__valid_INSDC_quals__2(self):
        ''' Test to evaluate example 2 of function `_valid_INSDC_quals` of 
            class `QualifierCheck`.
            This test evaluates the situation where invalid qualifiers are very 
            much present. '''
        lst_of_dcts = [
            {'allele': 'foobar', 'MyInvalidQual_1': 'foobar', 'anticodon': 'foobar'},
            {'MyInvalidQual_2': 'foobar', 'type_material': 'foobar', 'variety': 'foobar'}]
        with self.assertRaises(Exception):
            CkOps.QualifierCheck._valid_INSDC_quals(lst_of_dcts)


    def test_QualifierCheck__quality_of_qualifiers__1(self):
        ''' Test to evaluate example 2 of function `quality_of_qualifiers` of 
            class `QualifierCheck`.
            This test evaluates the situation where the input label is among the
            keys of the list of dictionaries. '''
        lst_of_dcts = [
            {'isolate': 'taxon_A', 'country': 'Ecuador'},
            {'isolate': 'taxon_B', 'country': 'Peru'}]
        label = 'isolate'
        self.assertTrue(CkOps.QualifierCheck(lst_of_dcts, label).quality_of_qualifiers())


    def test_QualifierCheck__quality_of_qualifiers__2(self):
        ''' Test to evaluate example 2 of function `quality_of_qualifiers` of 
            class `QualifierCheck`.
        This test evaluates the situation where the input label is NOT among
        the keys of the list of dictionaries. '''
        lst_of_dcts = [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                       {'isolate': 'taxon_B', 'country': 'Peru'}]
        label = 'sequence_name'
        with self.assertRaises(Exception):
            CkOps.QualifierCheck(lst_of_dcts, label).quality_of_qualifiers()


#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
