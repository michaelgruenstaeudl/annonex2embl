#!/usr/bin/env python
'''
Unit Tests for the classes of the module `CheckingOps`
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

import MyExceptions as ME
import CheckingOps as CkOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.02.18.1100'

#############
# DEBUGGING #
#############

import pdb
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
        ''' Test to evaluate example 2 of AnnoChecks.check

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
        ''' Test to evaluate example 3 of AnnoChecks.check

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

    def test_AnnoChecks_example_4(self):
        ''' Test to evaluate example 4 of AnnoChecks.check

        This test evaluates the situation where the input sequence only codes 
        for a single aminoacid. '''
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        from Bio.SeqFeature import FeatureLocation
        from Bio import SeqFeature
        
        extract = Seq("ATGTAA", generic_dna)
        loc = FeatureLocation(0, 6) # Stop_pos must be +1
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='CDS')
        record_id = 'foobar'
        with self.assertRaises(ME.MyException):
            CkOps.AnnoCheck(extract, feature, record_id).for_unittest()


class CheckCoordTestCases(unittest.TestCase):
    ''' Tests for class `CheckCoord` '''
    
    def test_CheckCoord_example_1(self):
        ''' Test to evaluate example 1 of CheckCoord().quality_of_qualifiers()

        This test evaluates the situation where the input label is among the
        keys of the list of dictionaries.
        '''
        
        label = 'isolate'
        lst_of_dcts = [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                       {'isolate': 'taxon_B', 'country': 'Peru'}]
        self.assertTrue(CkOps.CheckCoord().quality_of_qualifiers(lst_of_dcts,
                        label))

    def test_CheckCoord_example_2(self):
        ''' Test to evaluate example 2 of CheckCoord().quality_of_qualifiers()

        This test evaluates the situation where the input label is NOT among
        the keys of the list of dictionaries.
        '''
        
        label = 'sequence_name'
        lst_of_dcts = [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                       {'isolate': 'taxon_B', 'country': 'Peru'}]
        with self.assertRaises(ME.MyException):
            CkOps.CheckCoord().quality_of_qualifiers(lst_of_dcts, label)


class MetaCheckTestCases(unittest.TestCase):
    ''' Tests for class `MetaCheck` '''
    
    def test_MetaCheck_label_present_example_1(self):
        ''' Test to evaluate example 1 of MetaCheck.label_present

        This test evaluates the situation where the label is present in ALL 
        key lists.
        '''
        
        lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
            'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
            'baz': 'foobarbaz'}]
        label = 'foo'
        self.assertTrue(CkOps.MetaCheck(lst_of_dcts).label_present(label)) 
    
    def test_MetaCheck_label_present_example_2(self):
        ''' Test to evaluate example 2 of MetaCheck.label_present

        This test evaluates the situation where the label is not present in 
        EACH key list.
        '''
        
        lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
            'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
            'baz': 'foobarbaz'}]
        label = 'qux'
        with self.assertRaises(ME.MyException):
            CkOps.MetaCheck(lst_of_dcts).label_present(label)
    
    def test_MetaCheck_label_present_example_3(self):
        ''' Test to evaluate example 3 of MetaCheck.label_present

        This test evaluates the situation where the label is not present in 
        ANY key list.
        '''
        
        lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
            'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
            'baz': 'foobarbaz'}]
        label = 'norf'
        with self.assertRaises(ME.MyException):
            CkOps.MetaCheck(lst_of_dcts).label_present(label)

    def test_MetaCheck_valid_INSDC_quals_example_1(self):
        ''' Test to evaluate example 1 of MetaCheck.valid_INSDC_quals

        This test evaluates the situation where no invalid qualifiers are
        present.
        '''
        
        lst_of_dcts = [{'allele': 'foobar', 'altitude': 'foobar', 
            'anticodon': 'foobar'}, {'trans_splicing': 'foobar', 
            'type_material': 'foobar', 'variety': 'foobar'}]
        self.assertTrue(CkOps.MetaCheck(lst_of_dcts).valid_INSDC_quals())

    def test_MetaCheck_valid_INSDC_quals_example_2(self):
        ''' Test to evaluate example 2 of MetaCheck.valid_INSDC_quals

        This test evaluates the situation where invalid qualifiers are very 
        much present.
        '''
        
        lst_of_dcts = [{'allele': 'foobar', 'MyInvalidQual_1': 'foobar',
            'anticodon': 'foobar'}, {'MyInvalidQual_2': 'foobar',
            'type_material': 'foobar', 'variety': 'foobar'}]
        with self.assertRaises(ME.MyException):
            CkOps.MetaCheck(lst_of_dcts).valid_INSDC_quals()

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
