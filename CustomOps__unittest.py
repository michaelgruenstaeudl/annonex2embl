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
import CustomOps as CO
from CustomOps import MyException

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

class AnnoChecksTestCases(unittest.TestCase):
    ''' Tests for class `AnnoChecks` '''
    
    def test_AnnoChecks_example_1(self):
        ''' Test to evaluate example 1 of AnnoChecks.check

        This test evaluates the default behaviour of the function. '''
        
        extract = Seq("ATGGCCTAA", generic_dna)
        loc = FeatureLocation(0, 8)
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='cds')
        record_id = 'foobar'
        self.assertTrue(CO.AnnoChecks(extract, feature,
            record_id).for_unittest())    
    
    def test_AnnoChecks_example_2(self):
        ''' Test to evaluate example 2 of AnnoChecks.check

        This test evaluates the situation where an internal stop codon is 
        present in the input sequence. '''
        
        extract = Seq("ATGTAATAA", generic_dna)
        loc = FeatureLocation(0, 8)
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='cds')
        record_id = 'foobar'
        self.assertTrue(CO.AnnoChecks(extract, feature,
            record_id).for_unittest())

    def test_AnnoChecks_example_3(self):
        ''' Test to evaluate example 3 of AnnoChecks.check

        This test evaluates the situation where the input sequence does not 
        start with a Methionine. '''
        
        extract = Seq("AAGTAA", generic_dna)
        loc = FeatureLocation(0, 5)
        feature = SeqFeature.SeqFeature(loc, id='foobar', type='cds')
        record_id = 'foobar'
        with self.assertRaises(MyException):
            CO.AnnoChecks(extract, feature, record_id).for_unittest()


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
        self.assertTrue(CO.CheckCoord().quality_of_qualifiers(lst_of_dcts,
                        label))

    def test_CheckCoord_example_2(self):
        ''' Test to evaluate example 2 of CheckCoord().quality_of_qualifiers()

        This test evaluates the situation where the input label is NOT among
        the keys of the list of dictionaries.
        '''
        
        label = 'sequence_name'
        lst_of_dcts = [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                       {'isolate': 'taxon_B', 'country': 'Peru'}]
        with self.assertRaises(MyException):
            CO.CheckCoord().quality_of_qualifiers(lst_of_dcts, label)


class DegapButMaintainAnnoTestCases(unittest.TestCase):
    ''' Tests for class `DegapButMaintainAnno` '''

    def test_DegapButMaintainAnno_example_1(self):
        ''' Test to evaluate example 1 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains an internal gap.
        '''
        seq = "ATG-C"
        charsets = {"gene_1":[0,1],"gene_2":[2,3,4]}
        out_ideal = ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})

        out_actual = CO.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_2(self):
        ''' Test to evaluate example 2 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains start and end gaps.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = CO.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)
       
    def test_DegapButMaintainAnno_example_3(self):
        ''' Test to evaluate example 3 of DegapButMaintainAnno.degap

        This test evaluates the case where an entire gene is missing.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})

        out_actual = CO.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_4(self):
        ''' Test to evaluate example 4 of DegapButMaintainAnno.degap
        
        This test evaluates the case where genes with internal gaps are 
        overlapping.
        '''
        
        seq = "A--AT--T"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})

        out_actual = CO.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_5(self):
        ''' Test to evaluate example 5 of DegapButMaintainAnno.degap

        This test evaluates the case where genes with start and end gaps are 
        overlapping.
        '''
        seq = "AA----TT"
        charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
        out_ideal = ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = CO.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)

    def test_DegapButMaintainAnno_example_6(self):
        ''' Test to evaluate example 6 of DegapButMaintainAnno.degap

        This test evaluates the case where a gene contains start and end gaps 
        and the charset order is incorrect.
        '''
        seq = "AT----GC"
        charsets = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
        out_ideal = ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})

        out_actual = CO.DegapButMaintainAnno(seq, charsets).degap()
        self.assertTupleEqual(out_actual, out_ideal)    


class GenerateFeatLocTestCases(unittest.TestCase):
    ''' Tests for class `GenerateFeatLoc` '''

    def test_GenerateFeatLoc_example_1(self):
        ''' Test to evaluate example 1 of GenerateFeatLoc.exact

        This test evaluates an exact feature location.
        '''
        start_pos = 1
        stop_pos = 12

        out = CO.GenerateFeatLoc(start_pos, stop_pos).exact()
        self.assertIsInstance(out, Bio.SeqFeature.FeatureLocation)


class GenerateSeqFeatureTestCases(unittest.TestCase):
    ''' Tests for class `GenerateSeqFeature` '''

    def test_GenerateSeqFeature_source_feat_example_1(self):
        ''' Test to evaluate example 1 of GenerateSeqFeature().source_feat()

        This test evaluates the correct generation of the SeqFeature `source`.
        '''
        feat_len = 500
        quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
        transl_table = 11

        out = CO.GenerateSeqFeature().source_feat(feat_len, quals, 
                                                  transl_table)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)

    def test_GenerateSeqFeature_regular_feat_example_1(self):
        ''' Test to evaluate example 1 of GenerateSeqFeature().regular_feat()

        This test evaluates the correct generation of a regular SeqFeature.
        '''
        
        feat_name = 'psbI_CDS'
        feat_range = [2, 3, 4, 5]

        out = CO.GenerateSeqFeature().regular_feat(feat_name, feat_range)
        self.assertIsInstance(out, Bio.SeqFeature.SeqFeature)


class MetaChecksTestCases(unittest.TestCase):
    ''' Tests for class `MetaChecks` '''
    
    def test_MetaChecks_label_present_example_1(self):
        ''' Test to evaluate example 1 of MetaChecks.label_present

        This test evaluates the situation where the label is present in ALL 
        key lists.
        '''
        
        lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
            'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
            'baz': 'foobarbaz'}]
        label = 'foo'
        self.assertTrue(CO.MetaChecks(lst_of_dcts).label_present(label)) 
    
    def test_MetaChecks_label_present_example_2(self):
        ''' Test to evaluate example 2 of MetaChecks.label_present

        This test evaluates the situation where the label is not present in 
        EACH key list.
        '''
        
        lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
            'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
            'baz': 'foobarbaz'}]
        label = 'qux'
        with self.assertRaises(MyException):
            CO.MetaChecks(lst_of_dcts).label_present(label)
    
    def test_MetaChecks_label_present_example_3(self):
        ''' Test to evaluate example 3 of MetaChecks.label_present

        This test evaluates the situation where the label is not present in 
        ANY key list.
        '''
        
        lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
            'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
            'baz': 'foobarbaz'}]
        label = 'norf'
        with self.assertRaises(MyException):
            CO.MetaChecks(lst_of_dcts).label_present(label)

    def test_MetaChecks_valid_INSDC_quals_example_1(self):
        ''' Test to evaluate example 1 of MetaChecks.valid_INSDC_quals

        This test evaluates the situation where no invalid qualifiers are
        present.
        '''
        
        lst_of_dcts = [{'allele': 'foobar', 'altitude': 'foobar', 
            'anticodon': 'foobar'}, {'trans_splicing': 'foobar', 
            'type_material': 'foobar', 'variety': 'foobar'}]
        self.assertTrue(CO.MetaChecks(lst_of_dcts).valid_INSDC_quals())

    def test_MetaChecks_valid_INSDC_quals_example_2(self):
        ''' Test to evaluate example 2 of MetaChecks.valid_INSDC_quals

        This test evaluates the situation where invalid qualifiers are very 
        much present.
        '''
        
        lst_of_dcts = [{'allele': 'foobar', 'MyInvalidQual_1': 'foobar',
            'anticodon': 'foobar'}, {'MyInvalidQual_2': 'foobar',
            'type_material': 'foobar', 'variety': 'foobar'}]
        with self.assertRaises(MyException):
            CO.MetaChecks(lst_of_dcts).valid_INSDC_quals()
            

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
