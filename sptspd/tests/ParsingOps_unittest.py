#!/usr/bin/env python
'''
Unit Tests for the classes of the module `ParsingOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import unittest

import MyExceptions as ME
import ParsingOps as PrOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.02.24.1900'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

INSDC_feature_keys = ["assembly_gap", "C_region", "CDS", "centromere",
        "D-loop", "D_segment", "exon", "gap", "gene", "iDNA", "intron",
        "J_segment", "LTR", "mat_peptide", "misc_binding", "misc_difference",
        "misc_feature", "misc_recomb", "misc_RNA", "misc_structure",
        "mobile_element", "modified_base", "mRNA", "ncRNA", "N_region",
        "old_sequence", "operon", "oriT", "polyA_site", "precursor_RNA",
        "prim_transcript", "primer_bind", "protein_bind", "regulatory", 
        "repeat_region", "rep_origin", "rRNA", "S_region", "sig_peptide",
        "source", "stem_loop", "STS", "telomere", "tmRNA", "transit_peptide",
        "tRNA", "unsure", "V_region", "V_segment", "variation", "3'UTR",
        "5'UTR"]

###########
# CLASSES #
###########


class ParseCharsetNameTestCases(unittest.TestCase):
    ''' Tests to evaluate  class `ParseCharsetName` '''
    
    def test_1_ParseCharsetName(self):
        ''' This test evaluates the basic parsing ability of the class. '''
        
        charset_name = 'psbI_CDS'
        email_addr = 'mi.gruenstaeudl@gmail.com'
        handle = PrOps.ParseCharsetName(charset_name, email_addr).parse()      
        self.assertIsInstance(handle, tuple)
        self.assertIsInstance(handle[0], str)
        self.assertIsInstance(handle[1], str)
        self.assertTrue(handle[1] in INSDC_feature_keys)
    
    def test_2_ParseCharsetName(self):
        ''' This test evaluates the situation where a charset_name contains
        more than one charset_type. '''

        charset_name = 'psbI_rRNA_CDS' # Two feature keys present.
        email_addr = 'mi.gruenstaeudl@gmail.com'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()
    
    def test_3_ParseCharsetName(self):
        ''' This test evaluates the situation where a charset_name contains 
        more than one charset_sym. '''
        
        charset_name = 'psbI_matK_CDS' # Two gene symbols present.
        email_addr = 'mi.gruenstaeudl@gmail.com'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()
        
    def test_4_ParseCharsetName(self):
        ''' This test evaluates the situation where a charset_name contains an 
        unknown charset_sym. '''
        
        charset_name = 'xxxX_CDS'
        email_addr = 'mi.gruenstaeudl@gmail.com'
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()

    def test_5_ParseCharsetName(self):
        ''' This test evaluates the situation where a charset_name contains 
        only the charset_sym, not the charset_type. '''
        
        charset_name = 'matK'
        email_addr = 'mi.gruenstaeudl@gmail.com'        
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()
   
    def test_6_ParseCharsetName(self):
        ''' This test evaluates the situation where a charset_name contains
        only the charset_type, not the charset_sym. '''
        
        charset_name = 'CDS'
        email_addr = 'mi.gruenstaeudl@gmail.com'        
        with self.assertRaises(ME.MyException):
            PrOps.ParseCharsetName(charset_name, email_addr).parse()
   
    def test_7_ParseCharsetName(self):
        ''' This test evaluates the situation where the position of the 
        charset_sym and the charset_type have been inadvertedly switched. Such
        a switch should not matter. '''
        
        charset_name = 'CDS_psbI'
        email_addr = 'mi.gruenstaeudl@gmail.com'
        handle = PrOps.ParseCharsetName(charset_name, email_addr).parse()      
        self.assertIsInstance(handle, tuple)
        self.assertIsInstance(handle[0], str)
        self.assertIsInstance(handle[1], str)
        self.assertTrue(handle[1] in INSDC_feature_keys)

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
