#!/usr/bin/env python
'''
Unit Tests for the classes of the module `ParsingOps`
'''

#####################
# IMPORT OPERATIONS #
#####################

import unittest
import warnings
import sys, os
# Add specific directory to sys.path in order to import its modules
# Note: This relative importing is amateurish; why can I not replace it with 'import annonex2embl'?
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'annonex2embl'))

import GlobalVariables as GlobVars
import ParsingOps as PrOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2019 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2019.10.11.1900'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

###########
# CLASSES #
###########


class ParseCharsetNameTestCases(unittest.TestCase):
    ''' Tests to evaluate class `ParseCharsetName` '''
    def setUp(self):
        warnings.simplefilter('ignore')


    def test_ParseCharsetName__parse__1(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the basic parsing ability of the
            class. '''
        charset_name = 'psbI_CDS'
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        handle = PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()
        self.assertIsInstance(handle, tuple)
        self.assertIsInstance(handle[0], str)
        self.assertIsInstance(handle[1], str)
        self.assertTrue(handle[1] in GlobVars.nex2ena_valid_INSDC_featurekeys)


    def test_ParseCharsetName__parse__2(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name
            contains more than one charset_type. '''
        charset_name = 'psbI_rRNA_CDS' # Two feature keys present.
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        with self.assertRaises(Exception) as e: # Note the "as e" is important, as you otherwise receive exception messages during unittests
            PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()


    def test_ParseCharsetName__parse__3(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name contains
            more than one charset_sym. '''
        charset_name = 'psbI_matK_CDS' # Two gene symbols present.
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        with self.assertRaises(Exception) as e:
            PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()


    def test_ParseCharsetName__parse__4(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name contains an
            unknown charset_sym. '''
        charset_name = 'xxxX_CDS'
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        with self.assertRaises(Exception) as e:
            PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()


    def test_ParseCharsetName__parse__5(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name
            contains only the charset_sym, not the charset_type. '''
        charset_name = 'matK'
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        with self.assertRaises(Exception) as e:
            PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()


    def test_ParseCharsetName__parse__6(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where a charset_name
            contains only the charset_type, not the charset_sym. '''
        charset_name = 'CDS'
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        with self.assertRaises(Exception) as e:
            PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()


    def test_ParseCharsetName__parse__7(self):
        ''' This test evaluates the function `parse` of the class
            `ParseCharsetName`.
            This test evaluates the situation where the position of the
            charset_sym and the charset_type have been inadvertedly
            switched. Such a switch should not matter. '''
        charset_name = 'CDS_psbI'
        email_addr = 'your_email_here@yourmailserver.com'
        product_lookup = True
        handle = PrOps.ParseCharsetName(charset_name, email_addr, product_lookup).parse()
        self.assertIsInstance(handle, tuple)
        self.assertIsInstance(handle[0], str)
        self.assertIsInstance(handle[1], str)
        self.assertTrue(handle[1] in GlobVars.nex2ena_valid_INSDC_featurekeys)

class GetEntrezInfoTestCases(unittest.TestCase):
    ''' Tests to evaluate class `GetEntrezInfo` '''
    def setUp(self):
        warnings.simplefilter('ignore')


    def test_GetEntrezInfo__taxname_lookup_ena(self):
        ''' This test evaluates function `_taxname_lookup_ena` of class `GetEntrezInfo`.
            This test evaluates if the ENA taxonomy service responds correctly. '''
        taxon_name = 'Arabidopsis thaliana'
        email_addr = 'your_email_here@yourmailserver.com'
        handle = PrOps.GetEntrezInfo(email_addr)._taxname_lookup_ena(taxon_name)
        self.assertTrue(handle=="1")


    def test_GetEntrezInfo__does_taxon_exist__1(self):
        ''' This test evaluates function `does_taxon_exist` of class `GetEntrezInfo`.
            This test evaluates the case where a taxon name is used that has not been 
            registered via the ENA taxonomy service. '''
        taxon_name = "foo bar"
        email_addr = 'your_email_here@yourmailserver.com'
        handle = PrOps.GetEntrezInfo(email_addr).does_taxon_exist(taxon_name)
        self.assertFalse(handle)


    def test_GetEntrezInfo__does_taxon_exist__2(self):
        ''' This test evaluates function `does_taxon_exist` of class `GetEntrezInfo`.
            This test evaluates the case where a taxon name is used that has not been 
            registered via the ENA taxonomy service. '''
        taxon_name = 'Arabidopsis thaliana'
        email_addr = 'your_email_here@yourmailserver.com'
        handle = PrOps.GetEntrezInfo(email_addr).does_taxon_exist(taxon_name)
        self.assertTrue(handle)


#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
