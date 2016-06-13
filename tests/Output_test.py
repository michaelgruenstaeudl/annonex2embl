#!/usr/bin/env python
'''
Unit tests to compare actual and expected output
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

import MyExceptions as ME
import subprocess

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.06.10.1600'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

try:
    import inspect
    base_path = os.path.split(inspect.getfile(annonex2embl))[0] + '/'
except:
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

script_rel_path = 'scripts/annonex2embl.py'
script_abs_path = os.path.join(base_path, script_rel_path)

e_mail = 'mi.gruenstaeudl@gmail.com'

###########
# CLASSES #
###########


class OutputTestCases(unittest.TestCase):
    ''' Tests to evaluate miscellaneous operations'''
    
    def test_1_actual_vs_expected_output(self):
        ''' Assert that the actual and the expected output are 
        identical. If they are not, show their difference. '''
        
        infile_nex = 'Test1_ALIGNM.nex'
        infile_csv = 'Test1_META.csv'
        expected_outp = 'Test1_ALIGNM.embl'
        
        ## Name of this function
        my_name = sys._getframe().f_code.co_name

        infile_nex_rel_path = os.path.join('tests/data/input/', infile_nex)
        infile_csv_rel_path = os.path.join('tests/data/input/', infile_csv)
        actual_outp_rel_path = os.path.join('tests/data/temp/', my_name)
        expected_outp_rel_path = os.path.join('tests/data/output/', expected_outp)

        infile_nex_abs_path = os.path.join(base_path, infile_nex_rel_path)
        infile_csv_abs_path = os.path.join(base_path, infile_csv_rel_path)
        actual_outp_abs_path = os.path.join(base_path, actual_outp_rel_path)
        expected_outp_abs_path = os.path.join(base_path, expected_outp_rel_path)

        cmd_list = ['python2', script_abs_path,
                    '-n', infile_nex_abs_path,
                    '-c', infile_csv_abs_path,
                    '-o', actual_outp_abs_path,
                    '-e', e_mail]

        try:
            subprocess.check_output(' '.join(cmd_list), shell=True)
        except subprocess.CalledProcessError as e:
            print e.output

        expected_str = open(expected_outp_abs_path).read()
        ## Check if actual output exists
        if os.path.isfile(actual_outp_abs_path):
            actual_str = open(actual_outp_abs_path).read()
            
            ## Important: Remove actual output so that lines from 
            ## subsequent tests are not appended, rendering actual and 
            ## expected different!
            os.remove(actual_outp_abs_path)
            # (But keeping output can be helpful when generating new test files.)

        else:
            print 'annonex2embl TESTING ERROR: actual_str not found.'

        self.assertTrue(isinstance(expected_str, str),
                'Not a string: ' + expected_outp_abs_path)
        self.assertTrue(isinstance(actual_str, str),
                'Not a string: ' + actual_outp_abs_path)
        self.assertMultiLineEqual(expected_str, actual_str)
    
    
    def test_2_actual_vs_expected_output(self):
        ''' Assert that the actual and the expected output are 
        identical. If they are not, show their difference. '''
        
        infile_nex = 'Amaranths_ITS_subset_ALIGN.nex'
        infile_csv = 'Amaranths_ITS_subset_META.csv'
        expected_outp = 'Amaranths_ITS_subset_ALIGN.embl'
        
        ## Name of this function
        my_name = sys._getframe().f_code.co_name

        infile_nex_rel_path = os.path.join('tests/data/input/', infile_nex)
        infile_csv_rel_path = os.path.join('tests/data/input/', infile_csv)
        actual_outp_rel_path = os.path.join('tests/data/temp/', my_name)
        expected_outp_rel_path = os.path.join('tests/data/output/', expected_outp)

        infile_nex_abs_path = os.path.join(base_path, infile_nex_rel_path)
        infile_csv_abs_path = os.path.join(base_path, infile_csv_rel_path)
        actual_outp_abs_path = os.path.join(base_path, actual_outp_rel_path)
        expected_outp_abs_path = os.path.join(base_path, expected_outp_rel_path)

        cmd_list = ['python2', script_abs_path,
                    '-n', infile_nex_abs_path,
                    '-c', infile_csv_abs_path,
                    '-o', actual_outp_abs_path,
                    '-e', e_mail]

        try:
            subprocess.check_output(' '.join(cmd_list), shell=True)
        except subprocess.CalledProcessError as e:
            print e.output

        expected_str = open(expected_outp_abs_path).read()
        ## Check if actual output exists
        if os.path.isfile(actual_outp_abs_path):
            actual_str = open(actual_outp_abs_path).read()
            ## Important: Remove actual output so that lines from 
            ## subsequent tests are not appended, rendering actual and 
            ## expected different!
            os.remove(actual_outp_abs_path)
        else:
            print 'annonex2embl TESTING ERROR: actual_str not found.'

        self.assertTrue(isinstance(expected_str, str),
                'Not a string: ' + expected_outp_abs_path)
        self.assertTrue(isinstance(actual_str, str),
                'Not a string: ' + actual_outp_abs_path)
        self.assertMultiLineEqual(expected_str, actual_str)

#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
