#!/usr/bin/env python
'''
Unit tests to compare actual and expected output
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

import subprocess

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2019 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2019.10.10.1730'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

try:
    import inspect
    base_path = os.path.split(inspect.getfile(annonex2embl))[0] + '/'
except:
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # Note: Use dirname twice to make two steps downward

script = 'scripts/annonex2embl_launcher_CLI.py'
script_path = os.path.join(base_path, script)

metad_fn = 'Metadata.csv'
descr = '"a description here"'
e_mail = 'your_email_here@yourmailserver.com'
authors = '"LastName1 A.; LastName2 B."'

###########
# CLASSES #
###########

class OutputTestCases(unittest.TestCase):
    ''' Tests to evaluate if the actual and the expected output during
        the generation of EMBL flatfile are identical. If they are not, 
        show their difference. '''
    def setUp(self):
        warnings.simplefilter('ignore')


    def test_1_actual_vs_expected_output(self):
        infile_nex = 'TestData1.nex'
        expected_outp = 'TestData1.embl'
        my_name = sys._getframe().f_code.co_name  # Get name of this function
        infile_nex = os.path.join(base_path, 'examples/input/', infile_nex)
        infile_csv = os.path.join(base_path, 'examples/input/', metad_fn)
        actual_outp = os.path.join(base_path, 'examples/temp/', my_name+'.embl')
        expected_outp = os.path.join(base_path, 'examples/output/', expected_outp)
        cmd_list = ['python', script_path, '-n', infile_nex, '-c', infile_csv,
                    '-d', descr, '-e', e_mail, '-a', authors, '-o', actual_outp]
        try:
            subprocess.check_output(' '.join(cmd_list), shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)
        with open(expected_outp, 'r') as expected_hndl:
            expected_str = expected_hndl.read()
        if os.path.isfile(actual_outp):
            with open(actual_outp, 'r') as actual_hndl:
                actual_str = actual_hndl.read()
            ## Important: Remove actual output so that lines from 
            ## subsequent tests are not appended, rendering actual and 
            ## expected different!
#            os.remove(actual_outp)
            # (Although keeping output can be helpful when generating 
            # new test files.)
        else:
            print('annonex2embl TESTING ERROR: actual_str not found.')
        self.assertTrue(isinstance(expected_str, str), 'Not a string: ' + expected_outp)
        self.assertTrue(isinstance(actual_str, str), 'Not a string: ' + actual_outp)
        self.assertMultiLineEqual(expected_str, actual_str)
    
    
    def test_2_actual_vs_expected_output(self):
        infile_nex = 'TestData2.nex'
        expected_outp = 'TestData2.embl'
        my_name = sys._getframe().f_code.co_name  # Get name of this function
        infile_nex = os.path.join(base_path, 'examples/input/', infile_nex)
        infile_csv = os.path.join(base_path, 'examples/input/', metad_fn)
        actual_outp = os.path.join(base_path, 'examples/temp/', my_name+'.embl')
        expected_outp = os.path.join(base_path, 'examples/output/', expected_outp)
        cmd_list = ['python', script_path, '-n', infile_nex, '-c', infile_csv,
                    '-d', descr, '-e', e_mail, '-a', authors, '-o', actual_outp]
        try:
            subprocess.check_output(' '.join(cmd_list), shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)
        with open(expected_outp, 'r') as expected_hndl:
            expected_str = expected_hndl.read()
        if os.path.isfile(actual_outp):
            with open(actual_outp, 'r') as actual_hndl:
                actual_str = actual_hndl.read()
#            os.remove(actual_outp)
        else:
            print('annonex2embl TESTING ERROR: actual_str not found.')
        self.assertTrue(isinstance(expected_str, str), 'Not a string: ' + expected_outp)
        self.assertTrue(isinstance(actual_str, str), 'Not a string: ' + actual_outp)
        self.assertMultiLineEqual(expected_str, actual_str)
    

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
