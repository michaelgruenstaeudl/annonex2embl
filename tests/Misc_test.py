#!/usr/bin/env python
'''
Miscellaneous unit tests
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
__version__ = '2016.05.16.2000'

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
    
### CONTINUE HERE!

#    base_path = '/home/michael_science/git/michaelgruenstaeudl_annonex2embl/'
    base_path = 'test'
script_path = 'scripts/annonex2embl.py'

infile_nex_1_rel_path = 'tests/data/input/Amaranths_ITS_subset_ALIGN.nex'
infile_csv_1_rel_path = 'tests/data/input/Amaranths_ITS_subset_META.csv'
actual_outp_1_rel_path = 'tests/data/temp/test_1_output.embl'
expected_outp_1_rel_path = 'tests/data/output/Amaranths_ITS_subset_ALIGN.embl'


###########
# CLASSES #
###########


class MiscellaneousTestCases(unittest.TestCase):
    ''' Tests to evaluate miscellaneous operations'''
    
    def test_1_compare_actual_to_expected(self):
        ''' Assert that two multi-line strings are equal. If they 
        are not, show their difference. '''

        cmd_list = ['python2', ' ',
                    base_path, script_path, ' ',
                    '-n', ' ', base_path, infile_nex_1_rel_path, ' ',
                    '-c', ' ', base_path, infile_csv_1_rel_path, ' ',
                    '-o', ' ', base_path, actual_outp_1_rel_path, ' ',
                    '-e', ' ', 'mi.gruenstaeudl@gmail.com']
        try:
            subprocess.check_output(''.join(cmd_list), shell=True)
        except subprocess.CalledProcessError as e:
            print e.output
        # ---
        #stderr = subprocess.check_output(bash_cmd, stderr=subprocess.PIPE,
        #                           shell=True)
        #if len(stderr.split()) > 0:
        #    print stderr
        # ---

        expected_str = open(base_path + expected_outp_1_rel_path).read()
        ## Check if actual output exists
        if os.path.isfile(base_path + actual_outp_1_rel_path):
            actual_str = open(base_path + actual_outp_1_rel_path).read()
            ## Important: Remove actual output so that lines from different tests are not appended!
            os.remove(base_path + actual_outp_1_rel_path)
        else:
            print 'annonex2embl TESTING ERROR: actual_str not found.'

        self.assertTrue(isinstance(expected_str, str),
                'Not a string: ' + expected_outp_1_rel_path)
        self.assertTrue(isinstance(actual_str, str),
                'Not a string: ' + actual_outp_1_rel_path)
        self.assertMultiLineEqual(expected_str, actual_str)


#############
# FUNCTIONS #
#############

########
# MAIN #
########

if __name__ == '__main__':
    unittest.main()
