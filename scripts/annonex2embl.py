#!/usr/bin/env python2.7
'''
annonex2embl wrapper
'''

#####################
# IMPORT OPERATIONS #
#####################

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import Annonex2emblMain as AN2EMBLMain

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.01.22.1400'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

########
# TODO #
########

''' Include selection on topology of submission (linear [default] or circular) '''

############
# ARGPARSE #
############

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    
    # Required
    parser.add_argument('-n',
                        '--nexus',
                        help='absolute path to infile; infile in NEXUS format; Example: /path_to_input/test.nex',
                        default='/home/username/Desktop/test.nex',
                        required=True)

    parser.add_argument('-c',
                        '--csv',
                        help='absolute path to infile; infile in CSV format; Example: /path_to_input/test.csv',
                        default='/home/username/Desktop/test.csv',
                        required=True)

    parser.add_argument('-e',
                        '--email',
                        help='Your email address.',
                        default='my.username@gmail.com',
                        required=True)

    parser.add_argument('-o',
                        '--outfile',
                        help='absolute path to outfile; outfile in EMBL format; Example: /path_to_output/test.embl',
                        default='/home/username/Desktop/test.embl',
                        required=True)

    parser.add_argument('-u',
                        '--submmode',
                        help='a logical; enter into submission mode or not',
                        default='False',
                        required=True)

    # Optional
    parser.add_argument('-p',
                        '--topology',
                        help='`circular` or `linear`.', 
                        default='linear',
                        required=False)
                        
    parser.add_argument('-d',
                        '--taxdivision',
                        help='Any of the three letter codes specified in section 3.2 of the EMBL user manual.', 
                        default='PLN',
                        required=False)
    
    parser.add_argument('-f',
                        '--outformat',
                        help='Available arguments: embl, gb', 
                        default='embl',
                        required=False)

    parser.add_argument('-l',
                        '--collabel',
                        metavar='column specifying sequence names',
                        help='Name of column that specifies the sequence names.',
                        default='isolate',
                        required=False)

    parser.add_argument('-t',
                        '--transltable',
                        metavar='translation table',
                        help='Number of the translation table to translate coding regions with.'\
                        'For details, see: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi',
                        default='11',
                        required=False)

    parser.add_argument('-s',
                        '--seqversion',
                        metavar='sequence version',
                        help='An integer',
                        default='1',
                        required=False)

    parser.add_argument('--version', 
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

########
# MAIN #
########

    AN2EMBLMain.annonex2embl(   args.nexus,
                                args.csv,
                                args.email,
                                args.outfile,
                                args.submmode,
                                
                                args.topology,
                                args.taxdivision,
                                args.outformat,
                                args.collabel,
                                args.transltable,
                                args.seqversion )
