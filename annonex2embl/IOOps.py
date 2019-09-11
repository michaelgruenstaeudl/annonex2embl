#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import os

from csv import DictReader
from Bio.Nexus import Nexus
from Bio import SeqIO
from termcolor import colored

try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2019 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2019.09.11.1800'

#############
# DEBUGGING #
#############

import pdb
# pdb.set_trace()

###########
# CLASSES #
###########

class Inp:
    ''' This class contains functions to conduct miscellaneous input
        operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def extract_fn(self, in_path):
        ''' This function splits a the path from path+filename. '''
        path, fn = os.path.split(in_path)
        return fn

    def repl_fileend(self, fn, new_end):
        ''' This function replaces the file ending (e.g. ".csv") with a
        different ending. '''
        return fn[:fn.rfind('.')] + '.' + new_end


    def parse_csv_file(self, path_to_csv):
        ''' This function parses a csv file. '''
        try:
            reader = DictReader(open(path_to_csv, 'r'), delimiter=',',
                                quotechar='"', skipinitialspace=True)
            a_matrix = list(reader)
        except Exception as e:
            print(('\n annonex2embl ERROR: %s:\n %s' % (colored('Parsing of '
            '.csv-file unsuccessful', 'red'), e)))
            raise e
        return a_matrix

    def parse_nexus_file(self, path_to_nex):
        ''' This function parses a NEXUS file. '''
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_nex)
            charsets = aln.charsets
            matrix = aln.matrix
        except Nexus.NexusError as ne:
            raise ne
        except Exception as e:
            print(('\n annonex2embl ERROR: %s:\n %s' % (colored('Parsing of '
            '.nex-file unsuccessful', 'red'), e)))
            raise e
        return (charsets, matrix)


class Outp:
    ''' This class contains two functions for various output operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def write_SeqRecord(self, seq_name, seq_record, outp_handle, ENAstrict_bool):
        ''' This function writes a seqRecord in EMBL format. 
            Upon request (ENAstrict_bool), it masks the ID and AC
            lines as requested by ENA for submissions.
        Args:
            seq_name (str)
            seq_record (obj)
            outp_handle (obj)
            ENAstrict_bool(bool)
        Returns:
            currently nothing
        Raises:
            -
        '''

        SecRecord_handle = StringIO()
        try:
            SeqIO.write(seq_record, SecRecord_handle, 'embl')
        except Exception as e:
            print(('\n annonex2embl ERROR: %s: %s. Did not write to '
            'internal handle.\n %s' % (colored('Problem with sequence ', 'red'), seq_name, e)))
            raise e

        if ENAstrict_bool:
            SecRecord_handle_lines = SecRecord_handle.getvalue().splitlines()
            if SecRecord_handle_lines[0].split()[0] == 'ID':
                ID_line = SecRecord_handle_lines[0]
                ID_line_parts = ID_line.split('; ')
                if len(ID_line_parts) == 7:
                    ID_line_parts = ['XXX' if ID_line_parts.index(p) in
                                     [0, 1, 3, 4, 5, 6] else p for p in ID_line_parts]
                SecRecord_handle_lines[0] = 'ID   ' + '; '.join(ID_line_parts)
            if SecRecord_handle_lines[2].split()[0] == 'AC':
                SecRecord_handle_lines[2] = 'AC   XXX;'
            SecRecord_handle_new = '\n' + '\n'.join(SecRecord_handle_lines)
            SecRecord_handle.truncate(0)
            SecRecord_handle.write(SecRecord_handle_new)
        else:
            pass

        outp_handle.write(SecRecord_handle.getvalue())
        SecRecord_handle.close()


    def create_manifest(self, path_to_outfile, manifest_study, manifest_name, manifest_flatfile):
        ''' This function writes a manifest file. '''
        
        manifest_fn = ''.join(path_to_outfile.split('.')[:-1]) + '.manifest'
        manifest = open(manifest_fn, "w")
        manifest.write(("STUDY\t %s\n") % (manifest_study))
        manifest.write(("NAME\t %s\n") % (manifest_name))
        manifest.write(("FLATFILE\t %s\n") % (manifest_flatfile))
        manifest.close()
