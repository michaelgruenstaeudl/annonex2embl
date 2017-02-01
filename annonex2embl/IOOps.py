#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.01.24.1800'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

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
        import os
        path, fn =  os.path.split(in_path)
        return fn

    def repl_fileend(self, fn, new_end):
        ''' This function replaces the file ending (e.g. ".csv") with a 
        different ending. '''
        return fn[:fn.rfind('.')] + '.' + new_end

    def parse_csv_file(self, path_to_csv):
        ''' This function parses a csv file. '''
        from csv import DictReader
        try:
            reader = DictReader(open(path_to_csv, 'rb'), delimiter=',', 
                quotechar='"', skipinitialspace=True)
            a_matrix = list(reader)
        except:
            raise ME.MyException('Parsing of .csv-file unsuccessful.')
        return a_matrix    

    def parse_nexus_file(self, path_to_nex):
        ''' This function parses a NEXUS file. '''
        from Bio.Nexus import Nexus
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_nex)
            charsets = aln.charsets
            matrix = aln.matrix
        except:
            raise ME.MyException('Parsing of .nex-file unsuccessful.')
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

    def write_EntryUpload(self, seq_record, outp_handle, eusubm_bool):
        ''' This function writes a seqRecord in ENA format for a submission
            via Entry Upload. Upon request (eusubm_bool), it also masks the ID and AC 
            lines as requested by ENA for submissions.
        Args:
            seq_record (obj)
            outp_handle (obj)
            eusubm_bool(str)
        Returns:
            currently nothing
        Raises:
            -
        '''
        from StringIO import StringIO
        from Bio import SeqIO

        temp_handle = StringIO()
        try:
            SeqIO.write(seq_record, temp_handle, 'embl')
        except:
            raise ME.MyException('%s annonex2embl ERROR: Problem with \
            `%s`. Did not write to internal handle.' % ('\n', seq_name))
        
        if eusubm_bool:
            temp_handle_lines = temp_handle.getvalue().splitlines()
            if temp_handle_lines[0].split()[0] == 'ID':
                ID_line = temp_handle_lines[0]
                ID_line_parts = ID_line.split('; ')
                if len(ID_line_parts) == 7:
                    ID_line_parts = ['XXX' if ID_line_parts.index(p) in \
                        [0,1,3,4,5,6] else p for p in ID_line_parts]
                temp_handle_lines[0] = 'ID   ' + '; '.join(ID_line_parts)
            if temp_handle_lines[2].split()[0] == 'AC':
                temp_handle_lines[2] = 'AC   XXX;'
            temp_handle_new = '\n' + '\n'.join(temp_handle_lines)
            temp_handle.truncate(0)
            temp_handle.write(temp_handle_new)
        else:
            pass

        outp_handle.write(temp_handle.getvalue())
        temp_handle.close()
        
        #return something?


class ENAchecklist:
    ''' This class writes checklist in ENA format for a submission
        via ENA's checklist system.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass


    def matK_trnK(self, seq_record, counter, outp_handle):
        ''' This function writes a checklist in ENA format for a 
            submission via ENA's checklist system.
        Args:
            seq_record (obj)
            counter (int)
            outp_handle (obj)
            checklist_type (str)
        Returns:
            currently nothing
        Raises:
            -
        '''
        import Bio
        
        #ENTRYNUMBER
        entrynumber = str(counter+1) # enumerate counter starts counting at 0
            
        #ORGANISM_NAME
        organism_name = seq_record.name
            
        # trnK_intron
        trnK_intron = [f for f in seq_record.features \
            if f.id=='trnK' and f.type=='intron']
        try:
            trnK_intron = str(trnK_intron[0])
            trnK_intron_present = 'yes'
        except:
            trnK_intron_present = 'no'
            
        # matK
        matK_gene = [f for f in seq_record.features \
            if f.id=='matK' and f.type=='gene']
        try:
            matK_gene = matK_gene[0]
        except:
            try:
                matK_gene = [f for f in seq_record.features \
                    if f.id=='matK' and f.type=='CDS']
            except:
                raise ME.MyException('%s annonex2embl ERROR: Problem \
                    with `%s`. matK gene not found.' % ('\n', seq_name))
        #5'_CDS and 5'_PARTIAL
        fiveprime_cds = str(matK_gene.location.start.position)
        if type(matK_gene.location.start) == Bio.SeqFeature.ExactPosition:
            fiveprime_partial = 'no'
        if type(matK_gene.location.start) == Bio.SeqFeature.BeforePosition:
            fiveprime_partial = 'yes'
        #3'_CDS and 3'_PARTIAL
        threeprime_cds = str(matK_gene.location.end.position)
        if type(matK_gene.location.end) == Bio.SeqFeature.ExactPosition:
            threeprime_partial = 'no'
        if type(matK_gene.location.end) == Bio.SeqFeature.AfterPosition:
            threeprime_partial = 'yes'
        
        qualifiers = seq_record.features[0].qualifiers # source feature is always first in list
        #ISOLATE
        try:
            isolate = qualifiers['isolate']
        except:
            isolate = ''

        #SPEC_VOUCH
        try:
            spec_vouch = qualifiers['specimen_voucher']
        except:
            spec_vouch = ''

        #LOCALITY
        try:
            locality = qualifiers['country'] # tag 'locality' does not exist in INDSC, but 'country' does
        except:
            locality = ''

        #ECOTYPE
        try:
            ecotype = qualifiers['ecotype'] # tag 'locality' does not exist in INDSC, but 'country' does
        except:
            ecotype = ''

        #SEQUENCE
        sequence = str(seq_record.seq)
            
        out_list = [entrynumber,
                    organism_name,
                    fiveprime_cds,
                    threeprime_cds,
                    fiveprime_partial,
                    threeprime_partial,
                    trnK_intron_present,
                    isolate,
                    spec_vouch,
                    locality,
                    ecotype,
                    sequence
                   ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)
