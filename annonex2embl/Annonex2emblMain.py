#!/usr/bin/env python2.7

''' Main operations in annonex2embl '''

#####################
# IMPORT OPERATIONS #
#####################

from Bio import SeqIO
#from Bio.Alphabet import generic_dna
#from Bio.Seq import Seq
from Bio import SeqFeature
from collections import OrderedDict

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'annonex2embl'))

import MyExceptions as ME
import CheckingOps as CkOps
import DegappingOps as DgOps
import GenerationOps as GnOps
import ParsingOps as PrOps
import IOOps as IOOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.06.21.2200'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

start_codon = "ATG"
stop_codons = ["TAG", "TAA", "TGA"] # amber, ochre, opal

###########
# CLASSES #
###########

########
# TODO #
########
'''
TODO:
    (i) Include a function to check internet connectivity.

    (ii) Functionality for making ends fuzzy.

'''

#############
# FUNCTIONS #
#############

def annonex2embl(path_to_nex,
                 path_to_csv,
                 email_addr,
                 path_to_outfile,
                 out_format='embl',
                 seqname_col='isolate',
                 transl_table='11'):

# 1. Open outfile
    outp_handle = open(path_to_outfile, 'a')

# 2. Parse data from .nex-file
    try:
        charsets_withgaps, alignm = IOOps.Inp().parse_nexus_file(path_to_nex)
    except ME.MyException as e:
        sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))

# 3. Parse data from .csv-file
    try:
        qualifiers = IOOps.Inp().parse_csv_file(path_to_csv)
    except ME.MyException as e:
        sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))

# 4. Do quality checks on input data
    try:
        CkOps.CheckCoord().quality_of_qualifiers(qualifiers, seqname_col)
    except ME.MyException as e:
        sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))

# 5. Parse out feature key, obtain official gene name and gene product 
    charset_dict = {}
    for charset_name in charsets_withgaps.keys():
        try:
            charset_sym, charset_type, charset_product = PrOps.ParseCharsetName(
                charset_name, email_addr).parse()
        except ME.MyException as e:
            sys.exit('%s annonex2embl ERROR: %s' % ('\n', e))
        
        charset_dict[charset_name] = (charset_sym, charset_type,
            charset_product)

# 6. Create a full SeqRecord for each sequence of the alignment.
#    Work off the sequences alphabetically.
    sorted_seqnames = sorted(alignm.keys())
    for seq_name in sorted_seqnames:

# 6.1. Select current sequences and current qualifiers
        current_seq = alignm[seq_name]
        current_quals = [d for d in qualifiers\
            if d[seqname_col] == seq_name][0]

# 6.2. Generate the basic SeqRecord (i.e., without features or annotations)
        seq_record = GnOps.GenerateSeqRecord(current_seq,
            current_quals).base_record(seqname_col, charsets_withgaps)

# 6.3. Degap the sequence while maintaing correct annotations, which has to 
#      occur before (!) the SeqFeature 'source' is generated.
#      Note: Charsets are identical across all sequences.
        degap_handle = DgOps.DegapButMaintainAnno(seq_record.seq, charsets_withgaps)
        seq_record.seq, charsets_degapped = degap_handle.degap()

# 6.4. Generate SeqFeature 'source' and append to features list
        source_feature = GnOps.GenerateSeqFeature().source_feat(
            len(seq_record), current_quals, transl_table)
        seq_record.features.append(source_feature)

# 6.5. Populate the feature keys with the charset information
#      Note: Each charset represents a dictionary that must be added in 
#      full to the list "SeqRecord.features"
        for charset_name, charset_range in charsets_degapped.items():

# 6.5.1. Convert charset_range into Location Object
            location_object = GnOps.GenerateFeatLoc().make_location(charset_range)

# 6.5.2. Assign a gene product to a gene name
            charset_sym, charset_type, charset_product = charset_dict[charset_name]

# 6.5.3. Generate a regular SeqFeature and append to seq_record.features
#        Note: The position indices for the stop codon are truncated in 
#              this step.
            seq_feature = GnOps.GenerateSeqFeature().regular_feat(charset_sym,
                charset_type, location_object, charset_product)
            seq_record.features.append(seq_feature)

# 6.6. Sort all seq_record.features except the first one (which 
#      constitutes the source feature) by the start position
        sorted_features = sorted(seq_record.features[1:],
            key=lambda x: x.location.start.position)
        seq_record.features = [seq_record.features[0]] + sorted_features

# 6.7. Translate and check quality of translation
        removal_list = []
        for indx, feature in enumerate(seq_record.features):
            # Check if feature is a coding region
            if feature.type == 'CDS' or feature.type == 'gene':
                try:
                    feature = CkOps.CheckCoord().transl_and_quality_of_transl( \
                        seq_record, feature, transl_table)
                except ME.MyException as e:
                    print('%s annonex2embl WARNING: %s Feature `%s` '\
                        '(type: `%s`) of sequence `%s` is not saved to '\
                        'output.' % ('\n', e, feature.id, feature.type,
                        seq_record.id))
                    removal_list.append(indx)
        # TFL removes the objects in reverse order, as each removal
        # shifts the indices of subsequent objects to the left
        for indx in sorted(removal_list, reverse=True):
            seq_record.features.pop(indx)

# 6.8. Introduce fuzzy ends
        for feature in seq_record.features:
            # Check if feature is a coding region
            if feature.type == 'CDS' or feature.type == 'gene':
                # Note: Don't use "feature.extract(seq_record.seq)" in TFL,
                #       as stop codon was truncated from feature under 
                #       Step 6.5.3.
                coding_seq = ''.join([seq_record.seq[i] for i in charset_range])
                if not coding_seq.startswith(start_codon):
                    feature.location = GnOps.GenerateFeatLoc(
                        ).make_start_fuzzy(feature.location)
                if all([not coding_seq.endswith(c) for c in stop_codons]):
                    feature.location = GnOps.GenerateFeatLoc(
                        ).make_end_fuzzy(feature.location)

# 6.9. Write each completed record to file
        try:
            SeqIO.write(seq_record, outp_handle, out_format)
        except:
            sys.exit('%s annonex2embl ERROR: Problem with `%s`. Did not write to file.' % ('\n', seq_name))

# 7. Close outfile
    outp_handle.close()

