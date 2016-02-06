#!/usr/bin/env python
"""
Submission Preparation Tool for Sequences of Phylogenetic Datasets
"""
# CURRENT DESIGN:
# Charset-definitions in .nex-file must indicate annotation type (e.g. 'cds', 'gene', 'rrna', 'trna') in their names.
# 
# INPUT: .nex, .csv
# REQUIREMENTS:
# (i) One of the columns in the csv-file must contain the sequence names also found in the nex-file.

#####################
# IMPORT OPERATIONS #
#####################

from Bio import SeqFeature
from Bio import SeqIO
#from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus
#from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from csv import DictReader

import CustomOps as COps

import argparse

###############
# AUTHOR INFO #
###############

__author__ = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
__copyright__ = "Copyright (C) 2016 Michael Gruenstaeudl"
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic Datasets (SPTSPD)"
__version__ = "2016.02.06.2100"

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

#############
# FUNCTIONS #
#############

########
# MAIN #
########

def main(inFn_nex, inFn_csv, seqname_col_label, outformat):

# STEP 1
# Define output files; initialize output list
    outFn = inFn_nex[:inFn_nex.rfind('.')] + ".embl"
    records = []

# STEP 2
# Parse data from .nex-file.
    aln = Nexus.Nexus()
    aln.read(inFn_nex)
    charsets_full = aln.charsets
    alignment_full = aln.matrix

# STEP 3
# Parse data from .csv-file
    reader = DictReader(open(inFn_csv, 'rb'), delimiter=',', quotechar='"')
    qualifiers_full = list(reader)

# STEP 4
# Do quality checks

# i. Check if all of lines contain a field labelled by seqname_col_label
    if not all(seqname_col_label in dct.keys() for dct in qualifiers_full):
        raise ValueError('SPTSPD ERROR: csv-file does not contain a column labeled "%s"' % (seqname_col_label))
    
# IMPROVEMENTS NECESSARY:
# (ii) Check if sequence_names are also in .nex-file
# (iii) Compare other column titles to list of acceptable SeqFeatures

# STEP 5
# Create SeqRecords by providing sequence names and the sequences 
# themselves.
    for seq_name in alignment_full.keys():
        seq_record = SeqRecord(alignment_full[seq_name],
                               id=seq_name,
                               name=seq_name,
                               description=seq_name)
# IMPROVEMENTS NECESSARY:
# (i) check automatically if each "seq_name" in "qualifiers_full"
# (ii) automatically parse out name of taxon from from qualifiers_full


# STEP 6
# Degap the sequence while maintaing correct annotations; has to occur
# before (!) SeqFeature "source" is created.
        handle = COps.DegapButMaintainAnno(seq_record.seq, charsets_full)
        degapped_seq, degapped_charsets = handle.degap_2()
        seq_record.seq = degapped_seq

# STEP 7
# Create SeqFeature "source" for given seq_record; it is appended to 
# seq_record.features. 
# SeqFeature "source" is critical for submissions to EMBL or GenBank as it 
# contains all the relevant info on collection locality, herbarium voucher, etc.
        start_pos = SeqFeature.ExactPosition(1)
        end_pos = SeqFeature.ExactPosition(len(seq_record))
        seq_feature_Location = SeqFeature.FeatureLocation(start_pos, end_pos)
       
        try:
            seq_feature_Qualifiers = [q for q in qualifiers_full if q['sequence_name']==seq_name][0]
        except:
            raise ValueError('SPTSPD ERROR: Unable to generate SeqFeature "source"')
        source_feature = SeqFeature.SeqFeature(seq_feature_Location,
                                               type='source',
                                               qualifiers=seq_feature_Qualifiers)
        seq_record.features.append(source_feature)

# STEP 8
# Convert each charset (a dictionary) to a list element in the list SeqRecord.features
        for charset_name, charset_range in degapped_charsets.items():

# STEP 8.a
# Define the locations of the charsets
# Potential Improvements: MORE PRECISE POSITIONS OF CHARSETS (E.G. 
# AUTOMATIC IDENTIFICATION OF START CODON).
            #start_pos = SeqFeature.AfterPosition(charset_range[0])
            #end_pos = SeqFeature.BeforePosition(charset_range[-1])
            start_pos = SeqFeature.ExactPosition(charset_range[0])
            end_pos = SeqFeature.ExactPosition(charset_range[-1])
            seq_feature_Location = SeqFeature.FeatureLocation(start_pos, end_pos)

# STEP 8.b
# Define the annotation type
# Potential Improvements: AUTOMATICALLY IDENTIFY SEQFEATURE (E.G. SEARCH
# FOR TYPE IN DATABASE)
            anno_types = ['cds', 'gene', 'rrna', 'trna']
            keyw_present = [keyw for keyw in anno_types if keyw in charset_name.lower()]
            if keyw_present:
                type_info = keyw_present[0]
            else:
                type_info = 'misc_feature'
            seq_feature = SeqFeature.SeqFeature(seq_feature_Location,
                                                type=type_info,
                                                qualifiers={'note':charset_name})

# STEP 8.c
# Add to seq_record
            seq_record.features.append(seq_feature)

# STEP 9
# Save record to list "records"
        records.append(seq_record)

# STEP 10
# Export all seq_records as single file in embl-format

    outp_handle = open(outFn, 'w')
    SeqIO.write(records, outp_handle, outformat)
    outp_handle.close()


############
# ARGPARSE #
############

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument('-n', '--nexus',
                        help='/path_to_input/test.nex',
                        default='/home/michael_science/Desktop/test.nex',
                        required=True)
    parser.add_argument('-c', '--csv',
                        help='/path_to_input/test.csv',
                        default='/home/michael_science/Desktop/test.csv',
                        required=True)
    parser.add_argument('-f', '--outformat',
                        help='Available arguments: embl, genbank',
                        default='embl',
                        required=False)
    parser.add_argument('-l', '--label',
                        help='How the column specifying the sequence names is labelled',
                        default='sequence_name',
                        required=False)
    parser.add_argument('-V', '--version', 
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)
    args = parser.parse_args()


########
# MAIN #
########

main(args.nexus, args.csv, args.label, args.outformat)
