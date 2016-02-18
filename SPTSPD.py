#!/usr/bin/env python
"""
Submission Preparation Tool for Sequences of Phylogenetic Datasets
(SPTSPD)

EMBL Submission Preparation Tool for Sequences of Multiple Sequence Alignments
(EMBL-SPTSMSA)
"""
# CURRENT DESIGN:
# Charset-definitions in .nex-file must indicate annotation type (e.g. 'cds', 
# 'gene', 'rrna', 'trna') in their names.
# 
# INPUT: .nex, .csv
# REQUIREMENTS:
# (i) One of the columns in the csv-file must contain the sequence names 
#     also found in the nex-file.

#####################
# IMPORT OPERATIONS #
#####################


from Bio import SeqIO
#from Bio.Alphabet import generic_dna

#from Bio.Seq import Seq
#from Bio import SeqFeature


import CustomOps as CO

import argparse
import sys


###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'Submission Preparation Tool for Sequences of Phylogenetic '\
           'Datasets (SPTSPD)'
__version__ = '2016.02.18.1100'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################


###########
# CLASSES #
###########

class MyException(Exception):
    pass

"""
class GetGeneInfo:

    def get_gene_product:
        '''    
        >>> from Bio import Entrez
        >>> Entrez.email = "mi.gruenstaeudl@gmail.com"
        #>>> handle = Entrez.esearch(db="gene", retmax=10, term="psbI AND gene[FKEY]")
        >>> handle = Entrez.esearch(db="gene", retmax=10, term="psbI", retmode="text")
        >>> records = handle.read()
        #>>> record["IdList"]
        #['126789333', '37222967', '37222966', '37222965', ..., '61585492']
        
        
        >>> handle = Entrez.esummary(db="gene", id="30367")
        >>> record = Entrez.read(handle)
        >>> handle.close()
        >>> print(record[0]["Id"])
        30367
        >>> print(record[0]["Title"])

        
        >>> handle = Entrez.efetch(db="gene", id=records[0], rettype="gb", retmode="text")
        >>> print(handle.read())
        '''

    def get_gene_id:
        '''    
        '''
"""

#############
# FUNCTIONS #
#############

def extract_filename(in_path):
    ''' Split path from filename '''
    import os
    path, fn =  os.path.split(in_path)
    return fn


def replace_fileending(fn, new_end):
    ''' Replace file ending (e.g. ".csv") with different ending '''
    return fn[:fn.rfind('.')] + new_end


def parse_csv_file(path_to_csv):
    ''' This function parses the nexus file. '''
    from csv import DictReader
    
    reader = DictReader(open(path_to_csv, 'rb'), delimiter=',', quotechar='"',
        skipinitialspace=True)
    qualifiers_full = list(reader)
    return qualifiers_full


def parse_nexus_file(path_to_nex):
    ''' This function parses the nexus file. '''
    from Bio.Nexus import Nexus
    
    aln = Nexus.Nexus()
    aln.read(path_to_nex)
    charsets_full = aln.charsets
    alignment_full = aln.matrix
    return (charsets_full, alignment_full)


def check_quality_of_qualifiers(qualifiers_full, seqname_col_label):
    ''' This function conducts a series of quality checks on the qualifiers 
    list (a list of dictionaries).

    TODO:
        (i) Check if sequence_names are also in .nex-file
        (ii) Have all metadata conform to basic ASCII standards (not 
            extended ASCII)!
    '''
    
    qual_checks = CO.MetaChecks(qualifiers_full)
# i. Check if qualifier matrix (and, hence, each! entry) contains a column 
#    labelled by <seqname_col_label>
    try:
        qual_checks.label_present(seqname_col_label)
    except MyException as e:
        sys.exit('%s SPTSPD ERROR: %s' % ('\n', e))
# ii. Check if column names constitute valid INSDC feature table qualifiers 
#    (http://www.insdc.org/files/feature_table.html#7.3.1)
    try:
        qual_checks.valid_INSDC_quals()
    except MyException as e:
        sys.exit('%s SPTSPD ERROR: %s' % ('\n', e))

def transl_and_check_quality_of_transl(feature, transl_table, seq_record):
    ''' This function conducts a translation of a coding region and checks the
    quality of said translation.
    
    TODO:
        (i) Adjust code of AnnoChecks.check() so that the start position of a 
            subsequent feature is also adjusted.                
        (ii) Adjust the location position in code of AnnoChecks.check() so 
            that the stop codon is also included.
        (iii) SHOULD EXAMPLE 2 NOT RESULT IN A FEATURE LOCATION THAT ENDS 
            AT ExactPosition(5), I.E. AFTER THE STOP CODON ???
    '''

    if feature.type.lower() == 'cds': # Check if feature coding region
        extract = feature.extract(seq_record)
        try:
            transl, loc = CO.AnnoChecks(extract.seq, feature.location, 
                feature.type, seq_record.id, transl_table).check()
            feature.qualifiers["translation"] = transl
            feature.location = loc
        except MyException as e:
            sys.exit('%s SPTSPD ERROR: %s' % ('\n', e))
        return feature


########
# MAIN #
########

def main(path_to_nex, path_to_csv, email_addr, outformat, seqname_col_label,
         transl_table):

# STEP 01: Prepare output files
    in_fn = extract_filename(path_to_nex)
    out_fn = replace_fileending(in_fn, ".embl")  # Define output filename
    out_records = []                             # Initialize output list

# STEP 02: Parse data from .nex-file
    try:
        charsets_full, alignment_full = parse_nexus_file(path_to_nex)
    except:
        sys.exit('%s SPTSPD ERROR: %s' % ('\n',
            'Parsing of .nex-file unsuccessful.'))

# STEP 03: Parse data from .csv-file
    try:
        qualifiers_full = parse_csv_file(path_to_csv)
    except:
        sys.exit('%s SPTSPD ERROR: %s' % ('\n',
            'Parsing .csv-file unsuccessful.'))
            
# STEP 04: Do quality checks on input data
    check_quality_of_qualifiers(qualifiers_full, seqname_col_label)

# STEP 05: Create a full SeqRecord for each sequence of the alignment.
    for seq_name in alignment_full.keys():
        
# i. Select current sequences and current qualifiers
        current_seq = alignment_full[seq_name]
        current_quals = [d for d in qualifiers_full\
            if d[seqname_col_label] == seq_name][0]

# ii. Generate the basic SeqRecord (i.e., without features or annotations)
        seq_record = CO.GenerateSeqRecord(current_seq,
            current_quals).base_record(seqname_col_label, charsets_full)

# iii. Degap the sequence while maintaing correct annotations, which has to 
#      occur before (!) the SeqFeature 'source' is generated.
#      Note: Charsets are identical across all sequences.
        degap_handle = CO.DegapButMaintainAnno(seq_record.seq, charsets_full)
        seq_record.seq, degapped_charsets = degap_handle.degap()
            
# iv. Generate SeqFeature 'source' and append to features list
        source_feature = CO.GenerateSeqFeature().source_feat(len(seq_record),
            current_quals, transl_table)
        seq_record.features.append(source_feature)

# STEP 06: Populate the feature keys with the charset information
#          Note: Each charset represents a dictionary that must be added in 
#          full to the list "SeqRecord.features"
        for charset_name, charset_range in degapped_charsets.items():
            seq_feature = CO.GenerateSeqFeature().regular_feat(charset_name,
                                                               charset_range)
            seq_record.features.append(seq_feature)

# STEP 07: Translate and check quality of translation
        for feature in seq_record.features:
            feature = transl_and_check_quality_of_transl(feature,
                transl_table, seq_record)

# STEP 08: Save completed record to list "out_records"
        out_records.append(seq_record)

# STEP 09: Export all out_records as single file in embl-format
    outp_handle = open(out_fn, 'w')
    SeqIO.write(out_records, outp_handle, outformat)
    outp_handle.close()


############
# ARGPARSE #
############

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="  --  ".join([__author__,
        __copyright__, __info__, __version__]))
    parser.add_argument('-n', '--nexus', 
        help='/path_to_input/test.nex',
        default='/home/michael_science/Desktop/test.nex', required=True)
    parser.add_argument('-c', '--csv', 
        help='/path_to_input/test.csv',
        default='/home/michael_science/Desktop/test.csv', required=True)
    parser.add_argument('-e', '--email', 
        help='Your email address',
        default='mi.gruenstaeudl@gmail.com', required=True)
    parser.add_argument('-f', '--outformat', 
        help='Available arguments: embl, genbank', 
        default='embl', required=False)
    parser.add_argument('-l', '--label',
        help='Which xxx the column specifying the sequence names is labelled with.',
        default='isolate', required=False)
    parser.add_argument('-t', '--table',
        help='Which translation table coding regions shall be translate with.'\
        'For details, see: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi',
        default='11', required=False)
    parser.add_argument('-V', '--version', 
        help='Print version information and exit', 
        action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

# Include selection on topology of submission (linear [default] or circular)

########
# MAIN #
########

main(args.nexus, args.csv, args.email, args.outformat, args.label, args.table)
