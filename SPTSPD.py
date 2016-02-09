#!/usr/bin/env python
"""
Submission Preparation Tool for Sequences of Phylogenetic Datasets
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
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic"\
           " Datasets (SPTSPD)"
__version__ = "2016.02.06.2100"

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

# Possible feature table qualifiers as defined by the International Nucleotide 
# Sequence Database Collection (INSDC) 
# http://www.insdc.org/files/feature_table.html#7.3.1
valid_INSDC_qualifiers = ['allele','altitude','anticodon',
    'artificial_location','bio_material','bound_moiety','cell_line',
    'cell_type','chromosome','citation','clone','clone_lib','codon_start',
    'collected_by','collection_date','compare','country','cultivar',
    'culture_collection','db_xref','dev_stage','direction','EC_number',
    'ecotype','environmental_sample','estimated_length','exception',
    'experiment','focus','frequency','function','gap_type','gene',
    'gene_synonym','germline','haplogroup','haplotype','host',
    'identified_by','inference','isolate','isolation_source','lab_host',
    'lat_lon','linkage_evidence','locus_tag','macronuclear','map',
    'mating_type','mobile_element_type','mod_base','mol_type','ncRNA_class',
    'note','number','old_locus_tag','operon','organelle','organism','partial',
    'PCR_conditions','PCR_primers','phenotype','plasmid','pop_variant',
    'product','protein_id','proviral','pseudo','pseudogene','rearranged',
    'regulatory_class','replace','ribosomal_slippage','rpt_family','rpt_type',
    'rpt_unit_range','rpt_unit_seq','satellite','segment','serotype','serovar',
    'sex','specimen_voucher','standard_name','strain','sub_clone',
    'sub_species','sub_strain','tag_peptide','tissue_lib','tissue_type',
    'transgenic','translation','transl_except','transl_table','trans_splicing',
    'type_material','variety']

###########
# CLASSES #
###########

class MetaQualChecks:
    
    def __init__(self, qualifiers_full):
        self.quals = qualifiers_full
    
    def specific_label_present(self, label):
        ''' Check if all of lines contain a field labelled with <label>. '''

        if not all(label in dct.keys() for dct in self.quals):
            raise ValueError('SPTSPD ERROR: csv-file does not'\
            'contain a column labeled "%s"' % (label))
    
    def valid_INSDC_qualifiers(self):
        ''' Check if field labels are part of list "allowed_INSDC_qualifiers".
        
        Since all dictionaries have the same set of keys, it is sufficient to 
        check only the first dictionary. '''

        quals_present = self.quals[0].keys()
        not_valid = [q for q in quals_present if q not in\
        valid_INSDC_qualifiers]
        if not_valid:
            raise ValueError('SPTSPD ERROR: The following qualifiers are'\
            'invalid INSDC qualifiers:  "%s"' % (not_valid))


class GenerateSeqRecords:
    
    def __init__(self, alignment_full, qualifiers_full):
        self.alignm = alignment_full
        self.quals = qualifiers_full
    
    def generate_base_record(self, seq_name, seqname_col_label):
        ''' Generate a base SeqRecord '''
        
        current_seq = self.alignm[seq_name]
        try:
            current_qual = [d for d in self.quals\
                if d[seqname_col_label] == seq_name][0]
        except:
            current_qual = []
       
        seq_handle = current_seq
        id_handle = seq_name
        try:
            name_handle = current_qual['organism']
        except:
            name_handle = 'unknown organism'
        descr_handle = name_handle
        #try:
        #    descr_handle = current_qual['organism']
        #except:
        #    descr_handle = 'not specified'

        return SeqRecord(seq_handle, id=id_handle, name=name_handle, 
            description=descr_handle)  

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

########
# MAIN #
########

def main(path_to_nex, path_to_csv, seqname_col_label, outformat):

# STEP 01: Prepare output files
# i. Define output files
    in_fn = extract_filename(path_to_nex)
    out_fn = replace_fileending(in_fn, ".embl")
# ii. Initialize output list
    out_records = []

# STEP 02: Parse input data
# i. Parse data from .nex-file.
    aln = Nexus.Nexus()
    aln.read(path_to_nex)
    charsets_full = aln.charsets
    alignment_full = aln.matrix
# ii. Parse data from .csv-file
    reader = DictReader(open(path_to_csv, 'rb'), delimiter=',', quotechar='"',
        skipinitialspace=True)
    qualifiers_full = list(reader)

# STEP 03: Do quality checks on input data
    qual_checks = MetaQualChecks(qualifiers_full)
# i. Check if qualifier matrix (and, hence, each entry) contains a column 
#    labelled by <seqname_col_label>
    qual_checks.specific_label_present(seqname_col_label)
# ii. Check if column names constitute valid INSDC feature table qualifiers 
#    (http://www.insdc.org/files/feature_table.html#7.3.1)
    qual_checks.valid_INSDC_qualifiers()

# TO DO:
# (iii) Check if sequence_names are also in .nex-file
# (iv) Have all metadata conform to basic ASCII standards (not extended ASCII) !

# STEP 5
# Create SeqRecords for each sequence of the alignment.
    record_handle = GenerateSeqRecords(alignment_full, qualifiers_full)
    for seq_name in alignment_full.keys():
# i. Generate the basic SeqRecord (i.e., no features or annotations yet)
        seq_record = record_handle.generate_base_record(seq_name,
                                                        seqname_col_label)
# IMPROVEMENTS NECESSARY:
# (ii) include db_x in base_record


# STEP 6
# Degap the sequence while maintaing correct annotations; has to occur
# before (!) SeqFeature "source" is created.
        degap_handle = COps.DegapButMaintainAnno(seq_record.seq, charsets_full)
        seq_record.seq, degapped_charsets = degap_handle.degap_2()

# STEP 7
# Create SeqFeature "source" for given seq_record; it is appended to 
# seq_record.features. 
# SeqFeature "source" is critical for submissions to EMBL or GenBank as it 
# contains all the relevant info on collection locality, herbarium voucher, etc.
        start_pos = SeqFeature.ExactPosition(0)
        end_pos = SeqFeature.ExactPosition(len(seq_record))
        seq_feature_Location = SeqFeature.FeatureLocation(start_pos, end_pos)
       
        try:
            seq_feature_Qualifiers = [q for q in qualifiers_full\
                if q[seqname_col_label]==seq_name][0]
        except:
            raise ValueError('SPTSPD ERROR: Unable to generate SeqFeature'\
            '"source"')
        source_feature = SeqFeature.SeqFeature(seq_feature_Location,
            type='source', qualifiers=seq_feature_Qualifiers)
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
                type=type_info, qualifiers={'note':charset_name})

# STEP 8.c
# Add to seq_record
            seq_record.features.append(seq_feature)

# STEP 9
# Save record to list "out_records"
        out_records.append(seq_record)

# STEP 10
# Export all seq_out_records as single file in embl-format

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
    parser.add_argument('-f', '--outformat', 
        help='Available arguments: embl, genbank', 
        default='embl', required=False)
    parser.add_argument('-l', '--label',
        help='How the column specifying the sequence names is labelled',
        default='isolate', required=False)
    parser.add_argument('-V', '--version', 
        help='Print version information and exit', 
        action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

########
# MAIN #
########

main(args.nexus, args.csv, args.label, args.outformat)
