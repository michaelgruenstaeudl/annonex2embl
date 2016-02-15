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
from Bio.Nexus import Nexus
#from Bio.Seq import Seq
from Bio import SeqFeature
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
__version__ = "2016.02.09.1500"

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
    ''' Operations to evaluate the quality of metadata '''
    
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


class GenerateSeqRecord:
    ''' Operations to generate SeqRecords '''
        
    def __init__(self, current_seq, current_qual):
        self.seq = current_seq
        self.qual = current_qual
    
    def generate_base_record(self, seqname_col_label, charsets_full):
        ''' Generate a base SeqRecord '''
        from Bio.SeqRecord import SeqRecord
        
        gene_names = [k for k in charsets_full.keys()]
        gene_names = [gn.replace("_", " ") for gn in gene_names]
        gene_names_str = ' and '.join(gene_names)
        
        id_handle = self.qual[seqname_col_label]
        try:
            name_handle = self.qual['organism']
        except:
            name_handle = 'unknown organism'
        descr_handle = name_handle + ' ' + gene_names_str +' DNA.'

        return SeqRecord(self.seq, id=id_handle, name=name_handle, 
            description=descr_handle)




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

########
# MAIN #
########

def main(path_to_nex, path_to_csv, email_addr, outformat, seqname_col_label, transl_table):

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
# i. Check if qualifier matrix (and, hence, each! entry) contains a column 
#    labelled by <seqname_col_label>
    qual_checks.specific_label_present(seqname_col_label)
# ii. Check if column names constitute valid INSDC feature table qualifiers 
#    (http://www.insdc.org/files/feature_table.html#7.3.1)
    qual_checks.valid_INSDC_qualifiers()

# TO DO:
# (iii) Check if sequence_names are also in .nex-file
# (iv) Have all metadata conform to basic ASCII standards (not extended ASCII)!

# STEP 05: Create a full SeqRecord for each sequence of the alignment.
    for seq_name in alignment_full.keys():
# i. Select current sequences and current qualifiers
        try:
            current_seq = alignment_full[seq_name]
            current_qual = [d for d in qualifiers_full\
                if d[seqname_col_label] == seq_name][0]
        except:
            raise ValueError('SPTSPD ERROR: Sequence association between '\
                '%s and %s incorrect' % (path_to_nex, path_to_csv))
# ii. Generate the basic SeqRecord (i.e., without features or annotations)
        record_handle = GenerateSeqRecord(current_seq, current_qual)
        seq_record = record_handle.generate_base_record(seqname_col_label,
                                                        charsets_full)
# TO DO:
# (ii) include db_x in base_record

# iii. Degap the sequence while maintaing correct annotations
#     Note: Degapping has to occur before (!) the SeqFeature "source" is
#           generated.
#     Note: Charsets are identical across all sequences

        degap_handle = COps.DegapButMaintainAnno(seq_record.seq, charsets_full)
        seq_record.seq, degapped_charsets = degap_handle.degap_2()

# iv. Create a "source" feature for the seq_record.
#      Note: The SeqFeature "source" is critical for submissions to EMBL or 
#            GenBank, as it contains all the relevant info on collection 
#            locality, herbarium voucher, etc.
        feature_loc = COps.GenerateFeatureLocation(0, len(seq_record)).exact()
        source_feature = SeqFeature.SeqFeature(feature_loc, type='source',
            qualifiers=current_qual)

# v. Append translation table info to source qualifiers
        source_feature.qualifiers["transl_table"]=transl_table
            
# vi. Append to features list
        seq_record.features.append(source_feature)

# vii. Populate the feature keys with the charset information
#     Note: Each charset represents a dictionary that is added to the 
#           list "SeqRecord.features"
        for charset_name, charset_range in degapped_charsets.items():

# a. Define the locations of the charsets
            feature_loc = COps.GenerateFeatureLocation(charset_range[0],
                charset_range[-1]+1).exact()
# TO DO: 
# b. Include a greater number of possible feature location functions.
            #start_pos = SeqFeature.AfterPosition(charset_range[0])
            #end_pos = SeqFeature.BeforePosition(charset_range[-1])

# TO DO:
# c. AUTOMATICALLY IDENTIFY SEQFEATURE (E.G. SEARCH FOR TYPE IN DATABASE)

# vii. Define the annotation type
            anno_types = ['cds', 'gene', 'rrna', 'trna']
            keyw_present = [keyw for keyw in anno_types if keyw in charset_name.lower()]
            if keyw_present:
                type_info = keyw_present[0]
            else:
                type_info = 'misc_feature'
            seq_feature = SeqFeature.SeqFeature(feature_loc, type=type_info,
                qualifiers={'note':charset_name})

# ix. Append to features list
            seq_record.features.append(seq_feature)


# STEP 06: Perform translation and quality control on coding regions
        for feature in seq_record.features:
            if feature.type.lower() == 'cds': # Check if feature coding region
                extract = feature.extract(seq_record)
                transl, loc = COps.AnnoChecks(extract.seq, feature.location, 
                    feature.type, seq_record.id, transl_table).check()
                feature.qualifiers["translation"] = transl
                feature.location = loc

# STEP 07: Save completed record to list "out_records"
        out_records.append(seq_record)

# STEP 08: Export all out_records as single file in embl-format
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
