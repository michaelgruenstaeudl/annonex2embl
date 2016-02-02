#!/usr/bin/env python
"""
Submission Preparation Tool for Sequences of Phylogenetic Datasets
"""
# CURRENT DESIGN:
# Charset-definitions in .nex-file must indicate annotation type (e.g. 'cds', 'gene', 'rrna', 'trna') in their names.
# 
# INPUT: .nex, .csv

#####################
# IMPORT OPERATIONS #
#####################

from Bio import SeqFeature
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from csv import DictReader

import argparse
import sys

###############
# AUTHOR INFO #
###############

__author__ = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
__copyright__ = "Copyright (C) 2016 Michael Gruenstaeudl"
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic Datasets (SPTSPD)"
__version__ = "2016.02.02.2000"

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


class DegapButMaintainAnno:
    ''' This function degaps DNA sequences while maintaining annotations.
    
    This function removes dashes from strings while maintaining 
    annotations on these strings.
    
    Args:
        alignm (dict):  a dict with sequence names (str) as keys and 
                        sequences (str) as values; example: 
                        {seq_name: seq_record.seq}
        annot (dict):   a dictionary with sequence names (str) as keys 
                        and annotation info (dict) as values; example: 
                        {seq_name:charsets}
    '''
    
    def __init__(self, alignm, annot):
        self.alignm = alignm
        self.annot = annot
    
    def standard(self):
        ''' This function conducts the default operation.
        
        Returns:
            tupl.   The return consists of the degapped sequence and the 
                    corresponding annotations; example: 
                    {degapped_seq:adjusted_charsets}
        
        Raises:
            currently nothing
        
        Examples:
            >>> alignm = {"seq_1":"ATG-C"}
            >>> annot = {"seq_1":{"gene_1":[0,1],"gene_2":[2,3,4]}}
            >>> DegapButMaintainAnno(alignm, annot).standard()
            Out: ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})
            
        TODO:
            (i)     Why is the seq_name provided in the input. It is not used in the code.
            (ii)    The loop over self.alignm is unnecessary, since there is only a single key-value pair anyway.
            (iii)   Something is fishy, since the code generates negative numbers.
        '''
    
        for seq_name, seq in self.alignm.items():
            seq_out = ''
            gaps_in_gene = 0
            for s_name, index_list in self.annot[seq_name].items():
                gaps_in_list = 0
                for index, nucl in enumerate(seq):
                    if nucl == '-' and index in index_list:
                        index_list.remove(index)
                        gaps_in_list += 1
                    if nucl != '-' and index in index_list:
                        seq_out += nucl
                        index_list[index_list.index(index)] = index - gaps_in_list
                index_list = [i-gaps_in_gene for i in index_list]
                self.annot[seq_name][s_name] = index_list
                gaps_in_gene += gaps_in_list
            self.alignm[seq_name] = ''.join(seq_out)
            
            seq_degap = self.alignm.values()[0]
            charset_degap = self.annot.values()[0]
        return (seq_degap, charset_degap)

###########
# MODULES #
###########

def main(inFn_nex, inFn_csv, outformat):

# STEP 1
# Define the input files; initialize the output file.

# STEP 2
# Parse data from .nex-file.
    aln = Nexus.Nexus()
    aln.read(inFn_nex)
    charsets_full = aln.charsets
    alignment_full = aln.matrix

# STEP 3
# Parse data from .csv-file.
    reader = DictReader(open(inFn_csv, 'rb'), delimiter=',', quotechar='"')
    qualifiers_full = list(reader)

# STEP 4
# Create SeqRecords by providing sequence names and the sequences 
# themselves.
# Potential Improvements: an automated check if each "seq_name" is in
# "qualifiers_full"; also: automatically parse out taxon name
    for seq_name in alignment_full.keys():
        seq_record = SeqRecord(alignment_full[seq_name],
                               id=seq_name,
                               name=seq_name,
                               description=seq_name)

# STEP 5
# Degap the sequence while maintaing correct annotations; has to occur
# before (!) SeqFeature "source" is created.
        seq_degap, charset_degap = DegapButMaintainAnno({seq_name: seq_record.seq}, {seq_name:charsets_full}).standard()
        seq_record.seq = Seq(seq_degap, generic_dna)

# STEP 6
# Create SeqFeature "source" for given seq_record; is appended to 
# seq_record.features
        start_pos = SeqFeature.ExactPosition(1)
        end_pos = SeqFeature.ExactPosition(len(seq_record))
        seq_feature_Location = SeqFeature.FeatureLocation(start_pos, end_pos)
        try:
            seq_feature_Qualifiers = [q for q in qualifiers_full if q['sequence_name']==seq_name][0]
        except:
            sys.exit('Unable to generate SeqFeature "source"')
        source_feature = SeqFeature.SeqFeature(seq_feature_Location,
                                               type='source',
                                               qualifiers=seq_feature_Qualifiers)
        seq_record.features.append(source_feature)

# STEP 7
# Convert each charset (a dictionary) to a list element in the list SeqRecord.features
        for charset_name, charset_range in charset_degap.items():

# STEP 7.a
# Define the locations of the charsets
# Potential Improvements: MORE PRECISE POSITIONS OF CHARSETS (E.G. 
# AUTOMATIC IDENTIFICATION OF START CODON).
            #start_pos = SeqFeature.AfterPosition(charset_range[0])
            #end_pos = SeqFeature.BeforePosition(charset_range[-1])
            start_pos = SeqFeature.ExactPosition(charset_range[0])
            end_pos = SeqFeature.ExactPosition(charset_range[-1])
            seq_feature_Location = SeqFeature.FeatureLocation(start_pos, end_pos)

# STEP 7.b
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

# STEP 7.c
# Add to seq_record
            seq_record.features.append(seq_feature)

# STEP 8
# Export each seq_record in .gbf-format
        output_handle = open(seq_name + '.ena', 'w')
        SeqIO.write(seq_record, output_handle, outformat)
        output_handle.close()


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
    parser.add_argument('-V', '--version', 
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)
    args = parser.parse_args()


########
# MAIN #
########

main(args.nexus, args.csv, args.outformat)
