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
from copy import copy

###############
# AUTHOR INFO #
###############

__author__ = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
__copyright__ = "Copyright (C) 2016 Michael Gruenstaeudl"
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic Datasets (SPTSPD)"
__version__ = "2016.02.05.2000"

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

def DegapButMaintainAnno_2(seq, charsets):
    ''' This function degaps DNA sequences while maintaining annotations. Since 
    it works on overlapping charsets, this function replaces DegapButMaintainAnno_1.
    
    Source: http://stackoverflow.com/questions/35233714/maintaining-overlapping-annotations-while-removing-dashes-from-string
    
        Examples:
                  
        Example 4: # Overlapping genes with internal gaps
            >>> seq = "A--AT--T"
            >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
            >>> rewriteGene(seq, charsets)
            Out: ('AATTT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})
                
        Example 5: # Overlapping genes with start and end gaps
            >>> seq = "AA----TT"
            >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
            >>> rewriteGene(seq, charsets)
            Out: ('AATT', {'gene1': [0, 1], 'gene2': [1, 2]})
       
    '''
    
    annotations = copy(charsets)
    index = seq.find('-')
    while index > -1:
        for gene_name, indices in annotations.items():
            if index in indices:
                indices.remove(index)
            annotations[gene_name] = [e-1 if e > index else e for e in indices]
        seq = seq[:index] + seq[index+1:]
        index = seq.find('-')
    return seq, annotations


def intersection_exists(ranges_list):
    init_range = set(ranges_list[0])
    for r in ranges_list[1:]:
        if init_range.intersection(r):
            return True # Exits the entire function with 'True' (i.e., stops all of the loops)
    return False
    
def DegapButMaintainAnno_1(seq, charsets):
    ''' This function degaps DNA sequences while maintaining annotations.
    
    Specifically, this function removes dashes from strings while maintaining 
    annotations on these strings. In its current implementation, this function 
    only works if none of the charsets are overlapping; hence, the initial check.

    Args:
        seq (str):      a string that represents an aligned DNA sequence;
                        example: "ATG-C"
        charsets (dict):a dictionary with gene names (str) as keys and lists 
                        of nucleotide positions (list) as values; example: 
                        {"gene_1":[0,1],"gene_2":[2,3,4]}

    Returns:
        tupl.   The return consists of the degapped sequence and the 
                corresponding degapped charsets; example: 
                (degapped_seq, degapped_charsets)
    
    Raises:
        currently nothing
        
    Examples:
    
        Example 1: # Contains an internal gap
            >>> seq = "ATG-C"
            >>> annot = {"gene_1":[0,1],"gene_2":[2,3,4]}
            >>> DegapButMaintainAnno(seq, annot)
            Out: ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})
        
        Example 2: # Contains start and end gaps
            >>> seq = "AA----TT"
            >>> annot = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
            >>> DegapButMaintainAnno(seq, annot)
            Out: ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})
                
        Example 3: # Entire genes missing
            >>> seq = "AA----TT"
            >>> annot = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
            >>> DegapButMaintainAnno(seq, annot)
            Out: ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})
                
        Example 4: # Overlapping genes with internal gaps
            >>> seq = "A--AT--T"
            >>> annot = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
            >>> DegapButMaintainAnno(seq, annot)
            Out: ('AATTT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})
                
        Example 5: # Overlapping genes with start and end gaps
            >>> seq = "AA----TT"
            >>> annot = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
            >>> DegapButMaintainAnno(seq, annot)
            Out: ('AATT', {'gene1': [0, 1], 'gene2': [1, 2]})
        
        Example 6: # Contains start and end gaps; incorrect charset order
            >>> seq = "AT----GC"
            >>> annot = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
            >>> DegapButMaintainAnno(seq, annot)
            Out: ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})
        
    TODO:
        (i)   Error in example 4: Out: AATT[sic!]T
        (ii)  Error in example 5: Out: {'gene1': [0, 1], 'gene2': [1[sic!], 2]})
    
    Notes:
        (i) In its current implementation, this function only works if 
        none of the charsets are overlapping; hence, the initial check.
        (ii) Order of charset in charsets seems to be irrelevant.
        
    '''
    
    if intersection_exists(charsets.values()):
        raise ValueError('MY ERROR: Character sets are overlapping.')

    degapped_seq = ''
    degapped_charsets = {}
    gaps_cumulative = 0
    for gene_name, index_list in charsets.items():
        gaps_within_gene = 0
        for pos, nucl in enumerate(seq):
            if pos in index_list and nucl == '-':
                index_list.remove(pos)
                gaps_within_gene += 1
            if pos in index_list and nucl != '-':
                degapped_seq += nucl
                index_list[index_list.index(pos)] = pos - gaps_within_gene
        index_list = [i-gaps_cumulative for i in index_list]
        degapped_charsets[gene_name] = index_list
        gaps_cumulative += gaps_within_gene
    return (degapped_seq, degapped_charsets)


########
# MAIN #
########

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
        degapped_seq, degapped_charsets = DegapButMaintainAnno_2(seq_record.seq, charsets_full)
        seq_record.seq = Seq(degapped_seq, generic_dna)

# STEP 6
# Create SeqFeature "source" for given seq_record; is appended to 
# seq_record.features
        start_pos = SeqFeature.ExactPosition(1)
        end_pos = SeqFeature.ExactPosition(len(seq_record))
        seq_feature_Location = SeqFeature.FeatureLocation(start_pos, end_pos)
        try:
            seq_feature_Qualifiers = [q for q in qualifiers_full if q['sequence_name']==seq_name][0]
        except:
            raise ValueError('MY ERROR: Unable to generate SeqFeature "source"')
        source_feature = SeqFeature.SeqFeature(seq_feature_Location,
                                               type='source',
                                               qualifiers=seq_feature_Qualifiers)
        seq_record.features.append(source_feature)

# STEP 7
# Convert each charset (a dictionary) to a list element in the list SeqRecord.features
        for charset_name, charset_range in degapped_charsets.items():

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
