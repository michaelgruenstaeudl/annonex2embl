#!/usr/bin/env python
'''
Custom operations for EMBL submission preparation tool
'''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.06.13.1600'

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

class GenerateFeatLoc:
    ''' This class contains functions to generate feature locations.

    Args:
        charset_range (list): a list of index positions, example: [1,2,3,8,9 ...]
    Returns:
        FeatureLocation (obj):  A SeqFeature location object, either a
                                FeatureLocation or a CompoundLocation
    Raises:
        -
    '''

    def __init__(self, charset_range):
        self.csrange = charset_range

    @staticmethod
    def _exact(csrange):
        ''' An internal static function to generate an exact feature 
        location. '''
        from Bio.SeqFeature import ExactPosition, FeatureLocation
        start_pos = csrange[0]
        stop_pos = csrange[-1]+1
        start_exact = ExactPosition(start_pos)
        stop_exact = ExactPosition(stop_pos)
        return FeatureLocation(start_exact, stop_exact)

    @staticmethod
    def _extract_contiguous_subsets(compound_integer_range):
        ''' An internal static function to extract all contiguous
        integer ranges from compound integer range.
        
        Examples:
            Example 1:
            >>> compound_integer_range = [1,2,3,7,8,9]
            >>> _extract_contiguous_subsets(compound_integer_range)
            Out: [[1, 2, 3], [7, 8, 9]]
        '''
        from operator import itemgetter
        from itertools import groupby
        outlist = []
        for k, g in groupby(enumerate(compound_integer_range), lambda (i,x):i-x):
            outlist.append(map(itemgetter(1), g))
        return outlist

    def make_location(self):
        ''' This function goes through a decision tree and generates
        fitting feature locations.
            
        Examples:
            Example 1: # Default feature location
                >>> csrange = [1,2,3,7,8]
                >>> csrange
                Out: [1, 2, 3, 7, 8]
                >>> GenerateFeatLoc(csrange).make_location()
                Out: CompoundLocation([FeatureLocation(ExactPosition(1),
                ExactPosition(4)), FeatureLocation(ExactPosition(7),
                ExactPosition(9))], 'join')
        '''
        contiguous_ranges = GenerateFeatLoc._extract_contiguous_subsets(
            self.csrange)
        # Convert each contiguous range into an exact feature location
        for i,r in enumerate(contiguous_ranges):
            contiguous_ranges[i] = GenerateFeatLoc._exact(r)
        if len(contiguous_ranges) > 1:
            from Bio.SeqFeature import CompoundLocation
            return CompoundLocation(contiguous_ranges)
        else:
            return contiguous_ranges[0]


class GenerateSeqFeature:
    ''' This class contains functions to generate SeqFeatures. '''
        
    def __init__(self):
        pass
    
    def source_feat(self, full_len, quals, transl_table):
        ''' This function generates the SeqFeature `source` for a SeqRecord.

        The SeqFeature `source` is critical for submissions to EMBL or GenBank, 
        as it contains all the relevant info on collection locality, herbarium 
        voucher, etc. It also provides info on which translation table is used
        for subsequent CDS features.
            
        Args:
            full_len (int): the full length of the seq in question, example: 509
            quals (dict):   a dictionary of qualifiers; example: 
                            {'isolate': 'taxon_B', 'country': 'Ecuador'}
            transl_table (int): an integer; example: 11 (for bacterial code)
        Returns:
            SeqFeature (obj):   A SeqFeature object
        Raises:
            [currently nothing]
            
        Examples:
            Example 1: # Default evaluation
                >>> full_len = 509
                >>> quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
                >>> transl_table = 11
                >>> GenerateSeqFeature().source_feat(full_len, quals, transl_table)
                Out: SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(509)), type='source', id='source')
        '''
        
        from Bio import SeqFeature
        
        full_index = range(0, full_len)
        feature_loc = GenerateFeatLoc(full_index).make_location()
        source_feature = SeqFeature.SeqFeature(feature_loc, id='source',
            type='source', qualifiers=quals)
        source_feature.qualifiers["transl_table"]=transl_table
        return source_feature
    
    def regular_feat(self, feature_name, feature_type, feature_loc,
                     feature_product=None):
        ''' This function generates a regular SeqFeature for a SeqRecord.
            
        Args:
            feature_name (str):  usually a gene symbol; example: 'matK'
            feature_type (str):  an identifier as to the type of feature, 
                              example: 'intron'
            feature_loc (object): a SeqFeature object specifying a simple 
                                 or compund location on a DNA string
            feature_product (str): the product of the feature in question;
                                   example: 'maturase K'
        Returns:
            SeqFeature (obj):   A SeqFeature object
        Raises:
            -
    
        TODO: 
            (i) Automatically identify a SeqFeature (e.g. search for the type in
            a database)
    
        Examples:
            Example 1: #  Evaluates the correct generation of a regular, 
            non-coding SeqFeature.
                >>> feature_name = 'psbI'
                >>> feature_type = 'intron'
                >>> feature_loc = ... # To be completed
                Out: SeqFeature ... # To be completed
        '''
        from Bio import SeqFeature

        INSDC_feature_keys = ["assembly_gap", "C_region", "CDS", "centromere",
        "D-loop", "D_segment", "exon", "gap", "gene", "iDNA", "intron",
        "J_segment", "LTR", "mat_peptide", "misc_binding", "misc_difference",
        "misc_feature", "misc_recomb", "misc_RNA", "misc_structure",
        "mobile_element", "modified_base", "mRNA", "ncRNA", "N_region",
        "old_sequence", "operon", "oriT", "polyA_site", "precursor_RNA",
        "prim_transcript", "primer_bind", "protein_bind", "regulatory", 
        "repeat_region", "rep_origin", "rRNA", "S_region", "sig_peptide",
        "source", "stem_loop", "STS", "telomere", "tmRNA", "transit_peptide",
        "tRNA", "unsure", "V_region", "V_segment", "variation", "3'UTR",
        "5'UTR"]
        
        # a. Define the annotation type
        if feature_type not in INSDC_feature_keys:
            raise ME.MyException('%s nex2embl ERROR: Internal error: Name of '\
                'feature key not passed correctly.')
        # b. Generate qualifiers
        qualifiers={'note':feature_name}
        # c. Include product, if a coding feature
        if feature_product:
            if feature_type == 'CDS' or feature_type == 'gene':
                qualifiers['product'] = feature_product
        seq_feature = SeqFeature.SeqFeature(feature_loc,
            id=feature_name, type=feature_type, qualifiers=qualifiers)
        return seq_feature


class GenerateSeqRecord:
    ''' This class contains functions to generate SeqRecords.
    
    Args:
        current_seq (str): the DNA sequence of; example: 1
        
        current_qual (xxx): foobar; example: foobar
    
    Returns:
        SeqRecord (obj):   A SeqRecord object
    
    Raises:
        [specific to function]
    '''
        
    def __init__(self, current_seq, current_qual):
        self.seq = current_seq
        self.qual = current_qual
    
    def base_record(self, seqname_col_label, charsets_full):
        ''' This function generates a base SeqRecord (i.e., the foundation to
        subsequent SeqRecords).
            
        Args:
            seqname_col_label (str): the label of the .csv-file column that
                                     contains info on the sequence names; 
                                     example: "sequence_name"
        
            charsets_full (dict):    foobar; example: foobar
                
        
        TODO:
            (i) include info on linearity of molecule (i.e., linear or 
                circular)
            (ii) include db_x in base_record
            
        Examples:
        
            Example 1: # foobar
                >>> foobar
                >>> GenerateSeqRecord(foo, bar).generate_base_record(baz, qux)
                Out: ...
        '''

        from Bio.SeqRecord import SeqRecord
        
        gene_names = [k for k in charsets_full.keys()]
        gene_names = [gn.replace("_", " ") for gn in gene_names]
        gene_names_str = ' and '.join(gene_names)
        
        id_handle = self.qual[seqname_col_label]
        try:
            org_name = self.qual['organism']
        except:
            org_name = 'undetermined organism'
        descr_handle = org_name + ' ' + gene_names_str +' DNA.'

        return SeqRecord(self.seq, id=id_handle, name=org_name, 
            description=descr_handle)


#############
# FUNCTIONS #
#############

########
# MAIN #
########
