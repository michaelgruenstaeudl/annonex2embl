#!/usr/bin/env python
'''
Custom operations for EMBL submission preparation tool
'''

#####################
# IMPORT OPERATIONS #
#####################

import GlobalVariables as GlobVars
import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2017.01.31.1900'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()

###########
# CLASSES #
###########

class GenerateFeatLoc:
    ''' This class contains functions to generate or manipulate
        SeqFeature location objects.
    '''

    def __init__(self):
        pass

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
        '''
#        Examples:
#            Example 1:
#            >>> compound_integer_range = [1,2,3,7,8,9]
#            >>> _extract_contiguous_subsets(compound_integer_range)
#            Out: [[1, 2, 3], [7, 8, 9]]

        from operator import itemgetter
        from itertools import groupby
        outlist = []
        for k, g in groupby(enumerate(compound_integer_range), lambda (i,x):i-x):
            outlist.append(map(itemgetter(1), g))
        return outlist

    def make_location(self, charset_range):
        ''' This function goes through a decision tree and generates
            fitting feature locations.
        Args:
            charset_range (list): a list of index positions, example: [1,2,3,8,9 ...]
        Returns:
            FeatureLocation (obj):  A SeqFeature location object; either a
                                    FeatureLocation or a CompoundLocation
        Raises:
            -
        '''
        contiguous_ranges = GenerateFeatLoc._extract_contiguous_subsets(
            charset_range)
        # Convert each contiguous range into an exact feature location
        for i,r in enumerate(contiguous_ranges):
            contiguous_ranges[i] = GenerateFeatLoc._exact(r)
        if len(contiguous_ranges) > 1:
            from Bio.SeqFeature import CompoundLocation
            return CompoundLocation(contiguous_ranges)
        else:
            return contiguous_ranges[0]
    
    def make_start_fuzzy(self, location_object):
        ''' This function makes the start position of location 
            objects fuzzy.
        '''
        from Bio import SeqFeature
        if hasattr(location_object, 'parts'):
            if len(location_object.parts) == 1:
                new_start_pos = SeqFeature.BeforePosition(location_object.start)
                location_object = SeqFeature.FeatureLocation(new_start_pos,
                    location_object.end)
            if len(location_object.parts) > 1:
                new_start_pos = SeqFeature.BeforePosition(location_object.parts[0].start)
                location_object.parts[0] = SeqFeature.FeatureLocation(new_start_pos,
                    location_object.parts[0].end)
        return location_object
    
    def make_end_fuzzy(self, location_object):
        ''' This function makes the end position of location 
            objects fuzzy.
        '''
        
#        Examples:
#            Example 1:
#                >>> from Bio import SeqFeature
#                >>> start_pos = SeqFeature.ExactPosition(5)
#                >>> end_pos = SeqFeature.ExactPosition(9)
#                >>> location_object = SeqFeature.FeatureLocation(start_pos, end_pos)
#                >>> location_object
#                Out: FeatureLocation(ExactPosition(5), ExactPosition(9))
#                >>> new_loc = GenerateFeatLoc().make_end_fuzzy(location_object)
#                >>> new_loc
#                Out: FeatureLocation(ExactPosition(5), AfterPosition(9))
#
#            Example 2:
#                >>> from Bio import SeqFeature
#                >>> csrange = [1,2,3,7,8]
#                >>> location_object = GenerateFeatLoc().make_location(csrange)
#                >>> location_object
#                Out: CompoundLocation([FeatureLocation(ExactPosition(1),
#                ExactPosition(4)), FeatureLocation(ExactPosition(7),
#                ExactPosition(9))], 'join')
#                >>> new_loc = GenerateFeatLoc().make_end_fuzzy(location_object)
#                >>> new_loc
#                Out: CompoundLocation([FeatureLocation(ExactPosition(1), ExactPosition(4)), FeatureLocation(ExactPosition(7), AfterPosition(9))], 'join')
        
        from Bio import SeqFeature
        if hasattr(location_object, 'parts'):
            if len(location_object.parts) == 1:
                new_end_pos = SeqFeature.AfterPosition(location_object.end)
                location_object = SeqFeature.FeatureLocation(
                    location_object.start, new_end_pos)
            if len(location_object.parts) > 1:
                new_end_pos = SeqFeature.AfterPosition(location_object.parts[-1].end)
                location_object.parts[-1] = SeqFeature.FeatureLocation(
                    location_object.parts[-1].start, new_end_pos)
        return location_object


class GenerateSeqFeature:
    ''' This class contains functions to generate SeqFeatures. '''
        
    def __init__(self):
        pass
    
    def source_feat(self, full_len, quals, transl_table):
        ''' This function generates the SeqFeature `source` for a 
            SeqRecord. The SeqFeature `source` is critical for 
            submissions to EMBL or GenBank, as it contains all the 
            relevant info on collection locality, herbarium voucher, 
            etc. It also provides info on which translation table is 
            used for subsequent CDS features.
        Args:
            full_len (int): the full length of the seq in question;
                            example: 509
            quals (dict):   a dictionary of qualifiers; example: 
                            {'isolate': 'taxon_B', 'country': 'Ecuador'}
            transl_table (int): an integer; example: 11 (for bacterial code)
        Returns:
            SeqFeature (obj):   A SeqFeature object
        Raises:
            [currently nothing]
        '''
        from Bio import SeqFeature
        full_index = range(0, full_len)
        feature_loc = GenerateFeatLoc().make_location(full_index)
        source_feature = SeqFeature.SeqFeature(feature_loc, id='source',
            type='source', qualifiers=quals)
        # only if a CDS, should the trans_table be integrated
        source_feature.qualifiers["transl_table"] = transl_table
        return source_feature

    def regular_feat(self, feature_name, feature_type, feature_loc,
                     feature_product=None):
        ''' This function generates a regular SeqFeature for a SeqRecord.
        Args:
            feature_name (str):  usually a gene symbol; example: 'matK'
            feature_type (str):  an identifier as to the type of feature; 
                                 example: 'intron'
            feature_loc (object): a SeqFeature object specifying a simple 
                                  or compund location on a DNA string
            feature_product (str): the product of the feature in question;
                                   example: 'maturase K'
        Returns:
            SeqFeature (obj):   A SeqFeature object
        Raises:
            -
        '''
        from Bio import SeqFeature
        # 1. Define the annotation type
        if feature_type not in GlobVars.nex2ena_valid_INSDC_featurekeys:
            raise ME.MyException('%s nex2embl ERROR: Internal error: '\
                'Name of feature key not passed correctly.')
        # 2. Generate qualifiers
        qualifiers={'note':feature_name}
        # 3. Include product, if a coding feature
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
    
    def base_record(self, uniq_seqid, charsets_full):
        ''' This function generates a base SeqRecord (i.e., the foundation to
            subsequent SeqRecords).
        Args:
            uniq_seqid (str): the label of the .csv-file column that
                             contains info on the sequence names; 
                             example: "sequence_name"
            charsets_full (dict):    foobar; example: foobar
        
        TODO:
            (i) include info on linearity of molecule (i.e., linear or 
                circular)
            (ii) include db_x in base_record
        '''
        from Bio.SeqRecord import SeqRecord
        gene_names = [k for k in charsets_full.keys()]
        gene_names = [gn.replace("_", " ") for gn in gene_names]
        gene_names_str = ' and '.join(gene_names)
        ID_line = self.qual[uniq_seqid]
        try:
            org_name = self.qual['organism']
        except:
            org_name = 'undetermined organism'
        DE_line = ' '.join([org_name, gene_names_str, uniq_seqid])
        new_seqRecord = SeqRecord(self.seq, id=ID_line, name=org_name, 
            description=DE_line)
        return new_seqRecord


#############
# FUNCTIONS #
#############

########
# MAIN #
########
