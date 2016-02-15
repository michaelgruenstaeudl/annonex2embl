#!/usr/bin/env python
'''
Custom operations for EMBL submission preparation tool
'''

#####################
# IMPORT OPERATIONS #
#####################

from copy import copy

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation

###############
# AUTHOR INFO #
###############

__author__ = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
__copyright__ = "Copyright (C) 2016 Michael Gruenstaeudl"
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic Datasets (SPTSPD)"
__version__ = "2016.02.15.1900"

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

class AnnoChecks:
    ''' This class contains functions to evaluate the quality of an 
    annotation.
    
    Args:
        extract (obj):      a sequence object; example: Seq('ATGGAGTAA', 
                            IUPACAmbiguousDNA())

        location (obj):     a location object; example: FeatureLocation(
                            ExactPosition(0), ExactPosition(8))

        feature_type (str): a string detailing the type of the feature;
                            example: "cds"
        
        record_id (str):    a string deatiling the name of the sequence in 
                            question; example: "taxon_A"

        transl_table (int): an integer; example: 11 (for bacterial code)
    
    Returns:
        tupl.   The return consists of the translated sequence (a str) and the
                updated feature location (a location object); example: 
                (transl_out, feat_loc)
    
    Raises:
        -
    '''

    def __init__(self, extract, location, feature_type="foobar", 
                 record_id="foobar", transl_table=11):
        self.e = extract
        self.l = location
        self.f = feature_type
        self.i = record_id
        self.t = transl_table

    @staticmethod
    def _transl(extract, transl_table, to_stop=False, cds=False):
        ''' An internal static function to translate a coding region. '''
        transl = extract.translate(table=transl_table, to_stop=to_stop,
            cds=cds)
        return transl

    @staticmethod            
    def _check_protein_start(extract, transl_table):
        ''' An internal static function to translate a coding region and check
        if it starts with a methionine. '''            
        transl = extract.translate(table=transl_table)
        return transl.startswith("M")

    @staticmethod    
    def _adjust_feat_loc(location, with_internalStop, without_internalStop):
        ''' An internal static function to adjust the feature location if an
        internal stop codon were present. '''
        if len(without_internalStop) > len(with_internalStop):
            start_pos = location.start
            stop_pos = start_pos + (len(with_internalStop) * 3)
            feat_loc = GenerateFeatureLocation(start_pos, stop_pos).exact()
        if len(without_internalStop) == len(with_internalStop):
            feat_loc = location
        return feat_loc

    def check(self):
        ''' This function performs checks on a coding region.
        
        Specifically, the function tries to translate the coding region (CDS)
        directly, using the internal checker "cds=True". If a direct 
        translation fails, it confirms if the CDS starts with a methionine. If 
        the CDS does not start with a methionine, a ValueError is raised. If 
        the CDS does start with a methionine, translations are conducted with 
        and without regard to internal stop codons. The shorter of the two 
        translations is kept. The feature location is adjusted, where 
        necessary.
        
        Note:
            The asterisk indicating a stop codon is truncated under 
            _transl(to_stop=True) and must consequently be added again.
        
        Examples:
        
            Example 1: # Default behaviour
                >>> from Bio.Seq import Seq
                >>> from Bio.Alphabet import generic_dna
                >>> extract = Seq("ATGGCCTAA", generic_dna)
                >>> from Bio.SeqFeature import FeatureLocation
                >>> location = FeatureLocation(0, 8)
                >>> AnnoChecks(extract, location).check()
                Out: (Seq('MA*', ExtendedIUPACProtein()),
                     FeatureLocation(ExactPosition(0), ExactPosition(8)))
            
            Example 2: # Internal stop codon
                >>> from Bio.Seq import Seq
                >>> from Bio.Alphabet import generic_dna
                >>> extract = Seq("ATGTAATAA", generic_dna)
                >>> from Bio.SeqFeature import FeatureLocation
                >>> location = FeatureLocation(0, 8)
                >>> AnnoChecks(extract, location).check()
                Out: (Seq('M*', ExtendedIUPACProtein()),
                     FeatureLocation(ExactPosition(0), ExactPosition(3)))
            
            NOTE: SHOULD EXAMPLE 2 NOT RESULT IN A FEATURE LOCATION THAT ENDS 
                  AT ExactPosition(5), I.E. AFTER THE STOP CODON ???
             
            Example 3: # Does not start with a Methionine 
                >>> from Bio.Seq import Seq
                >>> from Bio.Alphabet import generic_dna
                >>> extract = Seq("AAGTAA", generic_dna)
                >>> from Bio.SeqFeature import FeatureLocation
                >>> location = FeatureLocation(0, 5)
                >>> AnnoChecks(extract, location).check()
                Out: ValueError: SPTSPD ERROR: Feature does not start with a 
                     Methionine.

        TODO:
            (i) Adjust code so that the start position of a subsequent feature
                is also adjusted.
                
            (ii) Adjust the location position so that the stop codon is also
                included.

        '''
        try:
            transl_out = AnnoChecks._transl(self.e, self.t, cds=True)
            feat_loc = self.l
        except:
            if not AnnoChecks._check_protein_start(self.e, self.t):
                raise ValueError('SPTSPD ERROR: Feature "%s" of '\
                    'sequence "%s" does not start with a Methionine.' 
                    % (self.f, self.i))
            else:
                try:
                    without_internalStop = AnnoChecks._transl(self.e, self.t)
                    with_internalStop = AnnoChecks._transl(self.e, self.t,
                        to_stop=True)
                    transl_out = with_internalStop
                    feat_loc = AnnoChecks._adjust_feat_loc(self.l, 
                        with_internalStop, without_internalStop)
                except:
                    raise ValueError('SPTSPD ERROR: Translation of feature '\
                    '"%s" of sequence "%s" not successful.' % (self.f, self.i))
        transl_out = transl_out + "*"
        return (transl_out, feat_loc)
    
    def for_unittest(self):
        transl_out, feat_loc = AnnoChecks(self.e, self.l, self.f, self.i,
            self.t).check()
        if isinstance(transl_out, Seq) and isinstance(feat_loc, 
            FeatureLocation):
            return True
        else:
            return False
        
        


class DegapButMaintainAnno:
    ''' This class contains functions to degap DNA sequences while maintaining 
    annotations. 
    
    Specifically, the functions remove dashes from strings while 
    maintaining annotations on these strings. Only some of the implementations 
    work if the charsets are overlapping.
    
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
    '''
    
    def __init__(self, seq, charsets):
        self.seq = seq
        self.charsets = charsets
        
    def _intersection_exists(ranges_list):
        init_range = set(ranges_list[0])
        for r in ranges_list[1:]:
            if init_range.intersection(r):
                return True # Exits the entire function with 'True' (i.e., stops all of the loops)
        return False

    def degap_1(self):
        '''
        In its current implementation, this function only works if none of the 
        charsets are overlapping; hence, the initial check.
    
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
        seq = self.seq
        charsets = self.charsets
        
        if DegapButMaintainAnno._intersection_exists(charsets.values()):
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

    def degap_2(self):
        ''' This function works on overlapping charsets and is preferable over 
        "degap_1".
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
        seq = self.seq
        charsets = self.charsets
        
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


class GenerateFeatureLocation:
    ''' Operations to generate feature locations
    '''

    def __init__(self, start_pos, stop_pos):
        self.start = start_pos
        self.stop = stop_pos

    def exact(self):
        ''' Generate an exact feature location
        '''
        
        from Bio import SeqFeature
        
        start_pos = SeqFeature.ExactPosition(self.start)
        end_pos = SeqFeature.ExactPosition(self.stop)
        return SeqFeature.FeatureLocation(start_pos, end_pos)


#############
# FUNCTIONS #
#############

########
# MAIN #
########
