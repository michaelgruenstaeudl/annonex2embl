#!/usr/bin/env python
'''
Custom operations for EMBL submission preparation tool
'''

#####################
# IMPORT OPERATIONS #
#####################

###############
# AUTHOR INFO #
###############

__author__ = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
__copyright__ = "Copyright (C) 2016 Michael Gruenstaeudl"
__info__ = "Submission Preparation Tool for Sequences of Phylogenetic Datasets (SPTSPD)"
__version__ = "2016.02.17.1900"

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

class MyException(Exception):
    pass

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
        MyException
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
        from Bio.Seq import Seq
        transl = extract.translate(table=transl_table, to_stop=to_stop,
            cds=cds)
        return transl

    @staticmethod            
    def _check_protein_start(extract, transl_table):
        ''' An internal static function to translate a coding region and check
        if it starts with a methionine. '''
        from Bio.Seq import Seq
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
             
            Example 3: # Does not start with a Methionine 
                >>> from Bio.Seq import Seq
                >>> from Bio.Alphabet import generic_dna
                >>> extract = Seq("AAGTAA", generic_dna)
                >>> from Bio.SeqFeature import FeatureLocation
                >>> location = FeatureLocation(0, 5)
                >>> AnnoChecks(extract, location).check()
                Out: ValueError: SPTSPD ERROR: Feature does not start with a 
                     Methionine.
        '''
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation

        try:
            transl_out = AnnoChecks._transl(self.e, self.t, cds=True)
            feat_loc = self.l
        except:
            if not AnnoChecks._check_protein_start(self.e, self.t):
                return MyException('Feature `%s` of sequence `%s` does not '\
                'start with a Methionine.' % (self.f, self.i))
            else:
                try:
                    without_internalStop = AnnoChecks._transl(self.e, self.t)
                    with_internalStop = AnnoChecks._transl(self.e, self.t,
                        to_stop=True)
                    transl_out = with_internalStop
                    feat_loc = AnnoChecks._adjust_feat_loc(self.l, 
                        with_internalStop, without_internalStop)
                except:
                    return MyException('Translation of feature `%s` of '\
                    'sequence `%s` unsuccessful.' % (self.f, self.i))
        transl_out = transl_out + "*"
        return (transl_out, feat_loc)
    
    def for_unittest(self):
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation

        try:
            transl_out, feat_loc = AnnoChecks(self.e, self.l, self.f, self.i,
                self.t).check()
            if isinstance(transl_out, Seq) and isinstance(feat_loc, 
                FeatureLocation):
                return True
            return False
        except ValueError:
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

    @staticmethod
    def _intersection_exists(ranges_list):
        init_range = set(ranges_list[0])
        for r in ranges_list[1:]:
            if init_range.intersection(r):
                return True # Exits the entire function with 'True' (i.e., stops all of the loops)
        return False
    
    def degap(self):
        ''' This function works on overlapping charsets and is preferable over 
        "degap_legacy".
        Source: http://stackoverflow.com/questions/35233714/
        maintaining-overlapping-annotations-while-removing-dashes-from-string
        
            Examples:
            
            Example 1: # Contains an internal gap
                >>> seq = "ATG-C"
                >>> charsets = {"gene_1":[0,1],"gene_2":[2,3,4]}
                >>> DegapButMaintainAnno(seq, annot).degap()
                Out: ('ATGC', {'gene_1': [0, 1], 'gene_2': [2, 3]})
            
            Example 2: # Contains start and end gaps
                >>> seq = "AA----TT"
                >>> charsets = {"gene1":[0,1,2,3], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, annot).degap()
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})
                    
            Example 3: # Entire genes missing
                >>> seq = "AA----TT"
                >>> charsets = {"gene1":[0,1,2], "gene2":[3,4], "gene3":[5,6,7]}
                >>> DegapButMaintainAnno(seq, annot).degap()
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [], 'gene3': [2, 3]})
                      
            Example 4: # Overlapping genes with internal gaps
                >>> seq = "A--AT--T"
                >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, charsets).degap()
                Out: ('AATT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})
                    
            Example 5: # Overlapping genes with start and end gaps
                >>> seq = "AA----TT"
                >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, charsets).degap()
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [2, 3]})
                        
            Example 6: # Contains start and end gaps; incorrect charset order
                >>> seq = "AT----GC"
                >>> charsets = {"gene2":[4,5,6,7], "gene1":[0,1,2,3]}
                >>> DegapButMaintainAnno(seq, annot).degap()
                Out: ('ATGC', {'gene1': [0, 1], 'gene2': [2, 3]})
        '''
        from copy import copy
        
        seq = self.seq
        charsets = self.charsets
        
        annotations = copy(charsets)
        index = seq.find('-')
        while index > -1:
            for gene_name, indices in annotations.items():
                if index in indices:
                    indices.remove(index)
                annotations[gene_name] = [e-1 if e > index else e \
                    for e in indices]
            seq = seq[:index] + seq[index+1:]
            index = seq.find('-')
        return seq, annotations
    
    """
    def degap_legacy(self):
        '''
        In its current implementation, this function only works if none of the 
        charsets are overlapping; hence, the initial check.
    
        Examples:
                   
            Example 4: # Overlapping genes with internal gaps
                >>> seq = "A--AT--T"
                >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, annot).degap_legacy()
                Out: ('AATTT', {'gene1': [0, 1, 2], 'gene2': [2, 3]})
                    
            Example 5: # Overlapping genes with start and end gaps
                >>> seq = "AA----TT"
                >>> charsets = {"gene1":[0,1,2,3,4], "gene2":[4,5,6,7]}
                >>> DegapButMaintainAnno(seq, annot).degap_legacy()
                Out: ('AATT', {'gene1': [0, 1], 'gene2': [1, 2]})
            
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
        
        if _intersection_exists(charsets.values()):
            raise MyException('SPTSPD ERROR: Character sets are overlapping.')
    
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
    """


class GenerateFeatureLocation:
    ''' This class contains functions to generate feature locations.
    
    Args:
        start_pos (int):    the start position of a feature; example: 1
        
        stop_pos (int):    the stop position of a feature; example: 12
    
    Returns:
        FeatureLocation (obj):   A FeatureLocation object
    
    Raises:
        -
    '''

    def __init__(self, start_pos, stop_pos):
        self.start = start_pos
        self.stop = stop_pos

    def exact(self):
        ''' This function generates an exact feature location.
            
        Examples:
                   
            Example 1: # Default feature location
                >>> start_pos = 1
                >>> stop_pos = 12
                >>> GenerateFeatureLocation(start_pos, stop_pos).exact()
                Out: FeatureLocation(ExactPosition(1), ExactPosition(12))
        '''
        from Bio import SeqFeature
        
        start_pos = SeqFeature.ExactPosition(self.start)
        end_pos = SeqFeature.ExactPosition(self.stop)
        return SeqFeature.FeatureLocation(start_pos, end_pos)


class MetaChecks:
    ''' This class contains functions to evaluate the quality of metadata.
    
    Args:
        lst_of_dcts (list): a list of dictionaries; example: 
                            [{'foo': 'foobarqux', 'bar': 'foobarqux', 
                              'qux': 'foobarqux'}, {'foo': 'foobarbaz', 
                              'bar': 'foobarbaz', 'baz': 'foobarbaz'}]    
    Returns:
        none
    
    Raises:
        MyException
    '''
    
    def __init__(self, lst_of_dcts):
        self.lst_of_dcts = lst_of_dcts
    
    def label_present(self, label):
        ''' Check if each (!) list of dictionary keys of a list of dictionaries
        encompass the element <label> at least once.
        
        Examples:

                >>> lst_of_dcts = [{'foo': 'foobarqux', 'bar': 'foobarqux', 
                'qux': 'foobarqux'}, {'foo': 'foobarbaz', 'bar': 'foobarbaz', 
                'baz': 'foobarbaz'}]
        
            Example 1: # Positive confirmation
                >>> label = 'foo'
                >>> MetaQualChecks(lst_of_dcts).label_present(label)
                Out: True
                
            Example 2: # Negative confirmation
                >>> label = 'qux'
                >>> MetaQualChecks(lst_of_dcts).label_present(label)
                Out: MyException: csv-file does not contain a column 
                labelled "qux"
        '''

        if not all(label in dct.keys() for dct in self.lst_of_dcts):
            return MyException('csv-file does not contain a column '\
                'labelled `%s`' % (label))
        return True
    
    def valid_INSDC_quals(self):
        ''' Check if field labels are part of list "allowed_INSDC_qualifiers".
        
        Since all dictionaries have the same set of keys, it is sufficient to 
        check only the first dictionary. '''
        
        # Valid feature table qualifiers as defined by the International
        # Nucleotide Sequence Database Collection (INSDC)
        # http://www.insdc.org/files/feature_table.html#7.3.1
        valid_INSDC_quals = ['allele','altitude','anticodon',
        'artificial_location','bio_material','bound_moiety','cell_line',
        'cell_type','chromosome','citation','clone','clone_lib','codon_start',
        'collected_by','collection_date','compare','country','cultivar',
        'culture_collection','db_xref','dev_stage','direction','EC_number',
        'ecotype','environmental_sample','estimated_length','exception',
        'experiment','focus','frequency','function','gap_type','gene',
        'gene_synonym','germline','haplogroup','haplotype','host',
        'identified_by','inference','isolate','isolation_source','lab_host',
        'lat_lon','linkage_evidence','locus_tag','macronuclear','map',
        'mating_type','mobile_element_type','mod_base','mol_type',
        'ncRNA_class','note','number','old_locus_tag','operon','organelle',
        'organism','partial','PCR_conditions','PCR_primers','phenotype',
        'plasmid','pop_variant','product','protein_id','proviral','pseudo',
        'pseudogene','rearranged','regulatory_class','replace',
        'ribosomal_slippage','rpt_family','rpt_type','rpt_unit_range',
        'rpt_unit_seq','satellite','segment','serotype','serovar','sex',
        'specimen_voucher','standard_name','strain','sub_clone','sub_species',
        'sub_strain','tag_peptide','tissue_lib','tissue_type','transgenic',
        'translation','transl_except','transl_table','trans_splicing',
        'type_material','variety']

        keys_present = self.lst_of_dcts[0].keys()
        not_valid = [k for k in keys_present if k not in \
            valid_INSDC_quals]
        if not_valid:
            return MyException('The following are invalid INSDC qualifiers: '\
                '`%s`' % (not_valid))
        return True

#############
# FUNCTIONS #
#############

########
# MAIN #
########
