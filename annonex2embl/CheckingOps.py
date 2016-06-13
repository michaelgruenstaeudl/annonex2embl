#!/usr/bin/env python
'''
Custom operations to check annotations
'''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME
import GenerationOps as GnOps

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016 Michael Gruenstaeudl'
__info__ = 'nex2embl'
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

class AnnoCheck:
    ''' This class contains functions to evaluate the quality of an 
    annotation.
    
    Args:
        extract (obj):      a sequence object; example: Seq('ATGGAGTAA', 
                            IUPACAmbiguousDNA())
        loc_object (obj):   a location object; example: FeatureLocation(
                            ExactPosition(0), ExactPosition(8))
        feature_type (str): a string detailing the type of the feature;
                            example: "CDS"
        record_id (str):    a string deatiling the name of the sequence in 
                            question; example: "taxon_A"
        transl_table (int): an integer; example: 11 (for bacterial code)
    Returns:
        tupl.   The return consists of the translated sequence (a str) and the
                updated feature location (a location object); example: 
                (transl_out, feat_loc)
    Raises:
        ME.MyException
    '''

    def __init__(self, extract, feature, record_id, transl_table=11):
        self.extract = extract
        self.feature = feature
        self.record_id = record_id
        self.transl_table = transl_table

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
    def _adjust_feat_loc(loc_object, with_internalStop, without_internalStop):
        ''' An internal static function to adjust the feature location if an
        internal stop codon were present. '''
        if len(without_internalStop) > len(with_internalStop):
            start_pos = loc_object.start
            stop_pos = start_pos + (len(with_internalStop) * 3)
            loc_range = range(start_pos, stop_pos)
            feat_loc = GnOps.GenerateFeatLoc(loc_range).exact()
        if len(without_internalStop) == len(with_internalStop):
            feat_loc = loc_object
        return feat_loc

    def check(self):
        ''' This function performs checks on a coding region.
        
        Specifically, the function tries to translate the coding region
        (CDS) directly, using the internal checker "cds=True". If a
        direct translation fails, it confirms if the CDS starts with a 
        methionine. If the CDS does not start with a methionine, a 
        ValueError is raised. If the CDS does start with a methionine, 
        translations are conducted with and without regard to internal 
        stop codons. The shorter of the two translations is kept. The 
        feature location is adjusted, where necessary.
        
        Note:
            The asterisk indicating a stop codon is truncated under 
            _transl(to_stop=True) and must consequently be added again
            (see line 137).
        '''
        
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation

        try:
            transl_out = AnnoCheck._transl(self.extract,
                self.transl_table, cds=True)
            feat_loc = self.feature.location
        except:
            if not AnnoCheck._check_protein_start(self.extract, 
                self.transl_table):
                raise ME.MyException('Feature "%s" of sequence "%s" does not '\
                    'start with a Methionine (ATG).' % (self.feature.id,
                                                        self.record_id))
            else:
                try:
                    without_internalStop = AnnoCheck._transl(self.extract,
                        self.transl_table)
                    with_internalStop = AnnoCheck._transl(self.extract,
                        self.transl_table, to_stop=True)
                    transl_out = with_internalStop
                    feat_loc = AnnoCheck._adjust_feat_loc(self.feature.location, 
                        with_internalStop, without_internalStop)
                except:
                    raise ME.MyException('Translation of feature `%s` of '\
                        'sequence `%s` unsuccessful.' % (self.feature.id,
                                                         self.record_id))
        transl_out = transl_out + "*"
        return (transl_out, feat_loc)
    
    def for_unittest(self):
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation

        try:
            transl_out, feat_loc = AnnoCheck(self.extract, self.feature,
                self.record_id, self.transl_table).check()
            if isinstance(transl_out, Seq) and isinstance(feat_loc, 
                FeatureLocation):
                return True
            return False
        #except ValueError: # Keep 'ValueError'; don't replace with 'ME.MyException'
        #    return False
        except ME.MyException as e:
            raise e


class CheckCoord:
    ''' This class contains functions to coordinate different checks. '''
        
    def __init__(self):
        pass

    def quality_of_qualifiers(self, lst_of_dcts, label):
        ''' This function conducts a series of quality checks on the qualifiers 
        list (a list of dictionaries).
        
        First (label_present), it checks if a qualifier matrix (and, hence, each 
        entry) contains a column labelled with <seqname_col_label>.
        Second (valid_INSDC_quals), it checks if column names constitute valid 
        INSDC feature table qualifiers.
            
        Args:
            label (str):  a string; example: 'isolate'
            lst_of_dcts (list): a list of dictionaries; example: 
                                [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                                 {'isolate': 'taxon_B', 'country': 'Peru'}] 
        Returns:
            True, unless exception
        Raises:
            passed exception

        Examples:
            Example 1: # Label is among keys of list of dicts
                >>> label = 'isolate'
                >>> lst_of_dcts = [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                                   {'isolate': 'taxon_B', 'country': 'Peru'}]
                >>> CheckCoord().quality_of_qualifiers(lst_of_dcts, label)
                Out: True

            Example 2: # Label is NOT among keys of list of dicts
                >>> label = 'sequence_name'
                >>> lst_of_dcts = [{'isolate': 'taxon_A', 'country': 'Ecuador'},
                                   {'isolate': 'taxon_B', 'country': 'Peru'}]
                >>> CheckCoord().quality_of_qualifiers(lst_of_dcts, label)
                Out: ME.MyException: csv-file does not contain a column labelled `sequence_name`
    
        TODO:
            (i) Check if sequence_names are also in .nex-file
            (ii) Have all metadata conform to basic ASCII standards (not 
                extended ASCII)!
        '''

        qual_checks = MetaCheck(lst_of_dcts)
        try:
            qual_checks.label_present(label)
        except ME.MyException as e:
            raise e
        try:
            qual_checks.valid_INSDC_quals()
        except ME.MyException as e:
            raise e
        return True
    
    def transl_and_quality_of_transl(self, seq_record, feature, transl_table):
        ''' This function conducts a translation of a coding region and checks 
        the quality of said translation.
            
        Args:
            seq_record (obj):   foobar; example: 'foobar'
            feature (obj):      foobar; example: 'foobar'
            transl_table (int): 

        Returns:
            True, unless exception
        Raises:
            feature

        Examples:
            Example 1: # 

            Example 2: # 

        TODO:
            (i) Adjust code of AnnoCheck.check() so that the start position of a 
                subsequent feature is also adjusted.                
            (ii) Adjust the location position in code of AnnoCheck.check() so 
                that the stop codon is also included.
            (iii) SHOULD EXAMPLE 2 NOT RESULT IN A FEATURE LOCATION THAT ENDS 
                AT ExactPosition(5), I.E. AFTER THE STOP CODON ???
        '''

        extract = feature.extract(seq_record)
        try:
            transl, loc = AnnoCheck(extract.seq, feature, seq_record.id, 
                                     transl_table).check()
            #pdb.set_trace()
            feature.qualifiers["translation"] = transl
            feature.location = loc
        except ME.MyException as e:
            raise e
        return feature


class MetaCheck:
    ''' This class contains functions to evaluate the quality of metadata.
    
    Args:
        lst_of_dcts (list): a list of dictionaries; example: 
                            [{'foo': 'foobarqux', 'bar': 'foobarqux', 
                              'qux': 'foobarqux'}, {'foo': 'foobarbaz', 
                              'bar': 'foobarbaz', 'baz': 'foobarbaz'}]    
    Returns:
        none
    
    Raises:
        ME.MyException
    '''
    
    def __init__(self, lst_of_dcts):
        self.lst_of_dcts = lst_of_dcts
    
    def label_present(self, label):
        ''' This function checks if each (!) list of dictionary keys of a 
        list of dictionaries encompass the element <label> at least once.
        
        Examples:
       
            Example 1: # Evaluates the situation where the label is  
                         present in ALL key lists.
                >>> lst_of_dcts = [{'foo': 'foobar', 'bar': 'foobar', 
                'qux': 'foobar'}, {'foo': 'foobar', 'bar': 'foobar', 
                'baz': 'foobar'}]
                >>> label = 'foo'
                >>> MetaCheck(lst_of_dcts).label_present(label)
                Out: True
                
            Example 2: # Evaluates the situation where the label is not 
                         present in EACH key list.
                >>> lst_of_dcts = [{'foo': 'foobar', 'bar': 'foobar', 
                'qux': 'foobar'}, {'foo': 'foobar', 'bar': 'foobar', 
                'baz': 'foobar'}]
                >>> label = 'qux'
                >>> MetaCheck(lst_of_dcts).label_present(label)
                Out: ME.MyException: csv-file does not contain a column labelled `qux`
                            
            Example 3: # Evaluates the situation where the label is not 
                         present in ANY key list.
                >>> lst_of_dcts = [{'foo': 'foobar', 'bar': 'foobar', 
                'qux': 'foobar'}, {'foo': 'foobar', 'bar': 'foobar', 
                'baz': 'foobar'}]
                >>> label = 'norf'
                >>> MetaCheck(lst_of_dcts).label_present(label)
                Out: ME.MyException: csv-file does not contain a column labelled `norf`
        '''

        if not all(label in dct.keys() for dct in self.lst_of_dcts):
            raise ME.MyException('csv-file does not contain a column '\
                'labelled `%s`' % (label))
        return True
    
    def valid_INSDC_quals(self):
        ''' This function checks if every (!) dictionary key in a list of 
        dictionaries is a valid INSDC qualifier.
        
        Examples:
        
            Example 1: # No invalid qualifiers present
                >>> lst_of_dcts = [{'allele': 'foobar', 'altitude': 'foobar', 
                'anticodon': 'foobar'}, {'trans_splicing': 'foobar', 
                'type_material': 'foobar', 'variety': 'foobar'}]
                >>> MetaCheck(lst_of_dcts).valid_INSDC_quals()
                Out: True
                
            Example 2: # Invalid qualifiers are very much present.
                >>> lst_of_dcts = [{'allele': 'foobar', 
                'MyInvalidQual_1': 'foobar', 'anticodon': 'foobar'}, 
                {'MyInvalidQual_2': 'foobar', 'type_material': 'foobar', 
                'variety': 'foobar'}]
                >>> MetaCheck(lst_of_dcts).valid_INSDC_quals()
                Out: ME.MyException: The following are invalid INSDC 
                qualifiers: `MyInvalidQual_1, MyInvalidQual_2`       
        '''
        
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
        
        from itertools import chain
        keys_present = list(chain.from_iterable([dct.keys() for dct in 
            self.lst_of_dcts]))
        not_valid = [k for k in keys_present if k not in \
            valid_INSDC_quals]
        if not_valid:
            raise ME.MyException('The following are invalid INSDC qualifiers: '\
                '`%s`' % (', '.join(not_valid)))
        return True
