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
    def _adjust_feat_loc(location, with_internalStop, without_internalStop):
        ''' An internal static function to adjust the feature location if an
        internal stop codon were present. '''
        if len(without_internalStop) > len(with_internalStop):
            start_pos = location.start
            stop_pos = start_pos + (len(with_internalStop) * 3)
            feat_loc = GenerateFeatLoc(start_pos, stop_pos).exact()
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
            transl_out = AnnoChecks._transl(self.extract, self.transl_table,
                cds=True)
            feat_loc = self.feature.location
        except:
            if not AnnoChecks._check_protein_start(self.extract, 
                self.transl_table):
                raise MyException('Feature "%s" of sequence "%s" does not '\
                    'start with a Methionine (ATG).' % (self.feature.id,
                                                        self.record_id))
            else:
                try:
                    without_internalStop = AnnoChecks._transl(self.extract,
                        self.transl_table)
                    with_internalStop = AnnoChecks._transl(self.extract,
                        self.transl_table, to_stop=True)
                    transl_out = with_internalStop
                    feat_loc = AnnoChecks._adjust_feat_loc(self.feature.location, 
                        with_internalStop, without_internalStop)
                except:
                    raise MyException('Translation of feature `%s` of '\
                        'sequence `%s` unsuccessful.' % (self.feature.id,
                                                         self.record_id))
        transl_out = transl_out + "*"
        return (transl_out, feat_loc)
    
    def for_unittest(self):
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation

        try:
            transl_out, feat_loc = AnnoChecks(self.extract, self.feature,
                self.record_id, self.transl_table).check()
            if isinstance(transl_out, Seq) and isinstance(feat_loc, 
                FeatureLocation):
                return True
            return False
        #except ValueError: # Keep 'ValueError'; don't replace with 'MyException'
        #    return False
        except MyException as e:
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
                Out: MyException: csv-file does not contain a column labelled `sequence_name`
    
        TODO:
            (i) Check if sequence_names are also in .nex-file
            (ii) Have all metadata conform to basic ASCII standards (not 
                extended ASCII)!
        '''

        qual_checks = MetaChecks(lst_of_dcts)
        try:
            qual_checks.label_present(label)
        except MyException as e:
            raise e
        try:
            qual_checks.valid_INSDC_quals()
        except MyException as e:
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
            (i) Adjust code of AnnoChecks.check() so that the start position of a 
                subsequent feature is also adjusted.                
            (ii) Adjust the location position in code of AnnoChecks.check() so 
                that the stop codon is also included.
            (iii) SHOULD EXAMPLE 2 NOT RESULT IN A FEATURE LOCATION THAT ENDS 
                AT ExactPosition(5), I.E. AFTER THE STOP CODON ???
        '''

        extract = feature.extract(seq_record)
        try:
            transl, loc = AnnoChecks(extract.seq, feature, seq_record.id, 
                                     transl_table).check()
            feature.qualifiers["translation"] = transl
            feature.location = loc
        except MyException as e:
            raise e
        return feature


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


class GenerateFeatLoc:
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
                >>> GenerateFeatLoc(start_pos, stop_pos).exact()
                Out: FeatureLocation(ExactPosition(1), ExactPosition(12))
        '''
        from Bio import SeqFeature
        
        start_pos = SeqFeature.ExactPosition(self.start)
        end_pos = SeqFeature.ExactPosition(self.stop)
        return SeqFeature.FeatureLocation(start_pos, end_pos)


class GenerateSeqFeature:
    ''' This class contains functions to generate SeqFeatures. '''
        
    def __init__(self):
        pass
    
    def source_feat(self, feat_len, quals, transl_table):
        ''' This function generates the SeqFeature `source` for a SeqRecord.

        The SeqFeature `source` is critical for submissions to EMBL or GenBank, 
        as it contains all the relevant info on collection locality, herbarium 
        voucher, etc. It also provides info on which translation table is used
        for subsequent CDS features.
            
        Args:
            feat_len (int): length of the feature; example: 500
            quals (dict):   a dictionary of qualifiers; example: 
                            {'isolate': 'taxon_B', 'country': 'Ecuador'}
            transl_table (int): an integer; example: 11 (for bacterial code)
        Returns:
            SeqFeature (obj):   A SeqFeature object
        Raises:
            [currently nothing]
            
        Examples:
            Example 1: # Default evaluation
                >>> feat_len = 500
                >>> quals = {'isolate': 'taxon_B', 'country': 'Ecuador'}
                >>> transl_table = 11
                >>> GenerateSeqFeature().source_feat(feat_len, quals, transl_table)
                Out: SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(500)), type='source')
        '''
        
        from Bio import SeqFeature
    
        feature_loc = GenerateFeatLoc(0, feat_len).exact()
        source_feature = SeqFeature.SeqFeature(feature_loc, id='source',
            type='source', qualifiers=quals)
        source_feature.qualifiers["transl_table"]=transl_table
        return source_feature
    
    def regular_feat(self, feat_name, feat_range):
        ''' This function generates a regular SeqFeature for a SeqRecord.
            
        Args:
            feat_name (str):  a string; example: 'psbI_CDS'
            feat_range (list): a list of indices; example: [2, 3, 4, 5]
        Returns:
            SeqFeature (obj):   A SeqFeature object
        Raises:
            -
    
        TODO: 
            (i) Include a greater number of possible feature location functions.
            #start_pos = SeqFeature.AfterPosition(feat_range[0])
            #end_pos = SeqFeature.BeforePosition(feat_range[-1])
            (ii) Automatically identify a SeqFeature (e.g. search for the type in
            a database)
            
        Examples:
            Example 1: # Default evaluation
                >>> feat_name = 'psbI_CDS'
                >>> feat_range = [2, 3, 4, 5]
                >>> GenerateSeqFeature().regular_feat(feat_name, feat_range)
                Out: SeqFeature(FeatureLocation(ExactPosition(2), ExactPosition(6)), type='cds')
        '''
        from Bio import SeqFeature
        
        # a. Define the locations of the charsets
        start_pos = feat_range[0]
        stop_pos = feat_range[-1]+1
        feature_loc = GenerateFeatLoc(start_pos, stop_pos).exact()
        # b. Define the annotation type
        anno_types = ['cds', 'gene', 'rrna', 'trna']
        kw_present = [kw for kw in anno_types if kw in feat_name.lower()]
        if kw_present:
            feat_type = kw_present[0]
        else:
            feat_type = 'misc_feature'
        seq_feature = SeqFeature.SeqFeature(feature_loc, id=feat_name,
            type=feat_type, qualifiers={'note':feat_name})
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
        ''' This function checks if each (!) list of dictionary keys of a 
        list of dictionaries encompass the element <label> at least once.
        
        Examples:
       
            Example 1: # Evaluates the situation where the label is  
                         present in ALL key lists.
                >>> lst_of_dcts = [{'foo': 'foobar', 'bar': 'foobar', 
                'qux': 'foobar'}, {'foo': 'foobar', 'bar': 'foobar', 
                'baz': 'foobar'}]
                >>> label = 'foo'
                >>> MetaChecks(lst_of_dcts).label_present(label)
                Out: True
                
            Example 2: # Evaluates the situation where the label is not 
                         present in EACH key list.
                >>> lst_of_dcts = [{'foo': 'foobar', 'bar': 'foobar', 
                'qux': 'foobar'}, {'foo': 'foobar', 'bar': 'foobar', 
                'baz': 'foobar'}]
                >>> label = 'qux'
                >>> MetaChecks(lst_of_dcts).label_present(label)
                Out: MyException: csv-file does not contain a column labelled `qux`
                            
            Example 3: # Evaluates the situation where the label is not 
                         present in ANY key list.
                >>> lst_of_dcts = [{'foo': 'foobar', 'bar': 'foobar', 
                'qux': 'foobar'}, {'foo': 'foobar', 'bar': 'foobar', 
                'baz': 'foobar'}]
                >>> label = 'norf'
                >>> MetaChecks(lst_of_dcts).label_present(label)
                Out: MyException: csv-file does not contain a column labelled `norf`
        '''

        if not all(label in dct.keys() for dct in self.lst_of_dcts):
            raise MyException('csv-file does not contain a column '\
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
                >>> MetaChecks(lst_of_dcts).valid_INSDC_quals()
                Out: True
                
            Example 2: # Invalid qualifiers are very much present.
                >>> lst_of_dcts = [{'allele': 'foobar', 
                'MyInvalidQual_1': 'foobar', 'anticodon': 'foobar'}, 
                {'MyInvalidQual_2': 'foobar', 'type_material': 'foobar', 
                'variety': 'foobar'}]
                >>> MetaChecks(lst_of_dcts).valid_INSDC_quals()
                Out: MyException: The following are invalid INSDC 
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
            raise MyException('The following are invalid INSDC qualifiers: '\
                '`%s`' % (', '.join(not_valid)))
        return True


#############
# FUNCTIONS #
#############

########
# MAIN #
########
