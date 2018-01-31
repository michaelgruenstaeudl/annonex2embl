#!/usr/bin/env python
'''
Custom operations to generate ENA checklists
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio
import MyExceptions as ME
import pdb

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2018.01.31.1600'

#############
# DEBUGGING #
#############

#pdb.set_trace()

###########
# CLASSES #
###########

class Writer:
    ''' This class writes a TSV spreadsheet for a submission via the 
    WEBIN checklist submission system.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def genomic_CDS(self, seq_record, counter, charset_sym, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via 
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            charset_sym (str)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''
        
        #ENTRYNUMBER
        entrynumber = str(counter+1) # enumerate counter starts counting at 0
        #ORGANISM_NAME
        organism_name = seq_record.name
        #ENV_SAMPLE
        env_sam = 'no'
        #GENE               # Symbol of the gene corresponding to a sequence region; example: RdRp, sigA, inv
        gene_symbol = "foo bar"
        #PRODUCT            # Name of the product associated with the feature; example: RNA dependent RNA polymerase, sigma factor A
        product_name = "foo bar"
        #TRANSLATION TABLE  # Translation table for this organism. Chose from a drop-down list; example: 1, 2, 3, 5, 11
        transl_table = "12345"

        # the gene
        the_gene = [f for f in seq_record.features \
            if f.type=='gene']
        try:
            the_gene = the_gene[0]
        except:
            try:
                the_gene = [f for f in seq_record.features \
                    if f.type=='CDS']
            except:
                raise ME.MyException('%s annonex2embl ERROR: Problem \
                    with `%s`. %s gene not found.' % ('\n', seq_name, 'The gene'))

        ## 5' CDS LOCATION and 5'_PARTIAL
            # 5' CDS LOCATION   # Start of the coding region relative to the submitted sequence. For a full length CDS this is the position of the first base of the start codon.
        fiveprime_cds = str(the_gene.location.start.position)
            # PARTIAL AT 5'? (yes/no)  # For an incomplete CDS with the start codon upstream of the submitted sequence.
        if type(the_gene.location.start) == Bio.SeqFeature.ExactPosition:
            fiveprime_partial = 'no'
        if type(the_gene.location.start) == Bio.SeqFeature.BeforePosition:
            fiveprime_partial = 'yes'
        ## 3' CDS LOCATION and 3'_PARTIAL
            # 3' CDS LOCATION # End of the coding region relative to the submitted sequence. For a full length CDS this is the position of the last base of the stop codon.
        threeprime_cds = str(the_gene.location.end.position)
            # PARTIAL AT 3'? (yes/no) # For an incomplete CDS with the stop codon downstream of the submitted sequence.
        if type(the_gene.location.end) == Bio.SeqFeature.ExactPosition:
            threeprime_partial = 'no'
        if type(the_gene.location.end) == Bio.SeqFeature.AfterPosition:
            threeprime_partial = 'yes'
        #READING FRAME  # Mandatory if your CDS is 5' partial as it defines the reading frame. Location of the first base of the first fully-encoded amino acid., Example: 1,2 or 3
        read_frame = "12345"

        qualifiers = seq_record.features[0].qualifiers # source feature is always first in list
        #ISOLATE
        try:
            isolate = qualifiers['isolate']
        except:
            isolate = ''
        #SPEC_VOUCH
        try:
            spec_vouch = qualifiers['specimen_voucher']
        except:
            spec_vouch = ''
        #LOCALITY
        try:
            country = qualifiers['country']
        except:
            country = ''
        #ECOTYPE
        try:
            ecotype = qualifiers['ecotype']
        except:
            ecotype = ''

        #SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    env_sam,
                    gene_symbol,
                    product_name,
                    transl_table,
                    fiveprime_cds,
                    threeprime_cds,
                    fiveprime_partial,
                    threeprime_partial,
                    read_frame,
                    isolate,
                    spec_vouch,
                    country,
                    ecotype,
                    sequence
                   ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)


    def trnK_matK(self, seq_record, counter, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via 
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''
        
        #ENTRYNUMBER
        entrynumber = str(counter+1) # enumerate counter starts counting at 0
        #ORGANISM_NAME
        organism_name = seq_record.name

        # trnK_intron
        trnK_intron = [f for f in seq_record.features \
            if f.id=='trnK' and f.type=='intron']
        try:
            trnK_intron = str(trnK_intron[0])
            trnK_intron_present = 'yes'
        except:
            trnK_intron_present = 'no'
        # matK
        matK_gene = [f for f in seq_record.features \
            if f.id=='matK' and f.type=='gene']
        try:
            matK_gene = matK_gene[0]
        except:
            try:
                matK_gene = [f for f in seq_record.features \
                    if f.id=='matK' and f.type=='CDS']
            except:
                raise ME.MyException('%s annonex2embl ERROR: Problem \
                    with `%s`. %s gene not found.' % ('\n', seq_name, 'matK'))

        ## 5'_CDS and 5'_PARTIAL
            # 5'_CDS: Start of the matK coding region relative to the submitted sequence. For a full length CDS this is the position of the first base of the start codon.
        fiveprime_cds = str(matK_gene.location.start.position)
            # 5'_PARTIAL: cds partial at 5'? (yes/no) For an incomplete CDS with the start codon upstream of the submitted sequence.
        if type(matK_gene.location.start) == Bio.SeqFeature.ExactPosition:
            fiveprime_partial = 'no'
        if type(matK_gene.location.start) == Bio.SeqFeature.BeforePosition:
            fiveprime_partial = 'yes'
        ## 3'_CDS and 3'_PARTIAL
            # 3'_CDS: End of the matK coding region relative to the submitted sequence. For a full length CDS this is the position of the last base of the stop codon.
        threeprime_cds = str(matK_gene.location.end.position)
            # 3'_PARTIAL: cds partial at 3'? (yes/no) For an incomplete CDS with the stop codon downstream of the submitted sequence.
        if type(matK_gene.location.end) == Bio.SeqFeature.ExactPosition:
            threeprime_partial = 'no'
        if type(matK_gene.location.end) == Bio.SeqFeature.AfterPosition:
            threeprime_partial = 'yes'
        
        qualifiers = seq_record.features[0].qualifiers # source feature is always first in list
        #ISOLATE
        try:
            isolate = qualifiers['isolate']
        except:
            isolate = ''
        #SPEC_VOUCH
        try:
            spec_vouch = qualifiers['specimen_voucher']
        except:
            spec_vouch = ''
        #LOCALITY
        try:
            country = qualifiers['country']
        except:
            country = ''
        #ECOTYPE
        try:
            ecotype = qualifiers['ecotype']
        except:
            ecotype = ''

        #SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    fiveprime_cds,
                    threeprime_cds,
                    fiveprime_partial,
                    threeprime_partial,
                    trnK_intron_present,
                    isolate,
                    spec_vouch,
                    country,
                    ecotype,
                    sequence
                   ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)


    def rRNA(self, seq_record, counter, charset_sym, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via 
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            charset_sym (str)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''
        
        #ENTRYNUMBER
        entrynumber = str(counter+1) # enumerate counter starts counting at 0
        #ORGANISM_NAME
        organism_name = seq_record.name
        #SEDIMENT
        sediment = charset_sym
        
        qualifiers = seq_record.features[0].qualifiers # source feature is always first in list
        #ISOLATE
        try:
            isolate = qualifiers['isolate']
        except:
            isolate = ''
        #ISOLATION_SOURCE
        try:
            isol_source = qualifiers['isolation_source']
        except:
            isol_source = ''
        #COUNTRY
        try:
            country = qualifiers['country']
        except:
            country = ''
        #ECOTYPE
        try:
            lat_lon = qualifiers['lat_lon']
        except:
            lat_lon = ''
        #COLLECTION_DATE
        try:
            collection_date = qualifiers['collection_date']
        except:
            collection_date = ''

        #SEQUENCE
        sequence = str(seq_record.seq)
            
        out_list = [entrynumber,
                    organism_name,
                    sediment,
                    isolate,
                    isol_source,
                    country,
                    lat_lon,
                    collection_date,
                    sequence
                   ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)


    def ITS(self, seq_record, counter, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via 
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''
        import Bio
        
        #ENTRYNUMBER
        entrynumber = str(counter+1) # enumerate counter starts counting at 0
            
        #ORGANISM_NAME
        organism_name = seq_record.name
       
        qualifiers = seq_record.features[0].qualifiers # source feature is always first in list
        #ISOLATE
        try:
            isolate = qualifiers['isolate']
        except:
            isolate = ''

        #ENV_SAMPLE
        env_sam = 'no'

        #COUNTRY
        try:
            country = qualifiers['country']
        except:
            country = ''

        #SPEC_VOUCH
        try:
            spec_vouch = qualifiers['specimen_voucher']
        except:
            spec_vouch = ''

        all_seqrec_features = [f.id for f in seq_record.features]
        # 18S
        if '18S' in all_seqrec_features:
            RNA_18S = 'partial'
        else:
            RNA_18S = 'no'
        # 26S
        if '26S' in all_seqrec_features:
            RNA_26S = 'partial'
        else:
            RNA_26S = 'no'
        # ITS1
        if '18S' in all_seqrec_features:
            ITS1_feat = 'complete'
        else:
            ITS1_feat = 'partial'
        # ITS2
        if '26S' in all_seqrec_features:
            ITS2_feat = 'complete'
        else:
            ITS2_feat = 'partial'
        # 58S
        if 'ITS1' in all_seqrec_features and 'ITS2' in all_seqrec_features:
            RNA_58S = 'complete'
        else:
            RNA_58S = 'partial'

        #SEQUENCE
        sequence = str(seq_record.seq)
            
        out_list = [entrynumber,
                    organism_name,
                    isolate,
                    env_sam,
                    country,
                    spec_vouch,
                    RNA_18S,
                    ITS1_feat,
                    RNA_58S,
                    ITS2_feat,
                    RNA_26S,
                    sequence
                   ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)


    def IGS(self, seq_record, counter, charset_sym, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via 
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            charset_sym (str)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''
        
        #ENTRYNUMBER
        entrynumber = str(counter+1) # enumerate counter starts counting at 0
        #ORGANISM_NAME
        organism_name = seq_record.name
        
        #ENV_SAMPLE
        env_sam = 'no'
        #GENE1
        gene1 = charset_sym[0:4]
        #G1PRESENT
        g1present = 'no' # TO BE IMPROVED
        #GENE2
        gene2 = charset_sym[4:8]
        #G2PRESENT
        g2present = 'no' # TO BE IMPROVED

        qualifiers = seq_record.features[0].qualifiers # source feature is always first in list
        #ISOLATE
        try:
            isolate = qualifiers['isolate']
        except:
            isolate = ''
        #SPEC_VOUCH
        try:
            spec_vouch = qualifiers['specimen_voucher']
        except:
            spec_vouch = ''
        #COUNTRY
        try:
            country = qualifiers['country']
        except:
            country = ''
        
        #SEQUENCE
        sequence = str(seq_record.seq)
        
        out_list = [entrynumber,
                    organism_name,
                    env_sam,
                    gene1,
                    g1present,
                    gene2,
                    g2present,
                    isolate,
                    spec_vouch,
                    country,
                    sequence
                   ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)
