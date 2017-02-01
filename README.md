# annonex2embl
### (future title: nex2enaflat)
Converts annotated DNA sequence alignments in NEXUS format to ENA flat files for submissions to [ENA] (http://www.ebi.ac.uk/ena) via Entry Upload (i.e., if no suitable ENA checklist is available).

#### CHANGELOG:
###### Version 0.2 (2017.01.22)
* Leading and trailing ambiguities per sequence are removed while maintaining correct annotations.
* Modifiers without content are no longer included in the output.
* Specification of taxonomic division and sequence version has been implemented.
* A submission mode that masks the ID line and the AC line has been implemented.

###### Version 0.3 (2017.01.25)
* Rudimentary checklist output has been implemented.
* Taxon names are compared to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) and corrected if not listed there.

###### Version 0.4 (2017.02.01)
* All qualifier values are formatted to consist of ASCII characters only.
* If a coding region is among the sequence features, the qualifier /\trans_table/ is added to the source feature.
