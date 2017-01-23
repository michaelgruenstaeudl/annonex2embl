# annonex2embl
Converts annotated DNA sequence alignments in NEXUS format to EMBL submission files.

#### CHANGELOG:

###### Version 0.2
* Leading and trailing ambiguities per sequence are removed while maintaining correct annotations.
* Modifiers without content are no longer included in the output.
* Specification of taxonomic division and sequence version has been implemented.
* A submission mode that masks the ID line and the AC line has been implemented.
