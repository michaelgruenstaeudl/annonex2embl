*annonex2embl*
===================

Converts an annotated DNA sequence alignment in NEXUS format to either 
(a) TSV spreadsheets for submission to [ENA](http://www.ebi.ac.uk/ena) via [Webin](https://www.ebi.ac.uk/ena/submit/sra/#home) checklist submissions, or
(b) ENA flat files for submission via Webin Entry Upload submissions, if no suitable checklist is available).


EXAMPLE USAGE
-------------

###### Example for EMBL mode

```
python2 /path_to_annonex2embl/scripts/annonex2embl.py
-n alignment.nex
-c metadata_table.csv
-d "description of alignment"
-e your_email_address@gmail.com
-o output_embl_file.embl
```

###### Example for checklist mode

```
python2 /path_to_annonex2embl/scripts/annonex2embl.py
-n alignment.nex
-c metadata_table.csv
-d "description of alignment"
-e your_email_address@gmail.com
-o output_checklist_file.csv
--clmode TRUE
--cltype trnK_matK
```


TO DO
-----
* Check at start of program if all sequence names have matching identifier names
* Have the CLMODE automatically add the colum names for the final checklists
* Have the CLMODE automatically add non-mandatory qualifiers as separate column


CHANGELOG
---------
###### Version 0.4.3 (2018.03.23)
* Improved formatting of Python code
###### Version 0.4.2 (2017.02.01)
* All qualifier values are formatted to consist of ASCII characters only.
* If a coding region is among the sequence features, the qualifier /\trans_table/ is added to the source feature.
* Implementation of customization of DE line.
###### Version 0.3 (2017.01.25)
* Rudimentary checklist output has been implemented.
* Taxon names are compared to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) and converted to informal names if not listed there.
###### Version 0.2 (2017.01.22)
* Leading and trailing ambiguities per sequence are removed while maintaining correct annotations.
* Modifiers without content are no longer included in the output.
* Specification of taxonomic division and sequence version has been implemented.
* A submission mode that masks the ID line and the AC line has been implemented.
