*annonex2embl*
===================

Converts an annotated DNA sequence alignment in NEXUS format to either 
(a) TSV spreadsheets for submission to [ENA](http://www.ebi.ac.uk/ena) via [Webin](https://www.ebi.ac.uk/ena/submit/sra/#home) checklist submissions, or
(b) ENA flat files for submission via Webin Entry Upload submissions, if no suitable checklist is available).


FILE PREPARATION
----------------
Annotations are specified via a SETS-BLOCK. Every gene and every exon charset must be accompanied by one CDS charset.
How do you add a SETS-BLOCK to a DATA-BLOCK in a NEXUS file?

###### Example for a SETS-BLOCK
```
BEGIN SETS;
CharSet trnK_intron = 0-928 2531-2813 2849-3152;
CharSet matK_gene = 929-2530;
CharSet matK_CDS = 929-2530;
CharSet trnK_exon = 2814-2848;
CharSet trnK_CDS = 2814-2848;
CharSet psbA_gene = 3153-3200;
CharSet psbA_CDS = 3153-3200;
END;
```


GENERAL USAGE
-------------

###### Example for EMBL mode

```
python2 $PWD/scripts/annonex2embl_CMD.py
-n examples/TestData1.nex
-c examples/TestData1.csv
-d "description of alignment"
-e your_email_address@gmail.com
-o examples/output_file.embl
```

###### Example for checklist mode

```
python2 $PWD/scripts/annonex2embl_CMD.py
-n examples/TestData1.nex
-c examples/TestData1.csv
-d "description of alignment"
-e your_email_address@gmail.com
-o examples/output_file.csv
--clmode TRUE
--cltype trnK_matK
```


TO DO
-----

###### 1. Add a function that does the following in order:
* reads and parses a bibtex file,
* extracts (a) the citation info and (b) the submitter references as required by EMBL, and 
* write the correctly formatted string-lines into the final file during post-processing.

###### 2. Add functions to read in a charset spec is forward or reverse and to adjust the info in the feature table.

###### 3. Add a function to compensate contraction of annotation due to identification of internal stop codons (see line 304 in Annonex2embl.py)
Since "CkOps.TranslCheck().transl_and_quality_of_transl()" shortens annotations to the first internal stop codon 
encountered, the subsequent intron or IGS needs to be extended towards 5' to compensate. This can be a general function without a priori info passed to it. The important aspect is that only the SUBSEQUENT feature (if it is an intron or an IGS!) can be extended; all other features cannot be extended and need to produce a warning.
Pseudocode:
```
Does a gap in the annotations exist?
  If yes, is the gap followed by an intron or an IGS?
    If yes, extend the intron or the IGS towards 5' to compensate.
    If no, print a warning.
  If no, continue without action.
```

###### 4. Implement various improvements of the checklist function
* 3.1. Have the CLMODE automatically add the colum names for the final checklists
* 3.2. Have the CLMODE automatically add non-mandatory qualifiers as separate column

###### 5. Write a GUI interface for input
* 4.1. The GUI should consist of just one Window, where all functions are immediately visible; the GUI should not have any dropdown-menus. In general, the simpler the interface, the better.

###### 6. Implement improvements of argparser (scripts/annonex2embl_CMD.py)
* 5.1. Currently, the "required" and "optional" parameters are not displayed currently when calling scripts/annonex2embl_CMD.py. It incorrectly says "optional parameters" for all.
* 5.2. Currently, --clmode requires "True" of "False" as parameters; how can I use it such that only the presence of --clmode indicates "True", whereas its abscence indicates "False"?


CHANGELOG
---------
###### Version 0.4.3 (2018.03.23)
* Improved formatting of Python code
* Checking if sequence names in NEX-file identical to sequence ids in csv-file
* Discarding charset_ranges that are empty
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
