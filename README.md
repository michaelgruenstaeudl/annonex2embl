*annonex2embl*
===================

Converts an annotated DNA sequence alignment in NEXUS format to an ENA flatfile for submission via an analysis XML (http://ena-docs.readthedocs.io/en/latest/prog_12.html#object-relationships).


INSTALLATION AND TESTING
------------------------
###### Installation
```
python2 setup.py install
```

###### Testing
```
python2 setup.py test
```


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

###### Example with supplied test data
```
SCRIPT=$PWD/scripts/annonex2embl_CMD.py
INF1=examples/TestData1.nex
INF2=examples/TestData1.csv
OUTF=examples/output_file.csv
DESCR="description of alignment"
EMAIL=mi.gruenstaeudl@gmail.com

python2 $SCRIPT -n $INF1 -c $INF2 -o $OUTF -d $DESCR -e $EMAIL
```


TO DO
-----

###### 1. Convert the code from Python2.7 to Python3.6.

###### 2. Add a function that does the following in order:
* reads and parses a bibtex file,
* extracts (a) the citation info and (b) the submitter references as required by EMBL, and 
* write the correctly formatted string-lines into the final file during post-processing.

###### 3. 
* The accession number shall be removed from the AC line ("AC   AC0663; SV 1; ..." --> "AC   XXX; SV 1; ...")
* The accession number shall be removed from the ID line ("ID   AC0663;" --> "ID   XXX;")

###### 4. Include additional functions 
* A function that converts missing sections of a sequence that are longer than 10 nucleotides into a "gap"-feature; that gap-feature has to have a mandatory qualifiers (/estimated_length=<integer>) and an optional qualifiers (/note="text").
Example error: 'Sequence contains a stretch of n characters between base 290 and 1.173 that is not represented with a "gap" feature (stretches of n greater than 0 gives a warning, greater than 10 gives an error). line: 2348-2373 of AC_trnLF_taxnamesIncl_adjusted_2018.04.03.1100.embl - AC_trnLF_taxnamesIncl_adjusted_2018.04.03.1100.embl'

What is needed is a function that return the start and stop of a poly-N strech in a string; see the following for a possible answer: https://stackoverflow.com/questions/25211905/determine-length-of-polypurine-tract

* Ensure that the sequence does not start and/or end with an 'n'
i.e., trim away any Ns from the start of the end of a sequence while adjusting all annotations

###### 5. Add functions to read in a charset spec is forward or reverse and to adjust the info in the feature table.

###### 6. Add a function to compensate contraction of annotation due to identification of internal stop codons (see line 304 in Annonex2embl.py)
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

###### 7. Implement improvements of argparser (scripts/annonex2embl_CMD.py)
* Currently, the "required" and "optional" parameters are not displayed currently when calling scripts/annonex2embl_CMD.py. It incorrectly says "optional parameters" for all.
* Currently, --taxcheck requires "True" of "False" as parameters; how can I use it such that only the presence of --taxcheck indicates "True", whereas its abscence indicates "False"?


DEVELOPMENT
-----------
###### Testing for development
To run the unittests outside of 'python setup.py test':
```
python -m unittest discover -s /home/michael_science/git/michaelgruenstaeudl_annonex2embl/tests -p "*_test.py"
```


CHANGELOG
---------
###### Version 0.4.4 (2018.03.29)
* Separate the annonex2embl from the embl2checklist function.
###### Version 0.4.3 (2018.03.23)
* Improved formatting of Python code
* Checking if sequence names in NEX-file identical to sequence ids in csv-file
* Discarding charset_ranges that are empty
* Source feature 'organelle' implemented
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
