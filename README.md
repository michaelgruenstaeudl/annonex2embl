*annonex2embl*
==============

Converts an annotated DNA multi-sequence alignment (in NEXUS format) to an EMBL flatfile for submission to [ENA](http://www.ebi.ac.uk/ena) via the [commandline submission system](https://ena-docs.readthedocs.io/en/latest/cli_05.html) of Webin.


## INSTALLATION AND TESTING
#### Installation
```
python2 setup.py install
```
#### Testing
```
python2 setup.py test
```

## FILE PREPARATION
Annotations are specified via a SETS-BLOCK. Every gene and every exon charset must be accompanied by one CDS charset.
How do you add a SETS-BLOCK to a DATA-BLOCK in a NEXUS file?
#### Example for a SETS-BLOCK
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

## GENERAL USAGE
#### Example with supplied test data
```
SCRPT=$PWD/scripts/annonex2embl_CMD.py
INPUT=examples/DNA_Alignment.nex
METAD=examples/Metadata.csv
DESCR="description of alignment"
EMAIL=your_email_here@gmail.com
AUTHR="Your_name_here"

python2 $SCRPT -n $INPUT -c $METAD -o ${INP%.nex*}.embl -d $DESCR -e $EMAIL -a $AUTHR
```

## TO DO
#### 1. Add "/codon_start=1" in CDS feature, if start and stop position of feature is uncertain (i.e., <100..>200)
#### 2. Remove all sequences that consist only of Ns (or ?s).
#### 3. Add a function to compensate contraction of annotation due to identification of internal stop codons (see line 304 in Annonex2embl.py)
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
#### 4. Add a function that does the following in order:
* reads and parses a bibtex file,
* extracts (a) the citation info and (b) the submitter references as required by EMBL, and
* write the correctly formatted string-lines into the final file during post-processing.
#### 5. Add functions to read in a charset spec is forward or reverse and to adjust the info in the feature table.
#### 6.
* The accession number shall be removed from the AC line ("AC   AC0663; SV 1; ..." --> "AC   XXX; SV 1; ...")
* The accession number shall be removed from the ID line ("ID   AC0663;" --> "ID   XXX;")
#### 7. Implement improvements of argparser (scripts/annonex2embl_CMD.py)
* Currently, the "required" and "optional" parameters are not displayed currently when calling scripts/annonex2embl_CMD.py. It incorrectly says "optional parameters" for all.
* Currently, --taxcheck requires "True" of "False" as parameters; how can I use it such that only the presence of --taxcheck indicates "True", whereas its abscence indicates "False"?
#### 8. Write GUI with similar to GUI of EMBL2checklists
#### 9. Integrate conversion of lat_long data via Canadensys API


## DEVELOPMENT
#### Testing for development
To run the unittests outside of 'python setup.py test':
```
python -m unittest discover -s /home/michael_science/git/michaelgruenstaeudl_annonex2embl/tests -p "*_test.py"
```

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

