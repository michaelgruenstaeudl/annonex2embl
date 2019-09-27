*annonex2embl*
==============

Converts an annotated DNA multi-sequence alignment (in [NEXUS](http://wiki.christophchamp.com/index.php?title=NEXUS_file_format) format) to an EMBL flatfile for submission to [ENA](http://www.ebi.ac.uk/ena) via the [Webin-CLI submission tool](https://ena-docs.readthedocs.io/en/latest/cli_05.html).


## INSTALLATION
```
python2 setup.py install  # Installation
python2 setup.py test     # Testing
```

## FILE PREPARATION
The annotations of a NEXUS file are specified via [SETS-block](http://hydrodictyon.eeb.uconn.edu/eebedia/index.php/Phylogenetics:_NEXUS_Format), which is located beneath a DATA-block and defines sets of characters in the DNA alignment. In such a SETS-block, every gene and every exon charset must be accompanied by one CDS charset. Other charsets can be defined unaccompanied.

#### Example of a complete SETS-BLOCK
```
BEGIN SETS;
CHARSET trnK_intron_forward = 0-928 2531-2813 2849-3152;
CHARSET matK_gene_forward = 929-2530;
CHARSET matK_CDS_forward = 929-2530;    # Accompanying the charset matK_gene
CHARSET trnK_exon_reverse = 2814-2848;
CHARSET trnK_CDS_reverse = 2814-2848;   # Accompanying the charset trnK_exon
CHARSET psbA_gene_forward = 3153-3200;
CHARSET psbA_CDS_forward = 3153-3200;   # Accompanying the charset psbA_gene
END;
```

#### Examples of DESCR variable
```
# ITS of nuclear ribosomal DNA
DESCR="18S rRNA gene (partial), ITS1, 5.8S rRNA gene, ITS2 and 28S rRNA gene (partial)"

# trnL-trnF intergenic spacer of plastid genome
DESCR="tRNA-Leu (trnL) gene, partial sequence; trnL-trnF intergenix spacer, complete sequence; and tRNA-Phe (trnF) gene"

# trnK-matK region of plastid genome
DESCR="tRNA-Lys (trnK) gene and intron, partial sequence; maturase K (matK) gene, complete cds; psbA gene, partial sequence"

# rpl16 intron of plastid genome
DESCR="rpl16 intron, partial sequence"
```

## USAGE
#### On Linux / MacOS
```
SCRPT=$PWD/scripts/annonex2embl_CMD.py
INPUT=examples/input/TestData1.nex
METAD=examples/input/TestData1.csv
OTPUT=examples/output/TestData1.embl
DESCR='description of alignment'  # Do not use double-quotes
EMAIL=your_email_here@yourmailserver.com
AUTHR='your name here'  # Do not use double-quotes
MNFTS=PRJEB00000
MNFTN=${DESCR//[^[:alnum:]]/_}

python3 $SCRPT -n $INPUT -c $METAD -d "$DESCR" -e $EMAIL -a "$AUTHR" -o $OTPUT --manifeststudy $MNFTS --manifestname $MNFTN  --productlookup True
```

#### On Windows
```
SET SCRPT=$PWD\scripts\annonex2embl_CMD.py
SET INPUT=examples\input\TestData1.nex
SET METAD=examples\input\TestData1.csv
SET OTPUT=examples\output\TestData1.embl
SET DESCR='description of alignment'
SET EMAIL=your_email_here@yourmailserver.com
SET AUTHR='your name here'
SET MNFTS=PRJEB00000
SET MNFTN=a_unique_description_here

python %SCRPT% -n %INPUT% -c %METAD% -d %DESCR% -e %EMAIL% -a %AUTHR% -o %OTPUT% --manifeststudy %MNFTS% --manifestname %MNFTN%  --productlookup True
```


## TO DO
* Replace section "POST-PROCESSING OF FILES" with code that reads in the file written to outp_handle and that edits the text string within Python (as opposed to calling 'sed' as done currently). That way, the code becomes platform independent. For example, the code can be loaded into a stringIO and edited via regular Python-string-functions (see function "write_SeqRecord" in IOOps.py as an example).
* After post-processing of the output file, write the as a gnu-zipped file. This is easy in Python:
```
import gzip
final_output = gzip.open('seqRecords.embl.gz', 'wb')
try:
    output.write(SeqRecords) ## See function 'write_SeqRecord' in IOOps.py for details
finally:
    output.close()
```
* Example files (./examples/input): Combine the alignments "fuzzy.nex" and "reverse.nex" into a single NEXUS file, while keeping the maximum sequence length of each sequence at 38 nucleotides (adjust the annotations to the new length accordingly)
* The taxonomy check (optional argument: --taxcheck) shall primarily be conducted against the ENA database, not the NCBI database. So far, the taxonomy check was conducted against the NCBI database (see the use of "Entrez.esearch" in the static function "_taxname_lookup" in ParsingOps.py). Instead, the taxonomy check shall be conducted as described in section "Fetch taxon by scientific name" on https://www.ebi.ac.uk/ena/browse/taxonomy-service. Hence, please comment out (don't replace!) the command "esearch_records=Entrez.esearch(db='taxonomy',term=query_term,retmax=retmax,retmod='xml')"
and replace it with something like:
"enaTaxonomy_records=urllib2.urlopen('http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/' + query_term).read()"
See this SO post for ideas: https://stackoverflow.com/questions/24124643/parse-xml-from-url-into-python-object  Probably you will also need to adjust the subsequent three or four commands of function "_taxname_lookup" to have the XML code parsed correctly.

<!---
NOT NECESSARY AT THIS POINT
* Currently, --taxcheck requires "True" of "False" as parameters; how can I use it such that only the presence of --taxcheck indicates "True", whereas its abscence indicates "False"?
* Implement improvements of argparser (scripts/annonex2embl_CMD.py): Currently, the "required" and "optional" parameters are not displayed when calling scripts/annonex2embl_CMD.py. It incorrectly says "optional parameters" for all.
* Add a function that (a) reads and parses a bibtex file, extracts the citation info as well as the submitter references as from that file, and write the correctly formatted string-lines into the EMBL output file during post-processing.
--->


## DEVELOPMENT
#### Testing for development
To run the unittests outside of 'python setup.py test':
```
python3 -m unittest discover -s tests -p "*_test.py"
```

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.
