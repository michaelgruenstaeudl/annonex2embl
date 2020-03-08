*annonex2embl*
==============

[![Build Status](https://travis-ci.com/michaelgruenstaeudl/annonex2embl.svg?branch=master)](https://travis-ci.com/michaelgruenstaeudl/annonex2embl)
[![PyPI status](https://img.shields.io/pypi/status/annonex2embl.svg)](https://pypi.python.org/pypi/annonex2embl/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/annonex2embl.svg)](https://pypi.python.org/pypi/annonex2embl/)
[![PyPI version shields.io](https://img.shields.io/pypi/v/annonex2embl.svg)](https://pypi.python.org/pypi/annonex2embl/)
[![PyPI license](https://img.shields.io/pypi/l/annonex2embl.svg)](https://pypi.python.org/pypi/annonex2embl/)

Converts an annotated DNA multi-sequence alignment (in [NEXUS](http://wiki.christophchamp.com/index.php?title=NEXUS_file_format) format) to an EMBL flatfile for submission to [ENA](http://www.ebi.ac.uk/ena) via the [Webin-CLI submission tool](https://ena-docs.readthedocs.io/en/latest/cli_05.html).


## INSTALLATION
To get the most recent stable version of *annonex2embl*, run:

    pip install annonex2embl

Or, alternatively, if you want to get the latest development version of *annonex2embl*, run:

    pip install git+https://github.com/michaelgruenstaeudl/annonex2embl.git


## INPUT, OUTPUT AND PREREQUISITES
* **Input**: an annotated DNA multiple sequence alignment in NEXUS format; and a comma-delimited (CSV) metadata table
* **Output**: a submission-ready, multi-record EMBL flatfile

#### Requirements / Input preparation
The annotations of a NEXUS file are specified via [SETS-block](http://hydrodictyon.eeb.uconn.edu/eebedia/index.php/Phylogenetics:_NEXUS_Format), which is located beneath a DATA-block and defines sets of characters in the DNA alignment. In such a SETS-block, every gene and every exon charset must be accompanied by one CDS charset. Other charsets can be defined unaccompanied.

#### Example of a complete SETS-BLOCK
```
BEGIN SETS;
CHARSET matK_gene_forward = 929-2530;
CHARSET matK_CDS_forward = 929-2530;
CHARSET trnK_intron_forward = 1-928 2531-2813;
END;
```

#### Examples of corresponding DESCR variable
```
DESCR="tRNA-Lys (trnK) intron, partial sequence; maturase K (matK) gene, complete sequence"
```

## EXAMPLE USAGE
#### On Linux / MacOS
```
SCRPT=$PWD/scripts/annonex2embl_launcher_CLI.py
INPUT=examples/input/TestData1.nex
METAD=examples/input/Metadata.csv
OTPUT=examples/temp/TestData1.embl
DESCR='description of alignment here'  # Do not use double-quotes
EMAIL=your_email_here@yourmailserver.com
AUTHR='your name here'  # Do not use double-quotes
MNFTS=PRJEB00000
MNFTD=${DESCR//[^[:alnum:]]/_}

python3 $SCRPT -n $INPUT -c $METAD -d "$DESCR" -e $EMAIL -a "$AUTHR" -o $OTPUT --productlookup --manifeststudy $MNFTS --manifestdescr $MNFTD --compress
```

#### On Windows
```
SET SCRPT=$PWD\scripts\annonex2embl_launcher_CLI.py
SET INPUT=examples\input\TestData1.nex
SET METAD=examples\input\Metadata.csv
SET OTPUT=examples\temp\TestData1.embl
SET DESCR='description of alignment here'
SET EMAIL=your_email_here@yourmailserver.com
SET AUTHR='your name here'
SET MNFTS=PRJEB00000
SET MNFTD=a_unique_description_here

python %SCRPT% -n %INPUT% -c %METAD% -d %DESCR% -e %EMAIL% -a %AUTHR% -o %OTPUT% --productlookup --manifeststudy %MNFTS% --manifestdescr %MNFTD% --compress
```

<!--
## TO DO
* Foo bar baz
* Foo bar baz
-->


<!--
## TESTING
    python3 -m unittest discover -s tests -p "*_test.py"
    python3 -m unittest discover -s tests -p "*_test.py" -v  # verbose version
    pytest  # on Linux only, if python-pytest installed via pip
-->

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.
