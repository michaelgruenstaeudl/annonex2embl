# annonex2embl
Converts annotated DNA sequence alignments in NEXUS format to EMBL submission files.

### Example Usage via Bash
```
SCRIPT=scripts/annonex2embl.py
INF1=tests/data/input/Amaranths_ITS_subset_ALIGN.nex
INF2=tests/data/input/Amaranths_ITS_subset_META.csv
OUTF=Amaranths_ITS_subset_ALIGN.embl
EMAIL=mi.gruenstaeudl@gmail.com

python2 $SCRIPT -n $INF1 -c $INF2 -o $OUTF -e $EMAIL
```
