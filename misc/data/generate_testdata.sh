#!/bin/bash

BASEDIR=/home/michael_science/git/michaelgruenstaeudl_annonex2embl/
SCRIPT=scripts/annonex2embl.py

## TEST 1
INF1=$BASEDIR/misc/data/input/Amaranths_ITS_subset_ALIGN.nex
INF2=$BASEDIR/misc/data/input/Amaranths_ITS_subset_META.csv
OUTF=Amaranths_ITS_subset_ALIGN.embl
TESTOUT=$BASEDIR/misc/data/output/Amaranths_ITS_subset_ALIGN.embl

python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -e mi.gruenstaeudl@gmail.com
#echo ""
#echo "Differences between current output and template:"
#diff $OUTF $TESTOUT
#rm $TESTOUT
#echo ""

## TEST 2
INF1=$BASEDIR/misc/data/input/Test1_ALIGNM.nex
INF2=$BASEDIR/misc/data/input/Test1_META.csv
OUTF=Test1_ALIGNM.embl
TESTOUT=$BASEDIR/misc/data/output/Test1_ALIGNM.embl

python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -e mi.gruenstaeudl@gmail.com
#echo ""
#echo "Differences between current output and template:"
#diff $OUTF $TESTOUT
#rm $TESTOUT
#echo ""
