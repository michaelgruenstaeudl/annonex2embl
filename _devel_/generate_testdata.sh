#!/bin/bash

BASEDIR=/home/michael_science/git/michaelgruenstaeudl_annonex2embl/
SCRIPT=scripts/annonex2embl.py

## TEST 1
INF1=$BASEDIR/tests/data/input/Amaranths_ITS_subset_ALIGN.nex
INF2=$BASEDIR/tests/data/input/Amaranths_ITS_subset_META.csv
OUTF=$BASEDIR/tests/data/temp/Amaranths_ITS_subset_ALIGN.embl
TEMPLATE=$BASEDIR/tests/data/output/Amaranths_ITS_subset_ALIGN.embl

touch $TEMPLATE
python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -d "nuclear ribosomal internal transcribed spacer" -e m.gruenstaeudl@fu-berlin.de --topol linear --taxdiv PLN
echo ""
echo "Differences between current output and template:"
diff $OUTF $TEMPLATE
echo ""

## TEST 2
INF1=$BASEDIR/tests/data/input/TestData_1.nex
INF2=$BASEDIR/tests/data/input/TestData_1.csv
OUTF=$BASEDIR/tests/data/temp/TestData_1.embl
TEMPLATE=$BASEDIR/tests/data/output/TestData_1.embl

touch $TEMPLATE
python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -d "chloroplast trnR-atpA intergenic spacer" -e m.gruenstaeudl@fu-berlin.de --topol linear --taxdiv PLN
echo ""
echo "Differences between current output and template:"
diff $OUTF $TEMPLATE
echo ""

## TEST 3
INF1=$BASEDIR/tests/data/input/TestData_2.nex
INF2=$BASEDIR/tests/data/input/TestData_2.csv
OUTF=$BASEDIR/tests/data/temp/TestData_2.embl
TEMPLATE=$BASEDIR/tests/data/output/TestData_2.embl

touch $TEMPLATE
python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -d "chloroplast trnR-atpA intergenic spacer" -e m.gruenstaeudl@fu-berlin.de --topol linear --taxdiv PLN
echo ""
echo "Differences between current output and template:"
diff $OUTF $TEMPLATE
echo ""

## TEST 4
INF1=$BASEDIR/tests/data/input/Pyrus_trnR_atpA.nex
INF2=$BASEDIR/tests/data/input/Pyrus_trnR_atpA.csv
OUTF=$BASEDIR/tests/data/temp/Pyrus_trnR_atpA.embl
TEMPLATE=$BASEDIR/tests/data/output/Pyrus_trnR_atpA.embl

touch $TEMPLATE
python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -d "chloroplast trnR-atpA intergenic spacer" -e m.gruenstaeudl@fu-berlin.de --topol linear --taxdiv PLN --linemask True
echo ""
echo "Differences between current output and template:"
diff $OUTF $TEMPLATE
echo ""

## TEST 5
INF1=$BASEDIR/tests/data/input/Pyrus_trnK_matK.nex
INF2=$BASEDIR/tests/data/input/Pyrus_trnK_matK.csv
OUTF=$BASEDIR/tests/data/temp/Pyrus_trnK_matK.embl
TEMPLATE=$BASEDIR/tests/data/output/Pyrus_trnK_matK.embl

touch $TEMPLATE
python2 $BASEDIR/$SCRIPT -n $INF1 -c $INF2 -o $OUTF -d "chloroplast matK-trnK intergenic spacer" -e m.gruenstaeudl@fu-berlin.de --topol linear --taxdiv PLN --taxcheck True --clmode True --cltype trnK_matK
echo ""
echo "Differences between current output and template:"
diff $OUTF $TEMPLATE
echo ""
