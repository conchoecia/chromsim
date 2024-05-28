#!/bin/sh


## +----------------------------------------------------+
## | This script gives a very basic example of how to   |
## | simulate a run using data from a .rbh file.        |
## | Compares FWMs of ALGs B1 with B2 and A1b with B3   |
## | in /Hydra vulgaris/ and /Rhopilema esculentum/.    |
## +----------------------------------------------------+

FILENAME="from_rbh"
OUTDIR="./chromsim_example/$FILENAME" # change this line if you want the output to be put somewhere else
MAIN="../inversion_sim/main.py" # change this line if you want to run the script from somewhere else
RBH="./HVU_RES_xy_reciprocal_best_hits.coloredby_BCnS_LGs.plotted.rbh" # change this line if you want to run the script from somewhere else or use a different .rbh file

echo
echo "executing $FILENAME.sh example script"
echo "====="

[ ! -f $MAIN ] && echo "main.py was not found. Make sure this script is executed from the repository's example_scripts directory, or change the path in the script." && exit 1

mkdir -p "$OUTDIR"

echo
echo "running simulation"
echo "====="

"$MAIN" -S -r "$RBH" -A "B1" -B "B2" -c "HVU7" -O "HVU" -M -f "$FILENAME""_1" -o "$OUTDIR"
"$MAIN" -S -r "$RBH" -A "B1" -B "B2" -c "RES4" -O "RES" -M -f "$FILENAME""_2" -o "$OUTDIR"
"$MAIN" -S -r "$RBH" -A "A1b" -B "B3" -c "RES13" -O "RES" -M -f "$FILENAME""_3" -o "$OUTDIR"
"$MAIN" -S -r "$RBH" -A "A1b" -B "B3" -c "HVU8" -O "HVU" -M -f "$FILENAME""_4" -o "$OUTDIR"

echo
echo "plotting results"
echo "===="

"$MAIN" -P -s "$OUTDIR/$FILENAME""_1.inv" -o "$OUTDIR"
"$MAIN" -P -s "$OUTDIR/$FILENAME""_2.inv" -o "$OUTDIR"
"$MAIN" -P -s "$OUTDIR/$FILENAME""_3.inv" -o "$OUTDIR"
"$MAIN" -P -s "$OUTDIR/$FILENAME""_4.inv" -o "$OUTDIR"
