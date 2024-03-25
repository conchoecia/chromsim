#!/bin/sh


## +----------------------------------------------------+
## | This script gives a very basic example of how to   |
## | simulate a run and plot the results.               |
## +----------------------------------------------------+

FILENAME="basic"
OUTDIR="./chromsim_example/$FILENAME" # change this line if you want the output to be put somewhere else
MAIN="../inversion_sim/main.py" # change this line if you want to run the script from somewhere else

echo
echo "executing $FILENAME.sh example script"
echo "====="

[ ! -f $MAIN ] && echo "main.py was not found. Make sure this script is executed from the repository's example_scripts directory, or change the path in the script." && exit 1

mkdir -p "$OUTDIR"

echo
echo "running simulation"
echo "====="

"$MAIN" -S -a 100 -b 95 -w 1 -o "$OUTDIR" -f "$FILENAME"

echo
echo "plotting results"
echo "====="

"$MAIN" -P -s "$OUTDIR/$FILENAME.inv" -o "$OUTDIR"
