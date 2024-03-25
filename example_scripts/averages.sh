#!/bin/sh


## +------------------------------------------------------+
## | This script runs multiple simulations and produces   |
## | their averages.                                      |
## | Simulates 10 runs for each combination of            |
## | (A,B)=(100,95),(50,40),(500,600) and w=1,2,5,10.     |
## +------------------------------------------------------+


FILENAME="averages"
OUTDIR="./chromsim_example/$FILENAME" # change this line if you want the output to be put somewhere else
MAIN="../inversion_sim/main.py" # change this line if you want to run the script from somewhere else

echo
echo "executing $FILENAME.sh example script"
echo "====="

[ ! -f $MAIN ] && echo "main.py was not found. Make sure this script is executed from the repository's example_scripts directory, or change the path in the script." && exit 1

mkdir -p "$OUTDIR"

echo
echo "running simulations"
echo "====="

for i in {0..9}; do
    "$MAIN" -S -a 100 -b 95 -w 1 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {10..19}; do
    "$MAIN" -S -a 100 -b 95 -w 2 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {20..29}; do
    "$MAIN" -S -a 100 -b 95 -w 5 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {30..39}; do
    "$MAIN" -S -a 100 -b 95 -w 10 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {40..49}; do
    "$MAIN" -S -a 50 -b 40 -w 1 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {50..59}; do
    "$MAIN" -S -a 50 -b 40 -w 2 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {60..69}; do
    "$MAIN" -S -a 50 -b 40 -w 5 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {70..79}; do
    "$MAIN" -S -a 50 -b 40 -w 10 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {80..89}; do
    "$MAIN" -S -a 500 -b 600 -w 1 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {90..99}; do
    "$MAIN" -S -a 500 -b 600 -w 2 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {100..109}; do
    "$MAIN" -S -a 500 -b 600 -w 5 -o "$OUTDIR" -f "$FILENAME""_""$i"
done
for i in {110..119}; do
    "$MAIN" -S -a 500 -b 600 -w 10 -o "$OUTDIR" -f "$FILENAME""_""$i"
done

echo
echo "collecting averages"
echo "====="

"$MAIN" -m -o "$OUTDIR" -f "$FILENAME"

echo
echo "plotting averages"
echo "====="

"$MAIN" -p -o "$OUTDIR" -f "$FILENAME"
