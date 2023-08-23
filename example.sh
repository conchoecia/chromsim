#!/bin/sh

if [ -z "$1" ]
then
    outdir=./
else
    outdir=$1
fi

echo simulation . . .
./inversion_sim/main.py -o $outdir -S -a 100 -b 200 -n example
echo plots . . .
./inversion_sim/main.py -o $outdir -P -s $outdir/example.inv
echo gif . . .
./inversion_sim/main.py -o $outdir -P -s $outdir/example.inv -G
