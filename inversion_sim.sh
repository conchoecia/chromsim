#!/bin/bash

#SBATCH --job-name=inversion_sim
#SBATCH --cpus-per-task=10
#SBATCH --mem=10000M
#SBATCH --time=0-00:05:00
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

module load python3/3.11
python3 inversion_sim.py $1 $2
