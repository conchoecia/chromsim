#!/bin/bash

#SBATCH --job-name=inversion_sim_window
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000M
#SBATCH --time=0-00:40:00
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --array=1-10

module load python3/3.11
path/to/main.py 400 500 -w $1
