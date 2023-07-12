#!/bin/bash

#SBATCH --job-name=inversion_sim_improved_convergence_array
#SBATCH --cpus-per-task=10
#SBATCH --mem=5000M
#SBATCH --time=0-02:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --array=1-10

config=./inversion_sim_array_args
Asize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
Bsize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

module load python3/3.11
python3 inversion_sim.py $Asize $Bsize
