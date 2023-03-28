#!/bin/bash
#
#SBATCH --job-name=tmp
#SBATCH --partition=long
#SBATCH --time 11
#SBATCH --mem-per-cpu=400
#SBATCH --ntasks=400
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stevewr@uw.edu

module load R
mpirun -np 1 Rscript my_parallel_code.R