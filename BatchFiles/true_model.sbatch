#!/bin/bash
#SBATCH --job-name cc_sims              # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition short             # Slurm partition to use
#SBATCH --ntasks 1                  # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-07:00                # Wall time limit in D-HH:MM
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000            # Memory limit for each tasks (in MB)
#SBATCH -o myscript_%j.out            # File to which STDOUT will be written
#SBATCH -e myscript_%j.err            # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=stevewr@uw.edu # Email to which notifications will be sent
module load R/4.0.0-foss-2020a
Rscript R/01c_full_model_estimation.R 
