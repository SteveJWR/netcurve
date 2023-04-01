#!/bin/bash
#SBATCH --job-name ls_geom              # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition medium             # Slurm partition to use
#SBATCH --ntasks 1                   # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --cpus-per-task 24
#SBATCH --time 8-00:00                # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=10             # Memory limit for each tasks (in MB)
#SBATCH -o myscript_%j.out            # File to which STDOUT will be written
#SBATCH -e myscript_%j.err            # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=stevewr@uw.edu # Email to which notifications will be sent
module load R #CVXR requires R version 4.0.0 on the cluster and a particular version of the C compiler
Rscript consistency_of_bootstrap_variance_cluster_script.R 