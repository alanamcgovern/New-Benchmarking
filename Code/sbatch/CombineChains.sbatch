#!/bin/bash
#SBATCH --job-name Combine      # Set a name for your job.
#SBATCH --partition short               # Slurm partition to use
#SBATCH --ntasks 1                      # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-00:10                  # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=1000              # Memory limit for each tasks (in MB)
#SBATCH --array=1-500
#SBATCH -o Combine.out               # File to which STDOUT will be written
#SBATCH -e Combine.err               # File to which STDERR will be written
#SBATCH --mail-type=ALL                 # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amcgov@uw.edu     # Email to which notifications will be sent

module load R

Rscript ExactBench_CombineChains.R $SLURM_ARRAY_TASK_ID 3 # data set number, simulation setting number
