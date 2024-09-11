#!/bin/bash
#SBATCH --job-name ExactBench      # Set a name for your job.
#SBATCH --partition medium               # Slurm partition to use
#SBATCH --ntasks 1                      # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 1-00:00                  # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=2000              # Memory limit for each tasks (in MB)
#SBATCH --array=1-500
#SBATCH -o ExactBench.out               # File to which STDOUT will be written
#SBATCH -e ExactBench.err               # File to which STDERR will be written
#SBATCH --mail-type=ALL                 # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amcgov@uw.edu     # Email to which notifications will be sent

module load R

for i in $(seq 4)
do
	Rscript run_ExactBenchSims.R $SLURM_ARRAY_TASK_ID $i 2000 3 # data set number, chain number, number of iterations, simulation setting number
done
