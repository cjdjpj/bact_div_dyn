#!/bin/bash
#SBATCH --job-name=beta_coalescent            # Job name
#SBATCH --mem=1G                              # Memory request
#SBATCH --ntasks=1                            # Number of tasks
#SBATCH --cpus-per-task=1                     # Number of CPU cores per task

srun python beta.py\
    -o "runs/r001"\
    -l 5000\
    -t 5\
    -n 1500\
    -m 0.0000006\
    -r 0.1\
    --continuous_genome\
