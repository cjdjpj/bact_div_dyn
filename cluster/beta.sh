#!/bin/bash
#SBATCH --job-name=beta_coalescent            # Job name
#SBATCH --mem=1G                              # Memory request
#SBATCH --ntasks=1                            # Number of tasks
#SBATCH --cpus-per-task=1                     # Number of CPU cores per task

srun python beta_mut_explicit.py\
    -o "runs/r9"\
    -l 3000\
    -t 3\
    -n 500\
    -m 0.00005\
    -r 0.003\
