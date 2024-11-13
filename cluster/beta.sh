#!/bin/bash
#SBATCH --job-name=beta_coalescent            # Job name
#SBATCH --mem=1G                              # Memory request
#SBATCH --ntasks=1                            # Number of tasks
#SBATCH --cpus-per-task=1                     # Number of CPU cores per task

srun python beta.py\
    --output "runs/r001"\
    --length 5000\
    --Ne 50000000\
    --track_length 5\
    --nsample 1500\
    --mu 0.0000000005\
    --r_m 0.1\
    --model "kingman"\
    --continuous_genome\
