#!/bin/bash
#
#SBATCH --job-name=batch_1
#SBATCH --output=output_slurm.%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-100

~/beast-2.7.6/bin/beast -working -overwrite -seed 1 ./$SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_ID.xml
