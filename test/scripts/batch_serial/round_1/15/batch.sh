#!/bin/bash
#
#SBATCH --job-name=serial_7
#SBATCH --output=output_slurm.%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-20

~/beast-2.7.6/bin/beast -working -overwrite -seed 1 ./$SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_ID.xml
