#!/bin/bash
#
#SBATCH --job-name=100sample_validation
#SBATCH --output=output_slurm.%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-15

~/beast-2.7.6/bin/beast -working -overwrite -seed 1 ./$SLURM_ARRAY_TASK_ID/100sample.xml
