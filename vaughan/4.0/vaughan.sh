#!/bin/bash
#
#SBATCH --job-name=beast2_job
#SBATCH --output=output_slurm.%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-10

~/beast-2.7.6/bin/beast -working -overwrite -seed 1 ./4.$SLURM_ARRAY_TASK_ID/testVaughan.xml
