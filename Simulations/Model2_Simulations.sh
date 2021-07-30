#!/bin/bash

#SBATCH -o Model2/Results/Model2_%a.Rout
#SBATCH --array=494
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 96:00:00

module load R/3.6

R CMD BATCH --vanilla Model2_Main.R  Model2/Results/Model2_${SLURM_ARRAY_TASK_ID}.Rout
