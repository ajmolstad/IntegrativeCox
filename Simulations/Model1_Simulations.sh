#!/bin/bash

#SBATCH -o Model1/Results/Model1_%a.Rout
#SBATCH --array=1-600
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 96:00:00

module load R/3.6

R CMD BATCH --vanilla Model1_Main.R  Model1/Results/Model1_${SLURM_ARRAY_TASK_ID}.Rout
