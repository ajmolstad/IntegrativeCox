#!/bin/bash

#SBATCH -o Model3/Results/Model3_%a.Rout
#SBATCH --array=1-500
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 96:00:00

module load R

R CMD BATCH --vanilla Model3_Main.R  Model3/Results/Model3_${SLURM_ARRAY_TASK_ID}.Rout
