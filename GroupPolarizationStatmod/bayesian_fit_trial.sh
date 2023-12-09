#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --output=bayes_fits.out
#SBATCH --partition=serc,normal

module load jags
module load R/4.2.0

R -e "source(\"experiments.R\"); singleBayesianFitTrial()"
