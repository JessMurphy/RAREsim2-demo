#!/bin/bash

#SBATCH --output=%x.out.%j
#SBATCH --error=%x.err.%j
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB

# 5. Plot results
singularity exec /storage/singularity/mixtures.sif Rscript ./5a_plot_results.R