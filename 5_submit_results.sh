#!/bin/bash

#SBATCH --output=%x.out.%j
#SBATCH --error=%x.err.%j
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB

container=/storage/singularity/mixtures.sif

# a. Plot power and type I error results
singularity exec "$container" Rscript ./5a_plot_results.R

# b. Plot allele frequency distribution (AFS) info for the target (gnomad) and simulated data
singularity exec "$container" Rscript ./5b_plot_afs.R
