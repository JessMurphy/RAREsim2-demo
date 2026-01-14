#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=100-40000:100

# Define the arrays (population and effective population size)
pop_list=(AFR EAS NFE SAS)
NE_list=(17469 14269 11418 14269)

# Calculate the population index
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 10000 ))  # 0-3 for populations

# Subset the population-specific variables
export pop=${pop_list[$pop_index]}
export NE=${NE_list[$pop_index]}
export DL=${DL_list[$pop_index]}
export end=$((SLURM_ARRAY_TASK_ID - 10000 * pop_index))

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop: $pop, Batch: $end"

# a. Run Hapgen2 to create the initial haplotypes for each population with an over-abundance of rare variants
singularity exec "$container" ./2a_run_hapgen2.sh