#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB
#SBATCH --array=100-40000:100

# exit pipeline if a step fails
set -e

# calculate necessary variables (population & batch) using the task array ID
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 10000 ))    # population index (0-3) for the current job
end=$(( $SLURM_ARRAY_TASK_ID - (10000 * $pop_index) )) # batch number (1-10) for the current job

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop Index: $pop_index, Batch: $end"

# a. Download and subset the data (only need to do this once overall)
#singularity exec "$container" ./1a_subset_data.sh "$SLURM_ARRAY_TASK_ID"

# b. Make the reference haplotype file and master legend file (only need to do once per population)
#singularity exec "$container" Rscript ./1b1_make_master_legend.R "$SLURM_ARRAY_TASK_ID"

# c. Subset the master legend file for each simulation replicate
singularity exec "$container" Rscript ./1c_subset_master_legend.R "$SLURM_ARRAY_TASK_ID"