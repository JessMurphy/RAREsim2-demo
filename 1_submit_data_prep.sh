#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB
#SBATCH --array=100-40000:100

# exit pipeline if a step fails
set -e

container=/storage/singularity/mixtures.sif

# calculate necessary variables (population & batch) using the task array ID
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 10000 ))    # population index (0-3) for the current job
end=$(( $SLURM_ARRAY_TASK_ID - (10000 * $pop_index) )) # batch number (1-10) for the current job

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop Index: $pop_index, Batch: $end"

# subset the master legend file for each simulation replicate
singularity exec "$container" Rscript ./1c_subset_master_legend.R "$SLURM_ARRAY_TASK_ID"