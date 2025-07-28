#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB
#SBATCH --array=100-40000:100

# exit pipeline if a step fails
set -e

# calculate necessary variables
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 10000 ))  # 0-3 for populations
end=$(( $SLURM_ARRAY_TASK_ID - (10000 * $pop_index) ))

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop Index: $pop_index, Batch: $end"

TIME0=$(date +%s)

# a. Download and subset the data (only need to do this once overall)
#singularity exec "$container" ./1a_subset_data.sh "$SLURM_ARRAY_TASK_ID"

TIME1=$(date +%s)
#echo "It took $(($TIME1 - $TIME0)) seconds to download and subset the datasets for all of the populations."

# b. Make the reference haplotype file and master legend file (only need to do once per population)
#singularity exec "$container" Rscript ./1b1_make_master_legend.R "$SLURM_ARRAY_TASK_ID"

TIME2=$(date +%s)
#echo "It took $(($TIME2 - $TIME1)) seconds to make a master legend file for one population."

TIME3=$(date +%s)

# c. Subset the master legend file
singularity exec "$container" Rscript ./1c_subset_master_legend.R "$SLURM_ARRAY_TASK_ID"

TIME4=$(date +%s)
echo "It took $(($TIME4 - $TIME3)) seconds to make 100 legend files for one population."