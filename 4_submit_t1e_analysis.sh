#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB
#SBATCH --array=100-40000:100

# exit pipeline if a step fails
set -e

source ~/.bashrc
conda activate

container=/storage/singularity/mixtures.sif

# Define the population array
pop_list=(AFR EAS NFE SAS) # population list

# Calculate the population index
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 1000 ))  # population index (0-3) for the current job

# Access and export the variables using the task ID
export pop=${pop_list[$pop_index]}                            # specific population for the current job
export end=$(( $SLURM_ARRAY_TASK_ID - (10000 * $pop_index) )) # batch number (1-10) for the current job

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop: $pop,  Batch: $end"

# a. Run RAREsim2 to create additional datasets for the type I error scenario
if [ "$end" -gt 1000 ]; then # the first 1,000 datasets were already created in Step 3

./4a_run_RAREsim2_t1e.sh > /dev/null 2>&1

fi

# b. Run the rare variant association methods for the type I error scenario
singularity exec "$container" Rscript ./4b_run_methods_t1e.R "$SLURM_ARRAY_TASK_ID"
