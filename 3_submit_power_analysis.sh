#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --mem-per-cpu=4096MB
#SBATCH --array=100-4000:100

# exit pipeline if a step fails
set -e

source ~/.bashrc
conda activate

# Define the population array
pop_list=(AFR EAS NFE SAS)

# Calculate the population index
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 1000 ))  # 0-3 for populations

# Access and export the variables using the task ID
export pop=${pop_list[$pop_index]}
export end=$(( $SLURM_ARRAY_TASK_ID - (1000 * $pop_index) ))

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop: $pop,  Batch: $end"

# a. Run RAREsim2 to create datasets for the same direction of effect and opposite directions of effect (50%/50%) scenarios
./3a_run_RAREsim2_power.sh > /dev/null 2>&1

# b. Run RAREsim2 to create datasets for the opposite directions of effect scenario with unequal amounts of risk/protective variants
./3b_run_RAREsim2_power_unequal.sh > /dev/null 2>&1

# c. Run the rare variant association methods for all of the scenarios
singularity exec /storage/singularity/mixtures.sif Rscript ./3c_run_methods_same_power.R "$SLURM_ARRAY_TASK_ID"
singularity exec /storage/singularity/mixtures.sif Rscript ./3c_run_methods_opp_power.R "$SLURM_ARRAY_TASK_ID"
singularity exec /storage/singularity/mixtures.sif Rscript ./3c_run_methods_opp_power_unequal.R "$SLURM_ARRAY_TASK_ID"
