#!/bin/bash

Nsim=10000 # number of individuals


# make initial result and dataset (sub)directories

mkdir -p ./results
mkdir -p ./datasets
mkdir -p ./datasets/Hapgen$((Nsim/1000))K
mkdir -p ./datasets/Hapgen$((Nsim/1000))K_pruned
mkdir -p ./datasets/Cases
mkdir -p ./datasets/Controls

pops=(AFR EAS NFE SAS) # population list

# make result and dataset subdirectories for each population

for pop in "${pops[@]}"; do

mkdir -p ./results/${pop}
mkdir -p ./datasets/Hapgen$((Nsim/1000))K/${pop}
mkdir -p ./datasets/Hapgen$((Nsim/1000))K_pruned/${pop}
mkdir -p ./datasets/Cases/${pop}
mkdir -p ./datasets/Controls/${pop}

# make result and dataset subdirectories within the population directories for each batch
# (10 batches of 1,000 simulation replicates for 10,000 total replicates)

for n in {1..10}; do

mkdir -p ./datasets/Hapgen$((Nsim/1000))K/${pop}/Round${n}
mkdir -p ./datasets/Hapgen$((Nsim/1000))K_pruned/${pop}/Round${n}
mkdir -p ./datasets/Cases/${pop}/Round${n}
mkdir -p ./datasets/Controls/${pop}/Round${n}

done

done


