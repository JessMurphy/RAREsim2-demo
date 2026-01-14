##!/bin/bash

# variables to be updated per simulation
Nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number

# input file paths
Hap=./input/1000G/1000G_chr19_block${b}_${pop}_ref.hap
Map=./input/1000G/genetic_map_chr${num}_combined_b37.txt

echo "Pop: $pop, Batch: $end"

# unzip the haplotype file initially
#if [ "$end" == "100" ]; then gunzip $Hap.gz; fi

start=$(($end-99)) # starting simulation replicate for the loop

for i in $(eval echo "{$start..$end}")
do

# determine n (batch number) from i (simulation replicate)
n=$(( (i - 1) / 1000 + 1 ))

# output file path
Output=./datasets/Hapgen$((Nsim/1000))K/${pop}/Round${n}/chr${num}.block${b}.${pop}.sim${i}

# save a copy of the legend file
cp $Output.legend $Output.copy.legend

# simulate with HAPGEN2
hapgen2 -h $Hap \
-m $Map \
-l $Output.copy.legend \
-o $Output.gz \
-n $Nsim 0 \
-dl $DL 1 1.0 1.0 \
-Ne $NE \
-no_gens_output

# remove unnecessary files produced by Hapgen
rm $Output.cases.*
rm $Output.controls.sample

echo $i

done

# rezip the haplotype file at the end
#if [ "$end" == "10000" ]; then gzip $Hap; fi

