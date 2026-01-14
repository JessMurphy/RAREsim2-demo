#!/bin/bash

# define variables
nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number
ncase=5000


start=$(($end-99))

# loop through the simulation replicates
for i in $(eval echo "{$start..$end}")
do 

echo "Simulation Replicate: ${i}"

# determine n (batch number) from i (simulation replicate)
n=$(( (i - 1) / 1000 + 1 ))

prefix=${pop}/Round${n}/chr${num}.block${b}.${pop}.sim${i}

# prune functional and synonymous variants down to 100% functional / synonymous and remove the rows of zeros (use the z flag)
python3 -m raresim sim \
    -m ./datasets/Hapgen$((nsim/1000))K/${prefix}.controls.haps.gz \
    --functional_bins ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_100_fun.txt \
    --synonymous_bins ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_100_syn.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K/${prefix}.copy.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.100fun.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.100fun.100syn.haps.gz \
    -z

# extract the type I error cases
python3 -m raresim extract \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.100fun.100syn.haps.gz \
    -o ./datasets/Cases/${prefix}.${ncase}.t1e.cases.100fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# move the controls
mv ./datasets/Cases/${prefix}.${ncase}.t1e.cases.100fun.100syn.haps-remainder.gz ./datasets/Controls/${prefix}.${ncase}.controls.same.100fun.100syn.haps-remainder.gz

done


# type I error cases: chr19.block37.${pop}.sim${rep}.${ncase}.t1e.cases.100fun.100syn.haps-sample.gz
# controls (same): chr19.block37.${pop}.sim${rep}.${ncase}.controls.same.100fun.100syn.haps-remainder.gz