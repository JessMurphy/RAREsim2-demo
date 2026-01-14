#!/bin/bash

# define variables
nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number
ncase=5000 #number of cases (also number of controls)
same=(160 140 120) #pcase necessary for same direction of effect scenario
opp=(130 120 110) #pcase neccessary for 50%/50% split of positive/negative directions of effect


# CREATE EXPECTED MAC BINS (need target data for all populations)

# combine all of the pcase percentages into one variable
pcase=("${same[@]}" "${opp[@]}")

for p in "${!pcase[@]}"; do

# calculate the expected number of variants for each pcase percentage
python3 -m raresim calc \
    --mac ./input/mac_bins/MAC_bins_${nsim}.txt \
    -o ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}.txt \
    -N $((2*$ncase)) \
    --nvar_target_data ./input/mac_bins/${pop}/chr19_block37_${pop}_nvar_target_data.txt \
    --afs_target_data ./input/mac_bins/${pop}/chr19_block37_${pop}_AFS_target_data.txt \
    --w_fun $(awk "BEGIN { printf \"%.1f\", ${pcase[$p]} / 100 }") \
    --reg_size 19.029

# remove the output files for the synonymous variants (only want to vary the percentage of functional variants)
rm ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}_syn.txt

done

# calculate the expected number of variants for pcase=100
python3 -m raresim calc \
    --mac ./input/mac_bins/MAC_bins_${nsim}.txt \
    -o ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_100.txt \
    -N $((2*$ncase)) \
    --nvar_target_data ./input/mac_bins/${pop}/chr19_block37_${pop}_nvar_target_data.txt \
    --afs_target_data ./input/mac_bins/${pop}/chr19_block37_${pop}_AFS_target_data.txt \
    --reg_size 19.029



start=$(($end-99))

# loop through the simulation replicates
for i in $(eval echo "{$start..$end}")
do 

echo "Simulation Replicate: ${i}"

# determine n (batch number) from i (simulation replicate)
n=$(( (i - 1) / 1000 + 1 ))

prefix=${pop}/Round${n}/chr${num}.block${b}.${pop}.sim${i}


# PRUNE / SUBSET HAPLOTYPES FOR THE SAME DIRECTION OF EFFECT SCENARIO

# loop through the indices of the same percentages
for s in "${!same[@]}"; do

if [[ $s -eq 0 ]]; then 

echo "Pruning down to ${same[$s]}% functional and 100% synonymous variants..."

# prune functional and synonymous variants down to same% functional / 100% synonymous and remove the rows of zeros (use the z flag)
python3 -m raresim sim \
    -m ./datasets/Hapgen$((nsim/1000))K/${prefix}.controls.haps.gz \
    --functional_bins ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${same[$s]}_fun.txt \
    --synonymous_bins ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_100_syn.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K/${prefix}.copy.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[$s]}fun.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$s]}fun.100syn.haps.gz \
    -z

else

echo "Pruning down from ${same[$((s-1))]} to ${same[$s]}% functional variants..."

# prune functional and synonymous variants down to same% functional / 100% synonymous and remove the rows of zeros (use the z flag)
python3 -m raresim sim \
    -m ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$((s-1))]}fun.100syn.haps.gz \
    --f_only ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${same[$s]}_fun.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[0]}fun.100syn.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[$s]}fun.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$s]}fun.100syn.haps.gz \

fi

# extract the power cases for the same direction of effect scenario
python3 -m raresim extract \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$s]}fun.100syn.haps.gz \
    -o ./datasets/Cases/${prefix}.${ncase}.power.cases.same.${same[$s]}fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# remove unnecessary files
rm ./datasets/Cases/${prefix}.${ncase}.power.cases.same.${same[$s]}fun.100syn.haps-remainder.gz

done


# PRUNE / SUBSET HAPLOTYPES FOR THE TYPE I ERROR SCENARIO

echo "Pruning down from ${same[$s]}% to 100% functional variants..."

# prune functional variants down to 100% in the entire sample 
# (the legend file is not output when using the z flag, only the pruned variants file)
python3 -m raresim sim \
    -m ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$s]}fun.100syn.haps.gz \
    --f_only ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_100_fun.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[0]}fun.100syn.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.100fun.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.t1e.100fun.100syn.haps.gz

# extract the cases and controls for the type I error scenario
python3 -m raresim extract \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.t1e.100fun.100syn.haps.gz \
    -o ./datasets/Cases/${prefix}.${ncase}.t1e.cases.100fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# move the controls for the same direction of effect and type I error scenario
mv ./datasets/Cases/${prefix}.${ncase}.t1e.cases.100fun.100syn.haps-remainder.gz ./datasets/Controls/${prefix}.${ncase}.controls.same.100fun.100syn.haps-remainder.gz


# PRUNE / SUBSET HAPLOTYPES FOR THE OPPOSITE DIRECTIONS OF EFFECT SCENARIO

# loop through the indices of the percentages for the opposite directions of effect scenario
for o in "${!opp[@]}"; do

echo "Pruning down from ${same[$o]}% to ${opp[$o]}% functional variants..."

# prune functional variants down from same% to opp% in the entire sample
python3 -m raresim sim \
    -m ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$o]}fun.100syn.haps.gz \
    --f_only ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${opp[$o]}_fun.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[0]}fun.100syn.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp[$o]}fun.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.${opp[$o]}fun.100syn.haps.gz

# record the pruned variants (make sure to account for the header or lack thereof)
awk -F'\t' -v OFS='\t' 'NR==FNR {if (FNR > 1) ids[$1]=1; next} FNR==1 {print $0, "protected"; next} {print $0, ($1 in ids) ? 1 : 0}' \
./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp[$o]}fun.100syn.legend-pruned-variants \
./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[0]}fun.100syn.legend > ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp[$o]}fun.100syn.protected.legend 

# extract the power cases for the opposite directions of effect scenario
python3 -m raresim extract \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.${opp[$o]}fun.100syn.haps.gz \
    -o ./datasets/Cases/${prefix}.${ncase}.power.cases.opp.${opp[$o]}fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

echo "Pruning down from ${same[$o]}% to ${opp[$o]}% functional variants, excluding protected variants..."

# prune functional variants down from same% to opp% in the entire sample again excluding the previously pruned variants
python3 -m raresim sim \
    -m ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$o]}fun.100syn.haps.gz \
    --f_only ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${opp[$o]}_fun.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp[$o]}fun.100syn.protected.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.opp.${opp[$o]}fun.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.protected.${opp[$o]}fun.100syn.haps.gz \
    --keep_protected

# extract the controls for the opposite direction of effect scenario
python3 -m raresim extract \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.protected.${opp[$o]}fun.100syn.haps.gz \
    -o ./datasets/Controls/${prefix}.${ncase}.controls.opp.${opp[$o]}fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# remove unnecessary files
rm ./datasets/Controls/${prefix}.${ncase}.controls.opp.${opp[$o]}fun.100syn.haps-sample.gz

done

done


# power cases (same): chr19.block37.${pop}.sim${rep}.${ncase}.power.cases.same.${same[$s]}fun.100syn.haps-sample.gz
# power cases (opp): chr19.block37.${pop}.sim${rep}.${ncase}.power.cases.opp.${opp[$o]}fun.100syn.haps-sample.gz
# controls (same): chr19.block37.${pop}.sim${rep}.${ncase}.controls.same.100fun.100syn.haps-remainder.gz
# controls (opp): chr19.block37.${pop}.sim${rep}.${ncase}.controls.opp.${opp[$o]}fun.100syn.haps-remainder.gz
# type I error cases: chr19.block37.${pop}.sim${rep}.${ncase}.t1e.cases.100fun.100syn.haps-sample.gz

