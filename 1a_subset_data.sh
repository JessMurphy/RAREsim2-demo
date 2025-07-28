#!/bin/bash

b=37 # block number on chromosome 19

# 1000G data: 1000GP_Phase3.tgz from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html

# unpack the tarball for the chromosome 19 1000 Genomes data
tar -xvzf ./input/1000G/1000GP_Phase3_chr19.tar.gz -C ./input/1000G/

# add row numbers, filter to the specified block, and filter to SNPs
zcat ./input/1000G/1000GP_Phase3_chr19.legend.gz | \
awk -v OFS=" " 'NR==1 {print "rownum", $0; next} {print NR-1, $0}' | \
awk -v OFS=" " 'NR==FNR {keep[$1]; next} (FNR==1 || $3 in keep)' ./input/positions/Block${b}_gencode_positions.txt - | \
awk -v OFS=" " 'NR==1 || $6 ~ /SNP/' > ./input/1000G/1000GP_Phase3_chr19_Block${b}.legend

# subset the .sample file to the NFE, AFR, EAS, and SAS samples
awk -F' ' 'NR==1 || ($3 == "EUR" && $2 != "FIN")' ./input/1000G/1000GP_Phase3.sample > ./input/1000G/1000GP_Phase3_NFE.sample
awk -F' ' 'NR==1 || ($3 == "AFR" && $2 != "ACB" && $2 != "ASW")' ./input/1000G/1000GP_Phase3.sample > ./input/1000G/1000GP_Phase3_AFR.sample
awk -F' ' 'NR==1 || ($3 == "SAS")' ./input/1000G/1000GP_Phase3.sample > ./input/1000G/1000GP_Phase3_SAS.sample
awk -F' ' 'NR==1 || ($3 == "EAS")' ./input/1000G/1000GP_Phase3.sample > ./input/1000G/1000GP_Phase3_EAS.sample

# get the row numbers from the legend file to keep
tail -n +2 ./input/1000G/1000GP_Phase3_chr19_Block${b}.legend | awk '{print $1}' > ./input/1000G/1000GP_Phase3_chr19_Block${b}_rows_to_keep.txt

# subset the haplotypes to the specific block
zcat ./input/1000G/1000GP_Phase3_chr19.hap.gz | awk 'NR==FNR {keep[$1]; next} (FNR in keep)' ./input/1000G/1000GP_Phase3_chr19_Block${b}_rows_to_keep.txt - | gzip > ./input/1000G/1000GP_Phase3_chr19_Block${b}.hap.gz


#############################################################

# gnomAD data (vcf): https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz 
#                    https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
  
# unpack the tarball for the chromosome 19 1000 Genomes data
tar -xvzf ./input/gnomad/gnomad_exomes_v2.1.1_chr19.tar.gz -C ./input/gnomad/

# get the starting and ending positions of the block
start=$(sort -n ./input/positions/Block${b}_gencode_positions.txt | head -n 1)
end=$(sort -n ./input/positions/Block${b}_gencode_positions.txt | tail -n 1)

# extract data for each population
pops=(afr eas nfe sas)
for pop in "${pops[@]}"; do

# filter gnomad data to block, only keep variants that passed all filters, remove indels, and get additional information for the population
vcftools --gzvcf ./input/gnomad/gnomad.exomes.r2.1.1.sites.19.vcf.bgz --out ./input/gnomad/gnomad.exomes.r2.1.1.sites.19.block${b}_${pop}.vcf \
--chr 19 --from-bp ${start} --to-bp ${end} --remove-indels --remove-filtered-all \
--get-INFO AF_${pop} --get-INFO AC_${pop} --get-INFO AN_${pop} --get-INFO nhomalt_${pop}

done
