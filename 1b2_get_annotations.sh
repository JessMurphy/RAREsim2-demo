#!/bin/bash

b=37 # block number on chromosome 19
pop=NFE # population

# get the starting and ending positions of the block
start=$(sort -n ./input/positions/Block${b}_gencode_positions.txt | head -n 1)
end=$(sort -n ./input/positions/Block${b}_gencode_positions.txt | tail -n 1)

# ANNOVAR (website): https://annovar.openbioinformatics.org/en/latest/
# ANNOVAR (wiki): https://davetang.org/wiki2/index.php?title=ANNOVAR

# download ANNOVAR
tar -xzf ./input/annovar.latest.tar.gz

# download whole genome FASTA files
./input/annovar/annotate_variation.pl -downdb -buildver hg19 seq ./input/annovar/humandb/hg19_seq/

# annotate all base pairs in block 37 with reference/alternate alleles
./input/annovar/convert2annovar.pl -format region -seqdir ./input/annovar/humandb/hg19_seq/ -out ./input/annovar/annovar.chr19.block${b}.txt chr19:${start}-${end}

# convert line endings of 1000G positions file from POS to UNIX format (if necessary)
#sed -i 's/\r//' ./input/1000G/1000G.chr19.block${b}.positions.txt

# select 1000G positions from ANNOVAR annotation
fgrep -wf ./input/positions/Block${b}_gencode_positions.txt ./input/annovar/annovar.chr19.block${b}.txt > ./input/annovar/annovar.chr19.block${b}.filtered.txt

#############################################################

# download gene annotation files (refGene database)
./input/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene ./input/annovar/humandb/

# perform functional annotation
./input/annovar/annotate_variation.pl -geneanno -buildver hg19 ./input/annovar/master.chr19.block${b}.${pop}.txt ./input/annovar/humandb/ -dbtype refGene