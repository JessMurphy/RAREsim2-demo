# Get command-line arguments
args = commandArgs(trailingOnly=TRUE)

id = as.numeric(args[[1]])
pop.index = ceiling(id/10000)
pop.list = c("AFR", "EAS", "NFE", "SAS")
pop = pop.list[pop.index]
print(paste0("Pop: ", pop))

# block number on chromosome 19
b = 37

# load libraries
library(dplyr)
library(tidyr)
library(DT, lib.loc="/home/math/murphjes/R/myLibs/")

# read in reference legend file (see original RAREsim paper for more details)
# original 1000G data downloaded from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html (1000GP_Phase3.tgz)
leg.ref = read.table(paste0("./input/1000G/", pop, "_Block", b, "_CDS_ref_added.legend"), header = TRUE) 
pos.1000G = leg.ref %>% filter(!grepl("Un_Known", id)) %>% select(position)
  
# write list of positions in legend file (necessary for annovar)
#write.table(leg.ref$position, paste0("./input/1000G/", pop, ".1000G.chr19.block", b, ".positions.txt"), 
#sep="\t", quote=F, row.names=F, col.names=F)
  
# code unknown positions with 0 for a0/a1 
leg.ref$a0 = ifelse(grepl("Un_Known", leg.ref$id), 0, leg.ref$a0)
leg.ref$a1 = ifelse(grepl("Un_Known", leg.ref$id), 0, leg.ref$a1)
  
# download gnomad data
# GnomAD (vcf): https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz 
# and https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
  
# filter gnomad data (command line) - update using min and max positions of leg.ref file
# vcftools --gzvcf ./input/gnomad/gnomad.exomes.r2.1.1.sites.19.vcf.bgz --out ./input/gnomad/gnomad.exomes.r2.1.1.sites.19.block37_NFE.vcf \
# --chr 19 --from-bp 14492336 --to-bp 14887568 --remove-indels --remove-filtered-all \
# --get-INFO AF_nfe --get-INFO AC_nfe --get-INFO AN_nfe --get-INFO nhomalt_nfe
  
# read in gnomad data
gnomad = read.table(paste0("./input/gnomad/gnomad.exomes.r2.1.1.sites.19.block", b, "_", pop, ".vcf.INFO"), sep='\t', header=T) %>% rename(position = POS)
names(gnomad)[5:8] = sapply(strsplit(names(gnomad)[5:8], "_", fixed=T), head, 1)
pos.gnomad = gnomad %>% filter(position %in% leg.ref$position) %>% distinct(position) # 4,447 positions
  
# merge reference legend with gnomad data 
combined = merge(leg.ref, gnomad, by="position", all.x=T) %>% arrange(position)
  
# subset the multiallelic SNVs
dups = combined %>% group_by(position) %>% filter(n()>1) 
  
# remove duplicates/triplicates from known multiallelic SNVs
dups.known = dups %>% filter(!grepl("Un_Known", id), a1==ALT) %>% mutate(prob = "1") 
  
# extract the positions of unknown multiallelic SNVs
dups.unknown = dups %>% filter(grepl("Un_Known", id)) 
dup.pos = levels(as.factor(dups.unknown$position)) 
  
out = c()
  
# loop through the unknown multiallelic SNVs to choose one if possible
for (i in dup.pos){
    
  temp = dups.unknown %>% filter(position==i)
    
  # choose the allele with the maximum allele count in the population
  result = temp %>% filter(AC==max(AC))
    
  # if only one max, prob=1; if two maxes, prob=0.5
  result$prob = ifelse(nrow(result)==1, "1", "0.5")
    
  # if three maxes, prob=. (will use transition/transversion probabilities)
  result$prob = ifelse(nrow(result)==3, ".", result$prob)
    
  out = rbind(out, result)
}
  
# subset the biallelic SNVs
singles = combined %>% group_by(position) %>% filter(n()==1) 
  
# set the probability of gnomad and known 1000G SNVs to 1, otherwise .
singles$prob = ifelse(is.na(singles$REF), ".", "1")
singles$prob = ifelse(!grepl("Un_Known", singles$id), "1", singles$prob)
  
# combine biallelic SNVS back with the known/unknown multiallelic SNVs
combined2 = union(singles, union(dups.known, out)) %>% arrange(position)
  
# ANNOVAR (website): https://annovar.openbioinformatics.org/en/latest/
# ANNOVAR (wiki): https://davetang.org/wiki2/index.php?title=ANNOVAR
  
# ANNOVAR steps (command line)
# tar -xzf ./input/annovar.latest.tar.gz
# ./input/annovar/annotate_variation.pl -downdb -buildver hg19 seq ./input/annovar/humandb/hg19_seq/
# ./input/annovar/convert2annovar.pl -format region -seqdir ./input/annovar/humandb/hg19_seq/ -out ./input/annovar/annovar.chr19.block37.txt chr19:14492336-14887568
# sed -i 's/\r//' ./1000G/1000G.chr19.block37.positions.txt
# fgrep -wf ./input/1000G/1000G.chr19.block37.positions.txt ./input/annovar/annovar.chr19.block37.txt > ./input/annovar/annovar.chr19.block37.filtered.txt

# read in ANNOVAR file
annovar = read.table(paste0("./input/annovar/annovar.chr19.block", b, ".filtered.txt"), sep="\t", header=F)
  
# unknown variants not in gnomad or 1000G
unknown = combined2 %>% filter(grepl("Un_Known", id), is.na(REF)) %>% distinct(position) # 14,537
  
# subset annovar to just the unknown alternate alleles
annovar.un = annovar %>% filter(V2 %in% unknown$position) 
  
# merge annovar variants with 1000G/gnomad variants
combined2 = merge(combined2, annovar.un, by.x="position", by.y="V2", all=T) %>% select(-V1, -V3)
  
# create master legend using the combined datasets
master = combined2
master$CHROM = "19"
  
# unknown positions
rows.un = which(grepl("Un_Known", master$id))
  
# loop through the unknown variants
for (i in rows.un){
    
  if (!is.na(master$REF[i])){
      
    # change a0/a1 from 1000G to the REF/ALT alleles from gnomad
    master[i, "a0"] = master[i, "REF"]
    master[i, "a1"] = master[i, "ALT"]
      
  } else if (!is.na(master$V4[i])){
      
    # change a0/a1 from 1000G to the REF/ALT alleles from annovar
    master[i, "a0"] = master[i, "V4"]
    master[i, "a1"] = master[i, "V5"]
  }
    
  # remove the "Un_Known" designation from the id
  id = substr(master[i, "id"], 1, nchar(master[i, "id"])-8)
    
  # rename the id based on the new REF/ALT alleles
  master[i, "id"] = paste0(id, master[i, "a0"], "_", master[i, "a1"])
}
  
master$AC = ifelse(is.na(master$AC), ".", master$AC)
  
# create file for annovar functional annotation
master2 = master %>% select(CHROM, START=position, END=position, REF=a0, ALT=a1, AC, prob)
  
#write.table(master2, paste0('./input/annovar/master.chr19.block', b, '.', pop, '.txt'), 
#row.names=F, col.names=F, quote=F, sep='\t') # necessary for annovar functional annotation
  
# ANNOVAR functional annotation steps (command line)
# ./input/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene ./input/annovar/humandb/
# ./input/annovar/annotate_variation.pl -geneanno -buildver hg19 ./input/annovar/master.chr19.block37.NFE.txt ./input/annovar/humandb/
  
# read in functional annotation files
anno = read.table(paste0("./input/annovar/master.chr19.block", b, ".", pop, ".txt.variant_function"), sep='\t') %>% 
  select(position2 = V4, InEx = V1, gene = V2) # all positions
anno.exo = read.table(paste0("./input/annovar/master.chr19.block", b, ".", pop, ".txt.exonic_variant_function"), sep='\t') %>% select(position2 = V5, fun = V2) # only exonic positions
  
# determine which positions are intronic
introns = which(anno$InEx=="intronic")
  
# add the functional annotation of the exons to the file with all positions
anno$fun = "."
anno[-introns, "fun"] = anno.exo$fun
  
# merge the functional annotations with the master legend file
leg.master = cbind(master, anno) %>% select(position, id, a0, a1, AC, prob, InEx, gene, fun)
write.table(leg.master, paste0('./input/chr19.block', b, '.', pop, '.master.legend'), row.names=F, col.names=F, quote=F, sep='\t')  