# Get command-line arguments
args = commandArgs(trailingOnly=TRUE)

# define variables
id = as.numeric(args[[1]])               # array task ID for the current job
pop.list = c("AFR", "EAS", "NFE", "SAS") # population list
pop = pop.list[id]                       # specific population for the current job
print(paste0("Pop: ", pop))

# block number on chromosome 19
b = 37

# load libraries
library(dplyr)
library(tidyr)
#library(DT)
library(data.table)
library(stringr)

# REFORMAT REFERENCE HAPLOTYPE AND LEGEND FILES

# read in the gencode coding positions for the specific block
code.pos = read.table(paste0("./input/positions/Block", b, "_gencode_positions.txt")) %>%
  rename(position=V1)

# read in the original sample file (the one used to generate the hap file)
sample.ids.all = read.table("./input/1000G/1000GP_Phase3.sample", header = TRUE)

# read in the sample file filtered by population
sample.ids.pop = read.table(paste0("./input/1000G/1000GP_Phase3_", pop, ".sample"), header = TRUE)

# find the column indices in the hap file
# (each sample has 2 haplotype columns, so the positions are 2*i-1 and 2*i)
id.indices = which(sample.ids.all$ID %in% sample.ids.pop$ID)
hap.cols = sort(c(2*id.indices-1, 2*id.indices))

# read in the haplotype file filtered to the specific block
hap.all = fread(paste0("./input/1000G/1000GP_Phase3_chr19_Block", b, ".hap.gz"))

# subset the columns
hap.pop = hap.all %>% select(paste0("V", hap.cols))

# read in the legend file filtered to the specific block
leg = read.table(paste0("./input/1000G/1000GP_Phase3_chr19_Block", b, ".legend"), header=T) %>% 
  select(-rownum) %>% mutate(AC=rowSums(hap.pop), row=1:nrow(hap.pop))

# for multiallelic positions, choose the variant with the largest AC
leg2 = leg %>% group_by(position) %>% filter(n()==1 | AC==max(AC)) %>%
  slice(1) %>% ungroup()

removed = setdiff(leg$row, leg2$row)

# remove the duplicated variants from the haplotype file
hap.pop2 = hap.pop[-removed,]

# add zeros back in to the legend and haplotype files for the coding positions

# merge the legend file with the gencode positions
leg.coding = merge(leg2, code.pos, by="position", all.y=T)

# create an empty haplotype matrix
hap.ref = as.data.frame(matrix(0, nrow=nrow(leg.coding), ncol=ncol(hap.pop2)))

# add the 1000G haplotypes back in
hap.ref[which(!is.na(leg.coding$id)),] = hap.pop2

# track positions explicitly so haplotype rows can be filtered to match the legend
hap.ref$position = leg.coding$position

# reformat the legend file
leg.ref = leg.coding %>% 
  mutate(hap.AC = rowSums(hap.ref), AC = ifelse(is.na(AC), 0, AC),
         id = if_else(str_detect(id, "^rs[0-9]+:"), str_replace(id, "^rs[0-9]+:", "19:"), id),
         id = case_when(is.na(id) ~ paste0("19:", position, "_Un_Known"),
                        str_detect(id, "^19:[0-9]+:[A-Z]+:[A-Z]+$") ~ {
                          parts = str_split(id, ":", simplify = TRUE)
                          paste0(parts[, 1], ":", parts[, 2], "_", parts[, 3], "_", parts[, 4])
                          }, TRUE ~ id)) %>%
  select(-c(TYPE:hap.AC))

# double check they were merged successfully
#which(leg.ref$AC!=leg.ref$hap.AC) #0

# CREATE MASTER LEGEND FILE

# subset the known 1000G positions
pos.1000G = leg.ref %>% filter(!grepl("Un_Known", id)) %>% select(position)

# read in gnomad data
gnomad = read.table(paste0("./input/gnomad/gnomad.exomes.r2.1.1.sites.19.block", b, "_", pop, ".vcf.INFO"), sep='\t', header=T) %>% rename(position = POS)
names(gnomad)[5:8] = sapply(strsplit(names(gnomad)[5:8], "_", fixed=T), head, 1)
pos.gnomad = gnomad %>% filter(position %in% leg.ref$position) %>% distinct(position) # 4,447 positions
  
# merge reference legend with gnomad data 
combined = merge(leg.ref, gnomad, by="position", all.x=T) %>% arrange(position)
  
# subset the multiallelic SNVs
dups = combined %>% group_by(position) %>% filter(n()>1) 
  
# remove duplicates/triplicates with the same alternate allele from known multiallelic SNVs
dups.known = dups %>% filter(!grepl("Un_Known", id), a1==ALT) %>% mutate(prob = "1") 

# remove duplicates/triplicates with different alternate alleles from known multiallelic SNVs
# (just use 1000G alternate allele)
dups.known2 = dups %>% filter(!grepl("Un_Known", id), !(position %in% dups.known$position)) %>%
  distinct(position, .keep_all=T) %>% mutate(prob = "1") 
  
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
combined2 = bind_rows(singles, dups.known, dups.known2, out) %>% arrange(position)
  
########## See 1b_get_annotations.sh file to run ANNOVAR ##########

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
  
########## See 1b_get_annotations.sh file again ##########

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

# remove intronic variants from the master legend
keep = is.na(leg.master$InEx) | leg.master$InEx != "intronic"
leg.master = leg.master[keep, ]

# keep reference haplotype rows for positions that remain in the legend
hap.ref = hap.ref %>% filter(position %in% unique(leg.master$position))

# order both outputs by position before writing
leg.master = leg.master %>% arrange(position)
hap.ref = hap.ref %>% arrange(position)

# save the filtered haplotypes for input into Hapgen2
write.table(hap.ref %>% select(-position), paste0("./input/1000G/1000G_chr19_block", b, "_", pop, "_ref.hap"), 
            quote=F, row.names=F, col.names=F)

write.table(leg.master, paste0('./input/chr19.block', b, '.', pop, '.master.legend'), row.names=F, col.names=F, quote=F, sep='\t')  
