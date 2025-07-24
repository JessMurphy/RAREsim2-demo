# Get command-line arguments
args = commandArgs(trailingOnly=TRUE)

id = as.numeric(args[[1]])
pop.index = ceiling(id/10000)
pop.list = c("AFR", "EAS", "NFE", "SAS")
pop = pop.list[pop.index]
end = id - 10000*(pop.index-1)
start = end - 99
Nsim = 10000

# block number on chromosome 19
b = 37

# Print confirmation
print(paste0("Running simulation with: Pop = ", pop, ", Batch = ", end))

library(dplyr) #lib.loc="/home/math/murphjes/R_libs/"

# read in reference legend file
leg.ref = read.table(paste0("./input/1000G/", pop, "_Block", b, "_CDS_ref_added.legend"), header = TRUE) 
  
# read in master legend file
master = read.table(paste0("./input/chr19.block", b, ".", pop, ".master.legend"), sep='\t')

colnames(master) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
master$alleles = paste0(master$a0, '/', master$a1)

# subset the file according to the number of alleles at each position
singles = master %>% filter(prob==1)
dups = master %>% filter(prob==0.5)
trips = master %>% filter(prob==".")

for (i in start:end){

  # determine the round number from i
  n = ceiling(i/1000)

  # get a list of all the  positions in trips
  trips.po = levels(droplevels(as.factor(trips$position)))

  # create a table of transitions/transversions
  trips.po1  =  as.data.frame(matrix(NA, nrow=length(trips.po), ncol=3))
  colnames(trips.po1) = c('position', 'draw', 'TiTv')
  trips.po1$position = trips.po
  trips.po1$draw = runif(nrow(trips.po1))
  trips.po1$TiTv[which(trips.po1$draw < 0.7396)] = 'transition'
  trips.po1$TiTv[which(trips.po1$draw >= 0.7396)] = 'transversion'
  #head(trips.po1)

  # subset the transitions
  ti = trips.po1 %>% filter(TiTv  == 'transition')
  ann.ti = trips %>% filter(position %in% ti$position, alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C')

  # remove all of the transitions
  ann.tv = trips %>% filter(!(position %in% ti$position), !(alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C'))

  # randomly pick an allele from the transversions
  ann.tv2 = ann.tv %>% group_by(position) %>% sample_n(size=1)

  # merge transitions and transversions
  trips2 = union(ann.ti, ann.tv2)

  # randomly pick an allele from the duplicates
  dups2 = dups %>% group_by(position) %>% sample_n(size=1)

  # merge all positions
  master2 = union(singles, union(dups2, trips2)) %>% select(id, position, a0, a1, AC, prob, exonic, gene, fun)

  master2$fun = ifelse(master2$fun=="synonymous SNV", "syn", "fun")
  master2$exonic[grepl("exonic", master2$exonic)] = "exonic"
  master2$gene[grepl("ZNF333", master2$gene)] = "ZNF333"
  
  # reorder the positions so they match the original reference legend file
  master3 = master2[match(leg.ref$position, master2$position),]
  
  # check that the positions were ordered correctly
  #identical(master3$position, leg.ref$position)

  # write output legend file
  out.name = paste0("./datasets/Hapgen10K/", pop, "/Round", n, "/chr19.block", b, ".", pop, ".sim", i, ".legend")
  write.table(master3, out.name, row.names=F, col.names=T, quote=F, sep='\t')

  print(i)
}

###########################################################################
  
# make input files for RAREsim2's expected_variants.py function (only need to do once per population)

if (end==100){
  
# read in the number of variants target data
nvar = read.table("./input/gnomad/Block", b, "_fun_syn_num_var.txt", header=T) %>% rename(n=downsample, Pop=pop)
     
# divide by the region size for per_kb
reg_size = 19.029 
nvar2 = nvar %>% filter(Pop==tolower(pop)) %>% mutate(fun_per_kb=obs_fun/reg_size, syn_per_kb=obs_syn/reg_size) %>% select(-block, -Pop)
  
# save number of variants target data
write.table(nvar2 %>% select(-c(obs_syn, obs_fun)), paste0("./input/mac_bins/chr19_block", b, "_", pop, "_nvar_target_data.txt"), row.names=F, quote=F, sep='\t')
  
# minor allele frequencies for the target data
Ntar = last(nvar$n)
tar_maf1 = round(0.01*(2*Ntar))
tar_maf0.5 = round(tar_maf1/2)
tar_maf0.25 = round(tar_maf0.5/2)
  
# minor allele count bins for the target data
if (Ntar > 3500){
    
  mac_tar = data.frame(Lower = c(1, 2, 3, 6, 11, 21, tar_maf0.5+1),
                       Upper = c(1, 2, 5, 10, 20, tar_maf0.5, tar_maf1))
    
} else {
  mac_tar = data.frame(Lower = c(1, 2, 3, 6, tar_maf0.25+1, tar_maf0.5+1),
                       Upper = c(1, 2, 5, tar_maf0.25, tar_maf0.5, tar_maf1))
}
  
# subset legend file to variants observed in gnomad (target data)
leg = master2 %>% filter(AC!='.', AC!=0, exonic=="exonic") # Megan just filtered by exonic
leg$AC = as.numeric(leg$AC)
  
# subset target data by functional status
leg_fun = leg %>% filter(fun=="fun")
leg_syn = leg %>% filter(fun=="syn")
  
# count the number of variants within each MAC bin (target data)
mac_fun_tar = mac_syn_tar = mac_tar
mac_fun_tar$count = mac_syn_tar$count = 0
for (k in 1:nrow(mac_tar)){
  mac_fun_tar$count[k] = sum(between(leg_fun$AC, mac_tar$Lower[k], mac_tar$Upper[k]))
  mac_syn_tar$count[k] = sum(between(leg_syn$AC, mac_tar$Lower[k], mac_tar$Upper[k]))
}
  
# convert the counts to proportions
mac_fun_tar2 = mac_fun_tar %>% mutate(Prop = count/nrow(leg_fun)) %>% select(-count)
mac_syn_tar2 = mac_syn_tar %>% mutate(Prop = count/nrow(leg_syn)) %>% select(-count)
  
# save allele frequency spectrum (AFS) target data
write.table(mac_tar %>% mutate(fun_prop=mac_fun_tar2$Prop, syn_prop=mac_syn_tar2$Prop), 
            paste0('./input/mac_bins/chr19_block', b, '_', pop, '_AFS_target_data.txt'), row.names = F, quote = F, sep = '\t')
  
# minor allele frequencies for the simulated data
sim_maf1 = round(0.01*(2*Nsim))
sim_maf0.5 = round(sim_maf1/2)
sim_maf0.25 = round(sim_maf0.5/2)
  
# minor allele count bins for the simulated data
if (Nsim > 3500){
  mac_sim <- data.frame(Lower = c(1, 2, 3, 6, 11, 21, sim_maf0.5+1),
                        Upper = c(1, 2, 5, 10, 20, sim_maf0.5, sim_maf1))
    
} else {
  mac_sim = data.frame(Lower = c(1, 2, 3, 6, sim_maf0.25+1, sim_maf0.5+1),
                       Upper = c(1, 2, 5, sim_maf0.25, sim_maf0.5, sim_maf1))
}
  
# save the mac bins for the simulated data
write.table(mac_sim, paste0("./input/mac_bins/MAC_bins_", format(Nsim, scientific=F), ".txt"), row.names=F, quote=F, sep='\t')

}