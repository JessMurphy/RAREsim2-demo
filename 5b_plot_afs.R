
# load libraries
library(tidyverse)

setwd("C:/Users/murphjes/Documents/Repositories/RAREsim2_demo")

# define variables
pop = c("AFR", "EAS", "NFE", "SAS") # population list
N = c(8128, 9197, 56885, 15308)     # gnomAD effective population sizes
maf = 0.01                          # minor allele frequency threshold to define rare variants

# list of genes within the median centimorgan block of chr19 & their starting/ending positions
genes = c("ADGRE5", "DDX39A", "PKN1", "PTGER1", "GIPC1", "DNAJB1", "TECR", "NDUFB7", "CLEC17A", "ADGRE3", "ZNF333", "ADGRE2")
start = c(14491313, 14519631, 14543865, 14583278, 14588572, 14625582, 14627897, 14676890, 14693896, 14729929, 14800613, 14843205)
end = c(14519537, 14530192, 14582679, 14586174, 14606944, 14640582, 14676792, 14682874, 14721969, 14800839, 14844558, 14889353)

sum.afs.gnomad = sum.afs.all = list()
for (i in 1:4){ # loop through the four populations
  
  # read in the master legend file
  master = read.table(paste0("./input/chr19.block37.", pop[i], ".master.legend"), sep='\t')
  colnames(master) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
  master$alleles = paste0(master$a0, '/', master$a1)
  
  ######################## GNOMAD AFS ########################
  # extracted from the master legend file, but could also use a legend file already subset from the master (for a single sim rep)
  # lines 26-29,36-69,96-109 of 1c_subset_master_legend.R are used below (functional annotations are required)
  
  # subset the file according to the number of alleles at each position
  singles = master %>% filter(prob==1)
  dups = master %>% filter(prob==0.5)
  trips = master %>% filter(prob==".")
  
  # get a list of all the  positions in trips
  trips.po = levels(droplevels(as.factor(trips$position)))
  
  # create a table of transitions/transversions
  trips.po1  =  as.data.frame(matrix(NA, nrow=length(trips.po), ncol=3))
  colnames(trips.po1) = c('position', 'draw', 'TiTv')
  trips.po1$position = trips.po
  trips.po1$draw = runif(nrow(trips.po1))
  trips.po1$TiTv[which(trips.po1$draw < 0.7396)] = 'transition'
  trips.po1$TiTv[which(trips.po1$draw >= 0.7396)] = 'transversion'
  
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
  master2 = union(singles, union(dups2, trips2)) %>% arrange(position) %>% select(id, position, a0, a1, AC, prob, exonic, gene, fun)
  
  # reformat some of the columns
  master2$fun = ifelse(master2$fun=="synonymous SNV", "syn", "fun")
  master2$exonic[grepl("exonic", master2$exonic)] = "exonic"
  master2$gene[grepl("ZNF333", master2$gene)] = "ZNF333"
  
  # filter to exonic, polymorphic variants & categorize variants as rare/common based on their AC
  master2.filtered = master2 %>% filter(AC!=".") %>% mutate(AC=as.numeric(AC)) %>% filter(AC>0, exonic=="exonic") %>%
    mutate(common=ifelse(AC/(2*N[i]) > maf & AC/(2*N[i]) < 1-maf, "common", "rare"))
  
  # minor allele frequencies for the gnomad target data
  tar_maf1 = round(0.01*(2*N[i]))
  tar_maf0.5 = round(tar_maf1/2)
  tar_maf0.25 = round(tar_maf0.5/2)
  
  # minor allele count bins for the gnomad target data
  if (N[i] > 3500){
    
    mac_tar = data.frame(Lower = c(1, 2, 3, 6, 11, 21, tar_maf0.5+1),
                         Upper = c(1, 2, 5, 10, 20, tar_maf0.5, tar_maf1))
    
  } else {
    mac_tar = data.frame(Lower = c(1, 2, 3, 6, tar_maf0.25+1, tar_maf0.5+1),
                         Upper = c(1, 2, 5, tar_maf0.25, tar_maf0.5, tar_maf1))
  }
  
  # extract the gnomad bin levels
  gnomad.bin.levels = c("1", "2", paste0(mac_tar$Lower[3:7], "-", mac_tar$Upper[3:7]), paste0(">", mac_tar$Upper[7]))
  bin.labels=c("1", "2", "3-5", "6-10", "11-20", "21-0.5%", "0.5%-1%", ">1%")
  
  # assign each variant to a MAC bin
  master2.bins = master2.filtered %>% mutate(Bin=cut(AC, breaks=c(mac_tar$Lower, mac_tar$Upper[7]+1, Inf), labels=gnomad.bin.levels, right=F))
  
  # count the number of variants within each MAC bin per gene
  sum.afs.gnomad[[i]] = master2.bins %>% filter(common=="rare", fun=="fun") %>% 
    group_by(gene, Bin) %>% summarize(n=n(), .groups = "drop") %>% complete(gene, Bin, fill=list(n=0)) %>% 
    group_by(gene) %>% mutate(prop=n/sum(n)) %>% ungroup() %>% 
    mutate(Bin = factor(Bin, levels=gnomad.bin.levels, labels=bin.labels), pop=pop[i])
  
  ######################## SIMULATED AFS ########################
  
  # extract the simulated MAC bin levels
  mac.bins = read.table("./input/mac_bins/MAC_bins_10000.txt", header=T)
  sim.bin.levels = c("1", "2", paste0(mac.bins$Lower[3:7], "-", mac.bins$Upper[3:7]), paste0(">", mac.bins$Upper[7]))
  
  # read in the AFS info for each simulation replicate
  sum.afs = read.table(paste0("./test/sum_afs_all_", pop[i], ".txt"), header=T) %>% mutate(files="old")
  
  # format the simulated AFS data
  sum.afs.all[[i]] = sum.afs %>% group_by(gene, rep) %>% mutate(prop=n/sum(n)) %>% ungroup() %>%
    mutate(Bin = factor(Bin, levels=sim.bin.levels, labels=bin.labels), pop=pop[i])
}

# combine the dataframes together across the populations & calculate
# the difference in proportions between gnomad and the simulated data
sum.afs.gnomad2 = do.call(rbind, sum.afs.gnomad)
sum.afs.all2 = do.call(rbind, sum.afs.all) %>% group_by(rep) %>% mutate(diff=prop-sum.afs.gnomad2$prop) %>% ungroup()

### PLOTS ###

# TECR AFS plot (Figure S4)
TECR.afs.diff.plot = ggplot() +
  geom_boxplot(data=sum.afs.all2 %>% filter(gene=="TECR"), aes(x=Bin, y=diff), fill="#56B4E9", position=position_dodge(width=1)) +
  geom_hline(yintercept=0, linetype=2, linewidth=1.1) +
  facet_wrap(~pop, nrow=2, ncol=2, scales="free_x") + ylim(-0.55, 0.55) +
  labs(y='Difference in Proportion of Functional Variants (Simulated - Target)', x='MAC bin', title="TECR (small)") +
  theme_bw(base_size=17)

ggsave(file = "./test/prop_diff_TECR.jpg", plot = TECR.afs.diff.plot, height = 8, width = 14, units = 'in')

# ADGRE3 AFS plot (Figure S5)
ADGRE3.afs.diff.plot = ggplot() +
  geom_boxplot(data=sum.afs.all2 %>% filter(gene=="ADGRE3"), aes(x=Bin, y=diff), fill="#56B4E9", position=position_dodge(width=1)) +
  geom_hline(yintercept=0, linetype=2, linewidth=1.1) +
  facet_wrap(~pop, nrow=2, ncol=2, scales="free_x") + ylim(-0.55, 0.55) +
  labs(y='Difference in Proportion of Functional Variants (Simulated - Target)', x='MAC bin', title="ADGRE3 (medium)") +
  theme_bw(base_size=17)

ggsave(file = "./test/prop_diff_ADGRE3.jpg", plot = ADGRE3.afs.diff.plot, height = 8, width = 14, units = 'in')

# ADGRE5 AFS plot (Figure S6)
ADGRE5.afs.diff.plot = ggplot() +
  geom_boxplot(data=sum.afs.all2 %>% filter(gene=="ADGRE5"), aes(x=Bin, y=diff), fill="#56B4E9", position=position_dodge(width=1)) +
  geom_hline(yintercept=0, linetype=2, linewidth=1.1) +
  facet_wrap(~pop, nrow=2, ncol=2, scales="free_x") + ylim(-0.55, 0.55) +
  labs(y='Difference in Proportion of Functional Variants (Simulated - Target)', x='MAC bin', title="ADGRE5 (large)") +
  theme_bw(base_size=17)

ggsave(file = "./test/prop_diff_ADGRE5.jpg", plot = ADGRE5.afs.diff.plot, height = 8, width = 14, units = 'in')

