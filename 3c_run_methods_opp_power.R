# Get command-line arguments
args = commandArgs(trailingOnly=TRUE)

# define variables
id = as.numeric(args[[1]])               # array task ID for the current job
pop.index = ceiling(id/1000)             # population index (1-4) for the current job
pop.list = c("AFR", "EAS", "NFE", "SAS") # population list
Pop = pop.list[pop.index]                # specific population for the current job
end = id - 1000*(pop.index-1)            # ending simulation replicate number for the loop
start = end-99                           # starting simulation replicate number for the loop

p.case = c(130, 120, 110)                # percentages for the expected number of functional variants in cases & controls
maf = 0.01                               # minor allele frequency threshold to define rare variants
Ncase = 5000                             # number of cases
Nsim = 10000                             # total simulated sample size (number of individuals)

# Print confirmation
print(paste0("Running simulation with: Pop = ", Pop, ", Batch = ", end))

# load libraries
library(dplyr)
library(tidyr)
library(SKAT)
library(data.table)

source("./methods_functions.R")

set.seed(0)

# loop through the pcase percentages
for (j in 1:length(p.case)){

# create empty vectors to store the p-values from each replicate
burden.genes.p = skato.genes.p = skat.genes.p = c()
out.genes.all.p = c()

# loop through the simulation replicates
for(i in start:end){
  
  # read in the legend file
  leg = read.table(paste0("./datasets/Hapgen", Nsim/1000, "K_pruned/", Pop, "/Round1/chr19.block37.", Pop, ".sim", i, ".", Nsim, ".", p.case[j], "fun.100syn.protected.legend"), header=T) %>%
    mutate(gene=ifelse(gene=="ZNF333;ZNF333(NM_001352243:exon9:UTR5)", "ZNF333", gene))
  leg$row = 1:nrow(leg)
  
  leg.pruned = read.table(paste0("./datasets/Hapgen", Nsim/1000, "K_pruned/", Pop, "/Round1/chr19.block37.", Pop, ".sim", i, ".", Nsim, ".opp.", p.case[j], "fun.100syn.legend-pruned-variants"), header=T)
  
  # identify pruned and protected variants for later sanity checks
  protected = which(leg$protected==1)
  pruned = which(leg$id %in% leg.pruned$id)
  
  # read in the haplotype files
  cases.hap = fread(paste0("./datasets/Cases/", Pop, "/Round1/chr19.block37.", Pop, ".sim", i, ".", Ncase, ".power.cases.opp.", p.case[j], "fun.100syn.haps-sample.gz"))
  cases.hap = as.data.frame(cases.hap)
  
  controls.hap = fread(paste0("./datasets/Controls/", Pop, "/Round1/chr19.block37.", Pop, ".sim", i, ".", Ncase, ".controls.opp.", p.case[j], "fun.100syn.haps-remainder.gz"))
  controls.hap = as.data.frame(controls.hap)
  
  # convert the haplotypes into genotypes
  cases.geno = make_geno(cases.hap)
  controls.geno = make_geno(controls.hap)
  
  # calculate the allele counts/frequencies
  cases.count = calc_allele_freqs(cases.geno, ncol(cases.geno))
  controls.count = calc_allele_freqs(controls.geno, ncol(controls.geno))
  
  # identify variants where the alternate and reference alleles would need to be flipped
  flip = which(cases.count$af >= 1-maf | controls.count$af >= 1-maf)
  
  print(paste0("There are ", length(flip), " variants that need to be flipped for the cases/controls."))
  
  if (length(flip) != 0) {
    
    flip.data = flip_data(leg, flip, cases.geno, cases.count, Ncase, cntrl="int", controls.geno, controls.count, Ncase)
    
    leg2 = flip.data[["leg"]]
    cases.geno2 = flip.data[["geno.case"]]
    controls.geno2 = flip.data[["geno.ic"]]
    cases.count2 = flip.data[["count.case"]]
    controls.count2 = flip.data[["count.ic"]]
    
    # double check all the variants have been flipped if needed
    flip2 = which(cases.count2$af >= 1-maf | controls.count2$af >= 1-maf)
    
    print(paste0("There are now ", length(flip2), " variants that need to be flipped."))
    
  } else {
    
    leg2 = leg
    cases.geno2 = cases.geno
    controls.geno2 = controls.geno
    cases.count2 = cases.count
    controls.count2 = controls.count
  }
  
  # double check there are no pruned variants in the controls
  cases.pruned = cases.count2[pruned,]
  controls.pruned = controls.count2[pruned,]
  
  print(paste0("There are ", sum(controls.pruned$ac), " pruned variants in the controls."))
  
  # double check there are no protected variants in the cases
  cases.protected = cases.count2[protected,]
  controls.protected = controls.count2[protected,]
  
  print(paste0("There are ", sum(cases.protected$ac), " protected variants in the cases."))
  
  # identify the common variants
  common = which((cases.count2$af > maf & cases.count2$af < 1-maf) | (controls.count2$af > maf & controls.count2$af < 1-maf))
  
  # create case/control phenotype matrices for SKAT
  pheno = rep(0, (ncol(cases.geno2) + ncol(controls.geno2))) 
  pheno[1:ncol(cases.geno2)] = 1
  
  # create combined genotype matrix
  geno.skat = cbind(leg2[, c("gene", "fun")], cases.geno2, controls.geno2)[-common,]
  colnames(geno.skat) = make.unique(colnames(geno.skat))   # make the column names of the genotype matrix unique
  geno.skat.fun = geno.skat %>% filter(fun=="fun")
  
  # create null model object for SKAT
  obj = SKAT_Null_Model(as.numeric(pheno) ~ 1, out_type="D") # D-dichotomous
  
  # create empty vectors to store the p-values from each gene
  burden.genes.p = skato.genes.p =  skat.genes.p = c()
  
  # call SKAT once per gene
  genes = levels(droplevels(as.factor(leg2$gene)))
  for(g in 1:length(genes)){ # loop through the genes
    
    # subset genotype matrix
    Z = geno.skat.fun %>% filter(gene==genes[g]) %>% select(-gene, -fun)
    
    # call the SKAT functions
    skat.genes = SKATBinary(t(Z), obj, method = 'SKAT')
    skato.genes = SKATBinary(t(Z), obj, method = 'SKATO')
    burden.genes = SKATBinary(t(Z), obj, method = 'Burden')
    
    # save the p-values
    burden.genes.p = c(burden.genes.p, burden.genes$p.value)
    skato.genes.p = c(skato.genes.p, skato.genes$p.value)
    skat.genes.p = c(skat.genes.p, skat.genes$p.value)
  }
  
  # combine the results from all three methods
  genes.all.p = as.data.frame(rbind(burden.genes.p, skato.genes.p, skat.genes.p))
  colnames(genes.all.p) = genes
  rownames(genes.all.p) = NULL
  genes.all.p$Method = c("Burden", "SKAT-O", "SKAT")
  
  out.genes.all.p = rbind(out.genes.all.p, genes.all.p)
  
  print(i)
}

# save the p-value results
write.table(out.genes.all.p, paste0("./results/", Pop, "/power_opp_pvalues_", Pop, "_", p.case[j], "_", end, ".txt"), quote=F, row.names=F, col.names=T)

# calculate the power for each method
power.genes.out = out.genes.all.p %>% group_by(Method) %>% group_map(~apply(.x, 2, FUN=my.power)) %>% setNames(unique(sort(out.genes.all.p$Method))) %>% bind_rows(.id="Method")

# save the power results
write.table(power.genes.out, paste0("./results/", Pop, "/power_opp_results_", Pop, "_", p.case[j], "_", end, ".txt"), quote=F, row.names=F, col.names=T)

}