
# function to convert the haplotypes into a genotype matrix
make_geno = function(hap){
  
  # create an empty genotype matrix
  geno = matrix(0, nrow(hap), ncol(hap)/2)
  
  # sum up the number of alleles in adjacent haplotypes (2 haplotypes per person)
  for (j in 1:(ncol(hap)/2)){
    geno[,j] = hap[,2*j] + hap[,2*j-1]
  }
  geno = as.data.frame(geno)
  
  return(geno)
}


# function to calculate the allele counts/frequencies
calc_allele_freqs = function(geno, n, sum=NULL, source=NULL) {
  
  #counts = data.frame(count = rowSums(geno)) %>%
  # mutate(mac = ifelse(count>n, 2*n-count, count)) %>%
  #  mutate(maf = mac/(2*n))

  if (!is.null(sum)) {
    counts = data.frame(ac = sum) %>% mutate(af = ac/(2*n))
  } else {
    counts = data.frame(ac = rowSums(geno)) %>% mutate(af = ac/(2*n))
  }
  
  if (!is.null(source)) {counts$source = source}
  
  return(counts)
}


# function to flip values for a specified file at variants with an AF >= 1-maf 
# (used in the flip_data function below)
flip_file = function(file_to_flip, flip, file_type, N=NULL) {
  
  # Make copy of original file
  file2 = file_to_flip
  
  if (file_type == "leg") {
    
    # Flip ref allele in file2 with alt allele in file
    file2[flip, "a0"] = file_to_flip[flip, "a1"]
    
    # Flip alt allele in leg2 with ref allele in leg
    file2[flip, "a1"] = file_to_flip[flip, "a0"]
    
    return(file2)
    
  } else if (file_type == "geno") {
    
    # Flip the alternate allele counts at the relevant variants
    file2[flip, ] <- 2-file_to_flip[flip, ]
    
    return(file2)
    
  } else if (file_type == "count") {
    
    # Flip the ACs at the variants that need to be flipped
    file2[flip, "ac"] = (2*N)-file_to_flip[flip, "ac"]
    
    # Flip the AFs at the variants that need to be flipped
    file2[flip, "af"] = 1-file_to_flip[flip, "af"]
    
    return(file2)
    
  } else {
    stop("ERROR: 'file_type' must be a string of either 'leg', 'geno', or 'count'")
  }
  
}


# function to flip the relevant datasets at variants where AF >= 1-maf for each possible scenario
flip_data = function(leg, flip, geno_case, count_case, Ncase, cntrl, geno_ic=NULL, count_ic=NULL, Nic=NULL, geno_cc=NULL, count_cc=NULL, Ncc=NULL) {
    
    # Create new leg file
    leg2 = flip_file(leg, flip, file_type="leg")
    
    # Create new case data files
    geno_case2 = flip_file(geno_case, flip, file_type="geno")
    count_case2 = flip_file(count_case, flip, file_type="count", N=Ncase)
    
    if (cntrl == "int") {
      
      # Update geno files
      geno_ic2 = flip_file(geno_ic, flip, file_type="geno")
      
      # Recalculate ac/af 
      count_ic2 = flip_file(count_ic, flip, file_type="count", N=Nic)
      
      # Return all changed files, note some may be NULL
      return(list(leg=leg2, geno.case=geno_case2, geno.ic=geno_ic2, count.case=count_case2, count.ic=count_ic2))
      
    } else if (cntrl == "ext") {
      
      # Update geno files
      geno_cc2 = flip_file(geno_cc, flip, file_type="geno")
      
      # Recalculate ac/af 
      count_cc2 = flip_file(count_cc, flip, file_type="count", N=Ncc)
      
      # Return all changed files, note some may be NULL
      return(list(leg=leg2, geno.case=geno_case2, geno.cc=geno_cc2, count.case=count_case2, count.cc=count_cc2))
      
    } else if (cntrl == "all") {
      
      # Update geno files
      geno_ic2 = flip_file(geno_ic, flip, file_type="geno")
      geno_cc2 = flip_file(geno_cc, flip, file_type="geno")
      
      # Recalculate ac/af 
      count_ic2 = flip_file(count_ic, flip, file_type="count")
      count_cc2 = flip_file(count_cc, flip, file_type="count")
      
      # Return all changed files, note some may be NULL
      return(list(leg=leg2, geno.case=geno_case2, geno.ic=geno_ic2, geno.cc=geno_cc2, count.case=count_case2, count.ic=count_ic2, count.cc=count_cc2))
    }
}


# define function for calculating power
my.power = function(values, alpha=0.05){
  values2 = values[!is.na(values)]
  sig = which(as.numeric(values2) <= alpha)
  out = length(sig)/length(values2)
  return(out)
}


# define function to add confidence intervals to power results
add_CIs = function(results, nsim){
  
  results$Lower = '.'
  results$Upper = '.'
  
  ## CI stack exchange:
  # https://stats.stackexchange.com/questions/82720/confidence-interval-around-binomial-estimate-of-0-or-1
  # paper: https://projecteuclid.org/journals/statistical-science/volume-16/issue-2/Interval-Estimation-for-a-Binomial-Proportion/10.1214/ss/1009213286.full?tab=ArticleFirstPage
  
  for(i in 1:nrow(results)){ # default level is 95% confidence
    results$Lower[i] = binom.confint(nsim*results$Value[i], nsim, method=c("wilson"), type="central")$lower
    results$Upper[i] = binom.confint(nsim*results$Value[i], nsim, method=c("wilson"), type="central")$upper
  }
  
  results$Lower = as.numeric(results$Lower)
  results$Upper = as.numeric(results$Upper)
  
  return(results)
}
                                