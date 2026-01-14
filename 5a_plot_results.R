
# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(binom, lib.loc="/home/math/murphjes/R/myLibs/")

source("./methods_functions.R")

# define variables
p.same = c(160, 140, 120)                # percentages for the expected number of functional variants in cases (same direction of effect)
p.opp = c(130, 120, 110)                 # percentages for the expected number of functional variants in cases & controls (opposite directions of effect)
p.opp1 = c(145, 130, 115)                # percentages for the expected number of functional variants in cases (opposite directions of unequal effect)
p.opp2 = c(115, 110, 105)                # percentages for the expected number of functional variants in controls (opposite directions of unequal effect)
pop.list = c("AFR", "EAS", "NFE", "SAS") # population list

# read in the results from each population
t1e.all = power.same.all = power.opp.all = power.opp.unequal.all = c()
for (pop in pop.list){
  
  # read in the type I error results for each batch of 100
  t1e.rounds = c()
  for (i in 1:100){
    
    t1e = read.table(paste0("./results/", pop, "/t1e_results_", pop, "_", i*100, ".txt"), header = TRUE) %>% mutate(Round=i)
    t1e.rounds = rbind(t1e.rounds, t1e)
  }
  
  t1e.all = rbind(t1e.all, t1e.rounds %>% group_by(Method) %>% summarize(across(ADGRE2:ZNF333, ~ sum(.)/100)) %>% mutate(Pop=pop))

  # read in the power results for each scenario and batch of 100
  for (p in 1:length(p.same)){

    power.same.rounds = power.opp.rounds = power.opp.unequal.rounds = c()
    for (i in 1:10){
    
       power.same = read.table(paste0("./results/", pop, "/power_same_results_", pop, "_", p.same[p], "_", i*100, ".txt"), header = TRUE) %>% mutate(Round=i)
       power.opp = read.table(paste0("./results/", pop, "/power_opp_results_", pop, "_", p.opp[p], "_", i*100, ".txt"), header = TRUE) %>% mutate(Round=i)
       power.opp.unequal = read.table(paste0("./results/", pop, "/power_opp_results_", pop, "_", p.opp1[p], "_", p.opp2[p], "_", i*100, ".txt"), header = TRUE) %>% mutate(Round=i)
    
       power.same.rounds = rbind(power.same.rounds, power.same)
       power.opp.rounds = rbind(power.opp.rounds, power.opp)
       power.opp.unequal.rounds = rbind(power.opp.unequal.rounds, power.opp.unequal)
    }

    power.same.all = rbind(power.same.all, power.same.rounds %>% group_by(Method) %>% summarize(across(ADGRE2:ZNF333, ~ sum(.)/10)) %>% mutate(Pop=pop, Scenario="Power (same)", pcase=p.same[p])) #/10
    power.opp.all = rbind(power.opp.all, power.opp.rounds %>% group_by(Method) %>% summarize(across(ADGRE2:ZNF333, ~ sum(.)/10)) %>% mutate(Pop=pop, Scenario="Power (opp)", pcase=p.opp[p])) #/10
    power.opp.unequal.all = rbind(power.opp.unequal.all, power.opp.unequal.rounds %>% group_by(Method) %>% summarize(across(ADGRE2:ZNF333, ~ sum(.)/10)) %>% mutate(Pop=pop, Scenario="Power (opp unequal)", pcase=paste0(p.opp1[p], ", ", p.opp2[p]))) #/10
  }
}
power.all = rbind(power.same.all, power.opp.all, power.opp.unequal.all)

# convert the results to long format
t1e.long = pivot_longer(t1e.all, ADGRE2:ZNF333, names_to="Gene", values_to="Value")
power.long = pivot_longer(power.all, ADGRE2:ZNF333, names_to="Gene", values_to="Value") 

# add confidence intervals
t1e.long2 = add_CIs(t1e.long, nsim=10000)
power.long2 = add_CIs(power.long, nsim=1000)

# colorblind friendly palette
colors.method = c("#56B4E9", "#0072B2", "#009E73")

# subset of genes to plot (big, medium, small)
genes.power = c("ADGRE5", "ADGRE3", "TECR") 

# type I error plot
t1e.plot = ggplot(t1e.long2 %>% filter(Gene %in% genes.power) %>% 
         mutate(Gene=factor(Gene, levels=c("TECR", "ADGRE3", "ADGRE5"), labels=c("TECR (small)", "ADGRE3 (medium)", "ADGRE5 (large)"))),
       aes(x=Pop, y=Value, color=Method)) +
  geom_point(size=2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1.1) +
  scale_y_continuous(limits=c(0.025,0.075)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
  scale_color_manual(values=colors.method) +
  facet_wrap(~Gene, nrow=1) +
  labs(y="Type I Error", x="Population", title="Type I Error Results") +
  theme_bw(base_size=17) #+ theme(legend.position="bottom")

ggsave(file = "./results/t1e_plot.jpg", plot = t1e.plot, height = 4, width = 14, units = 'in')

t1e.plot2 = ggplot(t1e.long2 %>% filter(Gene=="ADGRE3"),
                  aes(x=Pop, y=Value, color=Method)) +
  geom_point(size=2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1.1) +
  scale_y_continuous(limits=c(0.025,0.075)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
  scale_color_manual(values=colors.method) +
  labs(y="Type I Error", x="Population", title="Type I Error Results") +
  theme_bw(base_size=15) + theme(legend.position="bottom")

ggsave(file = "./results/t1e_plot2.jpg", plot = t1e.plot2, height = 4, width = 4.5, units = 'in')

# power (same direction of effect) plot
for (p in p.same){
  
  power.same.plot = ggplot(power.long2 %>% filter(Gene %in% genes.power, Scenario=="Power (same)", pcase==p) %>% 
                             mutate(Gene=factor(Gene, levels=c("TECR", "ADGRE3", "ADGRE5"), labels=c("TECR (small)", "ADGRE3 (medium)", "ADGRE5 (large)"))),
                           aes(x=Pop, y=Value, color=Method)) +
    geom_point(size=2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept=0.8, linetype=2, linewidth=1.1) +
    scale_y_continuous(limits=c(0,1)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
    scale_color_manual(values=colors.method) +
    facet_wrap(~Gene, nrow=1) +
    labs(y="Power", x="Population", title=paste0("Power Results (same direction of effect): ", p, "%")) +
    theme_bw(base_size=17) #+ theme(legend.position="top")
  
  ggsave(file = paste0("./results/power_same_plot_", p, ".jpg"), plot = power.same.plot, height = 6, width = 14, units = 'in')
}

power.same.plot2 = ggplot(power.long2 %>% filter(Gene=="ADGRE3", Scenario=="Power (same)", pcase==p.same[2]),
                   aes(x=Pop, y=Value, color=Method)) +
  geom_point(size=2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0.8, linetype=2, linewidth=1.1) +
  scale_y_continuous(limits=c(0,1)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
  scale_color_manual(values=colors.method) +
  labs(y="Power", x="Population", title=paste0("Power Results (same direction of effect): ", p.same[2], "%")) +
  theme_bw(base_size=15) + theme(legend.position="bottom")

ggsave(file = paste0("./results/power_same_plot2_", p.same[2], ".jpg"), plot = power.same.plot2, height = 6, width = 4.5, units = 'in')

# power (opposite direction of effect) plots
for (p in p.opp){
  
  power.opp.plot = ggplot(power.long2 %>% filter(Gene %in% genes.power, Scenario=="Power (opp)", pcase==p) %>% 
                            mutate(Gene=factor(Gene, levels=c("TECR", "ADGRE3", "ADGRE5"), labels=c("TECR (small)", "ADGRE3 (medium)", "ADGRE5 (large)"))),
                          aes(x=Pop, y=Value, color=Method)) +
    geom_point(size=2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept=0.8, linetype=2, linewidth=1.1) +
    scale_y_continuous(limits=c(0,1)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
    scale_color_manual(values=colors.method) +
    facet_wrap(~Gene, nrow=1) +
    labs(y="Power", x="Population", title=paste0("Power Results (opposite direction of effect): ", p, "%")) +
    theme_bw(base_size=17) #+ theme(legend.position="top")
  
  ggsave(file = paste0("./results/power_opp_plot_", p, ".jpg"), plot = power.opp.plot, height = 6, width = 14, units = 'in')
}

power.opp.plot2 = ggplot(power.long2 %>% filter(Gene=="ADGRE3", Scenario=="Power (opp)", pcase==p.opp[2]),
                          aes(x=Pop, y=Value, color=Method)) +
  geom_point(size=2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0.8, linetype=2, linewidth=1.1) +
  scale_y_continuous(limits=c(0,1)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
  scale_color_manual(values=colors.method) +
  labs(y="Power", x="Population", title=paste0("Power Results (opposite direction of effect): ", p.opp[2], "%")) +
  theme_bw(base_size=15) + theme(legend.position="bottom")

ggsave(file = paste0("./results/power_opp_plot2_", p.opp[2], ".jpg"), plot = power.opp.plot2, height = 6, width = 4.5, units = 'in')

# power (opposite direction of effect with unequal percentages) plots
for (p in 1:length(p.opp)){
  
  power.opp.unequal.plot = ggplot(power.long2 %>% filter(Gene %in% genes.power, Scenario=="Power (opp unequal)", pcase==paste0(p.opp1[p], ", ", p.opp2[p])) %>% 
                            mutate(Gene=factor(Gene, levels=c("TECR", "ADGRE3", "ADGRE5"), labels=c("TECR (small)", "ADGRE3 (medium)", "ADGRE5 (large)"))),
                          aes(x=Pop, y=Value, color=Method)) +
    geom_point(size=2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept=0.8, linetype=2, linewidth=1.1) +
    scale_y_continuous(limits=c(0,1)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(width = 0.5), linewidth=1.1) +
    scale_color_manual(values=colors.method) +
    facet_wrap(~Gene, nrow=1) +
    labs(y="Power", x="Population", title=paste0("Power Results (opposite direction of effect - unequal): ", p.opp1[p], "%, ", p.opp2[p], "%")) +
    theme_bw(base_size=17) #+ theme(legend.position="top")
  
  ggsave(file = paste0("./results/power_opp_plot_", p.opp1[p], "_", p.opp2[p], ".jpg"), plot = power.opp.unequal.plot, height = 6, width = 14, units = 'in')
}
