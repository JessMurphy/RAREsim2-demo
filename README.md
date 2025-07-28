# RAREsim2 Demonstration

This repository provides the workflow for the RAREsim2 demonstration detailed here (provide link to paper). 

## Computing Environment
We used a high-performance computing cluster with 2048 AMD cores and 16TB memory in 32 compute nodes. Each node has 2 AMD EPYC 7502 32 core processors for a total of 64 cores, 512GB DDR4 memory, and dual 960GB SSD. 

We ran the pipeline using a singularity container with R version 4.2.1, python version 3.10.6, and Hapgen2 version 2.2.0 installed. See below for a list of R packages that were used.

![R session info](R%20session%20info.png)

## Pipeline

### 0_make_directories.sh

Make directories to store the results and different datasets needed for each population.

### 1_submit_prep_data.sh

Submits the pipeline for preparing the input data for the simulations.

* **1a_subset_data.R**: unpacks the 1000 Genomes and gnomAD data for chromosome 19 and subsets them to the given block for each population
* **1b1_make_master_legend.R**: creates reference haplotypes for Hapgen2 based on the 1000G data and a master legend file with annotation information for all known and unknown coding positions
* **1b2_get_annotations.R**: obtains the reference and alternate alleles for all uknown coding positions and provides functional annotations for all positions (must be used alongside Step b1)
* **1c_subset_master_legend.R**: subsets the master legend file to produce a unique legend file for each simulation replicate with only one variant per position

The resulting datasets from steps 1a, 1b1, and 1b2 for Block 37 are provided in the Input directory.

### 2_submit_hapgen2.sh

Submits the script to run Hapgen2.

* **2a_run_hapgen2.sh**: produces simulated haplotypes with an over-abundance of rare variants (see below for the input parameters used for each population)

|**Population**	|**No. of<br>Haplotypes**|**Effective<br>Population Size**|
|:--------------|:----------------------|:------------------------------|
| AFR		| 1,008			| 17,469 			|
| EAS		| 1,008			| 14,269 			|
| NFE		| 808			| 11,418			|
| SAS		| 978			| 14,269			|	

Note: Heterozygote/homozygote risks of 1.00/1.00 were used along with a disease locus of 14499614 across the populations.

### 3_submit_power_analysis.sh

Submits the pipeline for running the power analysis.

* **3a_run_RAREsim2_power.sh**: runs RAREsim2 to prune the variants down for the same and opposite directions of effect power scenarios
* **3b_run_RAREsim2_power_unequal.sh**: runs RAREsim2 to prune the variants down for the opposite direction of effect power scenario with unequal proportions of risk and protective variants
* **3c_run_methods_opp_power.R**: runs the methods to calculate the power for the opposite direction of effect scenario
* **3c_run_methods_opp_power_unequal.R**: runs the methods to calculate the power for the opposite direction of effect scenario with unequal proportions of risk and protective variants
* **3c_run_methods_same_power.R**: runs the methods to calculate the power for the same direction of effect scenario

Note: SKATBinary was used in the R scripts.

### 4_submit_t1e_analysis.sh

* **4a_run_RAREsim2_t1e.sh**: runs RAREsim2 to prune the variants down for the type I error scenario with no effect
* **4b_run_methods_t1e.R**: runs the methods to calculate the type I error for no effect scenario

### 5_plot_results.R

## Run Time

|**Step**	|**Time (hh:mm:ss)**|**Notes** |
|:--------|:------------------|---------|
| 0_make_directories.sh | 00:00:00 | |
| 1a_subset_data.R | 00:09:15 | for all populations |
| 1b1_make_master_legend.R | 00:00:33 | per population |
| 1b2_get_annotations.sh | | per population |
| 1c_subset_master_legend.R | 00:01:03 | per batch of 100 replicates per population |
| 2a_run_hapgen2.sh | 01:20:41 | per batch of 100 replicates per population | 
| 3a_run_RAREsim2_power.sh | | |
| 3b_run_RAREsim2_power_unequal.sh | | |
| 3c_run_methods_*_power.R | | |
| 4a_run_RAREsim2_t1e.sh | | |
| 4b_run_methods_t1e.R | | |
| 5_plot_results.R | | |
|**Total**| | |

