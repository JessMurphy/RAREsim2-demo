# RAREsim2 Demonstration

This repository provides the workflow for the RAREsim2 demonstration detailed here (provide link to paper). 

## Computing Environment
We used a high-performance computing cluster with 2048 AMD cores and 16TB memory in 32 compute nodes. Each node has 2 AMD EPYC 7502 32 core processors for a total of 64 cores, 512GB DDR4 memory, and dual 960GB SSD. 

We ran the pipeline using a singularity container with R version 4.2.1, python version 3.10.6, and Hapgen2 version 2.2.0 installed.
* SKATBinary

![R session info](R%20session%20info.png)

## Pipeline

### 0_make_directories.sh

### 1_submit_legend_files.sh

* **1a_make_master_legend.R**:
* **1b_subset_master_legend.R**:

### 2_submit_hapgen2.sh

* **2a_run_hapgen2.sh**:

Hapgen2 requires reference haplotypes and an accompanying legend file as well as a recombination map and disease locus as input. 

*Table of Hapgen2 Input Parameters*

|**Population**	|**No. of<br>Haplotypes**|**Disease<br>Locus**	|**Effective<br>Population Size**|
|:--------------|:----------------------|:----------------------|:------------------------------|
| AFR		| 1,008			| 14627281 		| 17,469 			|
| EAS		| 1,008			| 14673368 		| 14,269 			|
| NFE		| 808			| 14705483		| 11,418			|
| SAS		| 978			| 14508902 		| 14,269			|	

Heterozygote/homozygote risks of 1.00/1.00 were used.

### 3_submit_power_analysis.sh

* **3a_run_RAREsim2_power.sh**:
* **3b_run_RAREsim2_power_unequal.sh**:
* **3c_run_methods_opp_power.R**:
* **3c_run_methods_opp_power_unequal.R**:
* **3c_run_methods_same_power.R**:

### 4_submit_t1e_analysis.sh

* **4a_run_RAREsim2_t1e.sh**:
* **4b_run_methods_t1e.R**:

### 5_plot_results.R

## Run Time

|**Step**	|**Time (hh:mm:ss)**|
|:--------------|:----------------------|
| 0_make_directories.sh | |
| 1a_make_master_legend.R | 00:00:33 |
| 1b_subset_master_legend.R | 00:01:03 |
| 2a_run_hapgen2.sh | 01:20:41 | |
| 3a_run_RAREsim2_power.sh | |
| 3b_run_RAREsim2_power_unequal.sh | |
| 3c_run_methods_*_power.R | |
| 4a_run_RAREsim2_t1e.sh | |
| 4b_run_methods_t1e.R | |
| 5_plot_results.R | |
|**Total**| |

