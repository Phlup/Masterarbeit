## Description
This is the code and data repository for a population genetics-based genotype simulator and ML-based phenotype prediction pipeline for plant breeding.
The pipeline is capable of effectively simulating the development of plant breeding populations across different mating designs and species-specific recombination parameters.
Given estimated marker effects in an additive phenotype model, it is capable of predicting phenotypic summary statistics based on genotypic data of potential parent crosses
and their simulated offspring. The method is validated on genotypic and phenotypic data of multiple studies on 25 maize NAM populations across multiple phenotypic traits
(see (1) for genotype data, (2) and (3) for phenotype data). The genotype simulator based on fixed breeding pedigrees is an extension of the population genetics tools
msprime (https://github.com/tskit-dev/msprime) and tskit (https://github.com/tskit-dev/tskit). The ML-based phenotype prediction uses CNNs, RandomForests, XGBoost and baseline models.
See below for an overview of the file structure, execution order of scripts to reproduce results and visualized data-flow in the pipeline.



## Project code and data structure:	
code/	  
	R_scripts/	  
		NAM_data_processing.R - Processes genotype and phenotype data from (1), (2) and (3)	  
		phenotype_prediction.R - phenotype prediction based on RandomForest, XGBoost and baseline models	  
		recomb_rates.R - Processes genetic map from (1) to get cM/Mb rates necessary for msprime	  
		rrBLUP_phenos_all.R - estimates rrBLUP marker effects and calculates additive traits and summary statistics for real and simulated populations	  
		stat_functions.R - provides functions for validation_stats.R and rrBLUP_phenos_all.R	  
		validation_stats.R - calculates statistics to validate genotype simulation	  
		xtables.R - generates tables for thesis	  
	data/ 	  
		Buckler_etal_2009_Science_flowering_time_data-090807/	  
			relevant data from (2)	  
		NAM_genotype_data/	  
			cleaned genotype data from (1)	  
		NAM_map_and_genos-121025/	  
			relevant data from (1)	  
		NAM_phenotype_data/	  
			cleaned phenotype data from (2) and (3)	  
		sim_data/	  
			cleaned data necessary for genotype simulation	  
		test_data/	  
			data for testing genotype simulation methods	  
		supp_pp.111.185033_185033Supplemental_Data_S1.xlsx relevant data from (3)	  
	-contains material from (1, NAM_map_and_genos-121025) and (2, Buckler_etal_2009_Science_flowering_time_data-090807).	  
	-cleaned genotype data for all 26 populations from (1) in NAM_genotype_data	   
	-data necessary for genotype simulation (reference allele, genmap, NAM_parent_genos) in sim_data	  
	plots/	  
		collects plots from ARGs, pedigrees, genotype distributions, traits (cumulative sum and summary statistics), 	  
		population genetics (LD-decay and rogers distance) and phenotype prediction scatter plots	  
	project_lib/	  
		genotype_simulation.py - pedigree generating functions and genotype simulation functions	  
		stat_functions.py - generates statistics and plots for genotype simulation	  
	sim_output/	  
		collects simulation output files for all real populations across recombination scenarios and new population simulations	  
	stats/	  
		pheno_prediction/	  
			estimated rrBLUP marker effects and calculated additive traits and summary stats	  
			results/	  
				phenotype prediction results for all models	  
		summary_stats/	  
			summary statistics for genotype simulation outputs	  
	test_notebooks/	  
		-various test .ipynb notebooks for developing methods	  
	sim_all_pop.py - simulates all real populations	  
	sim_sim_pops.py - simulates new populations	  
	gen_plots.py - simulates summary plots for real populations	  
	CNN_pheno_predict.py - phenotype prediction script based on CNNs	  
	genotype_sim.ipynb - demonstrates genotype simulation functionality	  
	pedigrees.ipynb - demonstrates pedigree generating functions functionality	  

## Bibliography
(1): Genetic Design and Statistical Power of Nested Association Mapping in Maize,
Yu et al. 2008, Genetics, Volume 178, Issue 1. https://doi.org/10.1534/genetics.107.074245

(2): The genetic architecture of maize flowering time, Buckler et al. 2009, 
Science. 2009 Aug 7;325(5941):714-8. doi: 10.1126/science.1174276.

(3): Genetic Architecture of Maize Kernel Composition in the Nested Association Mapping and Inbred Association Panels,
Cook et al. 2012, Plant Physiol. 2012 Feb; 158(2): 824–834. doi: 10.1104/pp.111.185033

## Reproduction of results:
.py:
Redo the conda environment using the enivornment.yaml file:
conda env create -f environment.yaml
Then use the interpreter of that environment to run .py or .ipynb files
.R:
Install packages found in .R scripts with install.packages("pkg")

#### execution order:
#get raw data 
NAM_data_processing.R +
recomb_rates.R +
### get genotype simulation
sim_all_pop.py +
sim_sim_pops.py +
### get genotype validation
gen_plots.py +
validation_stats.R + 
### get rrblup mrk effects and prep data for prediction
rrBLUP_phenos_all.R +
### run phenotype prediction
phenotype_prediction.R +
CNN_pheno_predict.py +
### notebooks showcasing pedigree functions and genotype simulation
pedigrees.ipynb
genotype_sim.ipynb

## Data flowchart:
![alt text](https://github.com/Phlup/Masterarbeit/blob/main/sim_predict_pipeline.jpg?raw=true)
