This is a short demonstration of using the population genetics tool msprime and
tskit to simulate entire genotypes of arbitrary pedigrees. Here, real data of
a cross between two recombinant inbred maize genotypes followed by 5 selfing
generations has been succesfully remodeled (1).

Redo the conda environment using the enivornment.yaml file:
conda env create -f environment.yaml

Then use the interpreter of that environment to run main.py or genotype_sim_demo.ipynb

Project structure:
data/ 
-contains material from (1, NAM_map_and_genos-121025) and (2, Buckler_etal_2009_Science_flowering_time_data-090807).
-cleaned genotype data for all 26 populations from (1) in NAM_genotype_data
-data necessary for genotype simulation (reference allele, genmap, NAM_parent_genos) in sim_data
plots/
-collects plots from genotype simulation
project_lib/
-contains .py files containing functions for genotype simulation and statistical evaluation of simulation
R_scripts/
-contains recomb_rates.R that shows how the genetic map from (1) was processed to get cM/Mb rates necessary for msprime
sim_output/
-collects simulation output files

The jupyter notebook genotype_sim_demo.ipynb demonstrates the genotype_simulation
function and data/models used by it.

main.py shows an example of simulating population 1 of NAM_genotype_data and using a summary plot to compare genotype distributions
simulation output is saved into sim_output/ and plot into plots/

(1): Genetic Design and Statistical Power of Nested Association Mapping in Maize,
Yu et al. 2008, Genetics, Volume 178, Issue 1. https://doi.org/10.1534/genetics.107.074245

(2): The genetic architecture of maize flowering time, Buckler et al. 2009, 
Science. 2009 Aug 7;325(5941):714-8. doi: 10.1126/science.1174276.

(3): Genetic Architecture of Maize Kernel Composition in the Nested Association Mapping and Inbred Association Panels,
Cook et al. 2012, Plant Physiol. 2012 Feb; 158(2): 824–834. doi: 10.1104/pp.111.185033


##execution order:
#get raw data 
NAM_data_processing.R +
recomb_rates.R +
#get genotype simulation
sim_all_pop.py +
sim_sim_pops.py +
#get genotype validation
gen_plots.py +
validation_stats.R + 
#get rrblup mrk effects and prep data for prediction
rrBLUP_phenos_all.R +
#run phenotype prediction
phenotype_prediction.R +
CNN_pheno_predict.py +
##notebooks showcasing pedigree functions and genotype simulation
pedigrees.ipynb
genotype_sim.ipynb

Data flowchart:
![alt text](https://github.com/Phlup/Masterarbeit/blob/main/sim_predict_pipeline.jpg?raw=true)
