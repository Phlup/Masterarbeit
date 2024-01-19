##execution order:
#get raw data (processes study data from ... and ...)
NAM_data_processing.R
recomb_rates.R
#get genotype simulation
sim_all_pop.py
sim_sim_pops.py
#get genotype validation
gen_plots.py
validation_stats.R
#get rrblup mrk effects and prep data for prediction
rrBLUP_phenos_all.R
#run phenotype prediction
phenotype_prediction.R
CNN_pheno_predict.py

##notebooks showcasing pedigree functions and genotype simulation
pedigrees.ipynb
genotype_sim.ipynb