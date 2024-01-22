##execution order:
#get raw data 
(processes study data from Yu et al. 2008 https://academic.oup.com/genetics/article/178/1/539/6062286?login=true,
Bukler et al. 2009 https://pubmed.ncbi.nlm.nih.gov/19661422/ and
Cook et al. 2012 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3271770/)
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