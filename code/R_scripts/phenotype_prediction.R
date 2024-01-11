library(randomForest)
library(openxlsx)


#predict simulated offspring phenotype from parental genotype features, best cross prediction

#load sim geno + marker effects
mrk_eff <- read.csv("../stats/pheno_prediction/rrBLUP_mrk_effects.csv")
#read populations
populations <- read.csv("../data/sim_data/populations.csv")

#partition genos into features with multiple markers
#name features eg 1:10, 501:510

#calc rolling correlation along features

#read parent genotypes
#random forest prediction for sim and real phenos
#loss: diff(diff(trait_mean_real_pop_i_rrblup, trait_mean_sim_pop_i_rrblup),
#           diff(trait_mean_real_pop_i_RF, trait_mean_sim_pop_i_RF)

#randomForest(y = real_pheno data = real_geno, importance = TRUE)

#predict(rf, sim_geno)

#mse train, predict
