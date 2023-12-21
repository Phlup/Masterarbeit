library(Matrix)
library(MASS)
library(crayon)
library(sommer)
library(openxlsx)

##train rrBLUP and gBLUP models to calc accuracy on real and sim datasets (MSE test/train traits) 
##and determine best cross of sim datasets (max traits with calc marker effects)
##calc accuracy of simulated best crosses based on real pheno max

#read in NAM phenotypes per individual and mean over pop+env
NAM_phenotypes <- read.csv("../data/NAM_phenotype_data/NAM_phenotypes.csv")

#read populations
populations <- read.csv("../data/sim_data/populations.csv")

#create output dfs
rrblup_preds <- data.frame("pop" = NULL, "ind" = NULL, "trait" = NULL, "pred" = NULL)
gblup_preds <- data.frame("pop" = NULL, "ind" = NULL, "trait" = NULL, "pred" = NULL)
pred_summary <- data.frame("population" = populations$pop, "trait_mean" = NA, "trait_max" = NA, "trait_min" = NA)
traits <- c("silk", "tassel", "oil", "protein", "starch")

#for now, focus on env Aurora NY 2006 (06A) and Columbia MO 2006 (06MO)
#for each env individually, no fixed effect calc
env <- "06A"
for(i in populations$pop){
  
  #load real and simulated genotypes (additive encoding)
  sim_add <- read.csv(paste("../sim_output/normal_rec/additive_encoding/add_",i,".csv", sep = ""))
  real_add <- read.csv(paste("../data/NAM_genotype_data/additive_encoding/pop_",i,"_add.csv", sep = ""))
  
  for(j in traits){
    #subset real and sim add and  phenotypes based on intersection of individual ids per env and pop
    phenos <- NAM_phenotypes[NAM_phenotypes$pop == i & NAM_phenotypes$env == env & NAM_phenotypes$trait == j,]
    IDs <- intersect(real_add$individual, phenos$individual)
    sim_add <- sim_add[IDs,]
    real_add <- real_add[real_add$individual %in% IDs,]
    phenos <- phenos[phenos$individual %in% IDs,]

    ##run rrBLUP
    rBLUP <- mmer(value ~ 1, random = ~vsr(list(real_add[-1])), rcov=~units, data = phenos, verbose = FALSE)
    
    #get marker effects
    effects <- rBLUP$U$`u:real_add`$value
    
    #calc sample values
    sample_values<-unlist(lapply(genotypes_to_predict, function(x){sum(mrk_effects*GT_wheat[x,])}))
    
    rrblup_preds[]
    
    pred_traits[pred_traits$pop == i, c(j_mean, j_max, j_min)] <- c(mean(sample_values), max(sample_values), min(sample_values)) 
  }
}

#load real mean values, calc mse between real and add per trait
MSE(pred_traits$trait_i, pheno_means$trait_i)

cor.test(pred_traits)

#max trait
max(pred_traits$trait_i)

#max over all traits
max(rowSums(pred_traits[traits]))

##model comparison similar to john et al comparison of classical and ml methods 2022
##-> similar model ranking? + overall accuracy + feature usage for prediction as three outcomes

##rrblup for marker effects -> sum up marker effects for total trait value and compare mean trait value of best
##performing cross with best performing real pop or do mse etc.

##gblup for genotype effects -> what do genotype effects imply
##xgboost feature importance: small number of always important features? ifnot, might not be usfeul
##cnn window approach: conserved/consistently important regions for features? ifnot, might not be useful
#-> can further information about breeding be gained from snp importance? if not, choose best performing model
#discussion