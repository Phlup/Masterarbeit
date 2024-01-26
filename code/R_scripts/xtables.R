library(xtable)
library(stringr)
##create latex compatible results tables

#validation stats tables
sumstats_full <- read.csv("../stats/summary_stats/sum_results.csv")
xtable(t(sumstats_full))

rows_sumstats <- c("population","p heterozygosity test", "Phi heterozygosity", "p genotype test",
                   "Cramers V genotype", "p KS-test alleles", "1D-WD alleles",
                   "p KS-test LD-decay", "1D-WD LD-decay", "p GC-dist. test",
                   "Cohens D GC-dist.", "p Rogers dist. test", "Cohens D Rogers dist.")
sumstats_norm <- read.csv("../stats/summary_stats/sum_stats_normal_rec.csv")
colnames(sumstats_norm) <- rows_sumstats
t_norm <- t(sumstats_norm)
xtable(t_norm[,1:13])
xtable(t_norm[,14:25])

sumstats_high <- read.csv("../stats/summary_stats/sum_stats_high_rec.csv")
colnames(sumstats_high) <- rows_sumstats
t_high <- t(sumstats_high)
xtable(t_high[,1:13])
xtable(t_high[,14:25])

sumstats_zero <- read.csv("../stats/summary_stats/sum_stats_zero_rec.csv")
colnames(sumstats_zero) <- rows_sumstats
t_zero <- t(sumstats_zero)
xtable(t_zero[,1:13])
xtable(t_zero[,14:25])

sumstats_mean <- read.csv("../stats/summary_stats/sum_stats_mean_rec.csv")
colnames(sumstats_mean) <- rows_sumstats
t_mean <- t(sumstats_mean)
xtable(t_mean[,1:13])
xtable(t_mean[,14:25])

##population tables
pop <- read.csv("../data/sim_data/populations.csv")
pop <- pop[,c("pop", "pedigree", "size")]
colnames(pop) <- c("Population", "Pedigree", "Size")
print(xtable(pop), include.rownames=FALSE)

sim_pops <- read.csv("../sim_output/parent_sim/sim_populations.csv")
sim_parents <- strsplit(sim_pops[, "name"], "_")
sim_pops$name <- paste("(",sapply(sim_parents, "[[", 1), "x", sapply(sim_parents, "[[", 2), ")S5", sep = "")
sim_pops$size <- "200"
colnames(sim_pops) <- c("Population", "Pedigree", "Size")
sim_pops <- rbind(sim_pops,rep(NA,3))

print(xtable(cbind(sim_pops[1:55,], sim_pops[56:110,], sim_pops[111:165,])), include.rownames = FALSE)

##sort rrblup summary trait table to show best performing populations.
real_summary <- read.csv("../stats/pheno_prediction/real_summary.csv")
sim_summary <- read.csv("../stats/pheno_prediction/sim_summary.csv")
real_mean <- real_summary[,c("pop", "trait", "trait_mean")]
real_mean$trait_mean <- round(real_mean$trait_mean,2)
real_mean_wide <- reshape(
  data = real_mean,
  direction = "wide",
  idvar = "pop",
  timevar = "trait",
  v.names = "trait_mean",
  drop = "pop_trait"
)
colnames(real_mean_wide) <- c("Population", "Silk", "Tassel", "Oil", "Protein", "Starch")

sim_mean <- sim_summary[,c("pop", "trait", "trait_mean")]
sim_mean$trait_mean <- round(sim_mean$trait_mean,2)
sim_mean_wide <- reshape(
  data = sim_mean,
  direction = "wide",
  idvar = "pop",
  timevar = "trait",
  v.names = "trait_mean",
  drop = "pop_trait"
)
colnames(sim_mean_wide) <- c("Population", "Silk", "Tassel", "Oil", "Protein", "Starch")

print(xtable(cbind(real_mean_wide, sim_mean_wide[-1])), include.rownames=FALSE)

real_95 <- real_summary[,c("pop", "trait", "trait_95_perc")]
real_95$trait_95_perc <- round(real_95$trait_95_perc,2)
real_95_wide <- reshape(
  data = real_95,
  direction = "wide",
  idvar = "pop",
  timevar = "trait",
  v.names = "trait_95_perc",
  drop = "pop_trait"
)
colnames(real_95_wide) <- c("Population", "Silk", "Tassel", "Oil", "Protein", "Starch")

sim_95 <- sim_summary[,c("pop", "trait", "trait_95_perc")]
sim_95$trait_95_perc <- round(sim_95$trait_95_perc,2)
sim_95_wide <- reshape(
  data = sim_95,
  direction = "wide",
  idvar = "pop",
  timevar = "trait",
  v.names = "trait_95_perc",
  drop = "pop_trait"
)
colnames(sim_95_wide) <- c("Population", "Silk", "Tassel", "Oil", "Protein", "Starch")

print(xtable(cbind(real_95_wide, sim_95_wide[-1])), include.rownames=FALSE)


pred_results_trees <- read.csv("../stats/pheno_prediction/results/pred_results_trees_BL.csv")
pred_results_CNN <- read.csv("../stats/pheno_prediction/results/pred_results_CNN.csv")

#trait mean
trees_mean <- pred_results_trees[pred_results_trees$task == "trait_mean",]
trees_mean$model <- ifelse(trees_mean$model == "baseline", "Baseline", trees_mean$model)
CNN_mean <- pred_results_CNN[pred_results_CNN$task == "trait_mean",]
mean_CNN_res <- CNN_mean[CNN_mean$model %in% c("p1_p2","p1_p2_rate", "corr"),]
trait_order <- c("Silk", "Tassel", "Oil", "Protein", "Starch")
mean_res <- rbind(trees_mean, mean_CNN_res)
mean_res$trait <- str_to_title(mean_res$trait)
mean_res$trait <- factor(mean_res$trait, levels = trait_order)
mean_res <- mean_res[order(mean_res[,"trait"]),]
mean_res$best_predict <- gsub("(?i)FALSE", "No", mean_res$best_predict)
mean_res$best_predict <- gsub("(?i)TRUE", "Yes", mean_res$best_predict)
mean_res$model <- ifelse(mean_res$model == "p1_p2", "CNN(p1,p2)", 
                         ifelse(mean_res$model == "p1_p2_rate", "CNN(p1,p2,rate)",
                                ifelse(mean_res$model == "corr", "CNN(corr)", mean_res$model)))
mean_res <- mean_res[,!colnames(mean_res) %in% "task"]
colnames(mean_res) <- c("Model", "Trait", "RMSE", "Corr. coef.", "p Corr.", "Best predict")
print(xtable(mean_res), include.rownames = FALSE)

mean_CNN_app <- CNN_mean[!(CNN_mean$model %in% c("p1_p2","p1_p2_rate", "corr")),]
mean_CNN_app$trait <- str_to_title(mean_CNN_app$trait)
mean_CNN_app$trait <- factor(mean_CNN_app$trait, levels = trait_order)
mean_CNN_app <- mean_CNN_app[order(mean_CNN_app[,"trait"]),]
mean_CNN_app$best_predict <- gsub("(?i)FALSE", "No", mean_CNN_app$best_predict)
mean_CNN_app$best_predict <- gsub("(?i)TRUE", "Yes", mean_CNN_app$best_predict)
mean_CNN_app$corr <- ifelse(is.na(mean_CNN_app$corr) ,"NA",mean_CNN_app$corr)
mean_CNN_app$model <- ifelse(mean_CNN_app$model == "p1_p2_cM", "CNN(p1,p2,cM)", 
                         ifelse(mean_CNN_app$model == "p1_p2_position", "CNN(p1,p2,position)",
                                ifelse(mean_CNN_app$model == "p1_p2_pos_diff", "CNN(p1,p2,interval)", 
                                       ifelse(mean_CNN_app$model == "p1_p2_pos_diff_rate", "CNN(p1,p2,interval,rate)", mean_CNN_app$model))))
mean_CNN_app <- mean_CNN_app[,!colnames(mean_CNN_app) %in% "task"]
colnames(mean_CNN_app) <- c("Model", "Trait", "RMSE", "Corr. coef.", "p Corr.", "Best predict")
print(xtable(mean_CNN_app), include.rownames = FALSE)

trees_95 <- pred_results_trees[pred_results_trees$task == "trait_95_perc",]
trees_95$model <- ifelse(trees_95$model == "baseline", "Baseline", trees_mean$model)
CNN_95 <- pred_results_CNN[pred_results_CNN$task == "trait_95_perc",]
CNN_95_res <- CNN_95[CNN_95$model %in% c("p1_p2","p1_p2_rate", "corr"),]
trait_order <- c("Silk", "Tassel", "Oil", "Protein", "Starch")
res_95 <- rbind(trees_95, CNN_95_res)
res_95$trait <- str_to_title(res_95$trait)
res_95$trait <- factor(res_95$trait, levels = trait_order)
res_95 <- res_95[order(res_95[,"trait"]),]
res_95$best_predict <- gsub("(?i)FALSE", "No", res_95$best_predict)
res_95$best_predict <- gsub("(?i)TRUE", "Yes", res_95$best_predict)
res_95$model <- ifelse(res_95$model == "p1_p2", "CNN(p1,p2)", 
                         ifelse(res_95$model == "p1_p2_rate", "CNN(p1,p2,rate)",
                                ifelse(res_95$model == "corr", "CNN(corr)", res_95$model)))
res_95 <- res_95[,!colnames(res_95) %in% "task"]
colnames(res_95) <- c("Model", "Trait", "RMSE", "Corr. coef.", "p Corr.", "Best predict")
print(xtable(res_95), include.rownames = FALSE)

CNN_95_app <- CNN_95[!(CNN_95$model %in% c("p1_p2","p1_p2_rate", "corr")),]
CNN_95_app$trait <- str_to_title(CNN_95_app$trait)
CNN_95_app$trait <- factor(CNN_95_app$trait, levels = trait_order)
CNN_95_app <- CNN_95_app[order(CNN_95_app[,"trait"]),]
CNN_95_app$best_predict <- gsub("(?i)FALSE", "No", CNN_95_app$best_predict)
CNN_95_app$best_predict <- gsub("(?i)TRUE", "Yes", CNN_95_app$best_predict)
CNN_95_app$corr <- ifelse(is.na(CNN_95_app$corr) ,"NA",CNN_95_app$corr)
CNN_95_app$model <- ifelse(CNN_95_app$model == "p1_p2_cM", "CNN(p1,p2,cM)", 
                             ifelse(CNN_95_app$model == "p1_p2_position", "CNN(p1,p2,position)",
                                    ifelse(CNN_95_app$model == "p1_p2_pos_diff", "CNN(p1,p2,interval)", 
                                           ifelse(CNN_95_app$model == "p1_p2_pos_diff_rate", "CNN(p1,p2,interval,rate)", CNN_95_app$model))))
CNN_95_app <- CNN_95_app[,!colnames(CNN_95_app) %in% "task"]
colnames(CNN_95_app) <- c("Model", "Trait", "RMSE", "Corr. coef.", "p Corr.", "Best predict")
print(xtable(CNN_95_app), include.rownames = FALSE)

