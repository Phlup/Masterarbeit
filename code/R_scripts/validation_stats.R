library(sommer)
library(stats)
library(transport)
library(effectsize)
library(adegenet)
source("stat_functions.R")

##generate summary statistics for genotype simulation across all populations and recombination scenarios

#Stats implemented:

#1. distributional differences
#Table of simulated alleles in diploid gene sequence format
#Heterozygous/reference/derived alleles
#chisq test between tables (hetero+all)
#ks test + wasserstein distance for dist diff
#gc content distribution across all individuals+sites

#2. additional popgen stats
#linkage disequilibrium decay (r2~d)
#ks test/density of r2 between sim and real encoding
#distance measure (e.g. euclidean/rogers distance)

#read genmap to constrain markers/calc ld decay etc.
genmap <- read.csv("../data/sim_data/B73_genmap.csv")
#read populations
populations <- read.csv("../data/sim_data/populations.csv")

#LD decay
#keep Locus/Marker, Position and LG/Chromosome
ld_map <- genmap[,c("Marker", "Map.cM.", "Chromosome")]
colnames(ld_map) <- c("Locus","Position","LG")
#vary over all recombination parameters
rec_param <- c("normal_rec", "high_rec", "zero_rec", "mean_rec")
for(i in rec_param){
  sum_stats <- data.frame("population" = populations$pop, "het_p" = NA, "het_phi" = NA, "geno_p" = NA,
                          "geno_cramersV" = NA, "ks_p" = NA, "w1d" = NA, "ld_ks_p" = NA, "ld_w1d" = NA,
                          "gc_p" = NA, "gc_cohens" = NA, "rogers_p" = NA, "rogers_cohens" = NA)
  #vary over all simulated pops
  for(j in populations$pop){
    #read real and sim pop
    sim_geno <- read.csv(paste("../sim_output/",i,"/geno_encoding/geno_",j,".csv", sep = ""))
    real_geno <- read.csv(paste("../data/NAM_genotype_data/geno_encoding/pop_",j,"_genos.csv", sep = ""))
    sim_add <- read.csv(paste("../sim_output/",i,"/additive_encoding/add_",j,".csv", sep = ""))
    real_add <- read.csv(paste("../data/NAM_genotype_data/additive_encoding/pop_",j,"_add.csv", sep = ""))
    
    #constrain to same markers
    sim_geno <- sim_geno[,genmap$Marker]
    real_geno <- real_geno[,genmap$Marker]
    sim_add <- sim_add[,genmap$Marker]
    real_add <- real_add[,genmap$Marker]

    ##1. dist stats
    #calculate # of het
    het_mat <- as.matrix(cbind(table(unlist(sim_add) == 0), table(unlist(real_add) == 0)))
    
    #calc # of all genos
    geno_names <- intersect(names(table(unlist(sim_geno))), names(table(unlist(real_geno))))
    geno_mat <- as.matrix(cbind(table(unlist(sim_geno))[geno_names],
                                table(unlist(real_geno))[geno_names]))
    
    #chisq.tests
    het_p <- chisq.test(het_mat)$p.value
    het_phi <- phi(het_mat)$phi_adjusted
    geno_p <- chisq.test(geno_mat)$p.value
    geno_cramersv <- cramers_v(geno_mat)$Cramers_v_adjusted
    #ks test
    ks_p <- ks.test(rowSums(sim_add), rowSums(real_add))$p.value
    #wasserstein dist
    w1d <- wasserstein1d(rowSums(sim_add), rowSums(real_add))
    
    #gc content
    haplo_mat_s <- cbind(apply(sim_geno,2,substr,1,1), apply(sim_geno,2,substr,2,2))
    gc_vec_s <- apply(haplo_mat_s, 1, GC_cont)
    haplo_mat_r <- cbind(apply(real_geno,2,substr,1,1), apply(real_geno,2,substr,2,2))
    gc_vec_r <- apply(haplo_mat_r, 1, GC_cont)
    
    gc_p <- t.test(gc_vec_s, gc_vec_r)$p.value
    gc_cohens <- cohens_d(gc_vec_s, gc_vec_r)$Cohens_d
    
    ##2. popgen stats
    #ld decay
    sim_ld <- LD.decay(sim_add, ld_map)
    sim_ld$all.LG <- sim_ld$all.LG[(sim_ld$all.LG$p < .001),]
    sim_r2_D <- sim_ld$all.LG[c("r2", "d")]
    
    real_ld <- LD.decay(real_add, ld_map)
    real_ld$all.LG <- real_ld$all.LG[(real_ld$all.LG$p < .001),]
    real_r2_D <- real_ld$all.LG[c("r2", "d")]

    #ks test + w1d on r2 measure
    ld_ks_p <- ks.test(sim_r2_D$r2, real_r2_D$r2)$p.value
    ld_w1d <- wasserstein1d(sim_r2_D$r2, real_r2_D$r2)
    #remove small number of "perfect" correlating sites for smoothing
    real_r2_D <- real_r2_D[real_r2_D$r2 < .999,]
    
    #rogers distance
    #adegenet does not like the original SNP names
    colnames(sim_add) <- c(1:974)
    genind <- df2genind(sim_add, sep = "")
    dist_sim <- as.numeric(dist.genpop(as.genpop(genind$tab), method = 4))
    colnames(real_add) <- c(1:974)
    genind <- df2genind(real_add, sep = "")
    dist_real <- as.numeric(dist.genpop(as.genpop(genind$tab), method = 4))
    
    #rogers dist test
    rogers_p <- t.test(dist_sim, dist_real)$p.value
    rogers_cohens <- cohens_d(dist_sim, dist_real)$Cohens_d
    
    plot_rogers_dist(real = dist_real, sim = dist_sim, 
                     out_path = paste("../plots/popgen_plots/rogers_dist/",i,"/rogers_dist_",j,".png", sep = ""))
    
    plot_ld_decay(real = real_r2_D, sim = sim_r2_D, ks_p = ifelse(ld_ks_p < 0.001, "<0.001", round(ld_ks_p,3)),
                  w1d = round(ld_w1d,2),
                  out_path = paste("../plots/popgen_plots/LD_decay/",i,"/ld_decay_",j,".png", sep = ""))
    
    #format sumstats
    p_vals <- c(het_p, geno_p, ks_p, ld_ks_p, gc_p, rogers_p)
    p_vals <- ifelse(p_vals < 0.001, "<0.001", round(p_vals,3))
    effects <- round(c(het_phi, geno_cramersv, w1d, ld_w1d, gc_cohens, rogers_cohens),2)
    
    stats <- c(p_vals[1], effects[1], p_vals[2], effects[2], p_vals[3], effects[3], p_vals[4], effects[4],
               p_vals[5], effects[5], p_vals[6], effects[6])
    sum_stats[sum_stats$population == j, c("het_p", "het_phi", "geno_p", "geno_cramersV", "ks_p", "w1d",
                                           "ld_ks_p", "ld_w1d", "gc_p", "gc_cohens", "rogers_p",
                                           "rogers_cohens")] <- stats
    
  }  
  print("finished")
  write.csv(sum_stats, paste("../stats/sum_stats_",i,".csv", sep = ""), row.names = FALSE)
}

#read in sum stats to compare rec scenarios
sum_results <- data.frame("rec_param" = c("normal_rec", "high_rec", "zero_rec", "mean_rec"),
                          "het_sig" = NA, "mean_sd_phi" = NA, "geno_sig" = NA, "mean_sd_cramersV" = NA,
                          "ks_sig" = NA, "mean_sd_w1d" = NA, "ld_ks_sig" = NA, "mean_sd_ld_w1d" = NA,
                          "gc_sig" = NA, "mean_sd_gc" = NA, "rogers_sig" = NA, "rogers_cohens" = NA)
for(i in rec_param){
  sum_stats <- read.csv(paste("../stats/sum_stats_",i,".csv", sep = ""))
  het_sig <- paste(table(sum_stats$het_p == "<0.001")[["TRUE"]],"/",
                   sum(table(sum_stats$het_p)), sep = "")
  het_phi <- paste("median: ", median(sum_stats$het_phi), ", IQR: ", IQR(sum_stats$het_phi),
                   ", mean: ", signif(mean(sum_stats$het_phi),3),
                   " ± ", signif(sd(sum_stats$het_phi),3), sep = "")
  geno_sig <- paste(table(sum_stats$geno_p == "<0.001")[["TRUE"]],"/",
                   sum(table(sum_stats$geno_p)), sep = "")
  geno_cramersV <- paste("median: ", median(sum_stats$geno_cramersV), ", IQR: ", IQR(sum_stats$geno_cramersV),
                   ", mean: ", signif(mean(sum_stats$geno_cramersV),3),
                   " ± ", signif(sd(sum_stats$geno_cramersV),3), sep = "")
  ks_sig <- paste(table(sum_stats$ks_p == "<0.001")[["TRUE"]],"/",
                  sum(table(sum_stats$ks_p)), sep = "")
  w1d <- paste("median: ", median(sum_stats$w1d), ", IQR: ", IQR(sum_stats$w1d),
                   ", mean: ", signif(mean(sum_stats$w1d),3),
                   " ± ", signif(sd(sum_stats$w1d),3), sep = "")
  ld_ks_sig <- paste(table(sum_stats$ld_ks_p == "<0.001")[["TRUE"]],"/",
                  sum(table(sum_stats$ld_ks_p)), sep = "")
  ld_w1d <- paste("median: ", median(sum_stats$ld_w1d), ", IQR: ", IQR(sum_stats$ld_w1d),
                  ", mean: ", signif(mean(sum_stats$ld_w1d),3),
                  " ± ", signif(sd(sum_stats$ld_w1d),3), sep = "")
  gc_sig <- paste(table(sum_stats$gc_p == "<0.001")[["TRUE"]],"/",
                     sum(table(sum_stats$gc_p)), sep = "")
  gc_cohens <- paste("median: ", median(sum_stats$gc_cohens), ", IQR: ", IQR(sum_stats$gc_cohens),
                  ", mean: ", signif(mean(sum_stats$gc_cohens),3),
                  " ± ", signif(sd(sum_stats$gc_cohens),3), sep = "")
  rogers_sig <- paste(table(sum_stats$rogers_p == "<0.001")[["TRUE"]],"/",
                  sum(table(sum_stats$rogers_p)), sep = "")
  rogers_cohens <- paste("median: ", median(sum_stats$rogers_cohens), ", IQR: ", IQR(sum_stats$rogers_cohens),
                     ", mean: ", signif(mean(sum_stats$rogers_cohens),3),
                     " ± ", signif(sd(sum_stats$rogers_cohens),3), sep = "")
  results <- c(het_sig, het_phi, geno_sig, geno_cramersV, ks_sig, w1d, ld_ks_sig, ld_w1d, gc_sig, gc_cohens,
               rogers_sig, rogers_cohens)
  sum_results[sum_results$rec_param == i, c("het_sig", "mean_sd_phi", "geno_sig", "mean_sd_cramersV",
                                            "ks_sig", "mean_sd_w1d", "ld_ks_sig", "mean_sd_ld_w1d",
                                            "gc_sig", "mean_sd_gc", "rogers_sig", "rogers_cohens")] <- results
}
write.csv(sum_results, "../stats/sum_results.csv", row.names = FALSE)
