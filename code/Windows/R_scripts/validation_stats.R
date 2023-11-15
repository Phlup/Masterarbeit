library(sommer)
library(stats)
source("nuc_diversity.R")

#Stats implemented to validate simulation
#Table of simulated alleles in diploid gene sequence format (python)
#Heterozygous/reference/derived alleles (python)
#Histogram/density plot of additive encoding based on reference allele (python)
#linkage disequilibrium decay (r2~d)
#ks test/density of r2 between sim and real encoding
#distance measure (e.g. euclidean/rogers distance)
#nucleotide diversity (tbd)
#others? -> measures mostly focused on popgenetic hypothesis, e.g. neutral vs non-neutral evolution tajimas D

map <- read.csv("../data/test_data/B73_genmap.csv")
#keep Locus/Marker, Position and LG/Chromosome
map <- map[,c("Marker", "Map.cM.", "Chromosome")]
colnames(map) <- c("Locus","Position","LG")

#data(DT_cpdata)
#CPgeno <- GT_cpdata
#CPgeno[1:5,1:5]
#
#mapCP <- MP_cpdata; head(mapCP)
#names(mapCP) <- c("Locus","Position","LG")
#
#mapCP <- mapCP[which(mapCP$LG <= 3),]
#
#res <- LD.decay(CPgeno, mapCP)
#names(res)
#
##constrain to significant LG with p < 0.001
#res$all.LG <- res$all.LG[which(res$all.LG$p < .001),]

plot_LD_decay <- function(est, title){
  with(est$all.LG, plot(r2~d,col=transp("cadetblue"),
                        xlim=c(0,55), ylim=c(0,1), 
                        pch=20,cex=0.5,ylab=expression(r^2),
                        xlab="Distance in cM",main=title)
  )
}

#LD decay with recombination
whole_add <- read.csv("../data/sim_output/whole_add.csv")
whole_add <- whole_add[-1]

res_rec <- LD.decay(whole_add, map)
res_rec$all.LG <- res_rec$all.LG[(res_rec$all.LG$p < .001),]
plot_LD_decay(res_rec, "R2~D in recombination simulation")

#LD decay without recombination
#no_add <- read.csv("../data/sim_output/no_add.csv")
#no_add <- no_add[-1]
#
#res_norec <- LD.decay(no_add, map)
#res_norec$all.LG <- res_norec$all.LG[(res_norec$all.LG$p < .001),]
#plot_LD_decay(res_norec)

#LD decay in real population
add_1 <- read.csv("../data/sim_output/add_1.csv")
add_1 <- add_1[-1]

res_pop1 <- LD.decay(add_1, map)
res_pop1$all.LG <- res_pop1$all.LG[(res_pop1$all.LG$p < .001),]
plot_LD_decay(res_pop1, "R2~D in real population")

#density plot of r2 of sim and real pop
plot(density(res_rec[["all.LG"]][["r2"]]))
plot(density(res_pop1[["all.LG"]][["r2"]]))


#ks test on simulated and real pop
ks.test(res_rec[["all.LG"]][["r2"]], res_pop1[["all.LG"]][["r2"]])
#large sample detects small divergence in distribution

#ks test on simulated and real pop
ks.test(sample(res_rec[["all.LG"]][["r2"]], 1000), sample(res_pop1[["all.LG"]][["r2"]]), 1000)
#smaller sample non-significant for alpha = 0.05
plot(density(sample(res_rec[["all.LG"]][["r2"]], 1000)))
plot(density(sample(res_pop1[["all.LG"]][["r2"]], 1000)))


#density of no rec and ks test between rec and no rec
#plot(density(res_norec[["all.LG"]][["r2"]]))
#ks.test(sample(res_rec[["all.LG"]][["r2"]], 100), sample(res_norec[["all.LG"]][["r2"]]), 100)


#distance measure (e.g. euclidian distance/rogers distance)
par(mfrow=c(1,2))
boxplot(dist(whole_add, method = "manhattan"), main = "Euclidean distance of additive encoding simulated population", ylim=c(200,1200))
boxplot(dist(add_1, method = "manhattan"), main = "Euclidean distance of additive encoding real population", ylim=c(200,1200))

#nucleotide diversity (overall and across sites)
rec_geno <- read.csv("../data/sim_output/whole_sim.csv")
rec_geno <- rec_geno[-c(1,2)]
rec_list <- apply(rec_geno, 1, function(row) paste(row, collapse = ""))
rec_list <- as.list(rec_list)

#sequence frequency?
duplicated(rec_geno)
#all dissimilar
seq_freq <- rep(1/194, 194)
rec_pd <- calculate_pairwise_differences(rec_list)
calculate_nucleotide_diversity(seq_freq, rec_pd, 194)

pop_1_geno <- read.csv("../data/test_data/pop_1_genos.csv")
pop_1_geno <- pop_1_geno[-1]

pop_1_list <- apply(pop_1_geno, 1, function(row) paste(row, collapse = ""))
pop_1_list <- as.list(pop_1_list)

#sequence frequency?
duplicated(pop_1_geno)
#all dissimilar
seq_freq <- rep(1/194, 194)
pop_1_pd <- calculate_pairwise_differences(pop_1_list)
calculate_nucleotide_diversity(seq_freq, pop_1_pd, 194)


