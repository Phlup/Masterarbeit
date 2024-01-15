library(AlphaSimR)

#Creating Founder Haplotypes
dir.create("/scratch/playground/federico/offspring_simulation/2_simulated_data/AlphaSim/Barley_DRR/B.1_simulated_genotypes")

pop_info<-read.csv("/scratch/playground/federico/offspring_simulation/1_real_data/Barley_DRR/A.0_pop.info.csv")
populations<-pop_info[,1]
populations <- populations[c(1,10,35)]

parental_geno<-read.csv("/scratch/playground/federico/offspring_simulation/1_real_data/Barley_DRR/A.1_cleaned_parental_data.csv")
row.names(parental_geno)<-parental_geno[,1]; parental_geno<-parental_geno[,-1]
parental_geno<-as.data.frame(t(parental_geno))
for (c in 1:ncol(parental_geno)){parental_geno[,c]<-as.numeric(as.character(parental_geno[,c]))}
row.names(parental_geno)<-gsub("W23829.803911", "W23829/803911", row.names(parental_geno))
row.names(parental_geno)<-gsub("Unumli.Arpa", "Unumli-Arpa", row.names(parental_geno))

# consensus_map<-read.csv("/scratch/playground/federico/offspring_simulation/1_real_data//Barley_DRR/sources/consensus_map.csv")
# consensus_map<-consensus_map[,-3]
# consensus_map$cM<-as.numeric(as.character(consensus_map$cM))

pop_maps<-readRDS("/scratch/playground/federico/offspring_simulation/1_real_data//Barley_DRR/sources/DRR_populations_genetic_maps.RDS")

#parameters
n_families<-150

#Recombination model
#The remainder (proportion 1 – p) come from a process following the chi-square model. The chi-square model has a single parameter, m, which is a non-negative integer.
#The chi-square model is a special case of the gamma model which has a positive parameter ν; 
#The chi-square model corresponds to the case that m = ν – 1 is an integer.
#So, m:10==v:11 and v:2.6==m:1.6

#p the proportion of crossovers coming from a non-interfering pathway, Default = 0.
Ps <- c(0, 0.125, 0.25, 0.5, 1) 
#Ps <- Ps[1]

#v the crossover interference parameter for a gamma model of recombination. 1 = NO interference (Haldane). Default = 2.6 (Kosambi)
#Vs<-c(0.000001,0.00001, 0.0001,0.001, 0.01, 0.1, 1, 10)  
Vs<-c(0.0001,0.001, 0.01, 0.1, 1, 10)  

for (P in Ps){ cat(P);cat(":")
  
if (P == 1){Vs<-0}

for (V in Vs){ cat(V); cat("-")

simulated_population_list<-list()

for (p in populations){ cat(p); cat("-")

#get parental geno
parents<-as.character(pop_info[which(pop_info[,1]==p),2:3])
pop_parental_geno<-parental_geno[parents,]

#take out NAs
to.delete<-unique(c(which(is.na(pop_parental_geno[1,])), which(is.na(pop_parental_geno[2,]))))
pop_parental_geno<-pop_parental_geno[,-to.delete]

#only use genetic map 
#genetic_map_pop<-consensus_map[which(consensus_map[,1]%in%colnames(pop_parental_geno)),]

#specific pop map
pop_map<-pop_maps[[p]][[1]][,2:4]
for (c in 2:7){pop_map<-rbind(pop_map, pop_maps[[p]][[c]][,2:4])}
pop_map<-as.data.frame(pop_map)
pop_map$cM<-as.numeric(as.character(pop_map$cM))
genetic_map_pop<-pop_map[which(pop_map[,1]%in%colnames(pop_parental_geno)),]

#G0 -It is important to consider that allele coding -1,0,1 is converted to 0,1,2. 
G0 <- importInbredGeno(geno = pop_parental_geno, genMap = genetic_map_pop)
SP <- SimParam$new(G0) #SimParam
SP$p <- P
SP$v <- V
SP$segSites
G0 <- newPop(G0, simParam = SP)
G0@id <- parents

#F1
biparental_cross <- matrix(parents, nrow=1, ncol=2)
F1 <- makeCross(G0, crossPlan = biparental_cross, simParam = SP, nProgeny = 1)
F1@id <- "F1"

#F2
selfing <- matrix(c("F1","F1"), nrow=n_families, ncol=2)
F2 <- makeCross(F1, crossPlan = selfing, simParam = SP, nProgeny = n_families)
F2@id <- paste("F2","-",1:n_families, sep = "")

#F3
selfing <- cbind(F2@id, F2@id)
F3 <- makeCross(F2, crossPlan = selfing, simParam = SP, nProgeny = 1)
F3@id <- paste("F3","-",1:n_families, sep = "")

#F4
selfing <- cbind(F3@id, F3@id)
F4 <- makeCross(F3, crossPlan = selfing, simParam = SP, nProgeny = 1)
F4@id <- paste("F4","-",1:n_families, sep = "")

#F5
selfing <- cbind(F4@id, F4@id)
F5 <- makeCross(F4, crossPlan = selfing, simParam = SP, nProgeny = 1)
F5@id <- paste("F5","-",1:n_families, sep = "")

#convert back
target_G_geno <- pullMarkerGeno(F5, markers = genetic_map_pop[,1])
row.names(target_G_geno)<-paste("RIL","-",1:n_families, sep = "")
#write.csv(target_G_geno, paste("/scratch/playground/federico/offspring_simulation/2_simulated_data/AlphaSim/Barley_DRR/B.1_simulated_genotypes/simulated_",p,".csv", sep = ""))
simulated_population_list[[p]]<-target_G_geno

}#p

saveRDS(simulated_population_list, paste("/scratch/playground/federico/offspring_simulation/2_simulated_data/AlphaSim/Barley_DRR/B.1_simulated_genotypes/simulated_pops_p=",P,"_v=",V,".RDS", sep = ""))

}#P
}#V
