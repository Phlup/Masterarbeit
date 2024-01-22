library(openxlsx)
library(data.table)

##clean genotype and phenotype data from studies on maize NAM population
#Genotypes from Yu et al. 2008 https://academic.oup.com/genetics/article/178/1/539/6062286?login=true,
#Phenotypes from Bukler et al. 2009 https://pubmed.ncbi.nlm.nih.gov/19661422/ and
#Cook et al. 2012 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3271770/

##Genotypes Yu et al. 2008
#read NAM genos extract parents and b73 allele
NAM_genos <- read.table("../data/NAM_map_and_genos-121025/NAM_SNP_genos_raw_20090921.txt",
                        row.names = 1)

#only parent genos
NAM_parent_genos <- as.data.frame(t(NAM_genos[,c(1:27)]))
colnames(NAM_parent_genos)[1] <- "RIL"

write.csv(NAM_parent_genos, "../data/sim_data/NAM_parent_genos.csv", row.names = FALSE)

#B73 allele
B73_allele <- data.table("SNP" = row.names(NAM_genos), "B73" = paste(substring(NAM_genos[,1], 1, 1), substring(NAM_genos[,1], 3, 3),
                                                                     sep = ""))
write.csv(B73_allele[-1], "../data/sim_data/B73_alleles.csv", row.names = FALSE)

#get population genotypes
populations <- read.xlsx("../data/NAM_map_and_genos-121025/NAM_populations.xlsx")

populations$name <- paste("B73_", populations$parent, sep = "")

genos_raw <- read.table("../data/NAM_map_and_genos-121025/NAM_SNP_genos_raw_20090921.txt",
                        row.names = 1)

genos_raw <- t(genos_raw)

genos_raw <- genos_raw[-c(1:27),]
start = 1
size = NA
for (i in 1:(length(genos_raw[,1])-1)){
  if(substring(genos_raw[i,1],0,4) != substring(genos_raw[i+1,1],0,4)){
    end = i
    pop = sub("^0+", "", substring(genos_raw[i,1],3,4))
    genos <- genos_raw[c(start:end),]
    genos[,1] <- sub("^0+", "", substr(genos[,1], nchar(genos[1,1])-2, nchar(genos[1,1])))
    colnames(genos)[1] <- "individual"
    write.csv(genos, paste("../data/NAM_genotype_data/geno_encoding/pop_", pop,
                                              "_genos.csv", sep = ""), row.names = FALSE)
    size <- c(size, length(start:end))
    start = i+1
  }
}
pop = substring(genos_raw[start,1],3,4)
genos <- genos_raw[c(start:i+1),]
genos[,1] <- sub("^0+", "", substr(genos[,1], nchar(genos[1,1])-2, nchar(genos[1,1])))
colnames(genos)[1] <- "individual"
write.csv(genos, paste("../data/NAM_genotype_data/geno_encoding/pop_", pop,
                                          "_genos.csv", sep = ""), row.names = FALSE)
size <- c(size, length(start:i+1))
populations$size <- size[-1]
write.csv(populations, "../data/sim_data/populations.csv", row.names = FALSE)

#NAM parent genos in additive encoding
populations <- read.csv("../data/sim_data/populations.csv")
NAM_parent_genos <- read.csv("../data/sim_data/NAM_parent_genos.csv")

populations <- rbind(populations, c(0, NA, "B73", NA, NA))

NAM_add <- NAM_parent_genos
NAM_add <- merge(populations[,c("pop", "parent")], NAM_add, by.x = "parent", by.y = "RIL", all = TRUE)
NAM_add <- NAM_add[!(is.na(NAM_add$pop)),]
hets <- apply(NAM_add[3:length(NAM_add[1,])], c(1,2), function(x){ifelse(substr(x,1,1)==substr(x,3,3),NA,0)})
homos <- t(apply(NAM_add[,3:length(NAM_add[1,])], 1, function(x){
  ifelse(NAM_add[1,3:length(NAM_add[1,])] == x, 1, -1)}))
homos[!is.na(hets)] <- hets[!is.na(hets)]
NAM_add[,3:length(NAM_add[1,])] <- homos

write.csv(NAM_add, "../data/sim_data/NAM_parent_add.csv", row.names = FALSE)

##Phenotype data

#extract target vectors for offspring populations and phenotype means
populations <- read.csv("../data/sim_data/populations.csv")

#Bukler et al. 2009
silk <- read.table("../data/Buckler_etal_2009_Science_flowering_time_data-090807/NAM_DaysToSilk.txt",
                   header = TRUE, fill = TRUE)

silk <- silk[,c("loc_name", "year", "pop", "accename", "value")]
silk <- silk[grepl("Z0", silk$accename),]
silk$value <- as.numeric(silk$value)
silk <- silk[order(silk$accename),]
silk$env <- paste(silk$loc_name, silk$year)
silk$env[silk$env == "Aurora, NY 2006"] <- "06A"
silk$env[silk$env == "Aurora, NY 2007"] <- "07A"
silk$env[silk$env == "Clayton, NC 2007"] <- "07CL"
silk$env[silk$env == "Columbia, MO 2006"] <- "06MO"
silk$env[silk$env == "Columbia, MO 2007"] <- "07MO"
silk$env[silk$env == "Urbana, IL 2006"] <- "06IL"
silk$env[silk$env == "Urbana, IL 2007"] <- "07IL"
silk <- silk[silk$pop %in% populations$pop,!(colnames(silk) %in% c("loc_name", "year"))]
silk$accename <- sub("^0+", "", substr(silk$accename, nchar(silk$accename)-2, nchar(silk$accename)))
silk <- silk[!(is.na(silk$value)),]
silk <- silk[!(duplicated(paste(silk$pop, silk$accename, silk$env))),]
silk$trait <- "silk"

tassel <- read.table("../data/Buckler_etal_2009_Science_flowering_time_data-090807/NAM_DaysToTassel.txt",
                     header = TRUE, fill = TRUE)

tassel <- tassel[,c("loc_name", "year", "pop", "accename", "value")]
tassel <- tassel[grepl("Z0", tassel$accename),]
tassel$value <- as.numeric(tassel$value)
tassel <- tassel[order(tassel$accename),]
tassel$env <- paste(tassel$loc_name, tassel$year)
tassel$env[tassel$env == "Aurora, NY 2006"] <- "06A"
tassel$env[tassel$env == "Aurora, NY 2007"] <- "07A"
tassel$env[tassel$env == "Clayton, NC 2007"] <- "07CL"
tassel$env[tassel$env == "Columbia, MO 2006"] <- "06MO"
tassel$env[tassel$env == "Columbia, MO 2007"] <- "07MO"
tassel$env[tassel$env == "Urbana, IL 2006"] <- "06IL"
tassel$env[tassel$env == "Urbana, IL 2007"] <- "07IL"
tassel <- tassel[tassel$pop %in% populations$pop,!(colnames(tassel) %in% c("loc_name", "year"))]
tassel$accename <- sub("^0+", "", substr(tassel$accename, nchar(tassel$accename)-2, nchar(tassel$accename)))
tassel <- tassel[!(is.na(tassel$value)),]
tassel <- tassel[!(duplicated(paste(tassel$pop, tassel$accename, tassel$env))),]
tassel$trait <- "tassel"

#Cook et al. 2012
spo <- read.xlsx("../data/supp_pp.111.185033_185033Supplemental_Data_S1.xlsx")
spo <- spo[spo$Rep == 1,]
spo <- spo[,c("pop", "Geno_code", "Starch_db", "Protein_db", "Oil_db", "env")]
spo <- spo[!(is.na(spo$Geno_code)),]
spo <- spo[order(spo$Geno_code),]
spo$accename <- sub("^0+", "", substr(spo$Geno_code, nchar(spo$Geno_code)-2, nchar(spo$Geno_code)))
spo <- spo[spo$pop %in% populations$pop,!(colnames(spo) %in% c("Geno_code"))]

spo_long <- reshape(
  spo,
  direction = "long",
  varying = list(c("Starch_db", "Protein_db", "Oil_db")),
  v.names = c("value"),
  timevar = "trait",
  times = c("starch", "protein", "oil"),
  idvar = c("pop", "env", "accename")
)
spo_long <- spo_long[colnames(tassel)]

NAM_phenotypes <- rbind(silk, tassel, spo_long)
colnames(NAM_phenotypes) <- c("pop", "individual", "value", "env", "trait")

mean_phenotypes <- do.call(data.frame,aggregate(value ~ pop + env + trait, data = NAM_phenotypes, 
                                       FUN = function(x) c(mean = mean(x), max = max(x), min = min(x), length = length(x))))
colnames(result) <- c("pop", "env", "trait", "trait_mean", "trait_max", "trait_min", "N") 

write.csv(NAM_phenotypes, "../data/NAM_phenotype_data/NAM_phenotypes.csv", row.names = FALSE)
write.csv(mean_phenotypes, "../data/NAM_phenotype_data/mean_phenotypes.csv", row.names = FALSE)


##----
##genotpyes encoded (unexpected encoding, see test_scripts/test_imputed_data.R)
##reduce nam offspring genotypes to subset
##load NAM genos imputed
#offspring_genos <- read.xlsx("../data/NAM_map_and_genos-121025/NAM_genos_imputed_20090807.xlsx")
#
##generate indices for 10 each of first 5 RILs 
#index <- rep(1:10) + rep(c(0,100,200,300,400), each = 10) + 2
#
##reduce offspring_genos
#genos_reduce <- offspring_genos[index,c(1:50)]
#colnames(genos_reduce) <- offspring_genos[2,c(1:50)]
#genos_reduce[,c(2:50)] <- sapply(genos_reduce[,c(2:50)], as.numeric)
#
#genos_reduce[genos_reduce == 1.5] <- 1.0
#
#genos_reduce[genos_reduce == 0.5] <- 1.0
#
#write.csv(genos_reduce, "../data/test_data/genos_reduce.csv", row.names = FALSE)
#
#
##get pop 1 genotypes
#pop_1_genos <- offspring_genos[c(3:196),]
#colnames(pop_1_genos) <- offspring_genos[2,]
#pop_1_genos[,c(2:1107)] <- sapply(pop_1_genos[,c(2:1107)], as.numeric)
##pop_1_genos[pop_1_genos == 1.5] <- 1.0
#pop_1_genos[pop_1_genos == 0.5] <- 1.0
#write.csv(pop_1_genos, "../data/test_data/pop_1_enc.csv", row.names = FALSE)
#
##get pop 2 genotypes
#pop_2_genos <- offspring_genos[c(197:392),]
#colnames(pop_2_genos) <- offspring_genos[2,]
#pop_2_genos[,c(2:1107)] <- sapply(pop_2_genos[,c(2:1107)], as.numeric)
#write.csv(pop_2_genos, "../data/test_data/pop_2_enc.csv", row.names = FALSE)
#
##get pop 3 genotypes
#pop_3_genos <- offspring_genos[c(393:582),]
#colnames(pop_3_genos) <- offspring_genos[2,]
#pop_3_genos[,c(2:1107)] <- sapply(pop_3_genos[,c(2:1107)], as.numeric)
#write.csv(pop_3_genos, "../data/test_data/pop_3_enc.csv", row.names = FALSE)
