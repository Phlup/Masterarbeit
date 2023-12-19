library(openxlsx)
library(data.table)

#construct melted tables of raw and imputed original data genotypes to test data quality
pop_1_enc <- read.csv("../data/test_data/pop_1_enc.csv")
intersect(colnames(NAM_parent_genos), colnames(pop_1_genos)) -> common_SNP
common_1 <- pop_1_genos[common_SNP]

NAM_genos <- read.table("../data/NAM_map_and_genos-121025/NAM_SNP_genos_raw_20090921.txt",
                        row.names = 1)

colnames(NAM_genos) <- NAM_genos[1,]

#reshape encoded geno into long format
melted_pop_1 <- reshape(common_1, idvar = "RIL", direction = "long", varying = list(colnames(common_1[-1])), 
                       times = colnames(common_1[-1]))

#reduce original to intersection of SNPs and reshape into long format
pop_1_genos <- NAM_genos[,(NAM_genos[c("SNP_NAME"),] %in% common_1$RIL),]
pop_1_genos <- as.data.frame(t(pop_1_genos))

melted_geno <- reshape(geno_1, idvar = "SNP_NAME", direction = "long", varying = list(colnames(geno_1[-1])),
        times = colnames(geno_1[-1]))

#merge by common SNPs and RIL name
melted_pop_1$full <- rownames(melted_pop_1)
melted_geno$full <- rownames(melted_geno)

compare <- merge(melted, melted_geno, by = "full")
compare <- compare[,c("full", "RIL", "time.x", "an1.5.x", "an1.5.y")]
colnames(compare) <- c("full", "RIL", "SNP", "encoding", "original")

#get only "heterozygous" encoded genotypes
het <- compare[compare$encoding == 1,]

#just 4224 of 14061 encoded heterozygotes are real heterozygotes!
table(substr(het$original, 1,1) == substr(het$original, 2,2))




