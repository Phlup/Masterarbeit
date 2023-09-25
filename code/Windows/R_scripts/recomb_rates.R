library(openxlsx)
library(data.table)

#calculate recombination rates (cM/Mb) for Maize NAM population (B73)
#read in cM and chr SNP data
NAM_map <- read.xlsx("../data/NAM_map_and_genos-121025/NAM_map_20080419_xlsx.xlsx")
SNP_pos <- read.xlsx("../data/NAM_map_and_genos-121025/NAM_1144SNPs_AGPv1_positions.xlsx")
SNP_pos <- SNP_pos[!((SNP_pos$Position == "NULL") | (SNP_pos$Position == "multiple")),]

#save SNP positions
SNP_pos <- SNP_pos[,c("Marker", "Chr", "Position")]
colnames(SNP_pos) <- c("Marker", "Chromosome", "Position(bp)")
write.csv(SNP_pos, "../data/test_data/NAM_SNP_positions.csv", row.names = FALSE)

#merge on SNP id
gen_map <- merge(NAM_map[,c("marker", "ch", "cumulative")], 
                 SNP_pos[,c("Marker", "Position(bp)")], by.x = "marker", by.y = "Marker")

gen_map$Position <- as.numeric(gen_map$Position)

#sort by chr then position
gen_map <- gen_map[order(gen_map$ch, gen_map$Position),]

gen_map <- data.table(gen_map)
gen_map[, order := ((shift(cumulative, n = -1) - cumulative) >= 0), by = ch]

#manually remove "misplaced SNPs"
gen_map <- gen_map[!(marker %in% c("PZA01072.1", "PZA00545.26", "PZA00963.3", "PZA01960.1",
                                   "zb7.2", "PHM1184.26", "PHM2438.28", "PZA03227.1"))]

#remove SNPs where cumulative cM order != SNP position order
gen_map <- gen_map[order == TRUE]

#(i+1-i)/(j+1-j) where i is cM and j is Mb, transform to Mb, group by chr
gen_map[, pos_diff := (shift(Position, n = -1) - Position), by = ch]
gen_map[, cM_diff := (shift(cumulative, n = -1) - cumulative), by = ch]
gen_map[, rate := cM_diff/(pos_diff/1000000)]

#impute rates for last SNPs per chr
gen_map[,rate := replace(rate, is.na(rate), median(rate, na.rm = TRUE)), by = ch]

gen_map <- gen_map[,c("ch", "Position", "rate", "cumulative")]
colnames(gen_map) <- c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
#gen_map$Chromosome <- as.character(gen_map$Chromosome)
#gen_map$Chromosome <- paste("chr", gen_map$Chromosome, sep = "")

#save chr1 and complete set
genmap_chr1 <- gen_map[Chromosome == 1]
write.csv(genmap_chr1, "../data/test_data/B73_genmap_chr1.csv", row.names = FALSE)

write.csv(gen_map, "../data/test_data/B73_genmap.csv", row.names = FALSE)


