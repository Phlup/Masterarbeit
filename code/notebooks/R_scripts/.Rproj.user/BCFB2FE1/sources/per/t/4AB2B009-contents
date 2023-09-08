library(openxlsx)

#genotpyes
#load NAM genos imputed
NAM_genos <- read.xlsx("../data/NAM_map_and_genos-121025/NAM_genos_imputed_20090807.xlsx")

#generate indices for 10 each of first 5 RILs 
index <- rep(1:10) + rep(c(0,100,200,300,400), each = 10) + 2

#reduce NAM_genos
genos_reduce <- NAM_genos[index,c(1:50)]
colnames(genos_reduce) <- NAM_genos[2,c(1:50)]
genos_reduce[,c(2:50)] <- sapply(genos_reduce[,c(2:50)], as.numeric)

genos_reduce[genos_reduce == 1.5] <- 1.0

genos_reduce[genos_reduce == 0.5] <- 1.0


write.csv(genos_reduce, "../data/test_data/genos_reduce.csv", row.names = FALSE)
