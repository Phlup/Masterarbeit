library(adegenet)

data(microsatt)
x <- as.genpop(microsatt$tab)

listDist <- lapply(1:5, function(i) cailliez(dist.genpop(obj,met=i)))
for(i in 1:5) {attr(listDist[[i]],"Labels") <- popNames(obj)}
listPco <- lapply(listDist, dudi.pco,scannf=FALSE)

par(mfrow=c(2,3))
for(i in 1:5) {scatter(listPco[[i]],sub=paste("Dist:", i))}

#allele counts
x@tab
#allele cardinality
nloc <- length(levels(x@loc.fac))
#allele factors/alleles per locus
loc.fac <- x@loc.fac
#calculate allele frequencies
X <- makefreq(x, missing = "mean", quiet = TRUE)
#population cardinality
nlig <- nrow(X)

x <- matrix(c(0,1,0,1,2,0,2,0,1,1,0,1), nrow = 3)

calculate_snp_frequency <- function(matrix) {
  num_individuals <- nrow(matrix)
  num_snps <- ncol(matrix)
  
  snp_frequencies <- matrix(0, nrow = num_snps, ncol = 3)
  
  for (j in 1:num_snps) {
    snp_counts <- table(matrix[, j])
    
    if (0 %in% names(snp_counts)) {
      snp_count_0 <- snp_counts[["0"]]
    } else {
      snp_count_0 <- 0
    }
    
    if (1 %in% names(snp_counts)) {
      snp_count_1 <- snp_counts[["1"]]
    } else {
      snp_count_1 <- 0
    }
    
    if (2 %in% names(snp_counts)) {
      snp_count_2 <- snp_counts[["2"]]
    } else {
      snp_count_2 <- 0
    }
    
    snp_frequency_0 <- snp_count_0 / num_individuals
    snp_frequency_1 <- snp_count_1 / num_individuals
    snp_frequency_2 <- snp_count_2 / num_individuals
    
    snp_frequencies[j, ] <- c(snp_frequency_0, snp_frequency_1, snp_frequency_2)
  }
  
  return(snp_frequencies)
}

example_matrix <- matrix(c(0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1), nrow = 3, ncol = 4, byrow = TRUE)

snp_frequencies <- calculate_snp_frequency(example_matrix)

print(snp_frequencies)

#neis distance
d <- X %*% t(X)
vec <- sqrt(diag(d))
d <- d/vec[col(d)]
d <- d/vec[row(d)]
d <- -log(d)
d <- as.dist(d)

#rogers distance
kX <- lapply(split(X, loc.fac[col(X)]), matrix, nrow = nlig)
dcano <- function(mat) {
  daux <- mat %*% t(mat)
  vec <- diag(daux)
  daux <- -2 * daux + vec[col(daux)] + vec[row(daux)]
  diag(daux) <- 0
  daux <- sqrt(0.5 * daux)
  return(daux)
}
d <- matrix(0, nlig, nlig)
for (i in 1:length(kX)) {
  d <- d + dcano(kX[[i]])
}
d <- d/length(kX)
d <- as.dist(d)

