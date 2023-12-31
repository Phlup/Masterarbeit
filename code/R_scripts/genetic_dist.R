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

