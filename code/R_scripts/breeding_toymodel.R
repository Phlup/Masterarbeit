# a simple toy computation for the effect of recombination on the distribution of traits in the offspring

#let a genome be a sequence of length 100 
#where each site is given a value n_i and the sum of the n_i is the genomes trait

num_sites = 100
num_genomes = 4

num_seeds = 500

set.seed(42)

genotype = matrix(rnorm(num_genomes * num_sites),num_genomes,num_sites)

genometraits = round(rowMeans(genotype) * num_sites,2)


pdf(paste0("parental_genomes.pdf"), width = 14)

par(mfrow = c(1,2))

plot(c(),xlim = c(0,num_sites), ylim = c(-3,3), ylab = "trait effect at genome position", xlab = "position along the genome", type = "p")
for(i in 1:num_genomes){
  points( genotype[i,], col = i+1, pch = i)  
}

plot(c(),xlim = c(0,num_sites), ylim = c(-15,15), ylab = "cumulative trait effect from position 1 to k", xlab = "position k")
for(i in 1:num_genomes){
  lines( cumsum(genotype[i,]), col = i+1)  
}
legend("topleft", legend = c(paste("genome", 1:4, "total trait effect:", genometraits)), fill = 1:4 + 1, cex = 1, xpd = TRUE, bty = "n")


dev.off()


#now breed them pairwise
parental_sets = combn(num_genomes,2)


# free recombination
childrens = list()
for (pind in 1:ncol(parental_sets)){
  #generate num_seeds children
  child = matrix(0,num_seeds,num_sites)
  for (j in 1:num_seeds){
    r1 = rbinom(num_sites, 1, 0.5) + 1
    child[j,] = genotype[parental_sets[,pind],] [cbind(r1,1:num_sites)]  
  }
  childrens[[pind]] = child
}


pdf(paste0("hybrids_recomb_high.pdf"), width = 14)
par(mfrow = c(2,3))

for (pind in 1:ncol(parental_sets)){
  
  main = paste("Trait of",num_seeds,"hybrids of genome", parental_sets[,pind][1], "and", parental_sets[,pind][2], "-- free recombination" )
  
  plot(c(),xlim = c(0,num_sites), ylim = c(-20,20), ylab = "cumulative trait up to pos k -- (total trait at pos 100)", xlab = "position k", main = main)
  for(i in 1:num_seeds){
    myrgb = col2rgb(i+1)
    lines( cumsum(childrens[[pind]][i,]), col = rgb(myrgb[1],myrgb[2],myrgb[3], max= 255, alpha = 75) )
  }
  lines( cumsum(genotype[parental_sets[,pind][1],]), col = parental_sets[,pind][1] + 1, lwd =8 )
  lines( cumsum(genotype[parental_sets[,pind][2],]), col = parental_sets[,pind][2] + 1, lwd =8 )
  
  # compute mean trait
  traits = rowSums(childrens[[pind]])
  
  text(x = 10, paste("var:", round(var(traits),2 )), adj = 0  )
  text(x = 12, paste("mean:", round(mean(traits),2 )), adj = 0  )
  text(x = 14, paste("max:", round(max(traits),2 )), adj = 0  )

}

dev.off()



# 10 recombination events on average
for (meanrecomb in c(1,5,10)){
  childrens = list()
  for (pind in 1:ncol(parental_sets)){
    #generate num_seeds children
    child = matrix(0,num_seeds,num_sites)
    for (j in 1:num_seeds){
      num_recombs = rpois(1,meanrecomb - 1) + 1
      positions = sort(sample(num_sites,num_recombs))
      r1 = rep( (seq( c(0,positions) ) + rbinom(1,1,0.5)) %% 2 + 1, c(positions,num_sites) - c(0,positions))
      child[j,] = genotype[parental_sets[,pind],] [cbind(r1,1:num_sites)]  
    }
    childrens[[pind]] = child
  }
  
  
  pdf(paste0("hybrids_recomb_", meanrecomb , ".pdf"), width = 14)
  
  par(mfrow = c(2,3))
  
  for (pind in 1:ncol(parental_sets)){
    
    main = paste("Trait of",num_seeds,"hybrids of genome", parental_sets[,pind][1], "and", parental_sets[,pind][2], "-- ~", meanrecomb, "unif recombs" )
    
    plot(c(),xlim = c(0,num_sites), ylim = c(-20,20), ylab = "cumulative trait up to pos k -- (total trait at pos 100)", xlab = "position k", main = main)
    for(i in 1:num_seeds){
      myrgb = col2rgb(i+1)
      lines( cumsum(childrens[[pind]][i,]), col = rgb(myrgb[1],myrgb[2],myrgb[3], max= 255, alpha = 75) )   
    }
    lines( cumsum(genotype[parental_sets[,pind][1],]), col = parental_sets[,pind][1] + 1, lwd =8 )
    lines( cumsum(genotype[parental_sets[,pind][2],]), col = parental_sets[,pind][2] + 1, lwd =8 )
    
    # compute mean trait
    traits = rowSums(childrens[[pind]])
    
    text(x = 10, paste("var:", round(var(traits),2 )), adj = 0  )
    text(x = 12, paste("mean:", round(mean(traits),2 )), adj = 0  )
    text(x = 14, paste("max:", round(max(traits),2 )), adj = 0  )
    
  }
  
  dev.off()
}






