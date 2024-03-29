library(ggplot2)
library(ggpubr)
#implemented here are functions used to generate the validation stats in validation_stats.R (helps declutter script)
#ld decay plot, rogers distance plot, GC content

#ld decay plot
plot_ld_decay <- function(real, sim, ks_p, w1d, pop, out_path){
  plot <- ggplot() +
    geom_point(data = real, aes(d, r2, color = "Real genotypes"), size = 1.5, alpha = 0.05) +
    geom_point(data = sim, aes(d, r2, color = "Simulated genotypes"), size = 1.5, alpha = 0.05) +
    geom_smooth(data = real, aes(d, r2), method = "loess", se = FALSE, color = "darkblue", span = 0.2) +
    geom_smooth(data = sim, aes(d, r2), method = "loess", se = FALSE, color = "darkorange", span = 0.2) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 1)) +
    labs(x = "Distance in cM", y = expression(r^2)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.box.background = element_rect(color = "grey", linewidth = 1),
          legend.background = element_rect(fill = "white", color = NA),
          legend.position = c(0.85,0.85),
          text = element_text(size = 18)) +
    scale_color_manual(name = NULL, values = c("Real genotypes" = "blue", "Simulated genotypes" = "orange"),
                       guide = guide_legend(override.aes = list(alpha = 0.7))) +
    annotate("text", x = 90, y = 0.75, label = paste("KS Test P-Value: ", ks_p),
             size = 4, color = "black") +
    annotate("text", x = 90, y = 0.70, label = paste("Wasserstein Distance: ", w1d),
             size = 4, color = "black")
  
  ggsave(out_path, plot, dpi = 300, width = 8, height = 6, units = "in", device = "png", bg = "white")
}

#rogers distance plot
plot_rogers_dist <- function(real, sim, out_path) {
  df <- data.frame(
    Group = c(rep("Real genotypes", length(real)), rep("Simulated genotypes", length(sim))),
    Value = c(real, sim)
  )
  plot <- ggboxplot(df, x = "Group", y = "Value", color = "Group", palette = c("blue", "orange")) +
    labs(x = "", y = "Rogers distance") + 
    theme(text = element_text(size=18))+
    stat_compare_means(method = "t.test", comparisons = list(c("Real genotypes", "Simulated genotypes")))

  ggsave(out_path, plot, dpi = 300, width = 8, height = 6, units = "in", device = "png", bg = "white")
}
#calculate gc content in vector of chars
GC_cont <- function(seq_vec) {
  GC_count <- sum(seq_vec %in% c("G", "C"))
  return(GC_count/length(seq_vec))
}

#build feature dataframe for rolling correlation of parental traits for i populations
calc_trait_cor <- function(traits, populations, window_size, genmap){
  trait_cors <- data.frame(NULL)
  trait_feat_num <- floor(length(genmap$Marker)/window_size)
  for(i in populations$pop){
    parents <- populations[populations$pop == i, c("parent_1", "parent_2")]
    pop_i_cors <- t(rollapply(t(as.data.frame(traits[traits$parent %in% parents, !colnames(traits) %in% c("pop", "parent")])),
                              width = window_size, function(x) cor(x[,1],x[,2]), by.column = FALSE, by = window_size))
    colnames(pop_i_cors) <- sapply(1:trait_feat_num, function(a) paste((a - 1) * window_size + 1, a * window_size, sep = ":"))
    pop_i_cors <- cbind(data.frame("pop" = i), pop_i_cors)
    trait_cors <- rbind(trait_cors, pop_i_cors)
  }
  return(trait_cors)
}
