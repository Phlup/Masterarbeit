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
    labs(title = paste("LD decay of real and simulated genotypes for population",pop), x = "Distance in cM", y = expression(r^2)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.box.background = element_rect(color = "grey", linewidth = 1),
          legend.background = element_rect(fill = "white", color = NA),
          legend.position = c(0.85,0.85)) +
    scale_color_manual(name = NULL, values = c("Real genotypes" = "blue", "Simulated genotypes" = "orange"),
                       guide = guide_legend(override.aes = list(alpha = 0.7))) +
    annotate("text", x = 90, y = 0.75, label = paste("KS Test P-Value: ", ks_p),
             size = 3, color = "black") +
    annotate("text", x = 90, y = 0.70, label = paste("Wasserstein Distance: ", w1d),
             size = 3, color = "black")
  
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
    stat_compare_means(method = "t.test", comparisons = list(c("Real genotypes", "Simulated genotypes")))

  ggsave(out_path, plot, dpi = 300, width = 8, height = 6, units = "in", device = "png", bg = "white")
}
#calculate gc content in vector of chars
GC_cont <- function(seq_vec) {
  GC_count <- sum(seq_vec %in% c("G", "C"))
  return(GC_count/length(seq_vec))
}