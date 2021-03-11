##############################################
##### Function to plot the panel dataset #####
##############################################
library(ggplot2)
library(data.table)
library(dplyr)

plot_df <- function(df, genes) {
  ## Factor the data so that it plots the x-axis in an order
  df$x_axis <- paste0(df$gene_name, " (n = " , df$total_exons, " exons; avg = ", df$avg_exon_size, ")")
  df$x_axis <- factor(df$x_axis, levels = unique(df$x_axis[order(df$total_nucleotides, decreasing = TRUE)]))
  
  ## Generate scatter plot --> Exons (y) per gene (x)
  p1 <- ggplot(data = df) + geom_point(aes(x = x_axis, y = exon_size, color = transcript_num)) + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 75, hjust = 1)) + 
    labs(x = "Gene", y = paste0("Exon Size"), color="# of Transcripts") + scale_y_log10() + 
    scale_color_gradient(low = "grey", high = "red")
  
  ## Show proportion of panel space taken up by each gene 
  p2 <- ggplot(data = df) + geom_line(aes(x = x_axis, y = proportion_total, group = 1), color = "Blue") + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + 
    labs(y = paste0("Proportion of", "\n", "Total (%)", "\n")) + 
    annotate(geom = "text", x = 40, y = 5, label = paste0("Total = ", sum(df$exon_size), " nucleotides")) + 
    annotate(geom = "text", x = 40, y = 3, label = paste0("Gene num = ", length(unique(df$gene_id)), "/", length(genes))) +
    geom_vline(xintercept = 9.5, linetype = 2, color = "red") 
  
  p <- ggarrange(p2, p1, nrow=2, ncol=1, common.legend = TRUE, legend = "right", heights = c(0.5, 2))
  return(p)
  
}