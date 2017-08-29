library(BuenColors)
library(ggplot2)
library(plotly)

df <- read.table("../../data/bulk/GWAS-Bulk/compare3.tsv", header = TRUE)               
df$weighted_pvalue <- pnorm(df$weighted_Zscore, lower.tail = FALSE)
df$chromVAR_pvalue <- pnorm(df$chromVAR_Zscore, lower.tail = FALSE)

p1 <- ggplot(df, aes(x = -log10(chromVAR_pvalue), y = -log10(weighted_pvalue), color = Celltype, value = Trait)) +
  pretty_plot() + geom_point() + scale_color_manual(values = ejc_color_maps) +
  coord_cartesian(xlim=c(0,7), ylim=c(0,7))

#ggplotly(p1)


cor(-log10(df$ldscore_pvalue), -log10(df$weighted_pvalue))
cor(-log10(df$ldscore_pvalue), -log10(df$chromVAR_pvalue))
cor(-log10(df$chromVAR_pvalue), -log10(df$weighted_pvalue))
