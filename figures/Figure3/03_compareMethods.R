library(BuenColors)
library(ggplot2)
library(plotly)
library(qvalues)

df <- read.table("../../data/bulk/GWAS-Bulk/compare3.tsv", header = TRUE)               
df$weighted_pvalue <- pnorm(df$weighted_Zscore, lower.tail = FALSE)
df$chromVAR_pvalue <- pnorm(df$chromVAR_Zscore, lower.tail = FALSE)

p1 <- ggplot(df[df$Trait=="MPV",], aes(x = -log10(weighted_pvalue), y = -log10(chromVAR_pvalue), color = Celltype, value = Trait)) +
  pretty_plot() + geom_point() + scale_color_manual(values = ejc_color_maps) +
  coord_cartesian(xlim=c(0,7), ylim=c(0,7))

ggplotly(p1)

cor(-log10(df$ldscore_pvalue), -log10(df$weighted_pvalue))
cor(-log10(df$ldscore_pvalue), -log10(df$chromVAR_pvalue))
cor(-log10(df$chromVAR_pvalue), -log10(df$weighted_pvalue))

df$outlier <- abs(-log10(df$ldscore_pvalue) - -log10(df$weighted_pvalue))

pvector <- df$weighted_pvalue
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector)))
plot(e,o,pch=20); abline(0,1,col="red")
pvector <- df$ldscore_pvalue
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector)))
points(e,o,pch=20,col="red")
pvector <- df$chromVAR_pvalue
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector)))
points(e,o,pch=20,col="purple")
abline(-log10(0.05/288),0,col="grey")

tab <- table(pred=cut(df$chromVAR_pvalue,c(0,0.05/288,1)),ldsr=cut(df$ldscore_pvalue,c(0,0.05/288,1)))
recall = tab[1,1] / sum(tab[,1]); recall # sensitivity / recall
specificity = tab[2,2] / sum(tab[,2]); specificity # specificity
accuracy = (tab[1,1] + tab[2,2]) / sum(tab); accuracy # accuracy
precision = (tab[1,1]) / sum(tab[1,]); precision # PPV / precision
f1 = 2 * precision * recall / (precision + recall); f1 



