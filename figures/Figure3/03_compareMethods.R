library(BuenColors)
library(ggplot2)
library(plotly)
library(qvalue)

df <- read.table("../../data/bulk/GWAS-Bulk/compare3.tsv", header = TRUE)               
df$weighted_pvalue <- pnorm(df$weighted_Zscore, lower.tail = FALSE)
df$chromVAR_pvalue <- pnorm(df$chromVAR_Zscore, lower.tail = FALSE)
gs <- readRDS("../Figure1/OR_Heme.rds")
df_gs <- as.data.frame(gs,stringsAsFactors=F)
names(df_gs) <- c("i","j","cell","trait","obs","perm","z")
df_gs$pvalue <- 2*pnorm(-abs(as.numeric(df_gs$z)))
df$goShifter_pvalue <- df_gs$pvalue

p1 <- ggplot(df[df$Trait=="MPV",], aes(x = -log10(weighted_pvalue), y = -log10(chromVAR_pvalue), color = Celltype, value = Trait)) +
  pretty_plot() + geom_point() + scale_color_manual(values = ejc_color_maps) +
  coord_cartesian(xlim=c(0,7), ylim=c(0,7))

ggplotly(p1)

cor(-log10(df$ldscore_pvalue), -log10(df$weighted_pvalue))
cor(-log10(df$ldscore_pvalue), -log10(df$chromVAR_pvalue))
cor(-log10(df$chromVAR_pvalue), -log10(df$weighted_pvalue))
cor(-log10(df$goShifter_pvalue), -log10(df$weighted_pvalue))
cor(-log10(df$goShifter_pvalue), -log10(df$ldscore_pvalue))

df$outlier <- abs(-log10(df$ldscore_pvalue) - -log10(df$weighted_pvalue))

pdf("../Figure3/methods_qq.pdf")
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
pvector <- OR.df$pvalue
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector)))
points(e,o,pch=20,col="green")
abline(-log10(0.05/288),0,col="grey")
dev.off()

tab <- table(pred=cut(df$weighted_pvalue,c(0,0.05/288,1)),ldsr=cut(df$ldscore_pvalue,c(0,0.05/288,1)))
recall = tab[1,1] / sum(tab[,1]); recall # sensitivity / recall
specificity = tab[2,2] / sum(tab[,2]); specificity # specificity
accuracy = (tab[1,1] + tab[2,2]) / sum(tab); accuracy # accuracy
precision = (tab[1,1]) / sum(tab[1,]); precision # PPV / precision
f1 = 2 * precision * recall / (precision + recall); f1 

df$uniq <- paste0(df$Celltype,"_",df$Trait)
df.ldsr <- df %>%
  dplyr::filter(ldscore_pvalue < 0.05/288) %>%
  .$uniq
df.weighted_pvalue <- df %>%
  dplyr::filter(weighted_pvalue < 0.05/288) %>%
  .$uniq
df.goShifter <- df %>%
  dplyr::filter(goShifter_pvalue < 0.05/288) %>%
  .$uniq
df_venn <- Venn(list(goShifter=df.goShifter,LDSR=df.ldsr,chromVAR=df.weighted_pvalue))  
pdf("../Figure3/methods_Venn.pdf")
  plot(df_venn)
dev.off()
