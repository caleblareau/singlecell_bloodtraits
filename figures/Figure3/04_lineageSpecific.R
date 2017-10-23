library(BuenColors)
library(ggplot2)

set.seed(14651)


# Define lineage specific attributes
Ery <- c("HSC", "MPP", "CMP", "MEP", "Ery")
Meg <- c("HSC", "MPP", "CMP", "MEP", "Mega")
Mye <- c("HSC", "MPP", "CMP", "LMPP", "GMP-A", "GMP-B", "GMP-C", "Mono", "mDC")
Lymph <- c("HSC", "MPP", "LMPP", "CLP", "NK", "pDC", "CD4", "CD8", "B")

make2df <- function(trait, celltypes){
  return(data.frame(Celltype = celltypes, Trait = rep(trait, length(celltypes))))
}

lineageSpecificDF <- rbind(
  make2df("BASO_COUNT", Mye),
  make2df("EO_COUNT", Mye),
  make2df("NEUTRO_COUNT", Mye),
  make2df("MONO_COUNT", Mye),
  make2df("WBC_COUNT", unique(c(Mye, Lymph))),
  make2df("LYMPH_COUNT", Lymph),
  make2df("PLT_COUNT", Meg),
  make2df("MPV", Meg),
  
  make2df("MEAN_RETIC_VOL", Ery),
  make2df("HCT", Ery),
  make2df("HGB", Ery),
  make2df("MCH", Ery),
  make2df("MCV", Ery),
  make2df("MCHC", Ery),
  make2df("RBC_COUNT", Ery),
  make2df("RETIC_COUNT", Ery)
)

# Return T/F whether rows in df1 are in df2
rowCheck <- function(df1, df2){
  xx <- apply(df1, 1, paste, collapse = "_")
  yy <- apply(df2, 1, paste, collapse = "_")
  return(xx %in% yy)
}

df <- read.table("../../data/bulk/GWAS-Bulk/compare3.tsv", header = TRUE)               
df$weighted_pvalue <- pnorm(df$weighted_Zscore, lower.tail = FALSE)
df$chromVAR_pvalue <- pnorm(df$chromVAR_Zscore, lower.tail = FALSE)
gs <- readRDS("../Figure1/OR_Heme.rds")
df_gs <- as.data.frame(gs,stringsAsFactors=F)
names(df_gs) <- c("i","j","cell","trait","obs","perm","z")
df_gs$pvalue <- 2*pnorm(-abs(as.numeric(df_gs$z)))
df$goShifter_pvalue <- df_gs$pvalue

df$hit_ldscore <- df$ldscore_pvalue < 0.05/288
df$hit_chromVAR <- df$chromVAR_pvalue < 0.05/288
df$hit_weighted <- df$weighted_pvalue < 0.05/288
df$hit_goShifter <- df$goShifter_pvalue < 0.05/288

df$lineageSpecific <- rowCheck(df[,c("Celltype", "Trait")], lineageSpecificDF)
odf <- df[,c(1,2,5:12)] 
#write.table(odf, file = "../../supplemental_tables/Bulk_Enrichments.tsv", 
#            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

permuted <- sapply(1:10000, function(i) sum(1:dim(df)[1] * sample(df$lineageSpecific, length(df$lineageSpecific))))
ldscore_ranksum <- sum(1:dim(df)[1]*df[order(df$ldscore_pvalue, decreasing = FALSE), "lineageSpecific"])
chromVAR_ranksum <- sum(1:dim(df)[1]*df[order(df$chromVAR_pvalue, decreasing = FALSE), "lineageSpecific"])
weighted_ranksum <- sum(1:dim(df)[1]*df[order(df$weighted_pvalue, decreasing = FALSE), "lineageSpecific"])
goShifter_ranksum <- sum(1:dim(df)[1]*df[order(df$goShifter_pvalue, decreasing = FALSE), "lineageSpecific"])

allstats <- data.frame(
  Yes = c(sum(df$lineageSpecific& df$hit_ldscore), sum(df$lineageSpecific& df$hit_chromVAR), sum(df$lineageSpecific& df$hit_weighted), sum(df$lineageSpecific& df$hit_goShifter)),
  No = c(sum(!df$lineageSpecific& df$hit_ldscore), sum(!df$lineageSpecific& df$hit_chromVAR), sum(!df$lineageSpecific& df$hit_weighted), sum(!df$lineageSpecific& df$hit_goShifter)),
  pval = c(pnorm((mean(permuted) - ldscore_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - chromVAR_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - weighted_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - goShifter_ranksum)/sd(permuted), lower.tail = FALSE)
  )
)
allstats$method <- c("LDscore", "chromVAR", "g-chromVAR", "goShifter")

ggplot(allstats[,c(3,4),drop=FALSE], aes(x = method, y = -1*log10(pval))) +
  geom_histogram(stat = "identity", color = "black", fill = c("blue", "red", "green", "pink")) +
  labs(x = "", y = "Lineage-Specific Enrichment -log10(p)") + coord_flip() +
  pretty_plot()

tfdf <- reshape2::melt(allstats[,c(1,2,4)], id.vars = "method")
tfdf$variable <- factor(as.character(tfdf$variable), c("No", "Yes"))
ggplot(tfdf, aes(x = method, y = value, fill = variable)) +
  geom_histogram(stat = "identity", color = "black") +
  labs(x = "", y = "Number of Significant Enrichments", fill = "Lineage Specific") +
  scale_fill_manual(values = c("firebrick","green4")) +
  pretty_plot() +  theme(legend.position = "bottom")


