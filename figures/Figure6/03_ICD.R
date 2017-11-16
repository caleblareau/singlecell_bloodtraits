library(data.table)
library(dplyr)
library(BuenColors)

icd <- data.frame(fread(paste0('zcat < ', "../../data/UKBB_BC_PP001/ICD10all_PP5_summaryTable.txt.gz")))

lapply(unique(icd$trait), function(trait){
  tdf <- icd[icd$trait == trait, ]
  odf <- tdf[tdf$pval < 0.05/(dim(tdf)[1]),]
  odf
}) %>% rbindlist() %>% as.data.frame() -> traitSignificant

df<- data.frame(
  rank = 1:length(table(unique(traitSignificant[,c("chr", "pos", "V2")])$V2)),
  counts = sort(table(unique(traitSignificant[,c("chr", "pos", "V2")])$V2), decreasing = TRUE)
)
df$label <- "Other"
keep <- c(1,2,6,7,10)
df$label[keep] <- gsub("Diagnoses - main ICD10: ", "", df$counts.Var1)[keep]

p1 <- ggplot(df, aes(x = rank, y= counts.Freq)) +
  geom_bar(width = 0.7, aes(fill = label), colour="black", stat = "identity", position = position_dodge(width=0.7)) +
  labs(fill = "ICD 10 Code: ", x = "Traits", y = "# PheWAS Variants Identified from PP > 0.5 Blood Variants")+
  pretty_plot() + scale_fill_manual(values = c("dodgerblue", "green4", "firebrick", "orange", "purple4", "black")) +
  theme(legend.position = "bottom")
cowplot::ggsave(p1, file = "ICD10_phewas_hightlights.pdf", height = 6, width = 8)
