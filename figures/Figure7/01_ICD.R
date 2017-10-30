library(data.table)
library(qqman)

icd <- fread(paste0('zcat < ', "data/ICD10all_PP5_summaryTable.txt.gz"))

pdf("qqPlot_ICD10.pdf", width = 6, height = 6)
qq(icd[["pval"]])
dev.off()