library(data.table)
library(qqman)

icd <- fread(paste0('zcat < ', "data/ICD10all_PP5_summaryTable.txt.gz"))

png("qqPlot_ICD10.png", width = 6, height = 6)
qq(icd[["pval"]])
dev.off()