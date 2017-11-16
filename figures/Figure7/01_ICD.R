library(data.table)
library(dplyr)

icd <- data.frame(fread(paste0('zcat < ', "data/ICD10all_PP5_summaryTable.txt.gz")))

lapply(unique(icd$trait), function(trait){
  tdf <- icd[icd$trait == trait, ]
  odf <- tdf[tdf$pval < 0.05/(dim(tdf)[1]),]
  odf
}) %>% rbindlist() %>% as.data.frame() -> traitSignificant