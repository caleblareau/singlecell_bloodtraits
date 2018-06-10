library(dplyr)
library(diffloop)
library(GenomicRanges)

regions <- bedToGRanges("../../data/revision/Pleiotropic_hits.txt")

traits <- gsub("_PP001.bed", "", list.files("../../data/UKBB_BC_PP001/", patter = ".bed$"))
lapply(traits, function(trait){
  t <- read.table(paste0("../../data/UKBB_BC_PP001/betas_added/",trait,"_PP001_betas.bed"))[,c("V1", "V2", "V3", "V5", "V8")]
  t <- t[t$V5 > 0.1, ]
  t$trait <- trait
  colnames(t) <- c("chr", "start", "stop", "PP", "Z", "trait")
  gr_t <- makeGRangesFromDataFrame(t)
  ov <- findOverlaps(gr_t, regions)
  return(t[queryHits(ov), ])
}) %>% rbindlist () %>% as.data.frame() -> SNPdf

SNPdf %>% arrange(chr, start)
