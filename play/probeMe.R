library(chromVAR)
library(chromVARxx)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(diffloop)
library(dplyr)

peaksdf <- fread("../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))

findPairsInEnhancer <-  function(file){
  gr <- bedToGRanges(file)  
  gr <- gr[gr$score > 0.5]
  kpeak <- which(countOverlaps(peaks, gr) > 1)
  varinfo <- lapply(kpeak, function(idx){
    df <- data.frame(gr)[subjectHits(findOverlaps(peaks[idx], gr)),c(1,2,3,7)]
    df2 <- (cbind(df[1,], df[2,]))
    colnames(df2) <- c("chr1", "start1", "end1", "PP1", "chr2", "start2", "end2", "PP2")
    df2
  }) %>% rbindlist() %>% as.data.frame()
  dfout <- cbind(data.frame(peaks)[kpeak,], data.frame(counts)[kpeak,], varinfo, file)
  return(dfout)
}

files <- c("../data/UKBB_BC_PP001/MCH_PP001.bed", "../data/UKBB_BC_PP001/MCV_PP001.bed", "../data/UKBB_BC_PP001/MEAN_RETIC_VOL_PP001.bed",
           "../data/UKBB_BC_PP001/MPV_PP001.bed", "../data/UKBB_BC_PP001/PLT_COUNT_PP001.bed", "../data/UKBB_BC_PP001/RBC_COUNT_PP001.bed",
           "../data/UKBB_BC_PP001/RETIC_COUNT_PP001.bed")

allDoubles <- lapply(files, findPairsInEnhancer)  %>% rbindlist() %>% as.data.frame()
write.table(allDoubles, file = "doubleVariantSingleEnhancer.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)


