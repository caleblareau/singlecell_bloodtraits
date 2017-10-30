library(data.table)
library(GenomicRanges)
library(chromVARxx)
library(stringr)
library(diffloop)
library(SummarizedExperiment)
library(magrittr)
library(BuenColors)
"%ni%" <- Negate("%in%")

# Import peak / mutation data / background peaks
peaksMutCount <- readRDS("data/peaks_by_muts.rds")
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
bg <- data.matrix(data.frame(fread(paste0("zcat < ", "../../data/bulk/ATAC/backgroundPeaks_500bp.txt.gz"))))
cpt <- sweep(peaksMutCount, 2, Matrix::colSums(peaksMutCount), FUN="/") * 1000

# Import FM variants
files <- list.files("../../data/UKBB_BC_PP001/",
                                         full.names = TRUE, pattern = "*.bed")
fileshorts <- list.files("../../data/UKBB_BC_PP001/",
                                         full.names = FALSE, pattern = "*.bed")

# Compute permuted peaksMutCount matrix; not exactly working yet
totalMat <- sapply(1:16, function(k){
  x <- data.frame(fread(files[k]))
  FM <- makeGRangesFromDataFrame(x[x$V5 > 0.1,], seqnames.field = "V1", start.field = "V2", end.field = "V3")
  ov_var <- findOverlaps(peaks, FM)
  Zs <- sapply(1:dim(peaksMutCount)[2], function(i){
    perm <- peaksMutCount[bg[1:length(peaks) %in% queryHits(ov_var), ],i]
    obs <- mean(peaksMutCount[1:length(peaks) %in% queryHits(ov_var),i])
    means <- mean(perm)
    Z <- (obs - mean(perm))/sd(perm)
    Z
  })
  Zs
})
totalMat
