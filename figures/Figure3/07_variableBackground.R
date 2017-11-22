set.seed(14651)

library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))

# Compute weighted deviation scores using gchromVAR for variable numbers of background peaks
excessive <- c(1:5,seq(10,100,5), seq(110,500,10))
lo <- lapply(excessive, function(bgn){
  bg <- getBackgroundPeaks(SE, bs = bgn)
  ukbb_wDEV <- computeWeightedDeviations(SE, ukbb, background_peaks = bg)
  zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]]))
  zscoreWeighted[,2] <- gsub("_PP001", "", zscoreWeighted[,2])
  zscoreWeighted$pvalue <- pnorm(zscoreWeighted$value, lower.tail = FALSE)
  saveRDS(zscoreWeighted, file = paste0("variableBG/variableBGpeak_", as.character(bgn), ".rds"))
  bgn
})


