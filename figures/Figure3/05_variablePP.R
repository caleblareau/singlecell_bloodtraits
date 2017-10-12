 library(chromVAR)
library(chromVARxx)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)
 library(BuenColors)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = ".bed"))
bgp <- getBackgroundPeaks(SE)

# Make variable PP thresholds
makePPthreshdf <- function(val){
  print(val)
  ukbbcp <- ukbb
  assays(ukbbcp)$weights[assays(ukbbcp)$weights < val] <- 0
  zw<- melt(t(assays( computeWeightedDeviations(SE, ukbbcp, background_peaks = bgp))[["z"]]))
  zw$PP <- val
  return(zw)
}

PPs <- c(1:10/1000, 2:9/10, 1)
listres <- lapply(PPs, makePPthreshdf)

zscoreWeighted <- rbindlist(listres)
saveRDS(zscoreWeighted, "../../data/supplement_rds/varyingPPcutoff.rds")


