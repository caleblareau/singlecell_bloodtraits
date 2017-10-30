library(data.table)
library(GenomicRanges)
library(dplyr)
library(diffloop)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)

# ---------------------------------------------
# Get background peaks
# ---------------------------------------------

# Import genomic coordinates
peaks_g <- makeGRangesFromDataFrame(fread("29August2017_EJCsamples_allReads_500bp.bed"),
                                    seqnames.field = "V1", start.field = "V2", end.field = "V3")

atac <- data.matrix(fread("29August2017_EJCsamples_allReads_500bp.counts.txt"))

SE <- SummarizedExperiment::SummarizedExperiment(
  rowRanges = peaks_g,
  colData = data.frame(names = colnames(atac)), 
  assays = list(counts = atac)
)
SE <- addGCBias(SE, BSgenome.Hsapiens.UCSC.hg19)
bg <- getBackgroundPeaks(SE)
write.table(bg, file = "backgroundPeaks_500bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
