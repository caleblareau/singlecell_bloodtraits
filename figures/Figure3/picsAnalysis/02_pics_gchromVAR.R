library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)
library(BuenColors)
library(dplyr)
library(matrixStats)
set.seed(14651)

# Create bulk Summarized Experiment
peaksdf <- fread("../../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
pics <- importBedScore(rowRanges(SE), list.files("pics_bedscore/", full.names = TRUE, pattern = "*.bed$"))
bg <- getBackgroundPeaks(SE)
dev <- computeWeightedDeviations(SE, pics, background_peaks = bg)
outdf <- melt(t(assays(dev)[["z"]]))
outdf$gchromVAR_pvalue <- pnorm(outdf$value, lower.tail = FALSE)
colnames(outdf) <- c("Celltype", "Trait", "Z", "pvalue")
fdf <- outdf[order(outdf$pvalue), ]
outdf[outdf$Trait == "Alzheimers_combined", ]

write.table(fdf, file = "pics_gchromVAR_dec30.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
