library(chromVAR)
library(chromVARxx)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

# Build (weighted) matches
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE))
assays(ukbb) <- list(weights = assays(ukbb)[["weights"]], matches = assays(ukbb)[["weights"]] > 0)

# Compute (weighted) deviation scores
ukbb_DEV <- computeDeviations(SE, ukbb)
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)

# Make tables and export
zscoreCHROMVAR <- melt(t(assays(ukbb_DEV)[["z"]]))
zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]]))
write.table(zscoreCHROMVAR, "../../data/bulk/GWAS-Bulk/bulkHeme_chromVAR_zscores.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(zscoreWeighted, "../../data/bulk/GWAS-Bulk/bulkHeme_weighted_zscores.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

