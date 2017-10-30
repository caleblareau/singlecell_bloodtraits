library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(dplyr)
library(chromVAR)
library(chromVARxx)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BuenColors)

SE <- readRDS("../../data/singlecell/scATAC/full_scHeme_QC.rds")

# Make one count summed over the cell type
typeCountVec <- table(colData(SE)$type)
mat <- sapply(names(typeCountVec), function(type) {
  cmat <- rowSums(assays(SE)[["counts"]][,which(colData(SE)$type == type)])
  cmat
})

sumCounts <- SummarizedExperiment::SummarizedExperiment(
  assays = list("counts" = mat),
  colData = S4Vectors::DataFrame(names = names(typeCountVec)),
  rowData = rowRanges(SE)
)

# Compute weighted deviation scores
sumCounts <- addGCBias(sumCounts, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(sumCounts), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))
ukbb_wDEV <- computeWeightedDeviations(sumCounts, ukbb)
zmat <- (t(assays(ukbb_wDEV)[["z"]]))
rownames(zmat) <- c("CLP", "CMP", "GMP", "GMP-A", "GMP-B", "GMP-C", "HSC", "LMPP", "MEP", "Mono", "MPP", "pDC")
zscoreWeighted <- melt(zmat)
zscoreWeighted[,2] <- gsub("_PP001", "", zscoreWeighted[,2])
sumSingleCells <- zscoreWeighted
colnames(sumSingleCells) <- c("Celltype", "Trait", "ZscoreSum")

# Import the true bulk deviations
wchromvar <- read.table("../../data/bulk/GWAS-Bulk/bulkHeme_weighted_zscores.txt", stringsAsFactors = FALSE)
colnames(wchromvar) <- c("Celltype", "Trait", "ZscoreBulk")

mergedf <- merge(sumSingleCells, wchromvar)

write.table(mergedf, file = "bulk_singleCellSum.gchromVAR.txt", row.names = FALSE,
            col.names = TRUE, sep = "\t", quote = FALSE)

cor(mergedf$ZscoreSum, mergedf$ZscoreBulk)


# Note: I tried log10p but the zscore makes more senese / looks better / is more correlated
p1 <- ggplot(shuf(mergedf), aes(x = ZscoreSum, y = ZscoreBulk, color = Celltype)) +
  labs(x = "Sum of Single Cells Z-score", y = "Bulk Z-score", color = "") +
  geom_point() +
  scale_color_manual(values = ejc_color_maps) +
  pretty_plot() + theme(legend.position = "bottom") + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + xlim(-6,10) + ylim(-6,10)
p1
ggsave(p1, file = "bulk_singleCellSum.gchromVAR.pdf", width = 4, height = 5)




