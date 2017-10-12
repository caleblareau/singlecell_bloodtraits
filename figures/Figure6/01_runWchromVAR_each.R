library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(chromVAR)
library(chromVARxx)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BuenColors)

set.seed(14651)

# Run weighted deviations
SE <- readRDS("../../data/singlecell/scATAC/full_scHeme_QC.rds")
pcaPlus <- data.table::fread("../../data/singlecell/scATAC/scHeme_meta.txt")

# Function to compute the g-chromVAR scores per each population separately
makeEachTable <- function(celltype, SE, pcaPlus){
  keepidx <- which(pcaPlus$type == celltype)
  SE <- filterPeaks(SE[,keepidx])
  pcaPlus <- pcaPlus[pcaPlus$type == celltype, ]
  SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
  ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001", pattern =  "*.bed", full.names = TRUE))
  ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
  zscoreRaw <- t(assays(ukbb_wDEV)[["z"]])
  
  # Make nice output tables
  bigdf <- data.frame(pcaPlus, data.matrix(zscoreRaw))
  colnames(bigdf) <- gsub("_PP001", "", colnames(bigdf))
  bigdf <- shuf(bigdf)
  write.table(bigdf, file = paste0("data/", celltype, "_weightedchromVAReach.txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
makeEachTable("HSC", SE, pcaPlus)
makeEachTable("MEP", SE, pcaPlus)
makeEachTable("CMP", SE, pcaPlus)
