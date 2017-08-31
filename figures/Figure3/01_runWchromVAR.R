library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(chromVAR)
library(chromVARxx)
library(BSgenome.Hsapiens.UCSC.hg19)
source("00_smooth.R")

SE <- readRDS("../../data/singlecell/scATAC/full_scHeme_QC.rds")
pcaPlus <- data.table::fread("../../data/singlecell/scATAC/scHeme_meta.txt")
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE))
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)

zscoreRaw <- t(assays(ukbbDEV)[["z"]])
colnames(zscoreRaw) <- paste0("raw_", colnames(zscoreRaw))
smoothCellxCell <- smoothPCAmat(pcaPlus[,-c(1,2)])
zscoreSmooth <- (smoothCellxCell %*% zscoreRaw)
colnames(zscoreSmooth) <- paste0("smooth_", colnames(zscoreSmooth))

bigdf <- data.frame(celltypes, pca, data.matrix(zscoreRaw), data.matrix(zscoreSmooth))
write.table(bigdf, file = "", row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")