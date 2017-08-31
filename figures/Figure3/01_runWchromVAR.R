library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(chromVAR)
library(chromVARxx)
library(BSgenome.Hsapiens.UCSC.hg19)
source("00_smooth.R")

SE <- readRDS("../../data/singlecell/scATAC/full_scHeme_QC.rds")
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE))
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)

