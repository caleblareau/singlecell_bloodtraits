library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(dplyr)
library(chromVAR)
library(chromVARxx)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg19)

SE <- readRDS("../singlecell/scATAC/full_scHeme_QC.rds")

# Make one synthetic counts
typeCountVec <- table(colData(SE)$type)
randCountsMat <- lapply(names(typeCountVec), function(type){
  cmat <- assays(SE)[["counts"]][,which(colData(SE)$type == type)]
  
  #Set up samples
  ns <- colSums(cmat) 
  nidx <- rep.int(1:length(ns), ns)
  
  # Filter peaks
  maxpeakidx <- dim(cmat)[1]
  keepPeaks <- which(rowSums(cmat) > 0)
  cmat <- cmat[keepPeaks, ]
  
  # Get random counts per each
  runifs <- runif(sum(ns))
  
  # Define break points
  breaks <- cut(runifs, c(0,cumsum(rowSums(cmat)/sum(rowSums(cmat)))))
  levels(breaks) <- as.character(keepPeaks)
  
  mat <- sparseMatrix(
    i = c(as.numeric(as.character(breaks)), maxpeakidx),
    j = c(nidx, 1), 
    x = c(rep(1, length(nidx)), 0)
  )
  mat
}) %>% do.call(what = cbind)

SEperm <- SummarizedExperiment(
  rowData = rowRanges(SE),
  assays = list(counts = randCountsMat), 
  colData = DataFrame(type = rep.int(names(typeCountVec), as.numeric(typeCountVec)))
)

SEperm <- filterPeaks(SEperm)
SEperm <- addGCBias(SEperm, genome = BSgenome.Hsapiens.UCSC.hg19)

ukbb <- importBedScore(rowRanges(SEperm), list.files("../UKBB_BC_PP001/", full.names = TRUE))
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
mat <- t(assays(ukbb_wDEV)[["z"]])
rownames(mat) <- rep.int(names(typeCountVec), as.numeric(typeCountVec))

pcalc <- function(vec){
  pchisq((length(vec) - 1) * var(vec), df = (length(vec) - 1), lower.tail = FALSE)
}

pvalMat <- sapply(1:16, function(i){
  tapply(zscoreRaw[,i],list(as.factor(rownames(mat))), pcalc)
})

colnames(pvalMat) <- colnames(zscoreRaw)
long <- reshape2::melt(pvalMat)
long <- long[order(long$value), ]
long$FDR <- p.adjust(long$value)

