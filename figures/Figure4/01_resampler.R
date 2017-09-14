library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(dplyr)
library(chromVAR)
library(chromVARxx)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg19)

args <- commandArgs(trailingOnly = TRUE)
simN <- args[1]

SE <- readRDS("../../data/singlecell/scATAC/full_scHeme_QC.rds")
# Make one synthetic counts
typeCountVec <- table(colData(SE)$type)

blappout <- lapply(1:100, function(i){

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
  
  names(typeCountVec) <- c("CLP", "CMP", "GMP", "GMP-A", "GMP-B", "GMP-C", "HSC", "LMPP", "MEP", "Mono", "MPP", "pDC")
  SEperm <- SummarizedExperiment(
    rowData = rowRanges(SE),
    assays = list(counts = randCountsMat), 
    colData = DataFrame(type = rep.int(names(typeCountVec), as.numeric(typeCountVec)))
  )
  
  SEperm <- filterPeaks(SEperm)
  SEperm <- addGCBias(SEperm, genome = BSgenome.Hsapiens.UCSC.hg19)
  
  ukbb <- importBedScore(rowRanges(SEperm), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))
  ukbb_wDEV <- computeWeightedDeviations(SEperm, ukbb)
  mat <- assays(ukbb_wDEV)[["z"]]
  colnames(mat) <- rep.int(names(typeCountVec), as.numeric(typeCountVec))
  rownames(mat) <- gsub("_PP001", "", rownames(mat))
  
  varMat <- sapply(1:16, function(i){
    tapply(mat[i,],list(as.factor(colnames(mat))), var)
  })
  colnames(varMat) <- rownames(mat)
  varMat
})

saveRDS(blappout, paste0("../../data/resample/resample.variance.", as.character(simN), ".rds"))
