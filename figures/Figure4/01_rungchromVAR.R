library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(chromVAR)
library(gchromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BuenColors)
library(reshape2)
library(data.table)
library(dplyr)

if(FALSE){
  
  # Run weighted deviations
  lapply(1:100, function(i){
    print(i)
    SE <- readRDS("../../data/singlecell/scATAC/full_scHeme_QC.rds")
    pcaPlus <- data.table::fread("../../data/singlecell/scATAC/scHeme_meta.txt")
    SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
    ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001", full.names = TRUE, pattern = ".bed$"))
    ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
    zscoreRaw <- t(assays(ukbb_wDEV)[["z"]])
    rownames(zscoreRaw) <- colData(SE)$name
    colnames(zscoreRaw) <- gsub("_PP001", "", colnames(zscoreRaw))
    mdf <- melt(zscoreRaw)
    mdf$iteration <- i
    return(mdf)
  }) %>% rbindlist() %>% as.data.frame() -> bdf
  
  saveRDS(bdf, file = "manyIteration.rds")
}

bdf <- readRDS("manyIteration.rds")
bdf %>% group_by(Var1, Var2) %>% summarize(zscore = mean(value))  %>% as.data.frame() -> meanDF
write.table(meanDF, file = "../../data/singlecell/scATAC/scGWASenrichments_average.tsv", sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

