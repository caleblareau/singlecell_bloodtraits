library(chromVAR)
library(chromVARxx)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)

# Create bulk Summarized Experiment
peaksdf <- fread("../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

# Return vector of peaks that contain a variant causal in both with opposite effect sizes
pleioMatch <- function(file1, file2){
  nnames <- c("chr", "start", "end", "region", "PP", "beta", "SE", "Z")
  # Import Bed and Make GRanges
  hitdf1 <- data.table::fread(as.character(file1))
  names(hitdf1) <- nnames 
  
  hitdf2 <- data.table::fread(as.character(file2))
  names(hitdf2) <- nnames 
  
  hitdf <- merge(hitdf1, hitdf2, by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"))
  hitdf <- hitdf[sign(hitdf$beta.x) != sign(hitdf$beta.y), ]
  if(dim(hitdf)[1] == 0) return(rep(0, length(peaks)))
  hitG <-makeGRangesFromDataFrame(hitdf, keep.extra.columns = TRUE)
  
  # Get overlaps
  ov <- findOverlaps(rowRanges(SE), hitG)
  ix <- Matrix(as.logical(as.integer((1:length(peaks) %in% unique(queryHits(ov))))))
  ix
}

files <- list.files("../data/UKBB_BC_PP001/betas_added",  full.names = TRUE)
pairs <- data.frame(t(combn(files, 2)), stringsAsFactors = FALSE)
pairs$short1 <- gsub("_PP001_betas.bed", "", gsub("../data/UKBB_BC_PP001/betas_added/", "", pairs[,1]))
pairs$short2 <- gsub("_PP001_betas.bed", "", gsub("../data/UKBB_BC_PP001/betas_added/", "", pairs[,2]))
pairs$traits <- paste0(pairs$short1, "_", pairs$short2)

ixveclist <- lapply(1:dim(pairs)[1], function(i){
  pleioMatch(pairs[i,1], pairs[i,2])
})

motifMatches<- do.call(cbind,ixveclist)
colSums(motifMatches)
idxkeep <- which(colSums(motifMatches) > 0)
motifMatches <- motifMatches[,idxkeep]
pairskeep <- pairs[idxkeep, ]

ixx <- SummarizedExperiment(assays = list(motifMatches = motifMatches),
                            rowRanges = rowRanges(SE),
                            colData = DataFrame(trait = pairskeep$traits))

# Compute deviation scores and save
ukbb_DEV <- computeDeviations(SE, ixx)

# Make tables and export
widemat <- t(assays(ukbb_DEV)[["z"]])
colnames(widemat) <-  pairskeep$traits
zscore_chromVAR <- melt(widemat)
zscore_chromVAR <- zscore_chromVAR[order(zscore_chromVAR$value, decreasing = TRUE),]
zscore_chromVAR[zscore_chromVAR$Var2 == "BASO_COUNT_RBC_COUNT",]
zscore_chromVAR[zscore_chromVAR$Var2 == "MPV_PLT_COUNT",]

#write.table(zscoreCHROMVAR, "../../data/bulk/GWAS-Bulk/bulkHeme_chromVAR_zscores.txt", 
#            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
