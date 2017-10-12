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
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE))

# Compute weighted deviation scores
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]]))
zscoreWeighted[,2] <- gsub("_PP001", "", zscoreWeighted[,2])

#write.table(zscoreWeighted, "../../data/bulk/GWAS-Bulk/bulkHeme_weighted_zscores.txt", 
#            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run regular chromVAR on genome-wide significant variants
ixveclist <- lapply(list.files("../../data/Significant_SNPs/", full.names = TRUE), function(file){
  print(file)
  # Import Bed and Make GRanges
  hitdf <- data.table::fread(paste0("zcat < ", file))
  names(hitdf) <- c("chr", "start", "end", "region", "PP")
  hitG <- addchr(makeGRangesFromDataFrame(hitdf, keep.extra.columns = TRUE))
  
  # Get overlaps
  ov <- findOverlaps(rowRanges(SE), hitG)
  ix <- Matrix(as.logical(as.integer((1:length(peaks) %in% unique(queryHits(ov))))))
  ix
})
motifMatches<- do.call(cbind,ixveclist)
traits <- gsub(".all.sigSNPs.txt.gz", "", list.files("../../data/Significant_SNPs/"))
ixx <- SummarizedExperiment(assays = list(motifMatches = motifMatches),
                            rowRanges = rowRanges(SE),
                            colData = DataFrame(trait = traits))

# Compute deviation scores and save
ukbb_DEV <- computeDeviations(SE, ixx)

# Make tables and export
widemat <- t(assays(ukbb_DEV)[["z"]])
colnames(widemat) <- traits
zscoreCHROMVAR <- melt(widemat)
write.table(zscoreCHROMVAR, "../../data/bulk/GWAS-Bulk/bulkHeme_chromVAR_zscores.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
