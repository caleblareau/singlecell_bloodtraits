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
counts_mix <- data.matrix(fread("../../data/bulk/simulatedATAC/MONO_RBC_mixed.counts.txt"))

mixnames <- c("MONO0_RBC100","MONO10_RBC90","MONO100_RBC0","MONO20_RBC80","MONO30_RBC70","MONO40_RBC60","MONO50_RBC50","MONO60_RBC40","MONO70_RBC30","MONO80_RBC20","MONO90_RBC10")

colnames(counts_mix) <- mixnames
SE <- SummarizedExperiment(assays = list(counts = cbind(counts, counts_mix)),
                               rowData = peaks, 
                               colData = DataFrame(names = c(colnames(counts), colnames(counts_mix))))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = ".bed"))

# Compute weighted deviation scores
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]]))
zscoreWeighted[,2] <- gsub("_PP001", "", zscoreWeighted[,2])
zscoreWeighted <- zscoreWeighted[zscoreWeighted$Var1 %in% mixnames, ]
zscoreWeighted[zscoreWeighted$Var2 == "RBC_COUNT", ]

#write.table(zscoreWeighted, "../../data/bulk/GWAS-Bulk/syntheticHeme_weighted_zscores.txt", 
#            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run regular chromVAR on genome-wide significant variants
ixveclist <- lapply(list.files("../../data/Significant_SNPs", full.names = TRUE), function(file){
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
zscoreCV <- melt(t(assays(ukbb_wDEV)[["z"]]))
zscoreCV[,2] <- gsub("_PP001", "", zscoreCV[,2])
zscoreCV <- zscoreCV[zscoreCV$Var1 %in% mixnames, ]
zscoreCV[zscoreCV$Var2 == "MONO_COUNT", ]



#write.table(zscoreCHROMVAR, "../../data/bulk/GWAS-Bulk/syntheticHeme_chromVAR_zscores.txt", 
#            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
