library(BuenColors)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)
library(matrixStats)
library(gchromVAR)
library(SummarizedExperiment)
library(Matrix)

colScale = function(x,center = TRUE, scale = TRUE){
  cm = colMeans(x, na.rm = TRUE)
  csd = colSds(x, center = cm)
  x = t( (t(x) - cm) / csd )
  return(x)
}

# Find peaks with FM variants
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
ukbb <- importBedScore(peaks, list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))

# Import counts and normalize
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
counts <- counts[,c("HSC", "MPP", "CMP", "MEP", "Ery")]
ATAC.cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

# Subset
boo <- rowMaxs(ATAC.cpm.log2) > 2 & rowSums(assays(ukbb)[["weights"]]) > 0
ATAC.cpm.log2 <- ATAC.cpm.log2[boo,]
peaksdf <- peaksdf[boo,]
peaks <- peaks[boo]
ATAC.cpm.log2.z <- t( colScale(t(ATAC.cpm.log2)) )

# Pull out the variable peaks in erythroid lineage
ery_countsZ <- ATAC.cpm.log2.z[,c("HSC", "MPP", "CMP", "MEP", "Ery")]

dim(ery_peaksdf)

pdf(file="plots/ery_cluster_peaks.pdf", width = 3, height = 8)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap(ery_countsZ, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 6),
        km = 8,show_heatmap_legend = FALSE,
        km_title = "C%i",
        name = "Peak\nAccessibility")
dev.off()


