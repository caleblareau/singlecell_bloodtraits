library(GenomicRanges)
library(data.table)
library(BuenColors)

# Choose Color Palette
THE_PALETTE <- jdb_palette("solar_rojos")

# Import data
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

# Import counts and normalize
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
ATAC.cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

cebpa <- makeGRangesFromDataFrame(data.frame(chr = "chr19", start = 33754548, end = 33754549))
kit <- makeGRangesFromDataFrame(data.frame(chr = "chr4", start = 55408999, end = 55409000))
bcl2 <- makeGRangesFromDataFrame(data.frame(chr = "chr18", start = 60920854, end = 60920855))
myc <- makeGRangesFromDataFrame(data.frame(chr = "chr8", start = 130604272, end = 130604273))
cebpe <- makeGRangesFromDataFrame(data.frame(chr = "chr14", start = 23589349, end = 23589350))


cebpa_counts <- ATAC.cpm.log2[subjectHits(findOverlaps(cebpa, peaks)),]
kit_counts <- ATAC.cpm.log2[subjectHits(findOverlaps(kit, peaks)),]
bcl2_counts <- ATAC.cpm.log2[subjectHits(findOverlaps(bcl2, peaks)),]
myc_counts <- ATAC.cpm.log2[subjectHits(findOverlaps(myc, peaks)),]
cebpe_counts <- ATAC.cpm.log2[subjectHits(findOverlaps(cebpe, peaks)),]

# Set up coordinates
cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "mDC", "Ery", "Mega"),
          x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
          y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
)

#---------------
# accessilibity plots per population
#---------------

makeAplot <- function(vec, plottrait){
  plotdf <- data.frame(cellCoordsDF, counts = vec[as.character(cellCoordsDF$CellLabel)] - min(vec))
  
  p1 <- ggplot(plotdf, aes(x = x, y = y, color = counts)) + 
    geom_point(size = 11) + pretty_plot() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3, color = "black") + 
    scale_color_gradientn(colors = THE_PALETTE) +
    scale_y_continuous(limits = c(0, 11)) + ggtitle(plottrait)
  
  cowplot::ggsave(p1, filename = paste0("accessibilityPlots/", plottrait, ".pdf"),
         height = 8, width = 10)
  return(plottrait)
}

makeAplot(cebpa_counts, "CEBPa")
makeAplot(kit_counts, "KIT")
makeAplot(bcl2_counts, "BCL2")
makeAplot(myc_counts, "MYC")
makeAplot(cebpe_counts, "CEBPE")