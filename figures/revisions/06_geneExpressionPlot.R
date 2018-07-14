library(GenomicRanges)
library(data.table)
library(BuenColors)

# Choose Color Palette
THE_PALETTE <- jdb_palette("solar_rojos")

# Import RNA
x <- read.table("../../data/bulk/RNA/16populations_RNAcounts.txt", header = TRUE)
rownames(x) <- make.unique(as.character(x$Genes))
x <- x[,2:17]

counts <- x/colSums(x) * 1000000
log2cpm <- data.matrix(log2(counts + 1))

gfi1b_counts <- log2cpm["GFI1B",]


# Set up coordinates
cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP.A", "GMP.B", "GMP.C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono",  "Ery"),
          x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     8), 
          y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,      2)
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
  
  cowplot::ggsave(p1, filename = paste0("expressionPlots/", plottrait, ".pdf"),
         height = 8, width = 10)
  return(plottrait)
}

makeAplot(gfi1b_counts, "GFI1B")

