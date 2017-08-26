library(ggplot2)
library(BuenColors)

# Choose Color Palette
THE_PALETTE <- jdb_palette("solar_rojos")

ldscore <- read.table("../../data/bulk/bulkHeme_baseline_plus1_pvalues.txt", stringsAsFactors = FALSE)

cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "mono", "mDC", "Ery", "Mega"),
          x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
          y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
)

#---------------
# LD Score plots
#---------------

makeLDScorePlot <- function(plottrait){
  plotdf <- merge(cellCoordsDF, ldscore[ldscore$trait == plottrait, ],
                  by.x = "CellLabel", by.y = "Category")
  
  p1 <- ggplot(plotdf, aes(x = x, y = y, color = -log10(pvalue))) + 
    geom_point(size = 11) + pretty_plot() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3) + 
    scale_color_gradientn(colors = THE_PALETTE) +
    scale_y_continuous(limits = c(0, 11)) + ggtitle(paste0("LDScore ", plottrait))
  
  ggsave(p1, filename = paste0("rawPDFsForIllustrator/LDscore/LDScore_", plottrait, ".pdf"),
         height = 8, width = 10)
  return(plottrait)
}
ldscoreout <- lapply(unique(ldscore$trait), makeLDScorePlot)
