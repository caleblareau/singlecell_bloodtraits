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
counts <- data.matrix(counts)
log2cpm <- data.matrix(log2(counts + 1))



# Set up coordinates
cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP.A", "GMP.B", "GMP.C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono",  "Ery"),
  x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     8), 
  y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,      2)
)

#---------------
# accessilibity plots per population
#---------------

makeAplot <- function(vec, gene){
  
  vec <- log2cpm[gene,]
  
  plotdf <- data.frame(cellCoordsDF, counts = vec[as.character(cellCoordsDF$CellLabel)] )
  
  p1 <- ggplot(plotdf, aes(x = x, y = y, color = counts)) + 
    geom_point(size = 11) + pretty_plot() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3, color = "black") + 
    scale_color_gradientn(colors = THE_PALETTE) +
    scale_y_continuous(limits = c(0, 11)) + ggtitle(gene)
  
  cowplot::ggsave(p1, filename = paste0("expressionPlots/", gene, ".pdf"),
                  height = 8, width = 10)
  return(gene)
}


makeAplot(gfi1b_counts, "USP8")
makeAplot(gfi1b_counts, "TUBB1")
makeAplot(gfi1b_counts, "HFE")

if(FALSE){
  df<- data.frame(counts = sort( counts["TMCC2",]), pop = names(sort( counts["TMCC2",])))
  p1 <- ggplot(df, aes(x = pop, y = counts))  + geom_bar(stat = "identity", fill = "lightgrey", color = "black") +
    pretty_plot(fontsize = 8) + L_border() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "TPM")
  
  cowplot::ggsave(p1, file = "expressionPlots/TMCC2-bar.pdf", width = 6.5, height = 2.75)
  
}



