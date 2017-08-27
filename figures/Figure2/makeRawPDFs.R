library(ggplot2)
library(BuenColors)

# Choose Color Palette
THE_PALETTE <- jdb_palette("solar_rojos")

# Import data
ldscore <- read.table("../../data/bulk/GWAS-Bulk/bulkHeme_baseline_plus1_pvalues.txt", stringsAsFactors = FALSE)
chromvar <- read.table("../../data/bulk/GWAS-Bulk/bulkHeme_chromVAR_zscores.txt", stringsAsFactors = FALSE)
wchromvar <- read.table("../../data/bulk/GWAS-Bulk/bulkHeme_weighted_zscores.txt", stringsAsFactors = FALSE)

# Set up coordinates
cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "mDC", "Ery", "Mega"),
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
  
  ggsave(p1, filename = paste0("rawPDFs/LDscore/LDScore_", plottrait, ".pdf"),
         height = 8, width = 10)
  return(plottrait)
}

#--------------------------
# chromVAR / weighted plots
#--------------------------

makeCVplot <- function(plottrait, method){
  
  if(method == "chromVAR"){
    df <- chromvar
  }else {
    df <- wchromvar
  }
  plotdf <- merge(cellCoordsDF, df[df$V2 == plottrait, ],
                  by.x = "CellLabel", by.y = "V1")
  plottrait <- gsub("_PP001", "", plottrait)
  
  p1 <- ggplot(plotdf, aes(x = x, y = y, color = V3)) + 
    geom_point(size = 11) + pretty_plot() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3) + 
    scale_color_gradientn(colors = THE_PALETTE, name = "Zscore") +
    scale_y_continuous(limits = c(0, 11)) + ggtitle(paste0(method, " ", plottrait))
  
  ggsave(p1, filename = paste0("rawPDFs/",method,"/",method,"_", plottrait, ".pdf"),
         height = 8, width = 10)
  return(plottrait)
}

chromvarout <- lapply(unique(chromvar$V2), makeCVplot, "chromVAR")
weightedout <- lapply(unique(wchromvar$V2), makeCVplot, "wChromVAR")
ldscoreout <- lapply(unique(ldscore$trait), makeLDScorePlot)

# Make combined table
stopifnot(all(wchromvar[,1] == chromvar[,1]))
stopifnot(all(wchromvar[,2] == chromvar[,2]))
outdf <- data.frame(
  Celltype = wchromvar[,1],
  Trait = gsub("_PP001", "", wchromvar[,2]),
  weightedchromVAR_Zscore =  wchromvar[,3],
  chromVAR_Zscore =  chromvar[,3]
)
mergedf <- merge(outdf, ldscore,
                 by.x = c("Celltype", "Trait"), by.y = c("Category", "trait"))
dim(mergedf)   
colnames(mergedf) <- c("Celltype", "Trait", "weighted_Zscore", "chromVAR_Zscore", "ldscore_pvalue")
write.table(mergedf, "../../data/bulk/GWAS-Bulk/compare3.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)                 


