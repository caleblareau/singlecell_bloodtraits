library(BuenColors)

makeHighLowPlot <- function(celltype, trait, X = "PC2", Y = "PC3"){
  df <- read.table(paste0("data/",celltype,"_weightedchromVAReach.txt"), sep = "\t", header = TRUE)
  plotdf <- data.frame(
    x = df[,X],
    y = df[,Y],
    color = df[,trait] > 0
  )
  
  ggplot(plotdf, aes(x, y, color =  color)) + geom_point(size = 4) +
    pretty_plot() + labs(x = X, y = Y, color = paste0(trait, " HIGH")) +
    theme(legend.position = "bottom") +
    scale_color_manual(values= c("dodgerblue4", "yellow")) + ggtitle(paste0(celltype, " ", trait))
  
}

makeHighLowPlot("CMP", "MONO_COUNT", "PC2", "PC3")
makeHighLowPlot("CMP", "MEAN_RETIC_VOL", "PC2", "PC3")
