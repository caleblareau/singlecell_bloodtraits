library(BuenColors)
library(ggplot2)
library(reshape2)
source("boostrapSplineFunctions.R")

# Import observed 
observed <- read.table("../../data/singlecell/scATAC/weightedSingleCellScores-shuffled.txt", sep = "\t", header = TRUE)
observeddf <- observed[,c(1,2,13:28)]
colnames(observeddf) <- gsub("raw_", "", colnames(observeddf))

# Load pseudotime DFs
ptDF <- lapply(list.files("../../data/pseudotime", pattern= ".txt", full.names = TRUE), function(file){
  read.table(file, header = FALSE, stringsAsFactors = FALSE)
})
names(ptDF) <- sapply(list.files("../../data/pseudotime", pattern= ".txt"), function(file) strsplit(file, split = "_")[[1]][2])
melt_ptDF <- melt(ptDF)[,c(1,3,4)]
names(melt_ptDF) <- c("name", "pseudotime", "trajectory")

# Merge data frames
mergedf <-  merge(melt_ptDF, observeddf, by = c("name"))
meltdf <- melt(mergedf, id.vars = c("name", "pseudotime", "trajectory", "type"))

# One plot for each trajectory for PLT
trait <- "PLT_COUNT"
lapply(c("ery", "lym", "myel"), function(traj){
  
  p1 <- ggplot(meltdf[meltdf$trajectory == traj & meltdf$variable == trait,], aes(pseudotime, value, color = type)) +
    geom_point() +
    stat_smooth(se = FALSE, inherit.aes = FALSE, aes(pseudotime, value), span = 0.6) +
    pretty_plot() + theme(legend.position="bottom") + 
    labs(list( x = "Pseudotime", y = "Zscore")) +
    scale_colour_manual(values = ejc_color_maps, name = "Cell Type") +
    theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill = "white") )
  ggsave(p1, file = paste0("PLT_pseudotime/",trait,".", traj,".pseudotime.pdf"), width = 5, height = 5)
  p1
})

