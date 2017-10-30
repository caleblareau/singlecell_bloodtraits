library(data.table)
library(BuenColors)
library(viridis)

dt <- fread("../../data/singlecell/scATAC/weightedSingleCellScores-shuffled.txt")

# Cell Labels plot setup
write.table(data.frame(dt[["colvec"]]), file = "raw3Dplot/cell_label/cellTypes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(data.frame(dt[,c(3,4,5)]), file = "raw3Dplot/cell_label/PC123.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Spit out raw scores
lapply(13:28, function(i){
  name <- colnames(dt)[i]
  name <- gsub("raw_", "", name)
  df <- data.frame(dt)[,c(3,4,5,i)]
  df[,4] <- scale(df[,4])
  df[df[,4] > 1.645, 4] <- 1.645
  df[df[,4] < - 1.645, 4] <- -1.645
  write.table(df, file = paste0("raw3Dplot/scores/dfs/raw/", name, ".txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
})

lapply(29:44, function(i){
  name <- colnames(dt)[i]
  name <- gsub("smooth_raw_", "", name)
  df <- data.frame(dt)[,c(3,4,5,i)]
  df[,4] <- scale(df[,4])
  df[df[,4] > 1.645, 4] <- 1.645
  df[df[,4] < - 1.645, 4] <- -1.645
  write.table(df, file = paste0("raw3Dplot/scores/dfs/smooth/", name, ".txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
})


# Make color scale for Mean Retic Volume
dt <- read.table("raw3Dplot/scores/dfs/raw/MEAN_RETIC_VOL.txt")
p1 <- ggplot(dt, aes(x = V2, y = V3, color = V4)) +
  geom_point() + coord_fixed() +
  scale_color_viridis() + pretty_plot()
ggsave(p1, file = "mrv_color.pdf", useDingbats=FALSE)
