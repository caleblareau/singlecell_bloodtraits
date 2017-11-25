library(data.table)
library(BuenColors)
library(viridis)

dt <- fread("../../data/singlecell/scATAC/scGWASenrichments_average.tsv")
coords <- fread("../../data/singlecell/scATAC/scHeme_meta.txt")

# Spit out raw scores
lapply(unique(dt[["V2"]]), function(name){
  df1 <- data.frame(coords)[,c(1,3,4,5)]
  df <- merge(df1, dt[dt$V2 == name, ], by.x = "name", by.y = "V1")
  df <- df[,c(2,3,4,6)]
  df[df[,4] > 1.645, 4] <- 1.645
  df[df[,4] < - 1.645, 4] <- -1.645
  write.table(shuf(df), file = paste0("raw3Dplot/scores/dfs/raw/", name, ".txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
})


# Make color scale for Mean Retic Volume
dt <- read.table("raw3Dplot/scores/dfs/raw/MEAN_RETIC_VOL.txt")
p1 <- ggplot(dt, aes(x = V2, y = V3, color = V4)) +
  geom_point() + coord_fixed() +
  scale_color_viridis() + pretty_plot()
ggsave(p1, file = "mrv_color.pdf", useDingbats=FALSE)
