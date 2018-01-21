library(BuenColors)
library(ggrepel)

df <- readRDS("../roadmap_dnase/deviationDataframe.rds")
df$Color <- (ifelse(df$celltype == "CD8_Memory_Primary_Cells",1,0) + 
  ifelse(df$celltype %in% c("HUES64_Cell_Line", "HUES6_Cell_Line", "HUES48_Cell_Line"),2,0) +
   ifelse(df$celltype == "IMR90_Cell_Line",3,0) ) * ifelse(df$P < 0.05/dim(df)[1],1,0) + 1
tv <-  c("Other", "CD8", "HESC", "IMR90")
df$ColorType <- factor(tv[df$Color], levels = tv)
df$CellTrait <- paste0(df$ColorType, " ", df$Trait)                       
p <- ggplot(df, aes(x = expected, y = logP, color = ColorType)) +
  geom_point() +
  scale_color_manual(values = c("black", "#001588", "green4", "purple4")) +
  geom_hline(yintercept = -log10(0.05/dim(df)[1]), linetype = 2) +
  geom_abline(slope=1, intercept=0) + pretty_plot() +
  labs(x = "Expected -log10(p)", y = "Observed g-chromVAR -log10(p)", color = "") +
  theme(legend.position = "bottom")
gglabeller_data <- p$data
gglabeller_data$gglabeller_labels <- df$CellTrait
gglabeller_data[c(9:848),'gglabeller_labels'] <- ''
p + geom_text_repel(data = gglabeller_data,mapping = aes(label = gglabeller_labels),
                    min.segment.length = unit(0.5, 'lines'), show.legend = FALSE,
                    box.padding = unit(0.25, 'lines'),point.padding = unit(1e-06, 'lines'))