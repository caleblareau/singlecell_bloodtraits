library(BuenColors)
library(ggrepel)
library(cowplot)

# ROADMAP

df <- readRDS("roadmap_deviationDataframe.rds")
df$Color <- (ifelse(df$celltype == "CD8_Memory_Primary_Cells",1,0) + 
  ifelse(df$celltype %in% c("HUES64_Cell_Line", "HUES6_Cell_Line", "HUES48_Cell_Line"),2,0) +
   ifelse(df$celltype == "IMR90_Cell_Line",3,0) )  + 1
tv <-  c("Other", "CD8", "HESC", "IMR90")
df$ColorType <- factor(tv[df$Color], levels = tv)
df$CellTrait <- paste0(df$ColorType, " ", df$Trait)      
df$FDR <- qvalue::qvalue(df$P)$qvalues
df$logFDR <- -1*log10(df$FDR)
df$Rank <- 1:dim(df)[1]

p1 <- ggplot(df, aes(x = Rank, y = logFDR, color = ColorType)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("black", "#001588", "firebrick", "green4")) +
  geom_hline(yintercept = -log10(0.01), linetype = 2) +
  pretty_plot() +
  labs(x = "Rank Sorted Enirchment", y = "-log10 FDR", color = "") +
  theme(legend.position = "bottom") + theme(legend.position = c(0.7, 0.7))


# PICS
df <- read.table("pics_gchromVAR.txt", header = TRUE)
df$FDR <- qvalue::qvalue(df$pvalue)$qvalues
df$logFDR <- -1*log10(df$FDR)
df$Rank <- 1:dim(df)[1]

p2 <- ggplot(df, aes(x = Rank, y = logFDR, color = Celltype)) +
  geom_point(size = 2) +
  scale_color_manual(values = ejc_color_maps) +
  geom_hline(yintercept = -log10(0.01), linetype = 2) +
  pretty_plot() +
  labs(x = "Rank Sorted Enirchment", y = "-log10 FDR", color = "") +
  theme(legend.position = "bottom") + theme(legend.position = c(0.5, 0.7)) +
  guides(color=guide_legend(title="", ncol = 3))


final_plot <- plot_grid(p1, p2, labels = c('a', 'b'), ncol = 2)

cowplot::ggsave(final_plot, filename = "../PDFs/FigureS7.pdf", height = 4, width = 7.5)
