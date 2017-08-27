library(ggplot2)
library(BuenColors)

cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "mono", "mDC", "Ery", "Mega"),
          x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
          y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
)

p1 <- ggplot(cellCoordsDF, aes(x = x, y = y, color = CellLabel)) + 
  geom_point(size = 11) + pretty_plot() +
  scale_color_manual(values = ejc_color_maps)+ 
  geom_text(aes(label=CellLabel),hjust=0.5, vjust=3) +
  scale_y_continuous(limits = c(0, 11))
ggsave(p1, file = "mainHierarchy.pdf")