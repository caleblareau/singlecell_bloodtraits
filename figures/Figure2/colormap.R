library(ggplot2)
library(BuenColors)

ldscore <- read.table("../../data/bulk/bulkHeme_baseline_plus1_pvalues.txt", stringsAsFactors = FALSE)

cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "mono", "mDC", "Ery", "Mega"),
          x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
          y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
)

ejc_color_maps <-c(
  "HSC" = "#00441B",
  "MPP" = "#46A040",
  "LMPP" = "#00AF99",
  "CMP" = "#FFC179",
  "CLP" = "#98D9E9",
  "MEP" = "#F6313E",
  "pDC" = "#CDA2DB",
  "mono" = "#FF5A00",
  "GMP-A" = "#FFCE00",
  "GMP-B" = "#FFA300",
  "GMP-C" =  "#FF7900",
  "Ery" = "#8F1336",
  "CD4" = "#0081C9",
  "CD8" = "#001588",
  "NK" = "#490C65",
  "B" = "#A65AC2",
  "Mega" = "#FF81AF",
  "mDC"= "#FFEC00")

df <- data.frame(names = names(ejc_color_maps), 
                 ejc_color_maps = unname(ejc_color_maps), stringsAsFactors = FALSE)

ggplot(df, aes(x = names, y = -1, fill = names)) +
  geom_bar(stat = "identity") + pretty_plot(fontsize = 20) + 
  scale_fill_manual(values = ejc_color_maps) + coord_flip() +
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Erik/Jacob/Caleb Lab Cell Schemes") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.border = element_blank())




ggplot(cellCoordsDF, aes(x = x, y = y, color = CellLabel)) + 
  geom_point(size = 11) + pretty_plot() +
  scale_color_manual(values = ejc_color_maps)