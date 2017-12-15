library(data.table)
library(BuenColors)

gene <- "AK3"

rna <- data.frame(fread(paste0("../../data/bulk/RNA/16populations_RNAcounts.txt")))
genes <- rna[,1]
counts <- data.matrix(rna[,-1])
cpm <- sweep(counts, 2, Matrix::colSums(counts), FUN="/") * 1000000
colnames(cpm) <- c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "Ery")

# Pull out gene and 
ix <- which(genes == gene)
plotdf <- data.frame(Name = colnames(cpm),
                     Expression = cpm[ix,], 
                     Expm = (cpm[ix,]/sum(cpm[ix,])) * 100)
plotdf <- plotdf[order(plotdf$Expm, decreasing = FALSE), ]

# plotdf$pos <- 0.5 * (cumsum(plotdf$Expm) + cumsum(c(0, plotdf$Expm[-length(plotdf$Expm)])))
# p <- ggplot(plotdf, aes(x=pos, y=Expression, width=Expm)) + 
#   geom_bar(aes(fill=Name), stat="identity",  color = "black") +
#   scale_fill_manual(values = ejc_color_maps) + pretty_plot() + coord_flip() +
#   theme(legend.position = "none") + ylab(paste0(gene, " Expression (TPM)")) +
#   xlab("") + theme(panel.grid = element_blank(), panel.border = element_blank()) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
plotdf$Name <- factor(as.character(plotdf$Name), levels = as.character(plotdf$Name))
p1 <- ggplot(plotdf, aes(x=Name, y=log2(Expression + 1))) + 
  geom_bar(aes(fill=Name), stat="identity",  color = "black") +
  scale_fill_manual(values = ejc_color_maps) + pretty_plot() + coord_flip() +
  theme(legend.position = "none") + ylab(paste0(gene, " Expression (log2 TPM)")) +
  xlab("") + theme(panel.grid = element_blank(), panel.border = element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


cowplot::ggsave(filename = paste0("plots/Expression_", gene, ".pdf"), p1, width = 5, height = 5)
