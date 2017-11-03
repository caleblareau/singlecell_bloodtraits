library(reshape2)
library(BuenColors)
"%ni%" <- Negate("%in%")

# Import data
observed <- read.table("../../data/singlecell/scATAC/weightedSingleCellScores-shuffled.txt", sep = "\t", header = TRUE)
observeddf <- observed[,c(1,2,13:28)]
colnames(observeddf) <- gsub("raw_", "", colnames(observeddf))

# Make p-values of variance of sc GWAS enrichments
pcalc <- function(vec){
  pchisq((length(vec) - 1) * var(vec), df = (length(vec) - 1), lower.tail = FALSE)
}

pvalMat <- sapply((1:16)+2, function(i){
  tapply(observeddf[,i],list(as.factor(observeddf[,2])), pcalc)
})

colnames(pvalMat) <- colnames(observeddf)[c((1:16)+2)]
long <- reshape2::melt(pvalMat)
long <- long[order(long$value), ]
long <- long[long$Var2 %ni% c("BASO_COUNT", "MCHC"),]
long$FDR <- p.adjust(long$value)
long$Rank <- 1:dim(long)[1]

p0<-ggplot(long, aes(x = Rank, y = -log10(FDR), color = Var1)) + 
  geom_point(size = 3) +
  geom_point(shape = 1,size = 3,colour = "black") +
  geom_text_repel(
    data = subset(long, Rank<=10),
    aes(label = paste(Var1,Var2,sep="/")),
    size = 2.5,
    nudge_y =0,
    nudge_x=0.6,
    point.padding=0.5,
    min.segment.length=0.25,
    segment.alpha=0.6,
    direction="both") + 
  labs(x = "Rank Sorted Enrichment", y = "- log10(FDR)", color = "") + pretty_plot() +
  theme(legend.position = "bottom") + scale_color_manual(values = ejc_color_maps) 

ggsave(p0, file = "varianceEnrichment_singleCell.pdf",width=6,height=6, useDingbats=F)
write.table(long, file = "varianceEnrichment_singlecell.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
