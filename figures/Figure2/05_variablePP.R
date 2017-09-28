 library(chromVAR)
library(chromVARxx)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = ".bed"))

# Make variable PP thresholds
makePPthreshdf <- function(val){
  ukbbcp <- ukbb
  assays(ukbbcp)$weights[assays(ukbbcp)$weights < val] <- 0
  zw<- melt(t(assays( computeWeightedDeviations(SE, ukbbcp))[["z"]]))
  zw$PP <- val
  return(zw)
}

PPs <- c(1:10/1000, 2:100/100)

zscoreWeighted <- rbindlist(lapply(head(PPs,2), makePPthreshdf))

ggplot(zscoreWeighted, aes(x = PP, y = -log10(pnorm(value,lower.tail = FALSE)))) +
  geom_line(aes(color = Var1)) +
  pretty_plot() + labs(x = "Posterior Probability", y = "g-chromVAR Enrichment (-log10 p)", color = "") +
   scale_color_manual(values = ejc_color_maps) + facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") 
#  geom_hline(yintercept = -log10(0.05 / (16*17)), linetype = 2) +

