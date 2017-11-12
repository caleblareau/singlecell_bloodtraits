library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)
library(BuenColors)

# For the bulk, import binarized peaks
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

fullnp <- list.files("../../data/bulk/ATAC/narrowpeaks", full.names = TRUE)
short <- gsub("_peaks.narrowPeak.gz", "", list.files("../../data/bulk/ATAC/narrowpeaks"))

sapply(1:length(short), function(i){
  dt <- data.frame(fread(paste0("zcat < ",fullnp[i])))
  g <- makeGRangesFromDataFrame(dt, seqnames = "V1", start.field = "V2", end.field = "V3")
  v <- 1:length(peaks) %in% subjectHits(findOverlaps(peaks, g))
  v
}) -> counts
dim(counts)
colnames(counts) <- short

SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
SE <- filterPeaks(SE)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))

# Compute weighted deviation scores
binPeakDev <- computeWeightedDeviations(SE, ukbb)
zscoreBinPeak <- melt(t(assays(binPeakDev)[["z"]]))
zscoreBinPeak[,2] <- gsub("_PP001", "", zscoreBinPeak[,2])

zscoreBinPeak$P <- pnorm(zscoreBinPeak$value, lower.tail = FALSE)
zscoreBinPeak$logP <- -log10(pnorm(zscoreBinPeak$value, lower.tail = FALSE))

ggplot(zscoreBinPeak, aes(x = Var1, y = logP)) +
  geom_bar(width = 1, aes(fill = Var1), colour="black",
           stat = "identity", position = position_dodge(width=1)) +
  pretty_plot() + labs(x = "", y = "Binarized Peak Enrichment (-log10 p)", fill = "") +
   scale_fill_manual(values = ejc_color_maps) + facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") +
  geom_hline(yintercept = -log10(0.05 / (16*17)), linetype = 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
