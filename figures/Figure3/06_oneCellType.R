library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(BuenColors)
library(Matrix)
library(magrittr)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = peaks, 
                               colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

# Loop over each cell type
lapply(1:dim(SE)[2], function(j){
  SEs <- SE[,j]
  SEs <- filterPeaks(SEs)
  ukbb <- importBedScore(rowRanges(SEs), list.files("../../data/UKBB_BC_PP001/",
                                                    full.names = TRUE, pattern = "*.bed$"))
  ukbb_wDEV <- computeWeightedDeviations(SEs, ukbb)
  zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]]))
  zscoreWeighted[,2] <- gsub("_PP001", "", zscoreWeighted[,2])
  zscoreWeighted
}) %>% rbindlist() %>% as.data.frame() -> oneCelltypeBG

saveRDS(oneCelltypeBG, file = "oneCelltypeBG.rds")

oneCelltype <- readRDS("oneCelltypeBG.rds")

oneCelltype$P <- pnorm(oneCelltype$value, lower.tail = FALSE)
oneCelltype$logP <- -log10(pnorm(oneCelltype$value, lower.tail = FALSE))

ggplot(oneCelltype, aes(x = Var1, y = logP)) +
  geom_bar(width = 1, aes(fill = Var1), colour="black",
           stat = "identity", position = position_dodge(width=1)) +
  pretty_plot() + labs(x = "", y = "One Celltype Enrichment (-log10 p)", fill = "") +
   scale_fill_manual(values = ejc_color_maps) + facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") +
  geom_hline(yintercept = -log10(0.05 / (16*17)), linetype = 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

