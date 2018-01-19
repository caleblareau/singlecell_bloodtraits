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
library(dplyr)
library(matrixStats)
library(cowplot)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))
bg <- getBackgroundPeaks(SE)

vec <- sort(as.numeric(assays(ukbb)[["weights"]]))

# Null simulation
n <- 100
weights <- Matrix(matrix(sample(vec, dim(SE)[1]*n, replace = TRUE), ncol = n))
colnames(weights) <-  paste0("X", 1:n)
simSE <- SummarizedExperiment(assays = list(weights = weights),
                              rowData = peaks, 
                              colData = DataFrame(names = paste0("X", 1:n)))
dev <- computeWeightedDeviations(SE, simSE, background_peaks = bg)
outdf <- melt(t(assays(dev)[["z"]]))

Nullsummdf <- outdf %>% group_by(Var1) %>%
  summarise( mean=mean(value), sd=sd(value)) %>%
  mutate(ymax = mean + sd, ymin = mean - sd)


# Non-null simulation
n <- 100 # how many traits
cs <- 0.9 # damping factor. Pushing this to 1 makes zscores really big. Much lower than 0.9 makes them small / not significant

cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
p <- pnorm((cpm-rowMeans(cpm))/(rowSds(as.matrix(cpm)))[row(cpm)]) *cs


simEnrich <- function(causalCelleType){
  
  i <- causalCelleType
  
  # Scale uniform observations such that the min is the CDF(Z) for chromatin
  unifR <- matrix(runif(dim(cpm)[1]*(n)), ncol = n)*(1-p[,i]) + p[,i]
  
  # Index the weights vector based on it's closest match
  vals <- Matrix(matrix(vec[round(unifR*length(vec))], ncol = n)) 
  
  # Down sample # of causal to what we'd normallly observe
  vdf <- data.frame(summary(vals))
  vdf$keep <- rbinom(dim(vdf)[1], 1,sum(weights != 0)/sum(vals != 0))
  vdf$newX <- vdf$x*vdf$keep
  vdfF <- vdf[vdf$newX > 0, ]
  valsDS <- sparseMatrix(i = c(vdfF$i, dim(vals)[1]), j = c(vdfF$j, dim(vals)[2]), x = c(vdfF$newX, 0))
  
  # Compute deviations
  colnames(valsDS) <- paste0("X", 1:n)
  simSEc2 <- SummarizedExperiment(assays = list(weights = valsDS),
                                  rowData = peaks, 
                                  colData = DataFrame(names = paste0("X", 1:n)))
  
  devC <- computeWeightedDeviations(SE, simSEc2, background_peaks = bg)
  outdf <- melt(t(assays(devC)[["z"]]))
  summdf <- outdf %>% group_by(Var1) %>%
    summarise( mean=mean(value), sd=sd(value)) %>%
    mutate(ymax = mean + sd, ymin = mean - sd)
  
  gp <- ggplot(summdf, aes(x = Var1, y = mean, color = Var1, ymax = ymax, ymin = ymin)) + geom_point() +
    scale_color_manual(values = ejc_color_maps) + pretty_plot() +
    labs(x = "", y = "Z-score") +  geom_errorbar(width = 0.2) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(legend.position = "none") +
    geom_point(data = data.frame(x = colnames(SE)[i] , y = max(summdf$ymax) + 0.3, stringsAsFactors = FALSE), 
               inherit.aes = FALSE, aes(x = x, y = y), shape = 8) + geom_hline(aes(yintercept=0), linetype = "dashed")
  
  return(gp)
}

# Simulated selected traits and make plot

BCellp <- simEnrich(1)
Eryp <- simEnrich(6)
Megap <- simEnrich(13)

CMPp <- simEnrich(5)
CLPp <- simEnrich(4)
HSCp <- simEnrich(10)

# Make  null plot and collate
NullPlot <- ggplot(Nullsummdf, aes(x = Var1, y = mean, color = Var1, ymax = ymax, ymin = ymin)) + geom_point() +
  scale_color_manual(values = ejc_color_maps) + pretty_plot() +
  labs(x = "", y = "g-chromVAR \n Enrichment Z-score") +  geom_errorbar(width = 0.2) +
  geom_hline(aes(yintercept=0), linetype = "dashed") + guides(color=guide_legend(title="Celltype", ncol = 3)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())



bottom_row <- plot_grid(BCellp, Eryp, Megap, CMPp, CLPp, HSCp, labels = c('b', 'c', 'd', 'e', 'f', 'g'), ncol = 2)
simFig <- plot_grid(NullPlot, bottom_row, labels = c('a', ''), ncol = 1, rel_heights = c(1, 3 ))
cowplot::ggsave(simFig, file = "SimulationFigure.pdf", height = 11, width = 8)

