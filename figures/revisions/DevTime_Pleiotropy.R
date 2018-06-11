# Libraries
library(BuenColors)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)
library(matrixStats)
library(gchromVAR)
library(SummarizedExperiment)
library(Matrix)
library(preprocessCore)

# Functions
colScale = function(x,center = TRUE, scale = TRUE){
  cm = colMeans(x, na.rm = TRUE)
  csd = colSds(x, center = cm)
  x = t( (t(x) - cm) / csd )
  return(x)
}

# Find peaks with FM variants
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
ukbb <- importBedScore(peaks, list.files("../../data/UKBB_BC_PP001/", full.names = TRUE, pattern = "*.bed$"))

# Import counts, filter, and normalize
counts.df <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
n = 0.9
counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- normalize.quantiles(as.matrix(counts.df))
keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))
counts.df <- counts.df[keep,]
ATAC.cpm <- round(sweep(counts.df, 2, colSums(counts.df), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

# Subset to only good peaks
peaks <- peaks[keep,]
ukbb <- ukbb[keep,]
ukbb.df <- data.frame(as.matrix(assays(ukbb)[["weights"]]))

# Subset to fine-mapped variants
ukbb.df.RBC <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(HCT_PP001 > 0.1 | HGB_PP001 > 0.1 | MCH_PP001 > 0.1, MCHC_PP001 > 0.1 | MCV_PP001 > 0.1 | MEAN_RETIC_VOL_PP001 > 0.1 | RBC_COUNT_PP001 > 0.1 | RETIC_COUNT_PP001 > 0.1)
ukbb.df.PLT <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(MPV_PP001 > 0.1 | PLT_COUNT_PP001 > 0.1)
ukbb.df.GRAN <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(BASO_COUNT_PP001 > 0.1 | EO_COUNT_PP001 > 0.1 | NEUTRO_COUNT_PP001 > 0.1)
ukbb.df.MONO <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(MONO_COUNT_PP001 > 0.1)
ukbb.df.LYMPH <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(LYMPH_COUNT_PP001 > 0.1)

# Identify pleiotropic variants
ukbb.df.PLEIT <- ukbb.df %>%
  rownames_to_column() %>%
  group_by(rowname) %>%
  mutate(ERY_PP = max(HCT_PP001, HGB_PP001, MCH_PP001, MCHC_PP001, MCV_PP001, MEAN_RETIC_VOL_PP001, RBC_COUNT_PP001, RETIC_COUNT_PP001),
         LYMPH_PP = max(LYMPH_COUNT_PP001),
         GRAN_PP = max(BASO_COUNT_PP001, EO_COUNT_PP001, NEUTRO_COUNT_PP001),
         PLT_PP = max(MPV_PP001, PLT_COUNT_PP001),
         MONO_PP = max(MONO_COUNT_PP001)) %>%
  mutate(ERY_PP_TF = (ERY_PP > 0.1),
         LYMPH_PP_TF = (LYMPH_PP > 0.1),
         GRAN_PP_TF = (GRAN_PP > 0.1),
         PLT_PP_TF = (PLT_PP > 0.1),
         MONO_PP_TF = (MONO_PP > 0.1)) %>%
  mutate(PLEIT = sum(ERY_PP_TF, LYMPH_PP_TF, GRAN_PP_TF, PLT_PP_TF, MONO_PP_TF)) %>%
  ungroup() %>%
  dplyr::filter(PLEIT >= 2)

# Combine 
ukbb.df.RBC$type <- "RBC"
ukbb.df.PLT$type <- "PLT"
ukbb.df.GRAN$type <- "GRAN"
ukbb.df.MONO$type <- "MONO"
ukbb.df.LYMPH$type <- "LYMPH"
ukbb.df.PLEIT$type <- "PLEIT"
ukbb.df.ALL <- do.call("rbind", list(ukbb.df.RBC[,c("rowname", "type")], ukbb.df.PLT[,c("rowname", "type")], ukbb.df.GRAN[,c("rowname", "type")], ukbb.df.MONO[,c("rowname", "type")], ukbb.df.LYMPH[,c("rowname", "type")], ukbb.df.PLEIT[,c("rowname", "type")]))
peaks.all <- peaks[as.numeric(ukbb.df.ALL$rowname),]
ATAC.cpm.log2.all <- ATAC.cpm.log2[as.numeric(ukbb.df.ALL$rowname),]

# Min / max scale
ATAC.cpm.log2.all.mm <- ATAC.cpm.log2.all / rowMax(ATAC.cpm.log2.all)
ATAC.cpm.log2.all.mm <- as.data.frame(ATAC.cpm.log2.all.mm)
ATAC.cpm.log2.all.mm$type <- ukbb.df.ALL$type
  
# Pull out the variable peaks in each lineage
ATAC.cpm.log2.all.mm.RBC <- ATAC.cpm.log2.all.mm[ATAC.cpm.log2.all.mm$type == "RBC", c("HSC", "MPP", "CMP", "MEP", "Ery")]
ATAC.cpm.log2.all.mm.PLT <- ATAC.cpm.log2.all.mm[ATAC.cpm.log2.all.mm$type == "PLT", c("HSC", "MPP", "CMP", "MEP", "Mega")]
ATAC.cpm.log2.all.mm.GRAN <- ATAC.cpm.log2.all.mm[ATAC.cpm.log2.all.mm$type == "GRAN", c("HSC", "MPP", "CMP", "GMP-A", "GMP-B", "GMP-C")]
ATAC.cpm.log2.all.mm.MONO <- ATAC.cpm.log2.all.mm[ATAC.cpm.log2.all.mm$type == "MONO", c("HSC", "MPP", "CMP", "GMP-A", "Mono")]
ATAC.cpm.log2.all.mm.LYMPH <- ATAC.cpm.log2.all.mm[ATAC.cpm.log2.all.mm$type == "LYMPH", c("HSC", "MPP", "LMPP", "CLP", "CD4", "CD8", "B", "NK")]
ATAC.cpm.log2.all.mm.PLEIT <- ATAC.cpm.log2.all.mm[ATAC.cpm.log2.all.mm$type == "PLEIT",]

# Plot stuff
pdf(file="plots/ery_cluster_peaks.pdf", width = 3, height = 8)  
par(cex.main=0.8,mar=c(1,1,1,1))
set.seed(123456)
km <- kmeans(ATAC.cpm.log2.all.mm.GRAN, centers = 5, nstart = 10000)
km.cluster <- factor(km$cluster)
#km.cluster <- factor(km$cluster, km$cluster, levels = c("2", "3", "1"))
hm <- Heatmap(ATAC.cpm.log2.all.mm.GRAN, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 6),
        split = km.cluster, show_heatmap_legend = FALSE,
        name = "Peak\nAccessibility")
hm
hm_ro <- row_order(hm)
hm_num <- unlist(lapply(hm_ro, length))
hm_num
dev.off()



