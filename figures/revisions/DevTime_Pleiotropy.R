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
library(tidyverse)
library(qvalue)
"%ni%" <- Negate("%in%")

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
peaks.gr <- peaks[keep,]
ukbb <- ukbb[keep,]
ukbb.df <- data.frame(as.matrix(assays(ukbb)[["weights"]]))

# Read in fine mapped variants
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)
CS.df$seqnames <- paste0("chr", CS.df$seqnames)
CS.gr <- GRanges(CS.df)

# UKBB exclusion list
exdf <- read.table("exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
CS.gr <- CS.gr[CS.gr$var %ni% exdf,]

# Subset to fine-mapped variants, exclude regions where sum(varPP) > 0.10 but no single variant has PP > 0.10
ukbb.df.RBC <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(HCT_PP001 > 0.1 | HGB_PP001 > 0.1 | MCH_PP001 > 0.1, MCHC_PP001 > 0.1 | MCV_PP001 > 0.1 | MEAN_RETIC_VOL_PP001 > 0.1 | RBC_COUNT_PP001 > 0.1 | RETIC_COUNT_PP001 > 0.1)
peaks.temp <- peaks[as.numeric(ukbb.df.RBC$rowname),]
peaks.keep <- unique(findOverlaps(peaks.temp, CS.gr[CS.gr$PP > 0.10 & CS.gr$trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")])@from)
ukbb.df.RBC <- ukbb.df.RBC[peaks.keep,]
ukbb.df.PLT <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(MPV_PP001 > 0.1 | PLT_COUNT_PP001 > 0.1)
peaks.temp <- peaks[as.numeric(ukbb.df.PLT$rowname),]
peaks.keep <- unique(findOverlaps(peaks.temp, CS.gr[CS.gr$PP > 0.10 & CS.gr$trait %in% c("MPV", "PLT_COUNT")])@from)
ukbb.df.PLT <- ukbb.df.PLT[peaks.keep,]
ukbb.df.GRAN <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(BASO_COUNT_PP001 > 0.1 | EO_COUNT_PP001 > 0.1 | NEUTRO_COUNT_PP001 > 0.1)
peaks.temp <- peaks[as.numeric(ukbb.df.GRAN$rowname),]
peaks.keep <- unique(findOverlaps(peaks.temp, CS.gr[CS.gr$PP > 0.10 & CS.gr$trait %in% c("BASO_COUNT", "EO_COUNT", "NEUTRO_COUNT")])@from)
ukbb.df.GRAN <- ukbb.df.GRAN[peaks.keep,]
ukbb.df.MONO <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(MONO_COUNT_PP001 > 0.1)
peaks.temp <- peaks[as.numeric(ukbb.df.MONO$rowname),]
peaks.keep <- unique(findOverlaps(peaks.temp, CS.gr[CS.gr$PP > 0.10 & CS.gr$trait %in% c("MONO_COUNT")])@from)
ukbb.df.MONO <- ukbb.df.MONO[peaks.keep,]
ukbb.df.LYMPH <- ukbb.df %>%
  rownames_to_column() %>%
  dplyr::filter(LYMPH_COUNT_PP001 > 0.1)
peaks.temp <- peaks[as.numeric(ukbb.df.LYMPH$rowname),]
peaks.keep <- unique(findOverlaps(peaks.temp, CS.gr[CS.gr$PP > 0.10 & CS.gr$trait %in% c("LYMPH_COUNT")])@from)
ukbb.df.LYMPH <- ukbb.df.LYMPH[peaks.keep,]

# Combine 
ukbb.df.RBC$type <- "RBC"
ukbb.df.PLT$type <- "PLT"
ukbb.df.GRAN$type <- "GRAN"
ukbb.df.MONO$type <- "MONO"
ukbb.df.LYMPH$type <- "LYMPH"
ukbb.df.ALL <- do.call("rbind", list(ukbb.df.RBC[,c("rowname", "type")], ukbb.df.PLT[,c("rowname", "type")], ukbb.df.GRAN[,c("rowname", "type")], ukbb.df.MONO[,c("rowname", "type")], ukbb.df.LYMPH[,c("rowname", "type")]))
peaks.all <- peaks[as.numeric(ukbb.df.ALL$rowname),]
ATAC.cpm.log2.all <- ATAC.cpm.log2[as.numeric(ukbb.df.ALL$rowname),]

# Min / max scale
ATAC.cpm.log2.all.mm <- ATAC.cpm.log2.all / rowMax(ATAC.cpm.log2.all)
ATAC.cpm.log2.all.mm <- as.data.frame(ATAC.cpm.log2.all.mm)
ATAC.cpm.log2.all.mm$type <- ukbb.df.ALL$type
  
# Pull out the variable peaks in each lineage
ATAC.cpm.log2.all.mm.RBC <- ATAC.cpm.log2.all.mm %>%
  dplyr::filter(type == "RBC") %>%
  dplyr::select("HSC", "MPP", "CMP", "MEP", "Ery")
ATAC.cpm.log2.all.mm.PLT <- ATAC.cpm.log2.all.mm %>%
  dplyr::filter(type == "PLT") %>%
  dplyr::select(c("HSC", "MPP", "CMP", "MEP", "Mega"))
ATAC.cpm.log2.all.mm.GRAN <- ATAC.cpm.log2.all.mm %>%
  dplyr::filter(type == "GRAN") %>%
  dplyr::select(c("HSC", "MPP", "CMP", "GMP-A", "GMP-B", "GMP-C"))
ATAC.cpm.log2.all.mm.MONO <- ATAC.cpm.log2.all.mm %>%
  dplyr::filter(type == "MONO") %>%
  dplyr::select(c("HSC", "MPP", "CMP", "GMP-A", "Mono"))
ATAC.cpm.log2.all.mm.LYMPH <- ATAC.cpm.log2.all.mm %>%
  dplyr::filter(type == "LYMPH") %>%
  dplyr::select(c("HSC", "MPP", "LMPP", "CLP", "CD4", "CD8", "B", "NK"))

# Load protein coding annotations
library(annotables)
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")

# Read in PCHIC
pchic <- readRDS("../../data/pchic/pchic.rds")
pchic.gr <- GRanges(pchic)
pchic.gr <- pchic.gr[pchic.gr$Gene %in% grch38.pc$symbol,]
# Read in PCHIC
pchic_cd34 <- readRDS("../../data/pchic/cd34_loops.rds")
pchic_cd34.gr <- GRanges(pchic_cd34)
pchic_cd34.gr <- pchic_cd34.gr[pchic_cd34.gr$gene %in% grch38.pc$symbol,]

# Read in enhancer gene correlations
pg.df <- fread(paste0("zcat < ", "../../data/bulk/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","j.start","j.end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>%
  dplyr::filter(qvalue < 0.001)

#' Interesect peaks, variants, peak to gene correlations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
peaks.df <- as.data.frame(peaks.gr)
peaks.df$peak <- paste0(peaks.df$seqnames,":",peaks.df$start,"-",peaks.df$end)
setkey(setDT(peaks.df), seqnames, start, end)
setkey(setDT(CS.df), seqnames, start, end)
setkey(setDT(pg.df), chrom, j.start, j.end)
peaks.CS.df <- foverlaps(peaks.df,CS.df, nomatch = 0)
setkey(setDT(peaks.CS.df), seqnames, start, end)
peaks.CS.pg.dt <- foverlaps(peaks.CS.df,pg.df, nomatch = 0)
temp <- peaks.CS.pg.dt %>% 
  dplyr::filter(PP > 0.1) %>%
  dplyr::select(trait,var,PP,region,peak,gene,cor,pvalue,B,CD4,CD8,CLP,CMP,Ery,GMP.A,GMP.B,GMP.C,HSC,LMPP,mDC,Mega,MEP,mono,MPP,NK,pDC)
CS.ATAC.cor <- peaks.CS.pg.dt %>%
  dplyr::filter(PP > 0.1) %>%
  dplyr::select(var, trait, PP, seqnames, start, end, gene, cor, qvalue)

# Merge with molecular mechanisms
# RBC ChIP + motif
peaks.RBC <- peaks.all[ukbb.df.ALL$type == "RBC",]
ALL.M_C.df.match.RBC.gr <- ALL.M_C.df.match %>%
  dplyr::filter(trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.RBC, ALL.M_C.df.match.RBC.gr)
ATAC.cpm.log2.all.mm.RBC.motifbreakannot <- data.frame("mc" = ifelse(seq(1:length(peaks.RBC)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.RBC.motif.peaks <- data.frame(cbind(as.data.frame(peaks.RBC[idx@from]), as.data.frame(ALL.M_C.df.match.RBC.gr[idx@to])))
look <- ATAC.cpm.log2.all.mm.RBC.motif.peaks %>%
  dplyr::select(SNP, trait, PP, antigen, celltype, geneSymbol, effect)
var_RBC <- c("chr4:55408875:A:T", "chr5:173287763:G:A", "chr6:41925159:G:A", "chr9:4852599:A:G", "chr16:87886545:C:T", "chr22:50964153:T:C", "chr22:50978262:G:A")
# RBC ChIP + motif proximity
ALL.MM_C.df.match.RBC.gr <- ALL.MM_C.df.match %>%
  dplyr::filter(trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.RBC, ALL.MM_C.df.match.RBC.gr)
ATAC.cpm.log2.all.mm.RBC.motifmatchannot <- data.frame("mc" = ifelse(seq(1:length(peaks.RBC)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.RBC.motifannot <- ifelse(ATAC.cpm.log2.all.mm.RBC.motifbreakannot == 1, 2, 0)
ATAC.cpm.log2.all.mm.RBC.motifannot[ATAC.cpm.log2.all.mm.RBC.motifmatchannot > ATAC.cpm.log2.all.mm.RBC.motifannot] <- 1
# RBC PCHIC
pchic.gr.RBC <- pchic.gr[pchic.gr$CellType == "Ery" & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT"),]
idx <- findOverlaps(peaks.RBC, pchic.gr.RBC)
ATAC.cpm.log2.all.mm.RBC.PCHIC1annot <- data.frame("mc" = ifelse(seq(1:length(peaks.RBC)) %in% unique(idx@from), 1, 0))
idx <- findOverlaps(peaks.RBC, pchic_cd34.gr)
ATAC.cpm.log2.all.mm.RBC.PCHIC2annot <- data.frame("mc" = ifelse(seq(1:length(peaks.RBC)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.RBC.PCHICannot <- apply(cbind(ATAC.cpm.log2.all.mm.RBC.PCHIC1annot, ATAC.cpm.log2.all.mm.RBC.PCHIC2annot), 1, max)
# RBC ATAC-RNA correlation
CS.ATAC.cor.RBC <- CS.ATAC.cor %>%
  dplyr::filter(trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.RBC, CS.ATAC.cor.RBC)
ATAC.cpm.log2.all.mm.RBC.ARCORannot <- data.frame("mc" = ifelse(seq(1:length(peaks.RBC)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.RBC.geneannot <- ATAC.cpm.log2.all.mm.RBC.PCHICannot + ATAC.cpm.log2.all.mm.RBC.ARCORannot

# Plot RBC
pdf(file="plots/ery_cluster_peaks.pdf", width = 7, height = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
set.seed(123456)
km <- kmeans(ATAC.cpm.log2.all.mm.RBC, centers = 5, nstart = 10000) # RBC = 5, LYMPH = 5, MONO = 4, GRAN = 5, PLT = 6
km.cluster <- factor(km$cluster, levels = c("3", "4", "1", "2", "5")) # RBC
hm <- Heatmap(ATAC.cpm.log2.all.mm.RBC, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
              cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 6),
              split = km.cluster, show_heatmap_legend = FALSE,
              name = "Peak\nAccessibility")
ha1 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.RBC.motifannot, col = list(mc = c("2" = jdb_palette("solar_basic")[1], "1" = jdb_palette("solar_basic")[5], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
ha2 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.RBC.geneannot, col = list(mc = c("2" = jdb_palette("ocean_green")[9], "1" = jdb_palette("ocean_green")[6], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
hm + ha1 + ha2
hm_ro <- row_order(hm)
hm_num <- unlist(lapply(hm_ro, length))
hm_num
dev.off()

# LYMPH ChIP + motif
peaks.LYMPH <- peaks.all[ukbb.df.ALL$type == "LYMPH",]
ALL.M_C.df.match.LYMPH.gr <- ALL.M_C.df.match %>%
  dplyr::filter(trait %in% c("LYMPH_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.LYMPH, ALL.M_C.df.match.LYMPH.gr)
ATAC.cpm.log2.all.mm.LYMPH.motifbreakannot <- data.frame("mc" = ifelse(seq(1:length(peaks.LYMPH)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.LYMPH.motif.peaks <- data.frame(cbind(as.data.frame(peaks.LYMPH[idx@from]), as.data.frame(ALL.M_C.df.match.LYMPH.gr[idx@to])))
look <- ATAC.cpm.log2.all.mm.LYMPH.motif.peaks %>%
  dplyr::select(SNP, trait, PP, antigen, celltype, geneSymbol, effect)
var_LYMPH <- c("chr2:112238447:G:A", "chr2:143886819:G:A", "chr4:38132999:A:G", "chr14:55368986:G:C", "chr16:30482494:T:C", "chr17:37912377:C:T")
# LYMPH ChIP + motif proximity
ALL.MM_C.df.match.LYMPH.gr <- ALL.MM_C.df.match %>%
  dplyr::filter(trait %in% c("LYMPH_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.LYMPH, ALL.MM_C.df.match.LYMPH.gr)
ATAC.cpm.log2.all.mm.LYMPH.motifmatchannot <- data.frame("mc" = ifelse(seq(1:length(peaks.LYMPH)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.LYMPH.motifannot <- ifelse(ATAC.cpm.log2.all.mm.LYMPH.motifbreakannot == 1, 2, 0)
ATAC.cpm.log2.all.mm.LYMPH.motifannot[ATAC.cpm.log2.all.mm.LYMPH.motifmatchannot > ATAC.cpm.log2.all.mm.LYMPH.motifannot] <- 1
# LYMPH PCHIC
pchic.gr.LYMPH <- pchic.gr[pchic.gr$CellType %in% c("nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB") & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% c("LYMPH_COUNT"),]
idx <- findOverlaps(peaks.LYMPH, pchic.gr.LYMPH)
ATAC.cpm.log2.all.mm.LYMPH.PCHIC1annot <- data.frame("mc" = ifelse(seq(1:length(peaks.LYMPH)) %in% unique(idx@from), 1, 0))
idx <- findOverlaps(peaks.LYMPH, pchic_cd34.gr)
ATAC.cpm.log2.all.mm.LYMPH.PCHIC2annot <- data.frame("mc" = ifelse(seq(1:length(peaks.LYMPH)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.LYMPH.PCHICannot <- apply(cbind(ATAC.cpm.log2.all.mm.LYMPH.PCHIC1annot, ATAC.cpm.log2.all.mm.LYMPH.PCHIC2annot), 1, max)
# LYMPH ATAC-RNA correlation
CS.ATAC.cor.LYMPH <- CS.ATAC.cor %>%
  dplyr::filter(trait %in% c("LYMPH_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.LYMPH, CS.ATAC.cor.LYMPH)
ATAC.cpm.log2.all.mm.LYMPH.ARCORannot <- data.frame("mc" = ifelse(seq(1:length(peaks.LYMPH)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.LYMPH.geneannot <- ATAC.cpm.log2.all.mm.LYMPH.PCHICannot + ATAC.cpm.log2.all.mm.LYMPH.ARCORannot

# Plot LYMPH
pdf(file="plots/lymph_cluster_peaks.pdf", width = 7, height = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
set.seed(123456)
km <- kmeans(ATAC.cpm.log2.all.mm.LYMPH, centers = 5, nstart = 10000) # LYMPH = 5, LYMPH = 5, MONO = 4, GRAN = 5, PLT = 6
km.cluster <- factor(km$cluster, levels = c("3", "4", "1", "5", "2")) # LYMPH
hm <- Heatmap(ATAC.cpm.log2.all.mm.LYMPH, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
              cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 6),
              split = km.cluster, show_heatmap_legend = FALSE,
              name = "Peak\nAccessibility")
ha1 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.LYMPH.motifannot, col = list(mc = c("2" = jdb_palette("solar_basic")[1], "1" = jdb_palette("solar_basic")[5], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
ha2 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.LYMPH.geneannot, col = list(mc = c("2" = jdb_palette("ocean_green")[9], "1" = jdb_palette("ocean_green")[6], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
hm + ha1 + ha2
hm_ro <- row_order(hm)
hm_num <- unlist(lapply(hm_ro, length))
hm_num
dev.off()

# MONO ChIP + motif
peaks.MONO <- peaks.all[ukbb.df.ALL$type == "MONO",]
ALL.M_C.df.match.MONO.gr <- ALL.M_C.df.match %>%
  dplyr::filter(trait %in% c("MONO_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.MONO, ALL.M_C.df.match.MONO.gr)
ATAC.cpm.log2.all.mm.MONO.motifbreakannot <- data.frame("mc" = ifelse(seq(1:length(peaks.MONO)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.MONO.motif.peaks <- data.frame(cbind(as.data.frame(peaks.MONO[idx@from]), as.data.frame(ALL.M_C.df.match.MONO.gr[idx@to])))
look <- ATAC.cpm.log2.all.mm.MONO.motif.peaks %>%
  dplyr::select(SNP, trait, PP, antigen, celltype, geneSymbol, effect)
#var_MONO <- c("")
# MONO ChIP + motif proximity
ALL.MM_C.df.match.MONO.gr <- ALL.MM_C.df.match %>%
  dplyr::filter(trait %in% c("MONO_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.MONO, ALL.MM_C.df.match.MONO.gr)
ATAC.cpm.log2.all.mm.MONO.motifmatchannot <- data.frame("mc" = ifelse(seq(1:length(peaks.MONO)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.MONO.motifannot <- ifelse(ATAC.cpm.log2.all.mm.MONO.motifbreakannot == 1, 2, 0)
ATAC.cpm.log2.all.mm.MONO.motifannot[ATAC.cpm.log2.all.mm.MONO.motifmatchannot > ATAC.cpm.log2.all.mm.MONO.motifannot] <- 1
# MONO PCHIC
pchic.gr.MONO <- pchic.gr[pchic.gr$CellType == c("Mon", "Mac0", "Mac1") & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% c("MONO_COUNT"),]
idx <- findOverlaps(peaks.MONO, pchic.gr.MONO)
ATAC.cpm.log2.all.mm.MONO.PCHIC1annot <- data.frame("mc" = ifelse(seq(1:length(peaks.MONO)) %in% unique(idx@from), 1, 0))
idx <- findOverlaps(peaks.MONO, pchic_cd34.gr)
ATAC.cpm.log2.all.mm.MONO.PCHIC2annot <- data.frame("mc" = ifelse(seq(1:length(peaks.MONO)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.MONO.PCHICannot <- apply(cbind(ATAC.cpm.log2.all.mm.MONO.PCHIC1annot, ATAC.cpm.log2.all.mm.MONO.PCHIC2annot), 1, max)
# MONO ATAC-RNA correlation
CS.ATAC.cor.MONO <- CS.ATAC.cor %>%
  dplyr::filter(trait %in% c("MONO_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.MONO, CS.ATAC.cor.MONO)
ATAC.cpm.log2.all.mm.MONO.ARCORannot <- data.frame("mc" = ifelse(seq(1:length(peaks.MONO)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.MONO.geneannot <- ATAC.cpm.log2.all.mm.MONO.PCHICannot + ATAC.cpm.log2.all.mm.MONO.ARCORannot

# Plot MONO
pdf(file="plots/MONO_cluster_peaks.pdf", width = 7, height = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
set.seed(123456)
km <- kmeans(ATAC.cpm.log2.all.mm.MONO, centers = 4, nstart = 10000) # MONO = 5, LYMPH = 5, MONO = 4, GRAN = 5, PLT = 6
km.cluster <- factor(km$cluster, levels = c("2", "4", "1", "3")) # MONO
hm <- Heatmap(ATAC.cpm.log2.all.mm.MONO, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
              cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 6),
              split = km.cluster, show_heatmap_legend = FALSE,
              name = "Peak\nAccessibility")
ha1 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.MONO.motifannot, col = list(mc = c("2" = jdb_palette("solar_basic")[1], "1" = jdb_palette("solar_basic")[5], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
ha2 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.MONO.geneannot, col = list(mc = c("2" = jdb_palette("ocean_green")[9], "1" = jdb_palette("ocean_green")[6], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
hm + ha1 + ha2
hm_ro <- row_order(hm)
hm_num <- unlist(lapply(hm_ro, length))
hm_num
dev.off()

# PLT ChIP + motif
peaks.PLT <- peaks.all[ukbb.df.ALL$type == "PLT",]
ALL.M_C.df.match.PLT.gr <- ALL.M_C.df.match %>%
  dplyr::filter(trait %in% c("MPV", "PLT_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.PLT, ALL.M_C.df.match.PLT.gr)
ATAC.cpm.log2.all.mm.PLT.motifbreakannot <- data.frame("mc" = ifelse(seq(1:length(peaks.PLT)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.PLT.motif.peaks <- data.frame(cbind(as.data.frame(peaks.PLT[idx@from]), as.data.frame(ALL.M_C.df.match.PLT.gr[idx@to])))
look <- ATAC.cpm.log2.all.mm.PLT.motif.peaks %>%
  dplyr::select(SNP, trait, PP, antigen, celltype, geneSymbol, effect)
var_PLT <- c("chr1:31241886:G:T", "chr1:247712303:T:C", "chr3:72397279:A:G", "chr9:4763491:G:A", "chr9:136128000:G:C", "chr11:113982321:T:A", "chr12:29435675:A:G", "chr13:95899607:A:G", "chr13:95899643:C:T", "chr16:85415838:T:C", "chr16:85450100:G:C", "chr22:50628937:G:A")
# PLT ChIP + motif proximity
ALL.MM_C.df.match.PLT.gr <- ALL.MM_C.df.match %>%
  dplyr::filter(trait %in% c("MPV", "PLT_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.PLT, ALL.MM_C.df.match.PLT.gr)
ATAC.cpm.log2.all.mm.PLT.motifmatchannot <- data.frame("mc" = ifelse(seq(1:length(peaks.PLT)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.PLT.motifannot <- ifelse(ATAC.cpm.log2.all.mm.PLT.motifbreakannot == 1, 2, 0)
ATAC.cpm.log2.all.mm.PLT.motifannot[ATAC.cpm.log2.all.mm.PLT.motifmatchannot > ATAC.cpm.log2.all.mm.PLT.motifannot] <- 1
# PLT PCHIC
pchic.gr.PLT <- pchic.gr[pchic.gr$CellType == "MK" & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% c("MPV", "PLT_COUNT"),]
idx <- findOverlaps(peaks.PLT, pchic.gr.PLT)
ATAC.cpm.log2.all.mm.PLT.PCHIC1annot <- data.frame("mc" = ifelse(seq(1:length(peaks.PLT)) %in% unique(idx@from), 1, 0))
idx <- findOverlaps(peaks.PLT, pchic_cd34.gr)
ATAC.cpm.log2.all.mm.PLT.PCHIC2annot <- data.frame("mc" = ifelse(seq(1:length(peaks.PLT)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.PLT.PCHICannot <- apply(cbind(ATAC.cpm.log2.all.mm.PLT.PCHIC1annot, ATAC.cpm.log2.all.mm.PLT.PCHIC2annot), 1, max)
# PLT ATAC-RNA correlation
CS.ATAC.cor.PLT <- CS.ATAC.cor %>%
  dplyr::filter(trait %in% c("MPV", "PLT_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.PLT, CS.ATAC.cor.PLT)
ATAC.cpm.log2.all.mm.PLT.ARCORannot <- data.frame("mc" = ifelse(seq(1:length(peaks.PLT)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.PLT.geneannot <- ATAC.cpm.log2.all.mm.PLT.PCHICannot + ATAC.cpm.log2.all.mm.PLT.ARCORannot

# Plot PLT
pdf(file="plots/plt_cluster_peaks.pdf", width = 7, height = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
set.seed(123456)
km <- kmeans(ATAC.cpm.log2.all.mm.PLT, centers = 5, nstart = 10000) # PLT = 5, LYMPH = 5, MONO = 4, GRAN = 5, PLT = 6
km.cluster <- factor(km$cluster, levels = c("2", "4", "5", "3", "1")) # PLT
hm <- Heatmap(ATAC.cpm.log2.all.mm.PLT, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
              cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 6),
              split = km.cluster, show_heatmap_legend = FALSE,
              name = "Peak\nAccessibility")
ha1 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.PLT.motifannot, col = list(mc = c("2" = jdb_palette("solar_basic")[1], "1" = jdb_palette("solar_basic")[5], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
ha2 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.PLT.geneannot, col = list(mc = c("2" = jdb_palette("ocean_green")[9], "1" = jdb_palette("ocean_green")[6], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
hm + ha1 + ha2
hm_ro <- row_order(hm)
hm_num <- unlist(lapply(hm_ro, length))
hm_num
dev.off()

# GRAN ChIP + motif
peaks.GRAN <- peaks.all[ukbb.df.ALL$type == "GRAN",]
ALL.M_C.df.match.GRAN.gr <- ALL.M_C.df.match %>%
  dplyr::filter(trait %in% c("BASO_COUNT", "EO_COUNT", "NEUTRO_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.GRAN, ALL.M_C.df.match.GRAN.gr)
ATAC.cpm.log2.all.mm.GRAN.motifbreakannot <- data.frame("mc" = ifelse(seq(1:length(peaks.GRAN)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.GRAN.motif.peaks <- data.frame(cbind(as.data.frame(peaks.GRAN[idx@from]), as.data.frame(ALL.M_C.df.match.GRAN.gr[idx@to])))
look <- ATAC.cpm.log2.all.mm.GRAN.motif.peaks %>%
  dplyr::select(SNP, trait, PP, antigen, celltype, geneSymbol, effect)
#var_GRAN <- 
# GRAN ChIP + motif proximity
ALL.MM_C.df.match.GRAN.gr <- ALL.MM_C.df.match %>%
  dplyr::filter(trait %in% c("BASO_COUNT", "EO_COUNT", "NEUTRO_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.GRAN, ALL.MM_C.df.match.GRAN.gr)
ATAC.cpm.log2.all.mm.GRAN.motifmatchannot <- data.frame("mc" = ifelse(seq(1:length(peaks.GRAN)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.GRAN.motifannot <- ifelse(ATAC.cpm.log2.all.mm.GRAN.motifbreakannot == 1, 2, 0)
ATAC.cpm.log2.all.mm.GRAN.motifannot[ATAC.cpm.log2.all.mm.GRAN.motifmatchannot > ATAC.cpm.log2.all.mm.GRAN.motifannot] <- 1
# GRAN PCHIC
pchic.gr.GRAN <- pchic.gr[pchic.gr$CellType == "Neu" & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% c("BASO_COUNT", "EO_COUNT", "NEUTRO_COUNT"),]
idx <- findOverlaps(peaks.GRAN, pchic.gr.GRAN)
ATAC.cpm.log2.all.mm.GRAN.PCHIC1annot <- data.frame("mc" = ifelse(seq(1:length(peaks.GRAN)) %in% unique(idx@from), 1, 0))
idx <- findOverlaps(peaks.GRAN, pchic_cd34.gr)
ATAC.cpm.log2.all.mm.GRAN.PCHIC2annot <- data.frame("mc" = ifelse(seq(1:length(peaks.GRAN)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.GRAN.PCHICannot <- apply(cbind(ATAC.cpm.log2.all.mm.GRAN.PCHIC1annot, ATAC.cpm.log2.all.mm.GRAN.PCHIC2annot), 1, max)
# GRAN ATAC-RNA correlation
CS.ATAC.cor.GRAN <- CS.ATAC.cor %>%
  dplyr::filter(trait %in% c("BASO_COUNT", "EO_COUNT", "NEUTRO_COUNT")) %>%
  GRanges()
idx <- findOverlaps(peaks.GRAN, CS.ATAC.cor.GRAN)
ATAC.cpm.log2.all.mm.GRAN.ARCORannot <- data.frame("mc" = ifelse(seq(1:length(peaks.GRAN)) %in% unique(idx@from), 1, 0))
ATAC.cpm.log2.all.mm.GRAN.geneannot <- ATAC.cpm.log2.all.mm.GRAN.PCHICannot + ATAC.cpm.log2.all.mm.GRAN.ARCORannot

# Plot GRAN
pdf(file="plots/GRAN_cluster_peaks.pdf", width = 7, height = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
set.seed(123456)
km <- kmeans(ATAC.cpm.log2.all.mm.GRAN, centers = 4, nstart = 10000) # GRAN = 5, LYMPH = 5, MONO = 4, GRAN = 5, GRAN = 6
km.cluster <- factor(km$cluster, levels = c("2", "3", "4", "1")) # GRAN
hm <- Heatmap(ATAC.cpm.log2.all.mm.GRAN, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
              cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 6),
              split = km.cluster, show_heatmap_legend = FALSE,
              name = "Peak\nAccessibility")
ha1 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.GRAN.motifannot, col = list(mc = c("2" = jdb_palette("solar_basic")[1], "1" = jdb_palette("solar_basic")[5], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
ha2 <- rowAnnotation(df = ATAC.cpm.log2.all.mm.GRAN.geneannot, col = list(mc = c("2" = jdb_palette("ocean_green")[9], "1" = jdb_palette("ocean_green")[6], "0" = "white")), width = unit(0.8, "cm"), show_legend = FALSE)
hm + ha1 + ha2
hm_ro <- row_order(hm)
hm_num <- unlist(lapply(hm_ro, length))
hm_num
dev.off()

# Target gene summaries
sum(table(ATAC.cpm.log2.all.mm.RBC.geneannot)[2:3]) + sum(table(ATAC.cpm.log2.all.mm.PLT.geneannot)[2:3]) + sum(table(ATAC.cpm.log2.all.mm.LYMPH.geneannot)[2:3]) + sum(table(ATAC.cpm.log2.all.mm.MONO.geneannot)[2:3]) + sum(table(ATAC.cpm.log2.all.mm.GRAN.geneannot)[2:3]) + sum(table(ATAC.cpm.log2.all.mm.RBC.geneannot)[2:3]) 
sum(table(ATAC.cpm.log2.all.mm.RBC.geneannot)[1:3]) + sum(table(ATAC.cpm.log2.all.mm.PLT.geneannot)[1:3]) + sum(table(ATAC.cpm.log2.all.mm.LYMPH.geneannot)[1:3]) + sum(table(ATAC.cpm.log2.all.mm.MONO.geneannot)[1:3]) + sum(table(ATAC.cpm.log2.all.mm.GRAN.geneannot)[1:3]) + sum(table(ATAC.cpm.log2.all.mm.RBC.geneannot)[1:3]) 

# Target gene examples
idx <- findOverlaps(peaks.RBC, pchic.gr.RBC[pchic.gr.RBC$PP > 0.5,])
pchic.gene <- unique(pchic.gr.RBC[idx@to]$Gene)
idx <- findOverlaps(peaks.RBC, CS.ATAC.cor.RBC[CS.ATAC.cor.RBC$PP > 0.5,])
cor.gene <- unique(CS.ATAC.cor.RBC[idx@to]$gene)
intersect(pchic.gene, cor.gene)

idx <- findOverlaps(peaks.PLT, pchic.gr.PLT[pchic.gr.PLT$PP > 0.1,])
pchic.gene <- unique(pchic.gr.PLT[idx@to]$Gene)
idx <- findOverlaps(peaks.PLT, CS.ATAC.cor.PLT[CS.ATAC.cor.PLT$PP > 0.1,])
cor.gene <- unique(CS.ATAC.cor.PLT[idx@to]$gene)
intersect(pchic.gene, cor.gene)
#cat(unique(c(pchic.gene, cor.gene)), " ")

idx <- findOverlaps(peaks.MONO, pchic.gr.MONO[pchic.gr.MONO$PP > 0.1,])
pchic.gene <- unique(pchic.gr.MONO[idx@to]$Gene)
idx <- findOverlaps(peaks.MONO, CS.ATAC.cor.MONO[CS.ATAC.cor.MONO$PP > 0.1,])
cor.gene <- unique(CS.ATAC.cor.MONO[idx@to]$gene)
intersect(pchic.gene, cor.gene)

idx <- findOverlaps(peaks.GRAN, pchic.gr.GRAN)
pchic.gene <- unique(pchic.gr.GRAN[idx@to]$Gene)
idx <- findOverlaps(peaks.GRAN, CS.ATAC.cor.GRAN)
cor.gene <- unique(CS.ATAC.cor.GRAN[idx@to]$gene)
intersect(pchic.gene, cor.gene)

idx <- findOverlaps(peaks.LYMPH, pchic.gr.LYMPH[pchic.gr.LYMPH$PP > 0.5,])
pchic.gene <- unique(pchic.gr.LYMPH[idx@to]$Gene)
idx <- findOverlaps(peaks.LYMPH, CS.ATAC.cor.LYMPH[CS.ATAC.cor.LYMPH$PP > 0.5,])
cor.gene <- unique(CS.ATAC.cor.LYMPH[idx@to]$gene)
intersect(pchic.gene, cor.gene)

# RBC PCHIC get the genes
pchic.gr.overlap <- findOverlaps(pchic.gr.RBC, peaks.RBC)
length(unique(pchic.gr.overlap@to))
temp2 <- pchic.gr.RBC[unique(pchic.gr.overlap@from)] %>%
  as.data.frame() %>% 
  group_by(variantPos) %>%
  dplyr::filter(Value == max(Value)) %>%
  dplyr::filter(Gene == min(Gene)) %>%
  .$Gene %>%
  unique() %>%
  sort()
  
