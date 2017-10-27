library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(diffloop)

# Import LD Scores 
ldscorefiles <- list.files("../data/eur_w_ld_chr/", pattern = "*.gz$", full.names = TRUE)

ldsc_g <- lapply(ldscorefiles, function(file){
  dt <- fread(paste0("zcat < ",file))
  dt 
}) %>% rbindlist() %>% as.data.frame() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "CHR", start.field = "BP", end.field = "BP")


# Import narrow peak files
narrowpeakfiles <- list.files("../data/bulk/ATAC/narrowpeaks", pattern = "*.gz$", full.names = TRUE)

list_np_g <- lapply(narrowpeakfiles, function(file){
  dg <- makeGRangesFromDataFrame(data.frame(fread(paste0("zcat < ",file), header = FALSE)), 
                                 seqnames.field = "V1", start.field = "V2", end.field = "V3") %>% rmchr()
  dg 
})

names(list_np_g) <- gsub("_peaks.narrowPeak.gz", "", list.files("../data/bulk/ATAC/narrowpeaks/", pattern = "*.gz$"))

# Extract cell-type specific LD scores
lapply(1:length(list_np_g), function(i){
  peaks <- list_np_g[[i]]
  snpsidx <- findOverlaps(peaks, ldsc_g) %>% subjectHits() %>% unique()
  ldscs <- mcols(ldsc_g)$L2[snpsidx]
  df <- data.frame(celltype = names(list_np_g)[i], scores = ldscs)
  return(df)
}) %>% rbindlist() %>% as.data.frame() -> SNPs_LDscores_celltypes

p1 <- ggplot(SNPs_LDscores_celltypes, aes(x = celltype, y = scores, fill = celltype)) +
  geom_boxplot(outlier.colour=NA) + scale_fill_manual(values = ejc_color_maps) +
  pretty_plot() +
    theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), legend.position = "none",
    axis.text.x = element_text(angle = 90)) +labs(x = NULL, y = "Variant-Level LD Scores", 
    fill = "Cell Type") + scale_y_continuous(limits = c(-3, 80))
ggsave(p1, file = "variantLevelLDscores.pdf")

anova(lm(scores ~ as.factor(celltype), SNPs_LDscores_celltypes))
