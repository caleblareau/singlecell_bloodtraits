library(data.table)
library(GenomicRanges)
library(dplyr)

# Compute per-gene correlations with ATAC peaks

# Determine TSS per gene
gtf <-rtracklayer::import("RNA/hg19_10X.gtf.gz")
gs <- tapply(start(gtf), list(mcols(gtf)@listData$gene_id), min)
ge <- tapply(start(gtf), list(mcols(gtf)@listData$gene_id), max)
chrs <- data.frame(chr = seqnames(gtf), gene = mcols(gtf)@listData$gene_name) %>% unique()
strand <- data.frame(strand =strand(gtf), gene = mcols(gtf)@listData$gene_name) %>% unique()
remove(gtf)
coords <- data.frame(chrs, gs, ge, strand); coords$gene <- rownames(coords)