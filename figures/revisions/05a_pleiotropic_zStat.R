library(gchromVAR)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(matrixStats)
library(data.table)

df <- read.table("pleiotropy/all.tsv", header = TRUE)
df <- df[df$class%in% c("accessible"),]
pleio_gr <- makeGRangesFromDataFrame(df)

# Import counts and normalize
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt")) 
ATAC.cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
ov <- findOverlaps(peaks, pleio_gr)

progenitors <- as.vector(ATAC.cpm[queryHits(ov),c("HSC", "CLP", "CMP", "GMP-A", "GMP-B", "GMP-C", "LMPP", "MEP", "MPP")])
committed <- as.vector(ATAC.cpm[queryHits(ov),c("B", "CD4", "CD8", "Ery", "mDC", "Mega", "Mono", "NK", "pDC")])

t.test(progenitors, committed)
