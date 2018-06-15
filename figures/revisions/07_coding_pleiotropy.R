library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(data.table)
library(annotables)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(BuenColors)
library(matrixStats)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
SNPdf <- read.table("pleiotropy/coding_lineages.tsv", sep = "\t", header = TRUE)

snp_g <- makeGRangesFromDataFrame(SNPdf, keep.extra.columns = TRUE)
loc <- locateVariants(snp_g, txdb, CodingVariants())

# Two merges one to annotate variants with genes; the second to get regular gene nmames
mdf <- merge(SNPdf, as.data.frame(loc))[,c("seqnames", "start", "end", "what", "GRAN", "LYMPH", "MONO", "PLT", "RBC", "GENEID")]
mdf <- mdf[!duplicated(mdf) & complete.cases(mdf),]
mdf2 <- merge(mdf, grch37[,c("entrez","symbol")], by.x = "GENEID", by.y = "entrez")


# Import RNA
x <- read.table("../../data/bulk/RNA/16populations_RNAcounts.txt", header = TRUE)
rownames(x) <- make.unique(as.character(x$Genes))
x <- x[,2:17]

counts <- x/colSums(x) * 1000000
log2cpm <- log2(counts + 1)

log2cpm_filt <- round(log2cpm[rownames(log2cpm)%in% unique(mdf2$symbol), ],2)
