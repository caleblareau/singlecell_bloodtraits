library(data.table)
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

"%ni%" <- Negate("%in%")

df <- data.frame(fread("SupplementalTable1.tsv"))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
coding <- predictCoding(gr, txdb, BSgenome.Hsapiens.UCSC.hg19, DNAStringSet(mcols(gr)$alt))
cdf <- as.data.frame(coding)
cdf <- unique(cdf[,colnames(cdf) %ni% c("width", "strand", "QUERYID", "Beta", "SE", "Z", "In_Gwas_catalog", "In_Heme_peak")])

write.table(cdf, file= "CodingTable.txt", sep = "\t", row.names = FALSE, col.names = TRUE)