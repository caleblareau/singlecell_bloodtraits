library(data.table)
library(GenomicRanges)
library(dplyr)
library(diffloop)

# ---------------------------------------------
# Compute per-gene correlations with ATAC peaks
# ---------------------------------------------

# Import genomic coordinates
peaks_g <- makeGRangesFromDataFrame(fread("ATAC/29August2017_EJCsamples_allReads_500bp.bed"),
                                    seqnames.field = "V1", start.field = "V2", end.field = "V3")
tss_g <- makeGRangesFromDataFrame(data.frame(fread("RNA/refGene_hg19_TSS.bed")), keep.extra.columns = TRUE,
                                  seqnames.field = "V1", start.field = "V2", end.field = "V2")

# Import other raw data
rna <- fread("RNA/16populations_RNAcounts.txt")
genes <- rna[["Genes"]]
rna <- data.matrix(rna[,Genes:=NULL])
atac <- data.matrix(fread("ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))[,colnames(rna)]

# Make it counts per million normalized
rna <- sweep(rna, 2, colSums(rna), FUN="/") * 1000000
atac <- sweep(atac, 2, colSums(atac), FUN="/") * 1000000

# And filter genes everywhere
keepGenes <- which(rowMeans(rna) > 1 & !is.na(match(genes, mcols(tss_g)$V5)))
length(keepGenes)
genes <- genes[keepGenes]
rna <- rna[keepGenes, ]
tss_g <- sort(padGRanges(tss_g[mcols(tss_g)[,3] %in% genes], pad = 1000000))

# Define a function to get pvalue from correlation given a sample size
rhoToP <- function(rho, n = 16){
  t <- abs(rho / sqrt((1-rho)^2/(n-2)))
  2*pt(t, n-1, lower=FALSE)
}

# Loop over genes to compute correlation, pvalue, etc. 
lappout <- lapply(genes, function(gene){
  atacidx <- unique(queryHits(findOverlaps(peaks_g, tss_g[which(mcols(tss_g)[,3] == gene)])))
  if(length(atacidx) > 1){
    rho <- cor(log2(rna[which(genes == gene)[1],]+1), t(log2(atac[atacidx,]+1)))[1,]
    return(data.frame(data.frame(peaks_g[atacidx])[,c(1,2,3)], gene = gene, rho = rho, p = rhoToP(rho)))
  } else {
    return(NULL)
  }
})

# Make final output
df <- rbindlist(Filter(Negate(is.null), setNames(lappout,seq_along(lappout))))
df$rho <- round(df$rho, digits = 3)
df$p <- formatC(df$p, format = "e", digits = 2)
write.table(df, file = "peakGeneCorrelation.tsv", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

<<<<<<< Updated upstream
# Make RDS
dt <- data.frame(fread(paste0("zcat < ","peakGeneCorrelation.tsv.gz")))
dt$V1 <- factor(dt$V1, levels = unique(dt$V1))
dt$V4 <- factor(dt$V4, levels = unique(dt$V4))
dt <- dt[dt$V6 < 0.01, ]
colnames(dt) <- c("chr", "start", "end", "gene", "correlation", "p-value")
saveRDS(dt, "peakGeneCorrelation.rds")
=======
# Read in output 
df <- fread(paste0("zcat < ", "../../data/bulk/peakGeneCorrelation.tsv.gz"))



>>>>>>> Stashed changes
