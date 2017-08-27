library(data.table)
library(GenomicRanges)
library(dplyr)
library(diffloop)

# ---------------------------------------------
# Compute per-gene correlations with ATAC peaks
# ---------------------------------------------

# Import genomic coordinates
peaks_g <- rmchr(makeGRangesFromDataFrame(fread("ATAC/26August2017_EJCsamples_allReads_250bp.bed"),
                                    seqnames.field = "V1", start.field = "V2", end.field = "V3"))
tss_g <- makeGRangesFromDataFrame(data.frame(fread("RNA/geneTSScoordaintes.txt")), keep.extra.columns = TRUE,
                                  seqnames.field = "V1", start.field = "V2", end.field = "V2")

# Import other raw data
rna <- fread("RNA/16populations_RNAcounts.txt")
genes <- rna[["Genes"]]
rna <- data.matrix(rna[,Genes:=NULL])
atac <- data.matrix(fread("ATAC/26August2017_EJCsamples_allReads_250bp.counts.txt"))[,colnames(rna)]

# Make it counts per million normalized
rna <- sweep(rna, 2, colSums(rna), FUN="/") * 1000000
atac <- sweep(atac, 2, colSums(atac), FUN="/") * 1000000

# And filter genes everywhere
keepGenes <- which(rowMeans(rna) > 1)
genes <- genes[keepGenes]
rna <- rna[keepGenes, ]
tss_g <- padGRanges(tss_g[mcols(tss_g)[,1] %in% genes], pad = 1000000)

# Define a function to get pvalue from correlation given a sample size
rhoToP <- function(rho, n = 16){
  t <- abs(rho / sqrt((1-rho)^2/(n-2)))
  2*pt(t, n-1, lower=FALSE)
}

# Loop over genes to compute correlation, pvalue, etc. 
lappout <- lapply(genes, function(gene){
  atacidx <- queryHits(findOverlaps(peaks_g, tss_g[which(mcols(tss_g)[,1] == gene)]))
  if(length(atacidx) > 1){
    rho <- cor(log2(rna[which(genes == gene)[1],]+1), t(log2(atac[atacidx,]+1)))[1,]
    return(data.frame(data.frame(peaks_g[atacidx])[,c(1,2,3)], gene = gene, rho = rho, p = rhoToP(rho)))
  } else {
    return(NULL)
  }
})

# Make final output
df <- rbindlist(Filter(Negate(is.null), setNames(lappout,seq_along(lappout))))
write.table(df, file = "peakGeneCorrelation.tsv", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


# Make TSS Data Frame
if(FALSE){
  
  # Determine TSS per gene
  gtf <-rtracklayer::import("RNA/hg19_10X.gtf.gz")
  gs <- tapply(start(gtf), list(mcols(gtf)@listData$gene_id), min)
  ge <- tapply(start(gtf), list(mcols(gtf)@listData$gene_id), max)
  chrs <- data.frame(chr = seqnames(gtf), gene = mcols(gtf)@listData$gene_id) %>% unique()
  strand <- data.frame(strand =strand(gtf), gene = mcols(gtf)@listData$gene_id) %>% unique()
  remove(gtf)
  coords <- data.frame(chrs, gs, ge, strand); coords$gene <- rownames(coords)
  
  # Import Common Gene Coordinates 
  cgc <- read.table("RNA/genes.tsv")
  stopifnot(all(coords$gene.1 == cgc$V1))
  
  coords$commonName <- cgc[,2]
  
  tssdf <- rbind(setNames(coords[coords$strand == "+", c(1,3,7)], c("chr", "start", "gene")), 
                 setNames(coords[coords$strand == "-", c(1,4,7)], c("chr", "start", "gene")))
  rownames(tssdf) <- NULL
  write.table(tssdf, file = "RNA/geneTSScoordaintes.txt",
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}


