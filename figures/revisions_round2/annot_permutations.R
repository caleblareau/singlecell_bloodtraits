library(data.table)
library(dplyr)
library(rtracklayer)
library(preprocessCore)
library(ggjoy)
library(BuenColors)
library(regioneR)
library(stringr)
library(qvalue)


#' Helper function for shifting
randomizeLocalRegions <- function(A, ...) {
  A <- toGRanges(A)
  rand <- ceiling(runif(1,-1500000,1500000))
  B <- GenomicRanges::shift(A,rand)
  return(B)
}

# Read in FINEMAP dataset
CS.gr <- readRDS("../../data/Finemap/UKBB_BC_v3_VEPannotations.rds")
CS.df <- as.data.frame(CS.gr)

#' Read in count matrix 
counts.df <- read.table("../../data/bulk/ATAC/26August2017_EJCsamples_allReads_250bp.counts.txt",header=T)
names(counts.df) <- c("B","CD4","CD8","CLP","CMP","Ery","GMP-A","GMP-B","GMP-C","HSC","LMPP","mDC","Mega","MEP","mono","MPP","NK","pDC")
# Remove "weak" peaks that aren't in top n% for at least one cell type
n=0.8
counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- normalize.quantiles(as.matrix(counts.df))
keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))

#' Read in consensus bed file of peaks and add count metadata
peaks.gr <- import("../../data/bulk/ATAC/26August2017_EJCsamples_allReads_250bp.bed",format="bed")
values(peaks.gr) <- counts.df
peaks.gr <- peaks.gr[keep,]

#' Read in other annotations
# ChIP dataset

# Motifs

# PChiC 
pchic.gr <- readRDS("../../data/pchic/pchic.rds")
# ATAC-RNA
pg.df <- fread(paste0("zcat < ", "../../data/bulk/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","j.start","j.end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>% dplyr::filter(qvalue < 0.001)
pg.gr <- makeGRangesFromDataFrame(pg.df, 
                                  seqnames = "chrom", start.field = "j.start", end.field = "j.end")

# Compile all annotations into a list
gr.list <- list(peaks.gr,pchic.gr,pg.gr)

# Run permutations on every annotation
enrichments <- NULL
perm <- NULL
trait <- unique(CS.gr@elementMetadata$trait)
for (i in 1:length(gr.list)) {
  perm[[i]] <- permTest(A=CS.gr, B=gr.list[[i]], 
                        ntimes=10000, 
                        alternative="auto", 
                        randomize.function=randomizeLocalRegions, 
                        evaluate.function=numOverlaps, 
                        force.parallel=FALSE, 
                        mc.set.seed=FALSE, 
                        mc.cores=8)
  enrichments <- rbind(OR,c(i,
                   perm[[i]]$numOverlaps$zscore,
                   pnorm(perm[[i]]$numOverlaps$zscore, lower.tail = FALSE)))
  print(i)
}

enrich.df <- as.data.frame(enrichments,stringsAsFactors=F)
names(enrich.df) <- c("i","zscore","pvalue")
enrich.df$i <- as.factor(OR.df$i)
enrich.df$qvalue <- qvalue::qvalue(enrich.df$pvalue)$qvalues
