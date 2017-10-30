#' ---
#' title: "Descriptive analyses of fine mapped variants"
#' author: "Jacob Ulirsch"
#' date: "July 31, 2017"
#' ---

library(data.table)
library(dplyr)
library(rtracklayer)
library(preprocessCore)
library(ggjoy)
library(BuenColors)
library(regioneR)
library(stringr)

#' Read in count matrix 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
counts.df <- read.table("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt",header=T)
names(counts.df) <- c("B","CD4","CD8","CLP","CMP","Ery","GMP-A","GMP-B","GMP-C","HSC","LMPP","mDC","Mega","MEP","mono","MPP","NK","pDC")
# Remove "weak" peaks that aren't in top n% for at least one cell type
n=0.8
counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- normalize.quantiles(as.matrix(counts.df))
keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))

#' Read in consensus bed file of peaks and add count metadata
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
peaks.gr <- import("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed",format="bed")
values(peaks.gr) <- counts.df
peaks.gr <- peaks.gr[keep,]

# Read in fine mapped variants
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)
CS.df$seqnames <- paste0("chr",CS.df$seqnames)

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
  dplyr::filter(PP > 0.5) %>%
  dplyr::select(trait,var,PP,region,peak,gene,cor,pvalue,B,CD4,CD8,CLP,CMP,Ery,GMP.A,GMP.B,GMP.C,HSC,LMPP,mDC,Mega,MEP,mono,MPP,NK,pDC)







CS.df$name <- gsub("_\\d*$","",CS.df$name)
setkey(setDT(peaks.df), seqnames, start, end)
setkey(setDT(CS.df), seqnames, start, end)
peaks.CS.df <- foverlaps(peaks.df,CS.df, nomatch = 0)
peaks.CS.tbl.df <- as.data.frame.matrix(table(peaks.CS.df$SNV,peaks.CS.df$name))
peaks.CS.tbl.df[peaks.CS.tbl.df >= 1] <- 1
peaks.CS.all.df <- peaks.CS.df %>% as.data.frame() %>% dplyr::select(-seqnames,-start,-end,-width,-strand,-i.start,-i.end,-i.width,-i.strand) %>% dplyr::select(SNV, region, peak, name, score, everything()) 
peaks.CS.df <- peaks.CS.all.df %>% group_by(SNV) %>% filter(row_number(SNV) == 1) %>% ungroup()
# In case of overlapping region remove duplicate SNVs
peaks.CS.annot.df <- peaks.CS.df[,1:5]
peaks.CS.data.df <- peaks.CS.df[,6:25]
#peaks.CS.data.df <- peaks.CS.df[,6:19]
peaks.CS.tbl.df <- peaks.CS.tbl.df[peaks.CS.annot.df$SNV,]

#peaks.CS.tbl.pleit.df <- data.frame(cbind(pleitropy=as.vector(rowSums(peaks.CS.tbl.df)),open=row.names(peaks.CS.tbl.df))

#' Min-max scale
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
rowmin.vec <- apply(peaks.CS.data.df,1,min)
rowmax.vec <- apply(peaks.CS.data.df,1,max)
peaks.CS.data.mm.df <- sweep(peaks.CS.data.df,1,rowmin.vec,"-")
peaks.CS.data.mm.df <- sweep(peaks.CS.data.mm.df,1,rowmax.vec - rowmin.vec,"/")
peaks.CS.data.mm.mat <- as.matrix(peaks.CS.data.mm.df)

#' Heatmap
#+ echo=FALSE, message=FALSE, warning=FALSE
main_heatmap(peaks.CS.data.mm.mat, y=peaks.CS.annot.df$SNV, name="High Confidence <br> Regulatory SNVs") %>%
  add_col_clustering() %>% 
  add_row_clustering() %>% 
  add_col_labels() %>%
  add_main_heatmap(as.matrix(peaks.CS.tbl.df), name="Blood Traits", colors=c("white",jdb_palette("Zissou")[1])) %>%
  add_col_clustering() %>%
  add_col_labels()







