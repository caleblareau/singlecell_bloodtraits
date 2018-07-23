#' ---
#' title: "gkm-SVM on hematopoiesis"
#' author: "Jacob Ulirsch"
#' date: "November 11, 2017"
#' ---

library(data.table)
library(dplyr)
library(rtracklayer)
library(preprocessCore)
library(ggjoy)
library(BuenColors)
library(regioneR)
library(stringr)
library(gkmSVM)
library(readr)
library("BSgenome")
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(tidyr)
library(ROCR)

#' Read in count matrix 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
counts.df <- read.table("../../data/bulk/ATAC/26August2017_EJCsamples_allReads_250bp.counts.txt",header=T)
names(counts.df) <- c("B","CD4","CD8","CLP","CMP","Ery","GMP-A","GMP-B","GMP-C","HSC","LMPP","mDC","Mega","MEP","mono","MPP","NK","pDC")
# Remove "weak" peaks that aren't in top n% for at least one cell type
n=0.9
counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- normalize.quantiles(as.matrix(counts.df))
keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))

#' Read in consensus bed file of peaks and add count metadata
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
peaks.gr <- import("../../data/bulk/ATAC/26August2017_EJCsamples_allReads_250bp.bed",format="bed")
values(peaks.gr) <- counts.df
peaks.gr <- peaks.gr[keep,]

#' Read in fine-mapped variants
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)
CS.df$seqnames <- paste0("chr",CS.df$seqnames)
CS.gr <- makeGRangesFromDataFrame(CS.df, keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)

# Function to find matching null sequences in genome and write out bed and fasta files
makeNull <- function(gr) {
  # Output files
  bedP <- paste0(paste0("bed/pos.bed"))
  bedN <- paste0(paste0("bed/neg.bed"))
  fastaP <- paste0(paste0("fasta/pos.fa"))
  fastaN <- paste0(paste0("fasta/neg.fa"))
  
  # Export bed file for specific region
  export(gr, bedP, format="bed")
  
  # Generate matched mull sequences
  genNullSeqs(bedP, xfold=1, nMaxTrials=25, outputPosFastaFN=fastaP, outputBedFN=bedN, outputNegFastaFN=fastaN)
}

# Function to make fasta files from an input bed file
makeFASTA <- function(x,gr,co) {
  # Output files
  bedP <- paste0(paste0("bed/", x, "_pos.bed"))
  fastaP <- paste0(paste0("fasta/", x, "_pos.fa"))
  null <- "null.txt"
  
  # Export bed file for specific region
  gr.sub <- gr[gr@elementMetadata[,x] > co]
  export(gr.sub, bedP, format="bed")
  
  # Make fasta
  genNullSeqs(bedP, xfold=1, nMaxTrials=1, batchsize=10, outputPosFastaFN=fastaP, outputBedFN=null, outputNegFastaFN=null)
}

# Where we actually do the things:
# Make negative control fasta file
makeNull(peaks.gr)
# Make positive controls for each cell type
cutoff <- mean(apply(counts.df,2,function(x) {quantile(x,n)}))
lapply(names(counts.df),function(x) {makeFASTA(x,peaks.gr,cutoff)})

# Plot precision recall curve for CV models and get AUPRC
PRC <- function(vals,labs) {
  library(ROCR)
  pred <- prediction(vals, labs)
  prc <- performance(pred, "prec", "rec")
  #plot(prc)
  return(prc)
}

AUPRC <- function(prc) {
  auprc <- sum(diff(prc@x.values[[1]][-1]) * (head(prc@y.values[[1]][-1],-1)+tail(prc@y.values[[1]][-1],-1)))/2
}

# Function to extract sequences for given variants
extract_oligos <- function(entry, genome, totlength=200, info) {
  chr.num <- trimws(as.character(entry['CHR']))
  bp <- as.integer(as.character(entry['BP']))
  ref <- as.character(entry['REF'])
  mut <- as.character(entry['ALT'])
  info <- as.character(entry['trait'])
  oldref <- as.character(entry['REF'])
  
  halfSW = (totlength-1)/2
  
  # Check if ref and alt are flipped; if so, flip them back
  testref <- getSeq(genome,paste0(chr.num),start=bp,end=bp+nchar(ref)-1,width=NA, as.character=TRUE, strand="+")
  if (testref != ref){
    #print("ref doesn't match") 
    if (testref == mut){
      #print("mut matches") 
      temp <- ref
      ref <- mut
      mut <- temp
    } else {
      print(paste(chr.num,bp,ref,mut))
      return(NULL)
    }
  }
  
  refoligos = getSeq(genome,paste0(chr.num),start=bp-halfSW,end=bp+nchar(ref)+(totlength-halfSW-1)-1,width=NA, as.character=TRUE, strand="+")
  mutoligos = paste0(substr(refoligos,1,halfSW+1),mut,substr(refoligos,halfSW+1+nchar(ref)+1,nchar(refoligos)))
  
  if (nchar(refoligos)>totlength){
    refoligos = trim_oligos(refoligos,totlength)
  }
  if (nchar(mutoligos)>totlength){
    mutoligos = trim_oligos(mutoligos,totlength)
  }
  
  output <- data.frame(variant=character(0),name.and.window=character(0),oligonucleotide=character(0),info=character(0),oldRef=character(0),stringsAsFactors=FALSE)
  output[1,]=c(paste0(paste(chr.num,bp,sep=":"),"_",paste(ref,mut,sep="_")),paste("Ref","Bot_1/2_Top_1/2",sep="_"),refoligos[[1]],info,oldref)
  output[2,]=c(paste0(paste(chr.num,bp,sep=":"),"_",paste(ref,mut,sep="_")),paste("Mut","Bot_1/2_Top_1/2",sep="_"),mutoligos[[1]],info,oldref)
  
  return(output)
}

# Trim the start and end of oligos by equal amounts to reach the desired total length
trim_oligos <- function(oligos, totlength=200){
  # trim one half SW
  output= substr(oligos[[1]],(nchar(oligos[[1]])-totlength)/2+1,nchar(oligos[[1]])-(nchar(oligos[[1]])-totlength)/2)
  return(list(output))
}

# Run 
genome <- BSgenome.Hsapiens.UCSC.hg19
totlength <- 50 
allVars <- CS.df %>%
  separate(var, c("CHR", "BP", "REF", "ALT"), ":|_")
allVars$BP <- as.numeric(allVars$BP)
allVars <- allVars[c('CHR','BP','REF','ALT',"trait")]
allVars$CHR <- paste0("chr", allVars$CHR)
last <- dim(allVars)[1]
constructs <- bind_rows(mclapply(X = 1:last, function(y) extract_oligos(allVars[y,], genome, totlength), mc.cores = 20))
saveRDS(constructs, "FM_fasta.rds")
ref <- constructs[grep("Ref",constructs$name.and.window),]
mut <- constructs[grep("Mut",constructs$name.and.window),]
write.fasta(as.list(ref$oligonucleotide),ref$variant, "FM_ref.fa")
write.fasta(as.list(mut$oligonucleotide),mut$variant, "FM_mut.fa")

# Identify deltaSVM output
files <- Sys.glob("../../data/gkm/*deltaSVM.txt")
cellnames <- gsub("../../data/gkm/", "", files)
cellnames <- gsub("_deltaSVM.txt", "", cellnames)

# Read in and create sparse matrix
gkm <- sapply(files, function(file) {
  data.frame(fread(file, select = 2, skip=0))[,1]
}) %>% data.matrix() %>% Matrix()

# Assign proper names
colnames(gkm) <- cellnames
row.names(gkm) <- data.frame(fread(files[1], select = 1, skip=0))[,1]
gkm <- gkm[!duplicated(row.names(gkm)),]
gkm.df <- data.frame(as.matrix(gkm))
row.names(gkm.df) <- gsub("chr", "", row.names(gkm.df))
gkm.df$var <- gsub("chr", "", row.names(gkm.df))

saveRDS(gkm, "gkm.rds")
cor(data.matrix(gkm[1:5,1:5]))

# Identify deltaSVM output
files <- Sys.glob("../../data/gkm/*.cvpred.txt")
cellnames <- gsub("../../data/gkm/", "", files)
cellnames <- gsub("_pos_cv.cvpred.txt", "", cellnames)

# Read in and create sparse matrix
gkm.cv <- lapply(files, function(file) {
  data.frame(fread(file, select = 2:3, skip=0))
})

# Plot PRC for each cell type
gkm.cv.plot <- lapply(gkm.cv, FUN = function(x) PRC(x$V2, x$V3))
gkm.cv.auprc <- unlist(lapply(gkm.cv.plot, FUN = function(x) AUPRC(x)))
mean(gkm.cv.auprc)

# Compare for variants
peaks.CS.df.gkm <- merge(peaks.CS.df, gkm.df, by = "var")
peaks.CS.df.gkm$PPbin <- cut(peaks.CS.df.gkm$PP, c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))
blah <- peaks.CS.df.gkm %>% dplyr::select(ends_with(".y"))
peaks.CS.df.gkm$max <- apply(blah, 1, function(x) {max(abs(x))})

peaks.CS.df.gkm.plot <- peaks.CS.df.gkm %>% 
  dplyr::filter(trait %in% "RBC_COUNT")
ggplot(peaks.CS.df.gkm.plot, aes(PPbin, max)) + geom_boxplot()






