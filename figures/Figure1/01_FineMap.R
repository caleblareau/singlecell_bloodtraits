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
counts.df <- read.table("../../data/bulk/ATAC/26August2017_EJCsamples_allReads_250bp.counts.txt",header=T)
names(counts.df) <- c("B","CD4","CD8","CLP","CMP","Ery","GMP-A","GMP-B","GMP-C","HSC","LMPP","mDC","Mega","MEP","mono","MPP","NK","pDC")
# Remove "weak" peaks that aren't in top n% for at least one cell type
n=0.8
counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- normalize.quantiles(as.matrix(counts.df))
keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))

#' Read in consensus bed file of peaks and add count metadata
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
peaks.gr <- import("../../data/bulk/ATAC/26August2017_EJCsamples_allReads_250bp.bed",format="bed")
values(peaks.gr) <- counts.df
peaks.gr <- peaks.gr[keep,]

#' Read in other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
coding.gr <- import("../../data/annotations/Coding_UCSC.bed",format="bed")
intron.gr <- import("../../data/annotations/Intron_UCSC.bed",format="bed")
promoter.gr <- import("../../data/annotations/Promoter_UCSC.fixed.bed",format="bed")
UTR.gr <- import("../../data/annotations/UTR_3_UCSC.bed",format="bed")

#' Read in heme narrow peaks
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
np_files <- list.files("../../data/bulk/ATAC/narrowpeaks/", full.names = TRUE)
np_files <- np_files[c(-7,-13,-21)]
heme.gr <- lapply(np_files, function(file){import(file, format = "BED", extraCols = extraCols_narrowPeak)})
names(heme.gr) <- c("B","CD4","CD8","CLP","CMP","Ery","GMP-A","GMP-B","GMP-C","HSC","LMPP","mDC","Mega","MEP","mono","MPP","NK","pDC")

#' Figure 1B
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)
CS.df.best <- CS.df %>%
  group_by(region,trait) %>%
  dplyr::summarize(PPmax=max(PP)) %>%
  as.data.frame()
CS.df.best$PPbin <- cut(CS.df.best$PPmax, c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))
CS.df.best.sum <- CS.df.best %>%
  group_by(trait,PPbin) %>%
  summarize(count=n())
ggplot(CS.df.best.sum,aes(y=count,x=trait,group=trait)) + 
  geom_bar(stat="identity",aes(fill=PPbin),position = position_stack(reverse = TRUE)) + 
  coord_flip() +
  theme_bw() +
  pretty_plot() +
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5]))

#' Figure 1C
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10
FM_exp <- readRDS("../../data/Finemap/allFINEMAPregions.rds")
ggplot(FM_exp,aes(x=expectedvalue,y=trait,fill=trait)) + 
  geom_joy(scale = 4) + 
  theme_joy() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_fill_cyclical(values = jdb_palette("GrandBudapest2")[c(3,4)])

#' Helper function for shifting
#+ echo=FALSE, message=FALSE, warning=FALSE
randomizeLocalRegions2 <- function(A, ...) {
  A <- toGRanges(A)
  rand <- ceiling(runif(1,-1500000,1500000))
  B <- GenomicRanges::shift(A,rand)
  return(B)
}

#' Figure 1E
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10
CS.df$seqnames <- paste0("chr",CS.df$seqnames)
CS.gr <- makeGRangesFromDataFrame(CS.df, keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)
OR <- NULL
perm <- NULL
trait <- unique(CS.gr@elementMetadata$trait)
gr.list <- list(peaks.gr,coding.gr,promoter.gr,intron.gr,UTR.gr)
bins = c(0.001,0.01,0.05,0.1,0.25,0.75,1)
for (i in seq(1,5,1)) {
  for (j in seq(1,6,1)) {
    CS.PP.gr <- subset(CS.gr,PP > bins[j] & PP <= bins[j+1]) 
    perm[[paste0(i,j)]] <- permTest(A=CS.PP.gr, B=gr.list[[i]], ntimes=10000, alternative="auto", randomize.function=randomizeLocalRegions2, evaluate.function=numOverlaps, force.parallel=FALSE, mc.set.seed=FALSE, mc.cores=4)
    OR <- rbind(OR,c(i,j,perm[[paste0(i,j)]]$numOverlaps$observed / length(CS.PP.gr), mean(perm[[paste0(i,j)]]$numOverlaps$permuted) / length(CS.PP.gr), perm[[paste0(i,j)]]$numOverlaps$zscore))
  print(c(i,j))
  }
}

#local <- localZScore(CS.PP.gr, gr.list[[2]], pt=perm[[paste0(2,6)]], window=1500000,step=1000, mc.set.seed=FALSE, mc.cores=32)

perm1 <- perm
OR1 <- OR

#' Figure 3A-C
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10
traits <- c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")
allregions <- lapply(traits, function(y) fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/LDtosentinel_alltraits/LDfiltered.0.8.",y,".txt")))
names(allregions) <- traits
allregions.df <- bind_rows(allregions)
allregions.df$trait <- rep(names(allregions), sapply(allregions, nrow))
allregions.df$cp <- paste0(allregions.df$trait,":","chr",allregions.df$chromosome,":",allregions.df$nonsentinel)
allregions.df$sp <- paste0(allregions.df$trait,":","chr",allregions.df$chromosome,":",allregions.df$sentinel)
CS.df$seqnames <- paste0("chr",CS.df$seqnames)
CS.df$cp <- paste0(CS.df$trait,":",CS.df$seqnames,":",CS.df$end)
LD.df <- CS.df %>% 
  dplyr::filter(cp %in% c(allregions.df$cp,allregions.df$sp))
LD.gr <- makeGRangesFromDataFrame(LD.df, keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)
OR <- NULL
perm <- NULL
trait <- unique(LD.gr@elementMetadata$trait)
gr.list <- heme.gr
for (i in seq(1,18,1)) {
  for (j in seq(8,8,1)) {
    LD.PP.gr <- LD.gr[LD.gr@elementMetadata$trait == trait[j],] 
    perm[[paste0(i,j)]] <- permTest(A=LD.PP.gr, B=gr.list[[i]], ntimes=10000, alternative="auto", randomize.function=randomizeLocalRegions2, evaluate.function=numOverlaps, force.parallel=FALSE, mc.set.seed=FALSE, mc.cores=4)
    OR <- rbind(OR,c(i,j,names(gr.list)[i],trait[j],perm[[paste0(i,j)]]$numOverlaps$observed / length(LD.PP.gr), mean(perm[[paste0(i,j)]]$numOverlaps$permuted) / length(LD.PP.gr), perm[[paste0(i,j)]]$numOverlaps$zscore))
    print(c(i,j))
  }
}

OR.df <- as.data.frame(OR,stringsAsFactors=F)
names(OR.df) <- c("i","j","cell","trait","obs","perm","z")
OR.df$i <- as.factor(OR.df$i)
OR.df$j <- as.factor(OR.df$j)
OR.df$OR <- as.numeric(OR.df$obs) / as.numeric(OR.df$perm)

ggplot(OR.df,aes(y=OR,x=i,group=j)) + 
  geom_bar(stat="identity",aes(fill=j),position = position_dodge(0.9)) +
  geom_hline(yintercept=1, colour="grey3", linetype =3) +
  theme_bw() +
  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5]))


