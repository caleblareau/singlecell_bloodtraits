#' ---
#' title: "Descriptive analyses of fine mapped variants"
#' author: "Jacob Ulirsch"
#' date: "July 31, 2017"
#' ---

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
    perm[[paste0(i,j)]] <- permTest(A=CS.PP.gr, B=gr.list[[i]], ntimes=10000, alternative="auto", randomize.function=randomizeLocalRegions, evaluate.function=numOverlaps, force.parallel=TRUE, mc.set.seed=FALSE, mc.cores=48)
    OR <- rbind(OR,c(i,j,perm[[paste0(i,j)]]$numOverlaps$observed / length(CS.PP.gr), mean(perm[[paste0(i,j)]]$numOverlaps$permuted) / length(CS.PP.gr), perm[[paste0(i,j)]]$numOverlaps$zscore))
  print(c(i,j))
  }
}

local <- localZScore(CS.PP.gr, gr.list[[2]], pt=perm[[paste0(2,6)]], window=1500000,step=1000, mc.set.seed=FALSE, mc.cores=32)

randomizeLocalRegions <- function(A, ...) {
  A <- toGRanges(A)
  rand <- ceiling(runif(1,-1500000,1500000))
  B <- shift(A,rand)
  return(B)
}

OR.df <- as.data.frame(OR)
names(OR.df) <- c("i","j","obs","perm","sd")
OR.df$i <- as.factor(OR.df$i)
OR.df$j <- as.factor(OR.df$j)
OR.df$OR <- OR.df$obs / OR.df$perm
OR.df$z <- (OR.df$obs - OR.df$perm) / (OR.df$sd / sqrt(1000))
OR.df <- OR.df %>% dplyr::select(i,j,OR)

ggplot(OR.df,aes(y=OR,x=i,group=j)) + 
  geom_bar(stat="identity",aes(fill=j),position = position_dodge(0.9)) +
  geom_hline(yintercept=1, colour="grey3", linetype =3) +
  theme_bw() +
  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5]))


