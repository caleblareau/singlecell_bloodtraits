library(data.table)
library(dplyr)
library(rtracklayer)
library(preprocessCore)
library(ggjoy)
library(BuenColors)
library(regioneR)
library(stringr)
library(qvalue)
library(annotables)
"%ni%" <- Negate("%in%")
set.seed(1026)

# Read input files ---------------------------------------
# Read in FINEMAP dataset
CS.gr <- readRDS("../../data/Finemap/UKBB_BC_v3_VEPannotations.rds")
CS.df <- as.data.frame(CS.gr)
trait <- unique(CS.gr@elementMetadata$trait)

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

# ChIP Atlas
chip_atlas.meta <- fread("../../data/chipatlas/experimentList.clean.tab", sep = "\t", header = F)
names(chip_atlas.meta) <- c("id", "genome", "class1", "antigen", "class2", "celltype")
chip_atlas.meta <- chip_atlas.meta %>%
  dplyr::filter(genome == "hg19", class1 == "TFs and others", class2 == "Blood", 
                antigen %ni% c("-", "5-hmC", "5-mC", "Cyclobutane pyrimidine dimers", "Epitope tags", "HIV Tat", "MethylCap", "pFM2")) 
# Read in bed file
chip_atlas.bed <- fread("zcat < /Volumes/broad_sankaranlab/ChIPAtlas/Oth.Bld.05.AllAg.AllCell.clean.bed.gz", sep = "\t", header = F)
names(chip_atlas.bed) <- c("seqnames", "start", "end", "id")
chip_atlas.bed <- chip_atlas.bed %>%
  dplyr::filter(id %in% chip_atlas.meta$id)
# Subset to 22 autosomes
setkey(setDT(chip_atlas.bed), id)
setkey(setDT(chip_atlas.meta), id)
chip_atlas.bed <- merge(chip_atlas.bed, chip_atlas.meta)
chip_atlas.bed.sub <- chip_atlas.bed %>%
  dplyr::filter(seqnames %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))
# Convert to GRanges since these are so much faster for overlaps
chip_atlas.bed.sub.gr <- GRanges(chip_atlas.bed.sub)

# MotifbreakR
ALL.motif <- readRDS("../../data/motifbreakR/alltraits.mbreaker.withPPs.rds")
exdf <- read.table("../revisions/exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
ALL.motif.sub <- ALL.motif %>%
  dplyr::filter(PP > 0.1, SNP %ni% gsub("_", ":", exdf)) %>%
  dplyr::select(seqnames, start, end, SNP, REF, ALT, geneSymbol, providerName, seqMatch, effect, PP, MAF, trait)
ALL.motif.sub.gr <- GRanges(ALL.motif.sub)

# PCHiC 
pchic.df <- readRDS("../../data/pchic/pchic.rds")
pchic.gr <- GRanges(pchic.df)
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")
pchic.gr <- pchic.gr[pchic.gr$Gene %in% grch38.pc$symbol,] 

# ATAC-RNA
pg.df <- fread(paste0("zcat < ", "../../data/bulk/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","j.start","j.end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>% dplyr::filter(qvalue < 0.001)
pg.gr <- makeGRangesFromDataFrame(pg.df, seqnames = "chrom", start.field = "j.start", end.field = "j.end")

# Perform permutations ----------------------------------------------------
#' Helper function for shifting
randomizeLocalRegions <- function(A, ...) {
  A <- toGRanges(A)
  rand <- ceiling(runif(1,-1500000,1500000))
  B <- GenomicRanges::shift(A,rand)
  return(B)
}

# PP filter for FINEMAP variants
CS.PP.gr <- CS.gr[CS.gr$PP > 0.10,]
# Compile all annotations into a list
gr.list <- list(peaks.gr,chip_atlas.bed.sub.gr,ALL.motif.sub.gr,pchic.gr,pg.gr)
names(gr.list) <- c("hemeATAC","ChIP","motifs","PCHiC","ATAC-RNA")

# Run permutations on each annotation
enrichments <- NULL
perm <- NULL
for (i in 1:length(gr.list)) {
  print(paste0("Starting annotation #",names(gr.list)[i]))
  perm[[i]] <- permTest(A=CS.PP.gr, B=gr.list[[i]], 
                        ntimes=1000, 
                        alternative="auto", 
                        randomize.function=randomizeLocalRegions, 
                        evaluate.function=numOverlaps, 
                        force.parallel=FALSE, 
                        mc.set.seed=T, 
                        mc.cores=2)
  enrichments <- rbind(enrichments,c(names(gr.list)[i],
                   perm[[i]]$numOverlaps$zscore,
                   pnorm(perm[[i]]$numOverlaps$zscore, lower.tail = FALSE)))
  print(paste0("Finished annotation #",names(gr.list)[i]))
}

enrich.df <- as.data.frame(enrichments,stringsAsFactors=F)
names(enrich.df) <- c("annot","zscore","pvalue")
enrich.df$annot <- as.factor(enrich.df$annot)
enrich.df$qvalue <- qvalue::qvalue(enrich.df$pvalue)$qvalues

# Plot enrichments
ggplot(enrich.df,aes(y=-log10(pvalue),x=annot)) + 
  geom_bar(stat="identity",fill=jdb_palette("Zissou")[1],position = position_dodge(0.9)) +
  geom_hline(yintercept= -log10(0.05 / length(gr.list)), colour="grey3", linetype =3) +
  theme_bw() +
  pretty_plot() + L_border()

