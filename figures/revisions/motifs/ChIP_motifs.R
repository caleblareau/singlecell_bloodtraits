# Load libraries
library(data.table)
library(dplyr)
library(stringr)
library(motifmatchr)
library(motifbreakR)
library(universalmotif)
library(TFBSTools)
library(GenomicRanges)
library(ggplot2)
library(BuenColors)
library(preprocessCore)
"%ni%" <- Negate("%in%")

#' Read in and munge ChIP-Atlas
# Read in meta data
chip_atlas.meta <- fread("../../data/chipatlas/experimentList.clean.tab", sep = "\t", header = F)
names(chip_atlas.meta) <- c("id", "genome", "class1", "antigen", "class2", "celltype")
chip_atlas.meta <- chip_atlas.meta %>%
  dplyr::filter(genome == "hg19", class1 == "TFs and others", class2 == "Blood", antigen %ni% c("-", "5-hmC", "5-mC", "Cyclobutane pyrimidine dimers", "Epitope tags", "HIV Tat", "MethylCap", "pFM2")) #%>%

# Read in bed file
chip_atlas.bed <- fread("zcat < ../../data/chipatlas/Oth.Bld.05.AllAg.AllCell.clean.bed.gz", sep = "\t", header = F)
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

#' Read in and munge fine mapped variants
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)
CS.df$seqnames <- paste0("chr",CS.df$seqnames)
CS.df$var <- paste0("chr",CS.df$var)
CS.df$var <- gsub("_", ":", CS.df$var)

# UKBB exclusion list
exdf <- read.table("exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
CS.gr <- CS.gr[CS.gr$var %ni% exdf,]

# Only PP > 10% variants
CS.df.sub <- CS.df %>%
  dplyr::filter(PP > 0.1)
highPP <- CS.df.sub %>% 
  .$var
CS.df.sub.gr <- GRanges(CS.df.sub)

# Read in and munge motifbreakR results
ALL.motif <- readRDS("../../data/motifbreakR/alltraits.mbreaker.withPPs.rds")
ALL.motif.sub <- ALL.motif %>%
  dplyr::filter(PP > 0.1, SNP %ni% gsub("_", ":", exdf)) %>%
  dplyr::select(seqnames, start, end, SNP, REF, ALT, geneSymbol, providerName, seqMatch, effect, PP, MAF, trait)
ALL.motif.sub.gr <- GRanges(ALL.motif.sub)

# Read in and munge motif/TF matching/similarity data
motif.cor <- readRDS("motifs/hocomoco_similarity.rds")
motif.cor <- rbind(motif.cor, data.frame(motif1 = levels(motif.cor$motif1), motif2 = levels(motif.cor$motif1), cor = 1))
motif.cor <- motif.cor %>%
  dplyr::filter(cor > 0.7)

# Function to make list of motif families
makeFam <- function(df, motif) {
    df.2 <- df %>%
      dplyr::filter(motif1 == motif | motif2 == motif)
    match <- c(as.character(df.2$motif1), as.character(df.2$motif2))
    match.gS <- hocomoco@elementMetadata[match,]$geneSymbol
    return(match.gS)
}

# Apply function
motif.cor.ls <- lapply(levels(motif.cor$motif1), function(x) {unique(makeFam(motif.cor, x))})
names(motif.cor.ls) <- gsub("Hsapiens-HOCOMOCO-", "", levels(motif.cor$motif1))

# Merge motifbreakR and ChIP-Atlas for high PP variants
idx <- findOverlaps(ALL.motif.sub.gr, chip_atlas.bed.sub.gr)
ALL.M_C.df <- data.frame(cbind(as.data.frame(ALL.motif.sub.gr[idx@from]), as.data.frame(chip_atlas.bed.sub.gr[idx@to])))
ALL.M_C.df.match <- ALL.M_C.df %>%
  rowwise() %>%
  mutate(keep = ifelse(antigen %in% motif.cor.ls[[providerName]], "yes", "no")) %>%
  dplyr::filter(keep == "yes") %>%
  dplyr::select(SNP, trait, PP, antigen, celltype, geneSymbol, providerName, seqMatch, effect, seqnames, start, end) %>% 
  arrange(SNP) %>%
  as.data.frame()
ALL.M_C.df.match$seqMatch <- gsub("[ ]*", "", ALL.M_C.df.match$seqMatch)

# Write table for supplement
write.table(ALL.M_C.df.match, "motifs/ALL.M_C.df.match.txt", sep = "\t", quote = F, row.names = F)

# JCU definition of lineage
lineage <- c("BASO_COUNT" = "GRAN", "EO_COUNT" = "GRAN", "HCT" = "RBC", "HGB" = "RBC", "LYMPH_COUNT" = "LYMPH", "MCH" = "RBC", "MCHC" = "RBC", "MCV" = "RBC", "MEAN_RETIC_VOL" = "RBC", "MONO_COUNT" = "MONO", "MPV" = "PLT", "NEUTRO_COUNT" = "GRAN", "PLT_COUNT" = "PLT", "RBC_COUNT" = "RBC", "RETIC_COUNT" = "RBC", "WBC_COUNT" = "LYMPH") 

# Unique total variants with at least one mechanisms
temp <- ALL.M_C.df.match %>% 
  group_by(SNP) %>%
  distinct(SNP, .keep_all = TRUE) 
temp %>% 
  ungroup() %>%
  summarize(count = n())

# Check overlaps with AC
idx <- findOverlaps(GRanges(temp), peaks)
length(table(idx@from))

# Unique variants per trait with at least one mechanism
ALL.M_C.df.match %>% 
  merge(., lineage, by.x = "trait", by.y = "row.names") %>%
  mutate(lineage = y) %>%
  group_by(SNP, lineage) %>%
  dplyr::select(SNP, lineage) %>%
  distinct() %>% 
  ungroup() %>%
  group_by(lineage) %>%
  summarize(count = n())

# Collapse to unique variant / antigen / lineage (trait) combinations
ALL.M_C.df.match.TFcount <- ALL.M_C.df.match %>% 
  merge(., lineage, by.x = "trait", by.y = "row.names") %>%
  mutate(lineage = y) %>%
  group_by(SNP, antigen, lineage) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(antigen, lineage) %>%
  summarize(count = sum(count > 0)) %>%
  arrange(desc(count)) %>%
  as.data.frame()

# Collapse to unique variant / antigen / lineage (trait) combinations
ALL.M_C.df.match.TFcountSum <- ALL.M_C.df.match %>% 
  merge(., lineage, by.x = "trait", by.y = "row.names") %>%
  mutate(lineage = y) %>%
  group_by(SNP, antigen) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(antigen) %>%
  summarize(count = sum(count > 0)) %>%
  arrange(desc(count)) %>%
  as.data.frame()
ALL.M_C.df.match.TFcountSum.vec <- as.vector(ALL.M_C.df.match.TFcountSum$count)
names(ALL.M_C.df.match.TFcountSum.vec) <- ALL.M_C.df.match.TFcountSum$antigen

# Long to wide for heatmap
ALL.M_C.df.match.TFcount.wide <- dcast(ALL.M_C.df.match.TFcount, antigen ~ lineage)
ALL.M_C.df.match.TFcount.wide[is.na(ALL.M_C.df.match.TFcount.wide)] <- 0
row.names(ALL.M_C.df.match.TFcount.wide) <- ALL.M_C.df.match.TFcount.wide$antigen
ALL.M_C.df.match.TFcount.wide <- ALL.M_C.df.match.TFcount.wide[,-1]
ALL.M_C.df.match.TFcount.wide <- as.matrix(ALL.M_C.df.match.TFcount.wide)
ALL.M_C.df.match.TFcount.wide <- ALL.M_C.df.match.TFcount.wide[names(sort(ALL.M_C.df.match.TFcountSum.vec)),]

# Plot ChIP-seq + motifBreakR overlap
hm <- Heatmap(ALL.M_C.df.match.TFcount.wide, col=as.character(jdb_palette("brewer_spectra", type="continuous")),
              cluster_rows = FALSE, cluster_columns = TRUE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              show_heatmap_legend = FALSE)
ha1 <- rowAnnotation(foo2 = row_anno_barplot(rev(ALL.M_C.df.match.TFcountSum.vec), axis = TRUE, width = unit(2, "cm")))
hm + ha1

# Read in heme ATAC-seq peaks
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

# Import counts, filter, and normalize
counts.df <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
n = 0.9
counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- normalize.quantiles(as.matrix(counts.df))
keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))

# Subset to only good peaks
peaks <- peaks[keep,]

# Calculate motifmatches for PP > 0.10 variants
data(hocomoco)
hocomoco.mm <- convert_motifs(hocomoco, "TFBSTools-PWMatrixList")
hocomoco.mm.name <- unlist(lapply(hocomoco.mm, function(x) {x@name}))
hocomoco.mm <- do.call(PFMatrixList, hocomoco.mm)
names(hocomoco.mm) <- hocomoco.mm.name
motif_ix <- motifmatchr::matchMotifs(hocomoco.mm, peaks, genome = "hg19", out = "positions") 
motif_ix.df <- as.data.frame(motif_ix)
motif_ix.df$start <- motif_ix.df$start - 20 
motif_ix.df$end <- motif_ix.df$end + 20
motif_ix.gr <- GRanges(motif_ix.df)
idx <- findOverlaps(CS.df.sub.gr, motif_ix.gr)
CS.df.sub.motif.df <- data.frame(cbind(as.data.frame(CS.df.sub.gr[idx@from]), as.data.frame(motif_ix.gr[idx@to])))
CS.df.sub.motif.gr <- GRanges(CS.df.sub.motif.df)

# Merge motifMatchR and ChIP-Atlas for high PP variants
idx <- findOverlaps(CS.df.sub.motif.gr, chip_atlas.bed.sub.gr)
ALL.MM_C.df <- data.frame(cbind(as.data.frame(CS.df.sub.motif.gr[idx@from]), as.data.frame(chip_atlas.bed.sub.gr[idx@to])))
ALL.MM_C.df.match <- ALL.MM_C.df %>%
  rowwise() %>%
  mutate(keep = ifelse(antigen %in% motif.cor.ls[[group_name]], "yes", "no")) %>%
  dplyr::filter(keep == "yes") %>%
  dplyr::select(var, trait, PP, antigen, celltype, group_name, seqnames, start, end) %>% 
  arrange(var) %>%
  as.data.frame()

# Write table for supplement
write.table(ALL.MM_C.df.match, "motifs/ALL.M_C.df.match.txt", sep = "\t", quote = F, row.names = F)

# Unique total variants with at least one mechanisms
temp <- ALL.MM_C.df.match %>% 
  group_by(var) %>%
  distinct(var, .keep_all = TRUE)
temp %>% 
  ungroup() %>%
  summarize(count = n())

# Check overlaps with AC
idx <- findOverlaps(GRanges(temp), peaks.all)
length(table(idx@from))

# Unique variants per trait with at least one mechanism
ALL.MM_C.df.match %>% 
  merge(., lineage, by.x = "trait", by.y = "row.names") %>%
  mutate(lineage = y) %>%
  group_by(var, lineage) %>%
  dplyr::select(var, lineage) %>%
  distinct() %>% 
  ungroup() %>%
  group_by(lineage) %>%
  summarize(count = n())

# Collapse to unique variant / antigen / lineage (trait) combinations
ALL.MM_C.df.match.TFcount <- ALL.MM_C.df.match %>% 
  merge(., lineage, by.x = "trait", by.y = "row.names") %>%
  mutate(lineage = y) %>%
  group_by(var, antigen, lineage) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(antigen, lineage) %>%
  summarize(count = sum(count > 0)) %>%
  arrange(desc(count)) %>%
  as.data.frame()

# Collapse to unique variant / antigen / lineage (trait) combinations
ALL.MM_C.df.match.TFcountSum <- ALL.MM_C.df.match %>% 
  merge(., lineage, by.x = "trait", by.y = "row.names") %>%
  mutate(lineage = y) %>%
  group_by(var, antigen) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(antigen) %>%
  summarize(count = sum(count > 0)) %>%
  arrange(desc(count)) %>%
  as.data.frame()
ALL.MM_C.df.match.TFcountSum.vec <- as.vector(ALL.MM_C.df.match.TFcountSum$count)
names(ALL.MM_C.df.match.TFcountSum.vec) <- ALL.MM_C.df.match.TFcountSum$antigen

# Long to wide for heatmap
ALL.MM_C.df.match.TFcount.wide <- dcast(ALL.MM_C.df.match.TFcount, antigen ~ lineage)
ALL.MM_C.df.match.TFcount.wide[is.na(ALL.MM_C.df.match.TFcount.wide)] <- 0
row.names(ALL.MM_C.df.match.TFcount.wide) <- ALL.MM_C.df.match.TFcount.wide$antigen
ALL.MM_C.df.match.TFcount.wide <- ALL.MM_C.df.match.TFcount.wide[,-1]
ALL.MM_C.df.match.TFcount.wide <- as.matrix(ALL.MM_C.df.match.TFcount.wide)
ALL.MM_C.df.match.TFcount.wide <- ALL.MM_C.df.match.TFcount.wide[names(sort(ALL.MM_C.df.match.TFcountSum.vec)),]

# Plot ChIP-seq + motifMatchR overlap
hm <- Heatmap(ALL.MM_C.df.match.TFcount.wide, col=as.character(jdb_palette("brewer_spectra", type="continuous")),
              cluster_rows = FALSE, cluster_columns = TRUE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              show_heatmap_legend = TRUE)
ha1 <- rowAnnotation(foo2 = row_anno_barplot(rev(ALL.MM_C.df.match.TFcountSum.vec), axis = TRUE, width = unit(2, "cm")))
hm + ha1

# Examples
results <- readRDS("../../data/motifbreakR/MEAN_RETIC_VOL_Motifbreakr_output_PP001.rds")
plotMB(results = results, rsid = "chr1:203275407:C:G", effect = "strong")
plotMB(results = results, rsid = "chr16:87886545:C:T", effect = "strong")
results <- readRDS("../../data/motifbreakR/MPV_Motifbreakr_output_PP001.rds")
plotMB(results = results, rsid = "chr10:97057370:C:T", effect = "strong")
results <- readRDS("../../data/motifbreakR/PLT_COUNT_Motifbreakr_output_PP001.rds")
plotMB(results = results[results$geneSymbol=="RUNX1",], rsid = "chr1:31241886:G:T", effect = "strong")
results <- readRDS("../..//LYMPH_COUNT_Motifbreakr_output_PP001.rds")
plotMB(results = results[results$geneSymbol=="RUNX1",], rsid = "chr2:143886819:G:A", effect = "strong")


chr1:203275157-203275657

  