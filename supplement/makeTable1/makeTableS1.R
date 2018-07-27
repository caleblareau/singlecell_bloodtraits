library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(diffloop)
"%ni%" <- Negate("%in%")

bdf <- lapply(list.files("../../data/UKBB_BC_PP001/betas_added/", full.names = TRUE), fread) %>% rbindlist %>% as.data.frame
splitdf <- str_split_fixed(bdf$V4, "-", 3)
splitdf2 <- str_split_fixed(splitdf[,2], "_", 3)

bbdf <- cbind(bdf, splitdf, splitdf2)
odf <- bbdf[,c(1,2,3,9,13,14,5,6,7,8)]
odf <- odf[!duplicated(odf),]

colnames(odf) <- c("chr", "start", "end", "trait", "ref", "alt", "PP", "Beta", "SE", "Z")

odf$var <- paste(odf$chr,odf$start,sep=":") %>% gsub("chr","",.)
odf$var <- paste(odf$var,odf$ref,odf$alt, sep="_")
odf <- odf %>% dplyr::select(var,everything())

# Exclude variants from exclude list
exdf <- read.table("../../figures/revisions/exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
odf <- odf[odf$var %ni% exdf,]

# Check for any missing variants
CS.df <- readRDS("../../data/Finemap/UKBB_BC_v3_VEPannotations.rds") %>% as.data.frame()
var1 <- CS.df$var %>% unique
var2 <- odf$var %>% unique 
setdiff(var1,var2)

g <- makeGRangesFromDataFrame(odf)

# Import NHGRI
gwas <- data.frame(fread(paste0("zcat < ", "gwas_catalog_v1.0-associations_e91_r2017-12-11.tsv.gz"), sep = "\t", header = TRUE))
gwas$CHR_POS <- as.numeric(gwas$CHR_POS )
gwas_g <- makeGRangesFromDataFrame(gwas[complete.cases(gwas),], seqnames.field = "CHR_ID", start.field = "CHR_POS", end.field = "CHR_POS")
ov1 <- findOverlaps(g, addchr(gwas_g))
odf$In_Gwas_catalog <- 1:length(g) %in% queryHits(ov1)

# Import heme ATAC peaks
heme <- data.frame(fread(paste0( "../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed"), sep = "\t", header = FALSE))
heme_g <- makeGRangesFromDataFrame(heme[complete.cases(heme),], seqnames.field = "V1", start.field = "V2", end.field = "V3")
ov2 <- findOverlaps(g, addchr(gwas_g))
odf$In_Heme_peak <- 1:length(g) %in% queryHits(ov2)

#write.table(odf, file = "SupplementalTable1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Redo rsIDs
odf$tomatch <- paste(odf$chr,odf$start,sep=":")
odf$rsID <- rsids[odf$tomatch]
odf <- odf %>% dplyr::select(-tomatch,-end)

# Load up previously made table
s1 <- readxl::read_xlsx("SupplementalTable1.xlsx")
s1$var <- paste(s1$chr,s1$start,sep=":") %>% gsub("chr","",.)
s1$var <- paste(s1$var,s1$ref,s1$alt, sep="_")
s1 <- s1[s1$var %ni% exdf,] %>% dplyr::select(-end)

# For variants with multiple Phylop, take mean
phylo_means <- s1 %>% group_by(var) %>%
  summarise(Phylop_100way = mean(Phylop_100way))
phylos <- phylo_means$Phylop_100way; names(phylos) <- phylo_means$var

# Redo PhyloP
odf$Phylop_100way <- phylos[odf$var]

write.table(odf, file = "SupplementalTable1_v2.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# List out incorrect rsID mappings
s1[duplicated(s1[,c("var","trait")]),c("var","rsID")] %>% nrow
