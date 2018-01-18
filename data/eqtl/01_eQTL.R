#' ---
#' title: "Match eQTLs in SMR with FM variants"
#' author: "Jacob Ulirsch"
#' date: "January 18, 2018"
#' ---

library(readr)
library(dplyr)

# Read in FM variants
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df$alleles <- gsub("^.*_", "", CS.df$var)
CS.df$var2 <- paste0(CS.df$seqnames, ":", CS.df$end)

# Read in esi files from SMR
westra <- read_delim("../../data/eqtl/westra_eqtl_hg19.esi.gz", delim = "\t", col_names = F)
mcrae <- read_delim("../../data/eqtl/mcrae_meQTL_hg19.esi.gz", delim = "\t", col_names = F)
jones <- read_delim("../../data/eqtl/jones_eqtl_hg19.esi.gz", delim = "\t", col_names = F)

# Merge with FM variants and get rsIDs for SMR
westra$var1 <- paste0(westra$X1, ":", westra$X4, "_", westra$X6, "_", westra$X5)
westra$var2 <- paste0(westra$X1, ":", westra$X4, "_", westra$X5, "_", westra$X6)
westra.fm <- westra %>%
  dplyr::filter(var1 %in% CS.df$var | var2 %in% CS.df$var)
mcrae$var1 <- paste0(mcrae$X1, ":", mcrae$X4, "_", mcrae$X6, "_", mcrae$X5)
mcrae$var2 <- paste0(mcrae$X1, ":", mcrae$X4, "_", mcrae$X5, "_", mcrae$X6)
mcrae.fm <- mcrae %>%
  dplyr::filter(var1 %in% CS.df$var | var2 %in% CS.df$var)
jones$var1 <- paste0(jones$X1, ":", jones$X4, "_", jones$X6, "_", jones$X5)
jones$var2 <- paste0(jones$X1, ":", jones$X4, "_", jones$X5, "_", jones$X6)
jones.fm <- jones %>%
  dplyr::filter(var1 %in% CS.df$var | var2 %in% CS.df$var)

# Write out rsIDs to query
write.table(westra.fm$X2, "../../data/eqtl/westra_FM_rsIDs.txt", quote = F, row.names = F, col.names = F)
write.table(mcrae.fm$X2, "../../data/eqtl/mcrae_FM_rsIDs.txt", quote = F, row.names = F, col.names = F)
write.table(jones.fm$X2, "../../data/eqtl/jones_FM_rsIDs.txt", quote = F, row.names = F, col.names = F)

# After running SMR to query results, read in files
westra.eqtl <- read_delim("../../data/eqtl/westra_FM_eqtl.txt.gz", delim = "\t", col_names = F)
mcrae.meqtl <- read_delim("../../data/eqtl/mcrae_FM_eqtl.txt.gz", delim = "\t", col_names = F)
jones.eqtl <- read_delim("../../data/eqtl/jones_FM_eqtl.txt.gz", delim = "\t", col_names = F)


