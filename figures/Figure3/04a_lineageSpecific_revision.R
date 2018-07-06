library(BuenColors)
library(ggplot2)
library(dplyr)
library(data.table)

set.seed(14651)


# Define lineage specific attributes
Ery <- c("HSC", "MPP", "CMP", "MEP", "Ery")
Meg <- c("HSC", "MPP", "CMP", "MEP", "Mega")
Mye <- c("HSC", "MPP", "CMP", "LMPP", "GMP-A", "GMP-B", "GMP-C", "Mono", "mDC")
Lymph <- c("HSC", "MPP", "LMPP", "CLP", "NK", "pDC", "CD4", "CD8", "B")

make2df <- function(trait, celltypes){
  return(data.frame(Celltype = celltypes, Trait = rep(trait, length(celltypes))))
}

lineageSpecificDF <- rbind(
  make2df("BASO_COUNT", Mye),
  make2df("EO_COUNT", Mye),
  make2df("NEUTRO_COUNT", Mye),
  make2df("MONO_COUNT", Mye),
  make2df("WBC_COUNT", unique(c(Mye, Lymph))),
  make2df("LYMPH_COUNT", Lymph),
  make2df("PLT_COUNT", Meg),
  make2df("MPV", Meg),
  
  make2df("MEAN_RETIC_VOL", Ery),
  make2df("HCT", Ery),
  make2df("HGB", Ery),
  make2df("MCH", Ery),
  make2df("MCV", Ery),
  make2df("MCHC", Ery),
  make2df("RBC_COUNT", Ery),
  make2df("RETIC_COUNT", Ery)
)

# Return T/F whether rows in df1 are in df2
rowCheck <- function(df1, df2){
  xx <- apply(df1, 1, paste, collapse = "_")
  yy <- apply(df2, 1, paste, collapse = "_")
  return(xx %in% yy)
}

df <- read.table("../../data/bulk/GWAS-Bulk/compare3.tsv", header = TRUE)               
df$gchromVAR_pvalue <- pnorm(df$weighted_Zscore, lower.tail = FALSE)
df$chromVAR_pvalue <- pnorm(df$chromVAR_Zscore, lower.tail = FALSE)
df <- df[,c(1,2,5:7)] # drop z score column 
gs <- readRDS("../Figure1/OR_Heme.rds")
df_gs <- as.data.frame(gs,stringsAsFactors=F)
names(df_gs) <- c("i","j","cell","trait","obs","perm","z")
df_gs$pvalue <- 2*pnorm(-abs(as.numeric(df_gs$z)))
# manually verified that df_gs and df have the same order of cell type and traits so we good
df$goShifter_pvalue <- df_gs$pvalue
ldscradjust <- read.table("../../data/bulk/GWAS-Bulk/bulkHeme_baseline_panheme.txt", header = TRUE)
colnames(ldscradjust) <- c("Celltype", "Trait", "panHemeLDSR_pvalue")
odf <- merge(df, ldscradjust)
df <- odf
df$lineageSpecific <- rowCheck(df[,c("Celltype", "Trait")], lineageSpecificDF)

dim(df)

# Custom function to import fgwas stuff
cells <- c("B","CD4","CD8","CLP","CMP","Ery","GMP.A","GMP.B","GMP.C","HSC","LMPP","mDC","Mega","MEP","Mono","MPP","NK","pDC")
traits <- c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

lapply(cells, function(cell){
  lapply(traits, function(trait){
    
    file <- paste0("../../data/fgwas/", trait, "_", cell, ".params")
    tab <- read.table(file, header = TRUE)
    sigma <- (as.numeric(as.character(tab[2,4])) - as.numeric(as.character(tab[2,2])))/(2*1.96)
    
    # if sigma is an NA, due to large confidence interval spanning -20; manually fix with approximation
    sigma <- ifelse(is.na(sigma), 5, sigma)
    
    estimate <- tab[2,3]
    
    cell <- ifelse(cell == "GMP.A", "GMP-A", cell)
    cell <- ifelse(cell == "GMP.B", "GMP-B", cell)
    cell <- ifelse(cell == "GMP.C", "GMP-C", cell)
    cell <- ifelse(cell == "mono", "Mono", cell)
    
    data.frame(cell, trait, z = estimate/sigma)
  }) %>% rbindlist() %>% data.frame() -> odf
  odf
}) %>% rbindlist() %>% data.frame() -> fgwas

gregor <- read.table("../../data/gregor/gregor_combinedenrichments.txt", header = TRUE)
gpa <- data.frame(rbindlist(lapply(list.files("../../data/GPA", full.names = TRUE), read.table, header = TRUE)))

# Fix names from GPA from cell types
gpa$cell <- gsub("GMP1low", "GMP-A", gpa$cell)
gpa$cell <- gsub("GMP2mid", "GMP-B", gpa$cell)
gpa$cell <- gsub("GMP3high", "GMP-C", gpa$cell)
gpa$cell <- gsub("Bcell", "B", gpa$cell)
gpa$cell <- gsub("MEGA", "Mega", gpa$cell)
gpa$cell <- gsub("Erythro", "Ery", gpa$cell)

colnames(gregor) <- c("Trait", "Celltype", "gregor_pvalue")
colnames(gpa) <- c("Celltype", "Trait", "GPA_chisq")
colnames(fgwas) <- c("Celltype", "Trait", "fgwas_z")


mdf1 <- merge(df, gpa)
mdf2 <- merge(mdf1, gregor)
df <- merge(mdf2, fgwas)

permuted <- sapply(1:10000, function(i) sum(1:dim(df)[1] * sample(df$lineageSpecific, length(df$lineageSpecific))))
gregor_ranksum <- sum(1:dim(df)[1]*df[order(df$gregor_pvalue, decreasing = FALSE), "lineageSpecific"])
pnorm((mean(permuted) - gregor_ranksum)/sd(permuted), lower.tail = FALSE)

GPA_ranksum <- sum(1:dim(df)[1]*df[order(df$GPA_chisq, decreasing = TRUE), "lineageSpecific"])
pnorm((mean(permuted) - GPA_ranksum)/sd(permuted), lower.tail = FALSE)

fgwas_ranksum <- sum(1:dim(df)[1]*df[order(df$fgwas_z, decreasing = TRUE), "lineageSpecific"])
pnorm((mean(permuted) - fgwas_ranksum)/sd(permuted), lower.tail = FALSE)


