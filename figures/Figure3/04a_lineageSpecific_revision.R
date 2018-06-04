library(BuenColors)
library(ggplot2)
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

gregor <- read.table("../../data/gregor/gregor_combinedenrichments.txt", header = TRUE)
gpa <- data.frame(rbindlist(lapply(list.files("../../data/GPA", full.names = TRUE), read.table, header = TRUE)))

gpa$cell <- gsub("GMP1low", "GMP-A", gpa$cell)
gpa$cell <- gsub("GMP2mid", "GMP-B", gpa$cell)
gpa$cell <- gsub("GMP3high", "GMP-C", gpa$cell)
gpa$cell <- gsub("Bcell", "B", gpa$cell)
gpa$cell <- gsub("MEGA", "Mega", gpa$cell)
gpa$cell <- gsub("Erythro", "Ery", gpa$cell)

colnames(gregor) <- c("Trait", "Celltype", "gregor_pvalue")
colnames(gpa) <- c("Celltype", "Trait", "GPA_chisq")

mdf1 <- merge(df, gpa)
mdf2 <- merge(mdf1, gregor)
df <- mdf2

permuted <- sapply(1:10000, function(i) sum(1:dim(df)[1] * sample(df$lineageSpecific, length(df$lineageSpecific))))
gregor_ranksum <- sum(1:dim(df)[1]*df[order(df$gregor_pvalue, decreasing = FALSE), "lineageSpecific"])
pnorm((mean(permuted) - gregor_ranksum)/sd(permuted), lower.tail = FALSE)

GPA_ranksum <- sum(1:dim(df)[1]*df[order(df$GPA_chisq, decreasing = TRUE), "lineageSpecific"])
pnorm((mean(permuted) - GPA_ranksum)/sd(permuted), lower.tail = FALSE)

