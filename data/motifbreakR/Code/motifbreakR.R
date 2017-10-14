#!/usr/bin/env Rscript

# Load Libraries
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(S4Vectors)
library(GenomicRanges)
library(matrixStats)
library(GO.db)
library(TFBSTools)
library("BiocParallel")
library(dplyr)

# First argument is path to input file without file type
args <- commandArgs(trailingOnly=TRUE)
snpList <- snps.from.file(file=args[1],search.genome=BSgenome.Hsapiens.UCSC.hg19,format="bed")
strand(snpList) <- "+"
trait <- args[2]


# Load processed ChromVar MotifList that can be used with motifBreaker
#pwmList <- readRDS("/broad/sankaranlab/ebao/MotifbreakR/Chromvar_motiflist.rds")
#data(motifbreakR_motif)
data(hocomoco)

results <- motifbreakR(snpList = snpList,
                       pwmList = hocomoco,
                       filterp = TRUE,
                       threshold = 1e-4,
                       show.neutral=FALSE,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())


saveRDS(results,paste0("/broad/sankaranlab/ebao/MotifbreakR/output/PP001_allvariants/",trait,"_Motifbreakr_output_PP001.rds"))
