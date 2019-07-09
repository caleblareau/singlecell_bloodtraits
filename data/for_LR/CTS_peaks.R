#' ---
#' title: "Peaks for LR"
#' author: "Jacob Ulirsch"
#' date: "July 9, 2019"
#' ---

library(data.table)
library(tidyr)
library(readr)
library(rtracklayer)
library(preprocessCore)

# Import peaks and counts
peaks.gr <- import("../bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed",format="bed")
counts.df <-  data.matrix(fread("../bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))
values(peaks.gr) <- counts.df

# Filter to cell-type specific peaks
filterPeaks <- function(df, gr, n) {
  df[1:dim(df)[1],1:dim(df)[2]] <- normalize.quantiles(as.matrix(df))
  out <- lapply(1:dim(df)[2], function(x) {peaks.gr[df[,x] > quantile(df[,x], n)]})
  names(out) <- colnames(df)
  return(out)
}

CTS_090.l <- filterPeaks(counts.df, peaks.gr, 0.9)
CTS_080.l <- filterPeaks(counts.df, peaks.gr, 0.8)

# Write out
lapply(1:length(CTS_090.l), function(x) {export(CTS_090[[x]], paste0("CTS_090/", names(CTS_090)[x], ".bed.gz"))})
lapply(1:length(CTS_080.l), function(x) {export(CTS_080[[x]], paste0("CTS_080/", names(CTS_080)[x], ".bed.gz"))})


