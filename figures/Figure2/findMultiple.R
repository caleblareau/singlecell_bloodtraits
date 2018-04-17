library(data.table)
library(GenomicRanges)
library(diffloop)
library(dplyr)
library(BuenColors)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

bedfiles <- list.files("../../data/UKBB_BC_PP001", full.names = TRUE, pattern = "*.bed")
counts <-  data.matrix(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt"))

# Import variants
df <- data.frame(rbindlist(lapply(bedfiles, read.table)))
df <- df[df$V5 > 0.5, c(1,2,3)]
df <- df[!duplicated(df),]
gr <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3")

ix <- subjectHits(findOverlaps(gr, peaks))[duplicated(subjectHits(findOverlaps(gr, peaks)))]
peaks_2cv <- peaks[ix]

# Function to pull out variants that localize in one peak
varsInRegion <- function(file, peak, cutoff_PP = 0.1){
  trait <- gsub(x = file, "../../data/UKBB_BC_PP001/", "") %>% gsub(pattern = "_PP001.bed", replace = "")
  gr <- sort(bedToGRanges(file))
  gr <- gr[mcols(gr)$score > cutoff_PP]
  ov <- findOverlaps(gr, peak)
  gr <- gr[queryHits(ov)]
  if(length(gr) > 0){
    return(data.frame(data.frame(gr)[,c(1,2,3,7)], trait))
  } else {
    return(data.frame(seqnames = NA, start = NA, end = NA, score = NA, trait = NA))
  }
}

lapply(1:length(peaks_2cv), function(i){
  lapply(bedfiles, varsInRegion, peaks_2cv[i]) %>% rbindlist() %>% as.data.frame() %>% filter(!is.na(seqnames)) -> odf
  odf$peakStart <- start(peaks_2cv[i])
  odf$peakEnd <- end(peaks_2cv[i])
  odf
}) %>% rbindlist() %>% as.data.frame() -> allVars

data.frame(peaks[ix], counts[ix, ])

