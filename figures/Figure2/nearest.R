
library(data.table)
library(GenomicRanges)
library(diffloop)
library(dplyr)
library(BuenColors)

# Create bulk Summarized Experiment
peaksdf <- fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

bedfiles <- list.files("../../data/UKBB_BC_PP001", full.names = TRUE, pattern = "*.bed")

# Define function to get the distance to nearest SNP
nearestDistance <- function(file, cutoff_PP = 0.001){
  trait <- gsub(x = file, "../../data/UKBB_BC_PP001/", "") %>% gsub(pattern = "_PP001.bed", replace = "")
  gr <- sort(bedToGRanges(file))
  gr <- gr[mcols(gr)$score > cutoff_PP]
  return(data.frame(distance = mcols(distanceToNearest(gr))$distance,
                    trait = trait, pp = cutoff_PP))
}

# Make plot of the distance density
lapply(c(0.001, 0.005, 0.01, 0.25, 0.5), function(pp) {
  lapply(bedfiles, nearestDistance, pp) %>% rbindlist() %>% as.data.frame()
}) %>% rbindlist() %>% as.data.frame() -> distdf

distdf$distance <- ifelse(distdf$distance > 10000, 10000, distdf$distance)
ggplot(distdf, aes(distance, fill = trait)) + facet_wrap(~as.factor(pp), scales = "free_y") +
  geom_histogram(binwidth = 100) + pretty_plot()

# Count occurences of variants in same ATAC peak
inSameATACPeak <- function(file, cutoff_PP = 0.001){
  trait <- gsub(x = file, "../../data/UKBB_BC_PP001/", "") %>% gsub(pattern = "_PP001.bed", replace = "")
  gr <- sort(bedToGRanges(file))
  gr <- gr[mcols(gr)$score > cutoff_PP]
  nPeaksMultiple <- length(subjectHits(findOverlaps(gr, peaks))[duplicated(subjectHits(findOverlaps(gr, peaks)))])
  return(data.frame(nPeaksMultiple, trait, cutoff_PP))
}

## Double Loop
lapply(c(0.001, 0.005, 0.01, 0.25, 0.5), function(pp) {
  lapply(bedfiles, inSameATACPeak, pp) %>% rbindlist() %>% as.data.frame()
}) %>% rbindlist() %>% as.data.frame() -> samePeakDF

ggplot(samePeakDF, aes(x = trait, y = nPeaksMultiple, color = as.factor(cutoff_PP))) + geom_point() + pretty_plot() + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(angle = 90)) +
  labs(x = "Trait", y = "Peaks with Multiple Variants", colour = "PP Cutoff") +
  scale_color_manual(values = jdb_palette("Zissou"))


