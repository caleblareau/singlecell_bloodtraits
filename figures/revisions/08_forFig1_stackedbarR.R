library(dplyr)
library(diffloop)
library(GenomicRanges)
library(data.table)
library(BuenColors)

"%ni%" <- Negate("%in%")

x01_coding_gr <- bedToGRanges("../../data/annotations/Coding_UCSC.bed")
x02_promoter_gr <- bedToGRanges("../../data/annotations/Promoter_UCSC.fixed.bed")
x03_utr_gr <- bedToGRanges("../../data/annotations/UTR_3_UCSC.bed")
x04_atac_gr <- bedToGRanges("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
x05_intron_gr <- bedToGRanges("../../data/annotations/Intron_UCSC.bed")
ex_var <- read.table("exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

traits <- gsub("_PP001.bed", "", list.files("../../data/UKBB_BC_PP001/", patter = ".bed$"))
lapply(traits, function(trait){
  
# Import and munge variants
  t <- read.table(paste0("../../data/UKBB_BC_PP001/betas_added/",trait,"_PP001_betas.bed"))[,c("V1", "V2", "V3", "V4","V5", "V8")]
  
  # Specify region
  reg_split <- stringr::str_split_fixed(as.character(t$V4), "-", 3)
  t$V4 <- reg_split[,3]
  
  # Filter out the exclude variants
  t <- t[reg_split[,2] %ni% ex_var, ]
  
  t <- t[t$V5 > 0.1, ]
  t$trait <- trait
  colnames(t) <- c("chr", "start", "stop", "region", "PP", "Z", "trait")
  gr_t <- unique(makeGRangesFromDataFrame(t, keep.extra.columns = TRUE))
  t <- data.frame(gr_t)
  # Do overlaps
  ov_1 <- findOverlaps(gr_t, x01_coding_gr)
  ov_2 <- findOverlaps(gr_t, x02_promoter_gr)
  ov_3 <- findOverlaps(gr_t, x03_utr_gr)
  ov_4 <- findOverlaps(gr_t, x04_atac_gr)
  ov_5 <- findOverlaps(gr_t, x05_intron_gr)
  
  # Classify each variant
  t$class <- ifelse(1:length(gr_t) %in% queryHits(ov_4), "coding",
                  ifelse(1:length(gr_t) %in% queryHits(ov_2), "promoter",
                         ifelse(1:length(gr_t) %in% queryHits(ov_3), "utr",
                                ifelse(1:length(gr_t) %in% queryHits(ov_1), "accessible",
                                       ifelse(1:length(gr_t) %in% queryHits(ov_5), "intron", "intergenic")))))

  t
  
}) %>% rbindlist () %>% as.data.frame() -> variant_class

variant_class %>% group_by(trait, class) %>% summarise(count = n()) %>% ungroup %>%
  group_by(trait) %>% mutate(freq = count / sum(count))-> groups_count

# Order for final presentation
groups_count$class <- factor(as.character(groups_count$class), levels = rev(c("coding", "promoter", "utr", "accessible", "intron", "intergenic")))
annotationColors <- jdb_palette("brewer_spectra")[c(1,3,4,6,7,9)]
names(annotationColors) <- c("coding", "promoter", "utr", "intergenic", "intron", "accessible")

p1 <- ggplot(groups_count, aes(trait, fill = class)) + 
  geom_bar(aes(y = freq), stat = "identity") + coord_flip() + pretty_plot(fontsize = 6) +
  L_border() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
  labs(x = "", y = "") + scale_fill_manual(values = annotationColors) +
  scale_y_continuous( expand = c(0, 0))

cowplot::ggsave(p1, file = "cumulateiveFraction.pdf", height = 1.85, width = 3.1)

