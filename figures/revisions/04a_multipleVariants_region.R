library(dplyr)
library(diffloop)
library(GenomicRanges)
library(data.table)
library(BuenColors)
library(cowplot)
"%ni%" <- Negate("%in%")

x01_coding_gr <- bedToGRanges("../../data/annotations/Coding_UCSC.bed")
x02_promoter_gr <- bedToGRanges("../../data/annotations/Promoter_UCSC.fixed.bed")
x03_utr_gr <- bedToGRanges("../../data/annotations/UTR_3_UCSC.bed")
x04_atac_gr <- bedToGRanges("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
x05_intron_gr <- bedToGRanges("../../data/annotations/Intron_UCSC.bed")

ex_var <- read.table("exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

traits <- gsub("_PP001.bed", "", list.files("../../data/UKBB_BC_PP001/", patter = ".bed$"))
lapply(traits, function(trait){
  
  print(trait)
  
  # Import and munge variants
  t <- read.table(paste0("../../data/UKBB_BC_PP001/betas_added/",trait,"_PP001_betas.bed"))[,c("V1", "V2", "V3", "V4","V5", "V8")]
  
  # Specify region
  reg_split <- stringr::str_split_fixed(as.character(t$V4), "-", 3)
  t$V4 <- reg_split[,3]
  
  # Filter out the exclude variants
  t <- t[reg_split[,2] %ni% ex_var, ]
  
  t <- t[t$V5 > 0.5, ]
  t$trait <- trait
  colnames(t) <- c("chr", "start", "stop", "region", "PP", "Z", "trait")
  gr_t <- unique(makeGRangesFromDataFrame(t, keep.extra.columns = TRUE))
  
  # Do overlaps
  ov_1 <- findOverlaps(gr_t, x01_coding_gr)
  ov_2 <- findOverlaps(gr_t, x02_promoter_gr)
  ov_3 <- findOverlaps(gr_t, x03_utr_gr)
  ov_4 <- findOverlaps(gr_t, x04_atac_gr)
  ov_5 <- findOverlaps(gr_t, x05_intron_gr)
  
  # Classify each variant
  class <- ifelse(1:length(gr_t) %in% queryHits(ov_1), "coding",
                  ifelse(1:length(gr_t) %in% queryHits(ov_2), "promoter",
                         ifelse(1:length(gr_t) %in% queryHits(ov_3), "utr",
                                ifelse(1:length(gr_t) %in% queryHits(ov_4), "accessible",
                                       ifelse(1:length(gr_t) %in% queryHits(ov_5), "intron", "intergenic")))))
  
  # Find nearest -- automatically handles chromosome
  df_raw <- data.frame(distanceToNearest(gr_t))
  
  lapply(1:dim(df_raw)[1], function(i){
    min_idx <- df_raw[i,1]
    max_idx <- df_raw[i,2]
    
    if( class[min_idx] > class[max_idx]){
      temp <- min_idx
      min_idx <- max_idx
      max_idx <- temp
    } else if(class[min_idx] == class[max_idx]){
      tempmin <- min(min_idx, max_idx)
      tempmax <- max(min_idx, max_idx)
      min_idx <- tempmin
      max_idx <- tempmax
    }
    
    region_vec <- mcols(gr_t)[,"region"]
    
    data.frame(Variant1 = as.character(as.character(gr_t)[min_idx]), Class1 = class[min_idx], PP1 = mcols(gr_t)[min_idx,"PP"],
               Variant2 = as.character(as.character(gr_t)[max_idx]), Class2 = class[max_idx], PP2 = mcols(gr_t)[max_idx,"PP"],
               trait = trait, distance = df_raw[i,3],
               regionA = region_vec[min_idx], regionB = region_vec[max_idx])
    
  }) %>% rbindlist() %>% data.frame() %>% distinct() %>% data.frame()-> out
  
  out
  
}) %>% rbindlist () %>% as.data.frame() -> nearestPairsDF_region

# Filter such that they are in the same region
nearestPairsDF_region %>% filter(as.character(regionA) == as.character(regionB)) %>% arrange(distance) -> PPvariantPairs

if(FALSE){
  rsid_vec <- readRDS("rsid_vector.rds")
  PPvariantPairs$Var1 <- stringr::str_split_fixed(PPvariantPairs$Variant1, "-", 2)[,1]
  PPvariantPairs$Var2 <- stringr::str_split_fixed(PPvariantPairs$Variant2, "-", 2)[,1]
  PPvariantPairs %>% mutate(rsID1 = rsid_vec[Var1], rsID2 = rsid_vec[Var2]) %>%
    dplyr::select(Var1, rsID1, Class1, PP1, 
                  Var2, rsID2, Class2, PP2, 
                  trait, distance) -> refined_region
  write.table(refined_region, file = "final_supplemental_tables/nearVariantsRegionPP50.tsv", sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

PPvariantPairs %>% dplyr::select(Variant1, Class1, Variant2, Class2) %>% distinct() %>%
  group_by(Class1, Class2) %>% summarise(count = n()) %>% data.frame() -> odf

# swap stuff around -- painful
odf[odf$Class1 == "accessible" & odf$Class2 == "promoter", c("Class1", "Class2")] <-
  odf[odf$Class1 == "accessible" & odf$Class2 == "promoter", c("Class2", "Class1")] 
odf[odf$Class1 == "accessible" & odf$Class2 == "coding", c("Class1", "Class2")] <-
  odf[odf$Class1 == "accessible" & odf$Class2 == "coding", c("Class2", "Class1")] 
odf[odf$Class1 == "intron" & odf$Class2 == "promoter", c("Class1", "Class2")] <-
  odf[odf$Class1 == "intron" & odf$Class2 == "promoter", c("Class2", "Class1")] 
odf[odf$Class1 == "intergenic" & odf$Class2 == "promoter", c("Class1", "Class2")] <-
  odf[odf$Class1 == "intergenic" & odf$Class2 == "promoter", c("Class2", "Class1")] 

# Swap around the UTR stuff
odf[odf$Class1 == "intergenic" & odf$Class2 == "intron", c("Class1", "Class2")] <-
  odf[odf$Class1 == "intergenic" & odf$Class2 == "intron", c("Class2", "Class1")] 

odf[odf$Class1 == "accessible" & odf$Class2 == "utr", c("Class1", "Class2")] <-
  odf[odf$Class1 == "accessible" & odf$Class2 == "utr", c("Class2", "Class1")] 

odf[odf$Class1 == "intron" & odf$Class2 == "utr", c("Class1", "Class2")] <-
  odf[odf$Class1 == "intron" & odf$Class2 == "utr", c("Class2", "Class1")] 

#odf <- rbind(odf, data.frame(Class1 = "utr", Class2 = "intergenic", count = 0))

# Order for final presentation
odf$Class1 <- factor(as.character(odf$Class1), levels = c("coding", "promoter", "utr", "accessible", "intron", "intergenic"))
odf$Class2 <- factor(as.character(odf$Class2), levels = rev(c("coding", "promoter", "utr", "accessible", "intron", "intergenic")))

p1 <- ggplot(odf, aes(Class1, Class2,fill = count)) + geom_tile( color = "black") +geom_text(aes(label = count), size = 2) +
  scale_fill_gradientn(colors = jdb_palette("brewer_red")[1:6]) +
  labs(x = "", y = "") + pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

cowplot::ggsave(p1, file = "nearestVariantOut/PP50_region_heatmap.pdf", height = 2, width = 2)
# write.table(PPvariantPairs, file = "nearestVariantOut/twoVar_region10kb_table.tsv",
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


