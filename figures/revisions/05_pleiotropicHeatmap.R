library(dplyr)
library(diffloop)
library(GenomicRanges)
library(data.table)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)

"%ni%" <- Negate("%in%")

# Import variants
traits <- gsub("_PP001.bed", "", list.files("../../data/UKBB_BC_PP001/", patter = ".bed$"))
count_traits <- c(traits[grep("COUNT", traits)])

# JCU definition of lineage
map <- c("BASO_COUNT" = "GRAN", "EO_COUNT" = "GRAN", "NEUTRO_COUNT" = "GRAN", "RBC_COUNT" = "RBC", "PLT_COUNT" = "PLT", 
         "MONO_COUNT" = "MONO", "LYMPH_COUNT" = "LYMPH")

# UKBB exclusion list
exdf <- stringr::str_split_fixed(read.table("exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1], "_", 4)

lapply(count_traits, function(trait){
  
  # Import and filter for PP > 0.1
  t <- read.table(paste0("../../data/UKBB_BC_PP001/betas_added/",trait,"_PP001_betas.bed"))[,c("V1", "V2", "V3", "V5", "V8")]
  t <- t[t$V5 > 0.1, ]
  
  # Filter for multi-region variants
  t %>% group_by(V1, V2, V3, V8) %>% summarise(V5 = max(V5)) %>% data.frame() -> t
  
  # Annotate with trait and lineage
  lineage <- map[trait]
  t$trait <- trait
  t$lineage <- lineage
  colnames(t) <- c("chr", "start", "stop", "Z",  "PP","trait", "lineage")
  t$UP <- t$Z > 0
  t$DOWN <- t$Z < 0
  
  t[complete.cases(t),]
}) %>% rbindlist () %>% as.data.frame() -> ukbb.df
ukbb.df$ID <- paste0(gsub("chr", "", ukbb.df$chr), ":", as.character(ukbb.df$start))
ukbb.df <- ukbb.df[ukbb.df$ID %ni% exdf[,1],]

# Identify pleiotropic variants
ukbb.df %>%
  group_by(chr, start, stop) %>%
  filter(n()>1) %>% arrange(desc(round(PP,2)),chr, start) %>%  summarise(totalUp = sum(UP), totalDown = sum(DOWN)) %>%
  mutate(switch = totalUp > 0 & totalDown > 0) %>% data.frame() -> summaryDF

dim(summaryDF)
sum(as.numeric(summaryDF$switch))

# Annotate with peaks or coding variants
peaks <- bedToGRanges("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
coding <- bedToGRanges("../../data/annotations/Coding_UCSC.bed")

ukbb.df[ukbb.df$start %in% summaryDF[summaryDF$switch,"start"],]  %>% group_by(chr, start, stop) %>% summarize(mean(PP)) %>%
  data.frame() %>% makeGRangesFromDataFrame() -> switch_gr

ukbb.df[ukbb.df$start %in% summaryDF[!summaryDF$switch,"start"],]  %>% group_by(chr, start, stop) %>% summarize(mean(PP)) %>%
  data.frame() %>% makeGRangesFromDataFrame() -> tune_gr


# Bin up pleitropy
x01_coding_gr <- bedToGRanges("../../data/annotations/Coding_UCSC.bed")
x02_promoter_gr <- bedToGRanges("../../data/annotations/Promoter_UCSC.fixed.bed")
x03_utr_gr <- bedToGRanges("../../data/annotations/UTR_3_UCSC.bed")
x04_atac_gr <- bedToGRanges("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
x05_intron_gr <- bedToGRanges("../../data/annotations/Intron_UCSC.bed")

assignVariantCategory <- function(gr_t){
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
  
  class
}


# Set up color annotation
annotationColors <- jdb_palette("brewer_spectra")[c(1,3,4,6,7,9)]
names(annotationColors) <- c("coding", "promoter", "utr", "intergenic", "intron", "accessible")

# Make Heatmap
tune_df <- data.frame(tune_gr)[,c(1,2,3)]; tune_df$class <- assignVariantCategory(tune_gr)
switch_df <- data.frame(switch_gr)[,c(1,2,3)]; switch_df$class <- assignVariantCategory(switch_gr)

# Function to annotate per variant what traits are implicated
makeBinaryPleiotropyMatrix <- function(df){
  lapply(1:dim(df)[1], function(i){
    
    # Define a mini function for togetherness
    together <- function(trait){
      ukbb.df %>% filter(chr == df[i,"seqnames"] & start == df[i,"start"]) -> subset
      return(as.numeric(trait %in% as.character(subset$lineage)))
    }
    
    data.frame(GRAN = together("GRAN"), 
               LYMPH = together("LYMPH"),
               MONO = together("MONO"),
               PLT = together("PLT"),
               RBC = together("RBC"))
  }) %>% rbindlist() %>% data.frame() -> mat
  return(mat)
}

tune_df <- cbind(tune_df, makeBinaryPleiotropyMatrix(tune_df))
switch_df <- cbind(switch_df, makeBinaryPleiotropyMatrix(switch_df))

tune_df %>% arrange(class, GRAN, LYMPH, MONO, PLT, RBC) -> tune_df
switch_df %>% arrange(class, GRAN, LYMPH, MONO, PLT, RBC) -> switch_df


if(FALSE){
  tune_df$Annotation <- "tune"
  switch_df$Annotation <- "switch"
  all_df <- rbind(tune_df, switch_df) 
  
  rsid_vec <- readRDS("rsid_vector.rds")
  all_df$Variant <- paste0(all_df$seqnames, ":", as.character(all_df$start))
  
  all_df %>% mutate(rsID = rsid_vec[Variant]) %>%
    dplyr::select(Variant, rsID, class, Annotation, GRAN, LYMPH, MONO, PLT, RBC) -> refined_pleio
  write.table(refined_pleio, file = "final_supplemental_tables/pleiotropic.tsv", sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
}

ha1 <- HeatmapAnnotation(df = data.frame(class = c(tune_df$class, switch_df$class)),
                         col = list(class = annotationColors)
)
mat <- cbind(t(data.matrix(tune_df[,c("GRAN", "LYMPH", "MONO", "PLT", "RBC")])),
             t(data.matrix(switch_df[,c("GRAN", "LYMPH", "MONO", "PLT", "RBC")])))
colnames(mat) <- paste0("lol", as.character(1:dim(mat)[2]))

pdf(file="plots/pleiotropy_heatmap.pdf", width = 6, height = 1.5)  
par(cex.main=0.8,mar=c(1,1,1,1))

Heatmap(mat, col=as.character(c("white", "black")),
        cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 0),
        bottom_annotation = ha1,
        name = "")
dev.off()



if(FALSE){
  
  countdf <- data.frame(table(assignVariantCategory(tune_gr)),
                        table(assignVariantCategory(switch_gr)))
  colnames(countdf) <- c("what", "tune", "what2", "switch")
  meltdf <- reshape2::melt(countdf[,c("what", "tune", "switch")], id.vars = "what")
  
  # Append
  meltdf$label_count <- paste0(as.character(meltdf$value))
  meltdf$percent <- meltdf$value/c(rep(sum( countdf[,"tune"]), 5),
                                   rep(sum( countdf[,"switch"]), 5)) * 100
  
  p1 <- ggplot(meltdf, aes(x = variable, y = percent, fill = what, label = label_count)) +
    geom_bar(stat = "identity", color = "black", alpha=0.5) +
    geom_text(size = 2, position = position_stack(vjust = 0.5)) +
    pretty_plot(fontsize = 8) + theme(legend.position = "bottom") +
    labs(x = "", y = "% of Variants", fill= "") +
    scale_fill_manual(values = annotationColors) + L_border()
  cowplot::ggsave(p1, file = "plots/pleiotropic_divy.pdf", width = 2, height = 2.5)
}
