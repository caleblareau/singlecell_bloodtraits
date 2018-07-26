library(dplyr)
library(diffloop)
library(GenomicRanges)
library(data.table)
library(BuenColors)
library(cowplot)
library(stringr)
library(plotly)
"%ni%" <- Negate("%in%")

# Import and normalize ATAC signal
x04_atac_gr <- bedToGRanges("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
counts <- data.matrix(data.frame(fread("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.counts.txt")))
ATAC.cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)

# Import variants
df <- read.table("final_supplemental_tables/nearVariantsRegionPP50.tsv", header = TRUE)
df %>% filter(Class1 %in% c( "accessible") & Class2 %in% c( "accessible")) %>% filter (distance < 1000000) -> df
vars1_gr <- makeGRangesFromDataFrame(data.frame(str_split_fixed(as.character(df$Var1), ":", 2)),
                                     seqnames.field = "X1", start.field = "X2", end.field = "X2")

vars2_gr <- makeGRangesFromDataFrame(data.frame(str_split_fixed(as.character(df$Var2), ":", 2)),
                                     seqnames.field = "X1", start.field = "X2", end.field = "X2")

ov1 <- findOverlaps(vars1_gr, x04_atac_gr)
ov2 <- findOverlaps(vars2_gr, x04_atac_gr)
ov_2val <- data.frame(hits1 = subjectHits(ov1), 
                      hits2 = subjectHits(ov2), 
                      variants = paste0(df$Var1[queryHits(ov1)], "--", df$Var2[queryHits(ov2)]),
                      trait = df$trait, distance = df$distance)
ov_2val %>% filter(hits1 != hits2) -> ov_2uniqueEnhancers

lapply(1:dim(ov_2uniqueEnhancers)[1], function(i){
  vals1 <- log2(ATAC.cpm[ov_2uniqueEnhancers[i,1],] + 1)
  vals2 <- log2(ATAC.cpm[ov_2uniqueEnhancers[i,2],] + 1)
  data.frame(vals1, vals2, cellType = names(vals2), variants = ov_2uniqueEnhancers[i,3],
             trait = ov_2uniqueEnhancers[i,4], distance = ov_2uniqueEnhancers[i,5])
  
}) %>% rbindlist() %>% data.frame() -> odf

odf$cellType2 <- ifelse(odf$cellType == "GMP.A", "GMP-A",
                        ifelse(odf$cellType == "GMP.B", "GMP-B",
                               ifelse(odf$cellType == "GMP.C", "GMP-C", as.character(odf$cellType))))
                        

ggplot(odf, aes(x = vals1, y = vals2, color = cellType2)) + geom_point() +
  scale_color_manual(values = ejc_color_maps) +
   geom_abline(intercept = 0, slope = 1, linetype = 2) + pretty_plot(fontsize = 7) + 
    theme(legend.position = "none") + L_border() + 
  scale_x_continuous(limits = c(0,6)) +
  scale_y_continuous(limits = c(0,6)) +
  labs(x = "Accessibility at variant 1", y = "Accessibility at variant 2") -> p1

cowplot::ggsave(p1, file = "nearestVariantOut/accessibility2plot.pdf", width = 2, height = 2)
