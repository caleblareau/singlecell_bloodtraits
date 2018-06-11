library(dplyr)
library(diffloop)
library(GenomicRanges)
library(data.table)
library(BuenColors)

# Import variants
traits <- gsub("_PP001.bed", "", list.files("../../data/UKBB_BC_PP001/", patter = ".bed$"))
count_traits <- traits[grep("COUNT", traits)]

# JCU definition
map <- c("BASO_COUNT" = "GRAN", "EO_COUNT" = "GRAN", "NEUTRO_COUNT" = "GRAN", "RBC_COUNT" = "RBC", "PLT_COUNT" = "PLT", 
         "MONO_COUNT" = "MONO", "LYMPH_COUNT" = "LYMPH")

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
  t[complete.cases(t),]
}) %>% rbindlist () %>% as.data.frame() -> ukbb.df
ukbb.df$UP <- ukbb.df$Z > 0
ukbb.df$DOWN <- ukbb.df$Z < 0

# Identify pleiotropic variants
ukbb.df %>%
  group_by(chr, start, stop) %>%
  filter(n()>1) %>% arrange(desc(round(PP,2)),chr, start) %>%  summarise(totalUp = sum(UP), totalDown = sum(DOWN)) %>%
  mutate(switch = totalUp > 0 & totalDown > 0) %>% data.frame() -> summaryDF

dim(summaryDF)
sum(as.numeric(summaryDF$switch))

peaks <- bedToGRanges("../../data/bulk/ATAC/29August2017_EJCsamples_allReads_500bp.bed")
ukbb.df[ukbb.df$start %in% summaryDF[summaryDF$switch,"start"],]  %>% group_by(chr, start, stop) %>% summarize(mean(PP)) %>%
  data.frame() %>% makeGRangesFromDataFrame() -> switch_gr

ukbb.df[ukbb.df$start %in% summaryDF[!summaryDF$switch,"start"],]  %>% group_by(chr, start, stop) %>% summarize(mean(PP)) %>%
  data.frame() %>% makeGRangesFromDataFrame() -> tune_gr

ov_switch <- findOverlaps(peaks, switch_gr)
ov_tune <- findOverlaps(peaks, tune_gr)


