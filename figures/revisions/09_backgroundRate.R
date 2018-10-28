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
  
  df <- t[complete.cases(t),]
  data.frame(trait , UP = sum(df$UP), DOWN = sum(df$DOWN))
}) %>% rbindlist () %>% as.data.frame() -> ukbb.df

