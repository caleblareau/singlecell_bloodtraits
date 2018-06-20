library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(BuenColors)

set.seed(14561)
"%ni%" <- Negate("%in%")

# Import gtf
if(FALSE){
  gtf <- import("../../data/revision/gencode.v19.annotation.gtf.gz")
  exdf <- stringr::str_split_fixed(read.table("exclude_list.txt", header = FALSE, stringsAsFactors = FALSE)[,1], "_", 4)
  
}

map <- c("BASO_COUNT" = "GRAN", "EO_COUNT" = "GRAN", "NEUTRO_COUNT" = "GRAN", "RBC_COUNT" = "RBC", "PLT_COUNT" = "PLT", 
         "MONO_COUNT" = "MONO", "LYMPH_COUNT" = "LYMPH", "MCHC" = "RBC", "MCH" = "RBC", "HGB" = "RBC", "HCT" = "RBC", 
         "MPV" = "PLT", "RETIC_COUNT" = "RBC", "MEAN_RETIC_VOL" = "RBC", "MCV" = "RBC")
length(map)

lapply(names(map), function(trait){
  
  # Import and filter for PP > 0.1
  t <- read.table(paste0("../../data/UKBB_BC_PP001/betas_added/",trait,"_PP001_betas.bed"))[,c("V1", "V2", "V4", "V5", "V8")]
  t <- t[t$V5 > 0.1, ]
  
  # Filter for multi-region variants
  t %>% group_by(V1, V2, V8) %>% summarise(long = head(V4,1), V5 = max(V5)) %>% data.frame() -> t
  
  # Annotate with trait and lineage
  lineage <- map[trait]
  t$trait <- trait
  t$lineage <- lineage
  colnames(t) <- c("chr", "start", "Z", "long",  "PP", "trait", "lineage")
  t$stop <- t$start + 1
  gr <- makeGRangesFromDataFrame(t)
  
  # Find coding variants
  ov <- findOverlaps(gr, gtf)
  data.frame(idx = queryHits(ov), gene = mcols(gtf)$gene_name[subjectHits(ov)]) %>% unique() %>% data.frame() -> geneDf
  geneDf <- geneDf[!grepl("^RP",geneDf$gene),]
  
  # Make them accessible
  vec <- as.character(geneDf$gene); names(vec) <- as.character(geneDf$idx)
  t$gene <- vec[as.character(1:dim(t)[1])]
  
  # Exclude problematic variant
  t <- t[complete.cases(t),]
  t$ID <- paste0(gsub("chr", "", t$chr), ":", as.character(t$start))
  ukbb.df <- t[t$ID %ni% exdf[,1],]
  ukbb.df
  
}) %>% rbindlist () %>% as.data.frame() -> ukbb.df

if(FALSE){
  v <- stringr::str_split_fixed(ukbb.df$long, "-", 3)[,2]
  x <- stringr::str_split_fixed(v, "_", 3)
  x2 <- stringr::str_split_fixed(x[,1], ":", 2)
  odf <- data.frame(chr = x2[,1], one = x2[,2], two = x2[,2], 
                    allele = paste0(x[,2],"/",x[,3]), a = 1)
  write.table(odf, file = "coding_forVEP.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

ukbb.df %>% group_by(chr, start) %>% summarize(gene = head(long,1)) -> sumDF

# Import RNA
x <- read.table("../../data/bulk/RNA/16populations_RNAcounts.txt", header = TRUE)
rownames(x) <- make.unique(as.character(x$Genes))
x <- data.matrix(x[,2:17])

# Min / Max Normalize
counts <- x/colSums(x) * 1000000
log2cpm <- log2(counts + 1)
log2cpm_minmax <- (((log2cpm) - matrixStats::rowMins(log2cpm))/matrixStats::rowMaxs(log2cpm))

# Cluster
merged_mega <- merge(ukbb.df, log2cpm_minmax, by.x = "gene", by.y = "row.names")
merged_mega <- merged_mega[complete.cases(merged_mega), ]
all.RBC <- unique(merged_mega[merged_mega$lineage == "RBC", c("HSC", "MPP", "CMP", "MEP", "Ery")])
all.PLT <- unique(merged_mega[merged_mega$lineage == "PLT", c("HSC", "MPP", "CMP", "MEP", "Ery")])
all.GRAN <- unique(merged_mega[merged_mega$lineage == "GRAN", c("HSC", "MPP", "CMP", "MEP", "Ery")])
all.MONO <- unique(merged_mega[merged_mega$lineage == "MONO", c("HSC", "MPP", "CMP", "MEP", "Ery")])
all.LYMPH <- unique(merged_mega[merged_mega$lineage == "LYMPH", c("HSC", "MPP", "CMP", "MEP", "Ery")])

dim(all.RBC)
dim(all.PLT)
dim(all.GRAN)
dim(all.MONO)
dim(all.LYMPH)

km <- kmeans(all.RBC, centers = 6, nstart = 1000)
km.cluster <- factor(km$cluster)
#km.cluster <- factor(km$cluster, km$cluster, levels = c("2", "3", "1"))
hm <- Heatmap(all.RBC, col=as.character(jdb_palette("solar_extra",type="continuous")),
              cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 6),
              split = km.cluster, show_heatmap_legend = FALSE,
              name = "Gene\nExpression")
hm
