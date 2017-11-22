library(data.table)
library(dplyr)
library(annotables)
library(dplyr)

gtex <- fread(paste0('zcat < ', "../../data/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz"))
pj <- stringr::str_split_fixed(gtex$variant_id, "_", n = 5)
gtex$chr <- paste0("chr", pj[,1])
gtex$bp <- as.numeric(pj[,2])
gtex$ref <- pj[,3]
gtex$alt <- pj[,4]

lapply(list.files("../../data/UKBB_BC_PP001/", pattern = "*.bed$", full.names = TRUE), function(f){
  dt <- fread(f)
  dt$trait <- f
  dt
}) %>% rbindlist() %>% as.data.frame() -> allVars


allVars$trait <- gsub("_PP001.bed", "", gsub("../../data/UKBB_BC_PP001//", "", allVars$trait))
mdf <- merge(allVars, gtex, by.x = c("V1", "V2"), by.y = c("chr", "bp"))
cleandf <- mdf[,c("V1", "V2", "trait", "V5", "gene_id", "tss_distance", "pval_nominal", "ref", "alt", "slope")]
var <- stringr::str_split_fixed(cleandf$gene_id, "[.]", n = 2)
cleandf$gene_id2 <- var[,1]
cleandf2 <- merge(cleandf, grch37, by.x = "gene_id2", by.y = "ensgene", all.x = TRUE)
df <- cleandf2[,c("V1", "V2", "trait", "V5", "gene_id", "tss_distance",
                  "pval_nominal", "symbol", "description", "ref", "alt", "slope")]
saveRDS(df, file = "bloodEQTLS_gtex.rds")

df <- readRDS("bloodEQTLS_gtex.rds")

dff <- df[df$V1 != "chr17" & df$V1 != "chr6" & df$V5 > 0.5, ]
sort(table(unique(dff[,c("V2", "symbol", "trait")])[,c(1)]), decreasing = TRUE) %>% head(40)
dff[dff$V2 == "34195572", ]


ss <- (df[df$V5 > 0.05 ,c("V2", "symbol", "tss_distance")])
a <- ss %>% group_by(symbol) %>% summarise(mint = min(tss_distance), maxt = max(tss_distance))
a$diff <- a$maxt - a$mint
head(a[order(a$diff, decreasing = TRUE),], 20)

head(sort(table(unique(ss)$symbol), decreasing = TRUE), 80)

df[abs(df$tss_distance) > 100000 & df$V5 > 0.5 & df[,1] != "chr17",1:8 ]


LD <- readRDS("../../data/Finemap/LD80.rds")
plt <- LD[["PLT_COUNT"]]
plt[(plt$position1 %in% df[df$V5 > 0.05 & df$symbol == "RHD","V2"]) | (plt$position2 %in% df[df$V5 > 0.05 & df$symbol == "RHD","V2"]), ]

rbc <- LD[["RBC_COUNT"]]
rbc[(rbc$position1 %in% df[df$V5 > 0.1 & df$symbol == "ABO","V2"]) | (rbc$position2 %in% df[df$V5 > 0.05 & df$symbol == "ABO","V2"]), ]
