library(data.table)

df <- readRDS("deviationDataframe.rds")
df$Var1 <- gsub("_PP001", "", df$Var1)

dfo <- df[,c(1,2,3,4,5)]
colnames(dfo) <- c("Trait", "Roadmap_ID", "gchromVAR_zscore", "P", "gchromVAR_logP")
dfo <- dfo[,c(1:4)]
dfo$FDR <- qvalue::qvalue(dfo$P)$qvalues

write.table(dfo, file = "SupplementalTable5.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
