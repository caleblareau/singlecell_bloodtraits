library(data.table)
library(BuenColors)

gene <- "AK3"

cors <- data.frame(fread(paste0("zcat < ", "../../data/bulk/peakGeneCorrelation.tsv.gz")))
