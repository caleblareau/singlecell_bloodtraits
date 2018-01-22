library(dplyr)
library(qqman)

df <- read.table("pics_gchromVAR.txt", header = TRUE)
head(df, 20)
bonf <- 0.05/dim(df)[1]
bonf
df$FDR <- qvalue::qvalue(df$pvalue)$qvalues
write.table(df, file = "pics_gchromVAR_jan22.txt", row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)
