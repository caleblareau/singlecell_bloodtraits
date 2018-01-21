library(dplyr)
library(qqman)

df <- read.table("pics_gchromVAR.txt", header = TRUE)
head(df, 20)
bonf <- 0.05/dim(df)[1]
bonf
df$FDR <- qvalue::qvalue(df$pvalue)$qvalues
