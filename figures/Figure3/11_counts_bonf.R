library(dplyr)

df <- readRDS("allEnrichments-df.rds")

# Equivalent statistics 
 p = 0.0001736111
 Z = 3.75
 chisq = 14.1

allstats <- data.frame(
  Yes = c(sum(df$lineageSpecific& df$ldscore_pvalue < p), sum(df$lineageSpecific& df$chromVAR_pvalue < p),
          sum(df$lineageSpecific& df$gchromVAR_pvalue < p), sum(df$lineageSpecific& df$panHemeLDSR_pvalue < p),
          sum(df$lineageSpecific& df$goShifter_pvalue < p), sum(df$lineageSpecific& df$GPA_chisq > chisq),
          sum(df$lineageSpecific& df$gregor_pvalue < p), sum(df$lineageSpecific& df$fgwas_z > Z)
          ),
  No = c(sum(!df$lineageSpecific& df$ldscore_pvalue < p), sum(!df$lineageSpecific& df$chromVAR_pvalue < p),
         sum(!df$lineageSpecific& df$gchromVAR_pvalue < p), sum(!df$lineageSpecific& df$panHemeLDSR_pvalue < p),
         sum(!df$lineageSpecific& df$goShifter_pvalue < p),  sum(!df$lineageSpecific& df$GPA_chisq > chisq),
         sum(!df$lineageSpecific& df$gregor_pvalue < p), sum(!df$lineageSpecific& df$fgwas_z > Z)
  )
)

rownames(allstats) <- c("LD score", "chromVAR", "gchromVAR", "adjLDS", "goShifter", "GPA", "GREGOR", "FGWAS")

allstats
