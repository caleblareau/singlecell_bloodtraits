library(BuenColors)
library(dplyr)
library(data.table)

importT <- function(i){
  t <- readRDS(paste0("variablePP/variablePPpeak_iteration_",as.character(i),".rds"))
  t$it <- i
  return(t)
}

df <- lapply(1:100, importT) %>% rbindlist()
df$Var2 <- gsub("_PP001", "", df$Var2)
df %>% group_by(Var1, Var2, PP) %>% summarize(zscore = mean(value)) -> meanDF

saveRDS(meanDF, "../../supplement/final_figures/FigureS4/variedPP_meanDF.rds")

#p1 <- ggplot(meanDF, aes(x = PP, y = -log10(pnorm(zscore,lower.tail = FALSE)))) +
p1 <- ggplot(meanDF, aes(x = PP, y = zscore)) +
  geom_line(aes(color = Var1)) + geom_point(aes(color = Var1)) +
  pretty_plot() + labs(x = "Posterior Probability Cutoff", y = "g-chromVAR Enrichment (-log10 p)", color = "") +
  facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") + scale_color_manual(values = ejc_color_maps) +
#  geom_hline(yintercept = -log10(0.05 / (16*18)), linetype = 2) 
  geom_hline(yintercept = 0, linetype = 2) 

ggsave(p1, filename = "variablePP.pdf")
