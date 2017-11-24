library(BuenColors)
library(dplyr)
library(data.table)

df <- readRDS("allIterations.rds")
df %>% group_by(Var1, Var2, npeaks) %>% summarize(zscore = mean(value)) -> meanDF

p1 <- ggplot(meanDF, aes(x = npeaks, y = -log10(pnorm(zscore,lower.tail = FALSE)))) +
  geom_line(aes(color = Var1)) + geom_point(aes(color = Var1)) +
  pretty_plot() + labs(x = "Number of Background Peaks", y = "g-chromVAR Enrichment (-log10 p)", color = "") +
  facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") + scale_color_manual(values = ejc_color_maps) +
  geom_hline(yintercept = -log10(0.05 / (16*18)), linetype = 2) 

ggsave(p1, filename = "npeaks.pdf")