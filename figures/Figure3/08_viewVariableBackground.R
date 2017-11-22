library(BuenColors)
library(magrittr)
library(data.table)

all <- c(5,seq(10,100,5), seq(110,250,10))
lodf <- lapply(all, function(bgn){
  df <- readRDS(paste0("variableBG/variableBGpeak_", as.character(bgn), ".rds"))
  df$BackgroundPeaksN <- bgn
  df
}) %>% rbindlist() %>% as.data.frame()


ggplot(lodf, aes(x = BackgroundPeaksN, y = -log10(pnorm(value,lower.tail = FALSE)))) +
  geom_line(aes(color = Var1)) + geom_point(aes(color = Var1)) +
  pretty_plot() + labs(x = "Number of ", y = "g-chromVAR Enrichment (-log10 p)", color = "") +
  facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") + scale_color_manual(values = ejc_color_maps) +
  geom_hline(yintercept = -log10(0.05 / (16*18)), linetype = 2) 
