library(BuenColors)
library(magrittr)
library(data.table)

all <- c(5,seq(10,100,5), seq(110,250,10))
makeDF1 <- function(n){
lodf <- lapply(all, function(bgn){
  fn <- paste0("variableBG/variableBGpeak_", as.character(bgn), "_", "iteration_", as.character(n), ".rds")
  df <- readRDS(fn)
  df$BackgroundPeaksN <- bgn
  df$iteration <- n
  df
}) %>% rbindlist() %>% as.data.frame()
return(lodf)
}

alldf <- lapply(1:100, makeDF1) %>% rbindlist %>% as.data.frame
saveRDS(alldf, file = "allIterationsAllPeaks.rds")



ggplot(lodf, aes(x = BackgroundPeaksN, y = -log10(pnorm(value,lower.tail = FALSE)))) +
  geom_line(aes(color = Var1)) + geom_point(aes(color = Var1)) +
  pretty_plot() + labs(x = "Number of ", y = "g-chromVAR Enrichment (-log10 p)", color = "") +
  facet_wrap(~Var2, scales = "free_y") +
  theme(legend.position="bottom") + scale_color_manual(values = ejc_color_maps) +
  geom_hline(yintercept = -log10(0.05 / (16*18)), linetype = 2) 
