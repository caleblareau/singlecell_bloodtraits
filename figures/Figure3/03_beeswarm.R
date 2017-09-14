library(BuenColors)
library(ggplot2)
library(ggbeeswarm)

# Import observed 
observed <- read.table("../../data/singlecell/scATAC/weightedSingleCellScores-shuffled.txt", sep = "\t", header = TRUE)
observeddf <- t(observed[,c(13:28)])
colnames(observeddf) <- observed[,2]
rownames(observeddf) <- gsub("raw_", "", rownames(observeddf))
observed_long <- reshape2::melt(observeddf)
observed_long$type <- "Observed"

#### MONO

kct <- c("HSC", "MPP", "CMP", "GMP-A", "GMP-B", "GMP-C", "Mono")
mono_long <- observed_long[observed_long$Var1 == "MONO_COUNT" & observed_long$Var2 %in% kct,]
mono_long$Var2 <- as.ordered(factor(mono_long$Var2, levels = kct))

ggplot(mono_long, aes(x = Var2, y = value, color = Var2, group = Var2)) +
  geom_quasirandom(dodge.width=1, alpha = 0.8) + pretty_plot() + ggtitle("Mono Count") + 
  stat_summary(fun.y=mean, geom="point", shape=95, size=20, color = "black") +
  scale_color_manual(values = ejc_color_maps) + theme(legend.position = "bottom") +
  labs(x = "", y = "Monocyte Count Enrichment", color = "")

#### RBC
trait <- "PLT_COUNT"
kct <- c("HSC", "MPP", "CMP", "MEP")
df_long <- observed_long[observed_long$Var1 ==trait  & observed_long$Var2 %in% kct,]
df_long$Var2 <- as.ordered(factor(df_long$Var2, levels = kct))

ggplot(df_long, aes(x = Var2, y = value, color = Var2, group = Var2)) +
  geom_quasirandom(dodge.width=1, alpha = 0.8) + pretty_plot() + ggtitle(trait) + 
  stat_summary(fun.y=mean, geom="point", shape=95, size=20, color = "black") +
  scale_color_manual(values = ejc_color_maps) + theme(legend.position = "bottom") +
  labs(x = "", y = paste0(trait, " Enrichment"), color = "")

#### Lymph
trait <- "WBC_COUNT"
kct <- c("HSC", "MPP", "LMPP", "CLP", "pDC")
df_long <- observed_long[observed_long$Var1 ==trait  & observed_long$Var2 %in% kct,]
df_long$Var2 <- as.ordered(factor(df_long$Var2, levels = kct))

ggplot(df_long, aes(x = Var2, y = value, color = Var2, group = Var2)) +
  geom_quasirandom(dodge.width=1, alpha = 0.8) + pretty_plot() + ggtitle(trait) + 
  stat_summary(fun.y=mean, geom="point", shape=95, size=20, color = "black") +
  scale_color_manual(values = ejc_color_maps) + theme(legend.position = "bottom") +
  labs(x = "", y = paste0(trait, " Enrichment"), color = "")

