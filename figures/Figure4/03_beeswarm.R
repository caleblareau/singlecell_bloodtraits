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

#### RBC
trait <- "RBC_COUNT"
kct <- c("HSC", "MPP", "CMP", "MEP")
df_long <- observed_long[observed_long$Var1 ==trait  & observed_long$Var2 %in% kct,]
df_long$Var2 <- as.ordered(factor(df_long$Var2, levels = kct))
df_long$value <- ifelse(df_long$value > 2.5, 2.5, df_long$value)

rbc_plot <- ggplot(df_long, aes(x = Var2, y = value, color = Var2, group = Var2)) +
  geom_quasirandom(dodge.width=1, alpha = 0.6) + pretty_plot() +
  stat_summary(fun.y=mean, geom="point", shape=95, size=20, color = "black") +
  scale_color_manual(values = ejc_color_maps) + theme(legend.position = "bottom") +
  labs(x = "", y = paste0("Red Blood Cell Count Enrichment"), color = "")

ggsave(rbc_plot, filename = "rawBeePlot/RBC.pdf")

#### Lymph
trait <- "LYMPH_COUNT"
kct <- c("HSC", "MPP", "LMPP", "CLP", "pDC")
lymph_long <- observed_long[observed_long$Var1 ==trait  & observed_long$Var2 %in% kct,]
lymph_long$Var2 <- as.ordered(factor(lymph_long$Var2, levels = kct))
lymph_long$value <- ifelse(lymph_long$value > 2.5, 2.5, lymph_long$value)

lymph_plot <- ggplot(lymph_long, aes(x = Var2, y = value, color = Var2, group = Var2)) +
  geom_quasirandom(dodge.width=1, alpha = 0.6) + pretty_plot() + 
  stat_summary(fun.y=mean, geom="point", shape=95, size=20, color = "black") +
  scale_color_manual(values = ejc_color_maps) + theme(legend.position = "bottom") +
  labs(x = "", y = paste0("Lymph Count Enrichment"), color = "")

ggsave(lymph_plot, filename = "rawBeePlot/LYMPH.pdf")

#### MONO
# Have to do some extra stuff here to get
# GMP variable colors but grouped by basic GMP

# Differential analysis
GMPA <- observed_long[observed_long$Var1 == "MONO_COUNT" & observed_long$Var2 == "GMP-A", "value"]
GMPB <- observed_long[observed_long$Var1 == "MONO_COUNT" & observed_long$Var2 == "GMP-B", "value"]
GMPC <- observed_long[observed_long$Var1 == "MONO_COUNT" & observed_long$Var2 == "GMP-C", "value"]
summary(GMPA)
summary(GMPB)
summary(GMPC)

t.test(GMPA, GMPC)

kct <- c("HSC", "MPP", "CMP", "GMP-A", "GMP-B", "GMP-C", "GMP", "Mono")
mono_long <- observed_long[observed_long$Var1 == "MONO_COUNT" & observed_long$Var2 %in% kct,]
mono_long$Var2 <- as.ordered(factor(mono_long$Var2, levels = kct))
mono_long$celltype <- ifelse(as.character(mono_long$Var2) %in% c("GMP-A", "GMP-B", "GMP-C", "GMP"), "GMP", as.character(mono_long$Var2))
mono_long$celltype <- as.ordered(factor(mono_long$celltype, levels = c("HSC", "MPP", "CMP", "GMP", "Mono")))
mono_long$value <- ifelse(mono_long$value > 2.5, 2.5, mono_long$value)

mono_plot <- ggplot(mono_long, aes(x = celltype, y = value, color = Var2, group = celltype)) +
  geom_quasirandom(dodge.width=1, alpha = 0.6) + pretty_plot() + 
  stat_summary(fun.y=mean, geom="point", shape=95, size=20, color = "black") +
  scale_color_manual(values = ejc_color_maps) + theme(legend.position = "bottom") +
  labs(x = "", y = "Monocyte Count Enrichment", color = "") +
  scale_y_continuous(limits = c(-2.5, 2.5))
ggsave(mono_plot, filename = "rawBeePlot/MONO.pdf")
