library(dplyr)

df <- readRDS("allEnrichments-df.rds")

permuted <- sapply(1:10000, function(i) sum(1:dim(df)[1] * sample(df$lineageSpecific, length(df$lineageSpecific))))
ldscore_ranksum <- sum(1:dim(df)[1]*df[order(df$ldscore_pvalue, decreasing = FALSE), "lineageSpecific"])
chromVAR_ranksum <- sum(1:dim(df)[1]*df[order(df$chromVAR_pvalue, decreasing = FALSE), "lineageSpecific"])
gchromVAR_ranksum <- sum(1:dim(df)[1]*df[order(df$gchromVAR_pvalue, decreasing = FALSE), "lineageSpecific"])
goShifter_ranksum <- sum(1:dim(df)[1]*df[order(df$goShifter_pvalue, decreasing = FALSE), "lineageSpecific"])
panHemeLDSR_ranksum <- sum(1:dim(df)[1]*df[order(df$panHemeLDSR_pvalue, decreasing = FALSE), "lineageSpecific"])

gregor_ranksum <- sum(1:dim(df)[1]*df[order(df$gregor_pvalue, decreasing = FALSE), "lineageSpecific"])
GPA_ranksum <- sum(1:dim(df)[1]*df[order(df$GPA_chisq, decreasing = TRUE), "lineageSpecific"])
fgwas_ranksum <- sum(1:dim(df)[1]*df[order(df$fgwas_z, decreasing = TRUE), "lineageSpecific"])



allstats <- data.frame(
  pval = c(pnorm((mean(permuted) - ldscore_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - chromVAR_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - gchromVAR_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - goShifter_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - panHemeLDSR_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - gregor_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - GPA_ranksum)/sd(permuted), lower.tail = FALSE),
           pnorm((mean(permuted) - fgwas_ranksum)/sd(permuted), lower.tail = FALSE)
  ))
allstats$method <- c("LDscore", "chromVAR", "gchromVAR", "goShifter", "adjLDscore", "GREGOR", "GPA", "fgwas")
allstats$log10p <- -1*log10(allstats$pval)

plot_df <- data.frame(val = dim(df)[1]:1)
plot_df$gchromVAR <- df %>% arrange(gchromVAR_pvalue) %>% pull(lineageSpecific)
plot_df$chromVAR <- df %>% arrange(chromVAR_pvalue) %>% pull(lineageSpecific)
plot_df$LDscore <- df %>% arrange(ldscore_pvalue) %>% pull(lineageSpecific)
plot_df$adjLDscore <- df %>% arrange(panHemeLDSR_pvalue) %>% pull(lineageSpecific)
plot_df$goShifter <- df %>% arrange(goShifter_pvalue) %>% pull(lineageSpecific)
plot_df$GPA <- df %>% arrange(desc(GPA_chisq)) %>% pull(lineageSpecific)
plot_df$GREGOR <- df %>% arrange(gregor_pvalue) %>% pull(lineageSpecific)
plot_df$fgwas <- df %>% arrange(desc(fgwas_z)) %>% pull(lineageSpecific)

melt_plot_df <- reshape2::melt(plot_df, id.var = "val")

order <- c("gchromVAR", "fgwas", "GREGOR", "adjLDscore", "LDscore", "goShifter", "GPA", "chromVAR")
allstats$method <- factor(as.character(allstats$method), levels =  order)

p1 <- ggplot(allstats, aes(x = method, y = log10p)) + 
  geom_bar(stat = "identity", fill = "grey", color = "black") + pretty_plot() + L_border() +
  labs(x = "", y = "-log10p") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

melt_plot_df$variable <- factor(as.character(melt_plot_df$variable), levels = order)
p2 <- ggplot(melt_plot_df, aes(x = variable, y = val, fill = value)) + geom_tile() +
  labs(x = "", y = "Relative Enrichment", fill = "Lineage Specific") +
  pretty_plot() + L_border() + scale_fill_manual(values = c("grey", "dodgerblue")) +
  theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(panel.grid = element_blank(), panel.border = element_blank())
  

cowplot::ggsave(cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1,1.5)),
                height = 5, width = 3, filename = "lineageEnrichments.pdf") 



