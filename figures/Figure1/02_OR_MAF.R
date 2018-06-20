library(dplyr)
library(BuenColors)
library(ggbeeswarm)

CS.df <- readRDS("allfinemap.0.001snps.correctorder.rds")
CS.df$seqnames <- paste0("chr",CS.df$CHR)
CS.df$MAF_fix <- ifelse(CS.df$MAF > 0.5, 1-CS.df$MAF, CS.df$MAF)
bins = c(0.001,0.01,0.05,0.1,0.25,0.75,1)
CS.df$PP_bin <- cut(x = CS.df$PP, bins)

CS.df %>% group_by(PP_bin) %>% dplyr::select(seqnames, POS, MAF_fix) %>% unique() %>% summarise(MEAN = mean(MAF_fix), SD = sd(MAF_fix))

ggplot(CS.df) + 
  geom_bar(stat="identity",aes(x=PP_bin, y = MAF_fix, fill=PP_bin),
                               position = position_dodge(0.9), fun.y = "mean") +
  pretty_plot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5]))


