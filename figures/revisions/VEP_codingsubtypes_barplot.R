library(tidyverse)
library(annotables)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(plotly)
library(BuenColors)
library(cowplot)
"%ni%" <- Negate("%in%")

setwd("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits/")

# Load Finemap CS ---------------------------------------------------------

CS.gr <- readRDS("data/Finemap/UKBB_BC_v3_VEPannotations.rds")

# Read VEP table ----------------------------------------------------------
vep <- data.frame(fread("figures/revisions/VEP/vep_output_everything.txt"))
colnames(vep)[1] <- "Uploaded_variation"

vep$Uploaded_variation <- gsub("/","_",vep$Uploaded_variation) %>% str_replace(.,"_",":") 

# Coding consequences
coding_consequences <- c("missense_variant","synonymous_variant","frameshift_variant",
                         "splice_acceptor_variant","splice_donor_variant","splice_region_variant",
                         "inframe_insertion","stop_gained","stop_retained_variant",
                         "start_lost","stop_lost","coding_sequence_variant","incomplete_terminal_codon_variant")

all_consequences <- CS.gr[CS.gr$PP > 0.1,"Consequence"] %>% unique() %>% mcols() %>% unlist()
coding_con <- all_consequences[all_consequences %in% coding_consequences] %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(coding_con)[1] <- "Categorie"
coding_con$Categorie <- factor(coding_con$Categorie, levels = coding_con$Categorie)

# Barplot
p <- ggplot(coding_con, aes(x = Categorie, y = Freq)) + 
  geom_bar(stat = "identity", color = "black", fill = "firebrick") + pretty_plot() +
  L_border() +
  labs(x = "", y = "Frequency") +
  geom_text(data=coding_con,aes(label=Freq),vjust=-0.5, size=3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cowplot::ggsave(p, file = "figures/revisions/plots/coding_subtypes_PP10.pdf", height = 5, width = 8)

