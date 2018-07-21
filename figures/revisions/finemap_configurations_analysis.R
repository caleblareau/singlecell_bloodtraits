library(tidyverse)
library(annotables)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(BuenColors)
library(cowplot)
"%ni%" <- Negate("%in%")

setwd("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits/")
# Load in top configs -----------------------------------------------------
traits=c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

top_configs <-lapply(traits, function(t) {
  df <- fread(paste0(t,"/output/topconfigs.100configs.",t,".txt"))
  df$trait <- t
  return(df)
})
all_configs <- bind_rows(top_configs) 

# Exclude all configurations that contain a variant in excluded list
exdf <- read.table("figures/revisions/exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

filtered_all_configs <- all_configs
for (i in 1:length(exdf)){
  filtered_all_configs <- filtered_all_configs[!grepl(exdf[i],filtered_all_configs$config),]
}

# Rerank the configurations after excluding configs
filtered_all_configs <- filtered_all_configs %>% arrange(region,trait, desc(config_prob)) %>%
  group_by(region,trait) %>%
  mutate(rank=rank(-config_prob,ties.method="first")) %>% as.data.frame()

# Stacked bar plots -------------------------------------------------------
plot_summed_configs <- function(num_configs=1,all_configs){
  all_configs_sum <- all_configs %>% subset(rank<=num_configs) %>% 
    group_by(region,trait) %>%
    dplyr::summarize(PPsum=sum(config_prob)) %>%
    as.data.frame()
  
  all_configs_sum$config_PP <- cut(all_configs_sum$PPsum, c(0,0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))
  allconfigs_toplot <- all_configs_sum  %>% 
    group_by(trait,config_PP) %>% dplyr::summarize(count=n())
  
  if (allconfigs_toplot$config_PP %>% unique %>% length == 7){
    colors <- c(jdb_palette("GrandBudapest2")[c(3:4)],jdb_palette("Zissou")[1:5])
  } else{
    colors <- c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5])
  }
  
  p <- ggplot(allconfigs_toplot,aes(y=count,x=trait,group=trait)) + 
    geom_bar(stat="identity",aes(fill=config_PP),position = position_stack(reverse = TRUE)) + 
    coord_flip() +
    theme_bw() +
    pretty_plot() +
    scale_fill_manual(values = colors)+
    labs(x="",y="number of regions")
  # cowplot::ggsave(p, file = paste0("figures/revisions/topconfigs_PPs/topconfig_PPs.",num_configs,".pdf"), 
  #                 height = 4, width = 7)
  return(p)
}

# Plot stacked barplot for top n configs
#setwd("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits")
plot_summed_configs(num_configs = 1,all_configs = filtered_all_configs)

# Figure out percentage of top N configs in each bin ----------------------
all_configs_sum <- filtered_all_configs %>% subset(rank<=1) %>% 
  group_by(region,trait) %>%
  dplyr::summarize(PPsum=sum(config_prob)) %>%
  as.data.frame()

all_configs_sum$config_PP <- cut(all_configs_sum$PPsum, c(0,0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))
all_configs_sum  %>% group_by(config_PP) %>% dplyr::summarize(count=n())

# Stratify by number of causal variants -----------------------------------
allregions <- readRDS("data/Finemap/allFINEMAPregions.rds")

all_configs_sum <- filtered_all_configs %>% subset(rank<=1) %>% 
  group_by(region,trait) %>%
  dplyr::summarize(PPsum=sum(config_prob)) %>%
  as.data.frame()
all_configs_sum$config_PP <- cut(all_configs_sum$PPsum, c(0,0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))

num_causal_merged <- all_configs_sum %>% left_join(.,allregions[,c("region","trait","expectedvalue")],by=c("region","trait"))

colors <- c(jdb_palette("GrandBudapest2")[c(3:4)],jdb_palette("Zissou")[1:5])
p <- ggplot(num_causal_merged,aes(x=config_PP,group=config_PP,y=expectedvalue,fill=config_PP)) + 
  scale_fill_manual(values=colors) +
  geom_boxplot() +
  pretty_plot() + 
  L_border() + 
  labs(x = "Top configuration PP", y="# causal variants")
cowplot::ggsave(p, file = paste0("figures/revisions/topconfigs_PPs/configPP_vs_numcausalvariants_boxplot.pdf"),
                height = 4, width = 7)


# Supplementary Table for Top Configurations ------------------------------

CS.gr <- readRDS("data/Finemap/UKBB_BC_v3_VEPannotations.rds")
CS.df <- as.data.frame(CS.gr)
CS.df <- rename(CS.df,top_individual_variant_PP=PP)
CS.df$region <- gsub("region","",CS.df$region) %>% as.numeric()

CS.df.topvariantonly <- CS.df %>% dplyr::select(var,top_individual_variant_PP,region,trait) %>% 
  group_by(region,trait) %>% filter(top_individual_variant_PP == max(top_individual_variant_PP)) %>%
    filter(1:n() == 1) %>% as.data.frame()

all_top_configs <- filtered_all_configs[filtered_all_configs$rank == 1,]

topPP_merged <- left_join(all_top_configs, 
                          CS.df.topvariantonly[,c("var","trait","region","top_individual_variant_PP")],
                          by=c("trait","region")) %>% dplyr::select(-rank) %>% unique %>% arrange(trait) 

write.table(topPP_merged, file = "/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/Revisions/Supp_Tables/SupplementalTable_configs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)  


