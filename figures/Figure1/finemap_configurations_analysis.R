library(tidyverse)
library(annotables)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(plotly)
library(BuenColors)
"%ni%" <- Negate("%in%")

setwd("/Volumes/broad_sankaranlab/ebao/FINEMAP/")
# Load in top configs -----------------------------------------------------
traits=c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

top_configs <-lapply(traits, function(t) {
  df <- fread(paste0(t,"/output/topconfigs.100configs.",t,".txt"))
  df$trait <- t
  return(df)
})
all_configs <- bind_rows(top_configs) 


# Stacked bar plots -------------------------------------------------------
plot_summed_configs <- function(num_configs=1,all_configs){
  all_configs_sum <- all_configs %>% subset(rank<=num_configs) %>% 
    group_by(region,trait) %>%
    dplyr::summarize(PPsum=sum(config_prob)) %>%
    as.data.frame()
  
  all_configs_sum$PPbin <- cut(all_configs_sum$PPsum, c(0,0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))
  allconfigs_toplot <- all_configs_sum  %>% 
    group_by(trait,PPbin) %>% dplyr::summarize(count=n())
  
  if (allconfigs_toplot$PPbin %>% unique %>% length == 7){
    colors <- c(jdb_palette("GrandBudapest2")[c(3:4)],jdb_palette("Zissou")[1:5])
  } else{
    colors <- c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5])
  }
  
  ggplot(allconfigs_toplot,aes(y=count,x=trait,group=trait)) + 
    geom_bar(stat="identity",aes(fill=PPbin),position = position_stack(reverse = TRUE)) + 
    coord_flip() +
    theme_bw() +
    pretty_plot() +
    scale_fill_manual(values = colors)+
    labs(x="",y="number of regions")
}

# Plot stacked barplot for top n configs
plot_summed_configs(num_configs = 1,all_configs = all_configs)