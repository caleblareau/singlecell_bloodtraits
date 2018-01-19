library(corrplot)
library(BuenColors)
library(dplyr)
library(reshape2)
library(cowplot)
library(data.table)
library(plotly)
library(stringr)
library(ggrepel)
library(scales)

# Narrow-Sense Heritability Estimated Obtained from LDSC
h2 <-fread("../../data/LDSC_heritability/trait_h2.txt")
colnames(h2) <- c("trait","h_obs","h_obs_se")

# Read in enrichments for FM PP001 variants
reorderedtraits=traits=c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

ctHeme <- vector("list",length=16)
ct_enrichment <- vector("list",length=length(ctHeme))
celltypes <- paste0(traits,"_PP001L2")

FM_h2 <- NULL
FM_enrichment <- NULL
dir="../../data/LDSC_heritability/FM_PP001/"
for (i in 1:length(ctHeme)){
  ctHeme[[i]]<-fread(paste0(dir,traits[i],".",celltypes[i],".results"))
  FM_h2[i] <- ctHeme[[i]]$Prop._h2[1]
  FM_enrichment[i] <- ctHeme[[i]]$Enrichment[1]
}

h2$fm_herit <- FM_h2*h2$h_obs
h2$fm_h2 <- FM_h2
h2$trait <- factor(as.character(h2$trait), rev(unique(as.character(h2$trait))))

legendtitle <- "fine-mapped heritability"
p0 <- ggplot(data=h2, aes(x=trait,y=fm_herit)) + 
  geom_bar(aes(y=h_obs), fill="light gray",stat="identity",position = position_stack(reverse = TRUE)) +
  geom_bar(aes(y=fm_herit,fill=legendtitle),stat="identity") +
  geom_text(aes(label=percent(FM_h2)), hjust=-0.25, color="black", size=3.25) +
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5])) +
  ggtitle("Finemapped Trait Heritabilities") + theme_bw() + 
  coord_flip() + 
  labs(fill="") + 
  scale_y_continuous(expand = c(0.01, 0.0)) +
  theme(plot.subtitle = element_text(vjust = 1),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(.98, 0.99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical",
        legend.key.size = unit(0.15, "in"))  +
  geom_errorbar(aes(ymin=h_obs-h_obs_se, ymax=h_obs+h_obs_se),width=.2,position=position_dodge(.9),color="black") +
  ylab("Narrow Sense Heritability") + xlab("")

phenocors <-read.table("../../data/phenotypeCorrelations/raw_phenotypes.txt",header=TRUE,row.names=1)
dissimilarity <- 1 - cor(data.matrix(phenocors))
distance <- as.dist(dissimilarity) 
reorderedtraits=traits=rownames(phenocors)[hclust(distance)$order]

phenocors <- data.matrix(phenocors[reorderedtraits,reorderedtraits])
melt_pheno <- reshape2::melt(phenocors)
melt_pheno$Var1 <- ordered(as.character(melt_pheno$Var1), ordered(traits))
melt_pheno$Var2 <- ordered(as.character(melt_pheno$Var2), ordered(rev((traits))))

p1 <- ggplot(data = melt_pheno, aes(x=Var1, y=Var2, fill=value)) + pretty_plot() + 
  geom_tile() + scale_fill_gradientn(colors = jdb_palette("solar_flare"), limits = c(-1,1)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(angle = 90)) +labs(x = NULL, y = NULL, fill = "Pearson") +
  theme(legend.position="bottom")  + ggtitle("Phenotype Correlation")

### Constrained intercept genetic correlations
constrained_gencors <- vector(mode = "list", length = length(traits))

# Read unix files into R lists
constrained_gencors <- lapply(traits, function(y) read.table(paste0("../../data/phenotypeCorrelations/ldscore/",y,"_UK10K_constrained.gcsummary.txt"),header=TRUE,row.names=NULL))

make.gc.matrix <- function(list_of_gcs,traits){
  allgencors <- bind_rows(list_of_gcs)
  allgencors$rg <- as.numeric(as.character(allgencors$rg))
  
  gc.matrix <- matrix(ncol = 16, nrow = 16, dimnames=list(traits,traits))
  pv.matrix <- matrix(ncol = 16, nrow = 16, dimnames=list(traits,traits))
  
  for (i in traits) {
    for (j in traits) {
      if (i == j){
        gc.matrix[i,j]=1
      } else {
        gc.matrix[i,j] <- allgencors[allgencors$p1 == i & allgencors$p2 == j,'rg']
        if (gc.matrix[i,j]>1) gc.matrix[i,j]=1
        pv.matrix[i,j]<- allgencors[allgencors$p1 == i & allgencors$p2 == j,'p']
      }
    }
  }
  return(list(gc.matrix,pv.matrix))
}

# Plot LD Score correlation
ldscorecor <- make.gc.matrix(constrained_gencors,traits)[[1]]
melt_phenoLD <- reshape2::melt(ldscorecor)
melt_phenoLD$Var1 <- ordered(as.character(melt_phenoLD$Var1), ordered(traits))
melt_phenoLD$Var2 <- ordered(as.character(melt_phenoLD$Var2), ordered(rev((traits))))

p2 <- ggplot(data = melt_phenoLD, aes(x=Var1, y=Var2, fill=value)) + pretty_plot() + 
  geom_tile() + scale_fill_gradientn(colors = jdb_palette("solar_flare"), limits = c(-1,1)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(angle = 90)) +labs(x = NULL, y = NULL, fill = "LDScore Correlation  ") +
  theme(legend.position="bottom")  + ggtitle("Genetic Correlation  ")


bottom_row <- plot_grid(p1, p2, labels = c('b', 'c'), align = 'h', rel_widths = c(1, 1))
final_plot <- plot_grid(p0, bottom_row, labels = c('a', ''), rel_heights = c(1,1.35), ncol = 1)

cowplot::ggsave(final_plot, filename = "PDFs/FigureS1.pdf", width = 9, height = 10)

