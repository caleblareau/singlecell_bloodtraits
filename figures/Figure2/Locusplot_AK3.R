library(data.table)
library(BuenColors)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggrepel)
library(gglabeller)
library(ggplot2)
library(GGally)
library(dplyr)
library(stringr)

sz <- 1.5
width=3
locustheme <-  theme(plot.title = element_text(size=sz*4,hjust = 0.50,face="bold"),
                     text=element_text(size=sz*4),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     legend.position = c(.05, .9),
                     legend.justification = c("left", "top"),
                     legend.box.just = "left",
                     legend.margin = margin(1, 1, 1, 1),
                     legend.direction = "horizontal",
                     legend.key.size = unit(width/30, "in"),
                     legend.text = element_text(size=sz*2),
                     legend.title = element_text(face="bold",size=sz*2))

###################################################################################################
# AK3 standard locus plot
rsid <- "rs12005199"
trait <- "PLT_COUNT"

locuszoom <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/",rsid,"/final_locuslist.",trait,".txt"))

# Tag the two causal variants identified from finemapping
locuszoom$sentinel <- ifelse(locuszoom$SNP == "rs12005199"| 
                               locuszoom$SNP =="rs409950" |
                               locuszoom$SNP =="", "yes", "no")

p1 <- ggplot(locuszoom,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(fill=RSQR),shape=21,size=sz)  +
  pretty_plot()+
  scale_y_continuous(expand = c(0.05, 0))+
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1],name="R2") +
  guides(fill=guide_colorbar(title.vjust=0.75))+
  locustheme+
  labs(x="Position on Chromosome 6 (Mb)",y="-log10(P-value)") + 
  geom_point(data=subset(locuszoom,sentinel=="yes"),
             aes(x=POS/(10^6),y=-log10(PVAL)),
             fill="yellow",shape=21,size=sz)+
  geom_text_repel(data = subset(locuszoom, sentinel=="yes"),
    aes(label = paste(SNP,RSQR,sep=" - ")),
    size = sz,
    force=TRUE,
    nudge_x = 0.75) 

dir="/AK3_PLT_COUNT/"
ggsave(paste0(rsid,"_locusplot.pdf"), plot = p1, device = NULL, path = dir,
       scale = 1, width = width, height=width*1/2, units = "in",
       dpi = 300, limitsize = TRUE)

###################################################################################################
####Conditional Locus plot (conditioning on rs12005199)
rsid <- "rs12005199"
trait <- "PLT_COUNT"
conditional <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/",rsid,"/conditional_sumstats.",rsid,".",trait,".txt")) %>% arrange(PVAL)
rsqre <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/",rsid,"/locuszoom_ld.",trait,".txt"))
conditional <- merge(conditional, rsqre, by.x="SNP", by.y="snp1")
conditional$sentinel <- ifelse(conditional$rsquare ==1 | 
                                 conditional$SNP =="rs409950", "yes", "no")

conditional<-ggplot(conditional,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(fill=RSQR),shape=21,size=sz)  +
  pretty_plot()+
  scale_y_continuous(expand = c(0.05, 0))+
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1],name="R2") +
  guides(fill=guide_colorbar(title.vjust=0.75))+
  locustheme+
  labs(x="Position on Chromosome 9 (Mb)",y="-log10(P-value)") + 
  geom_point(data=subset(conditional,sentinel=="yes"),
             aes(x=POS/(10^6),y=-log10(PVAL)),
             fill="yellow",shape=21,size=sz)+
  geom_label_repel(data = subset(conditional, sentinel=="yes"),
                   aes(label = paste(SNP,RSQR,sep=" - ")),
                   size = sz,
                   force=TRUE,
                   nudge_x = 0.75,
                   nudge_y=0) 

ggsave("rs12005199_conditional_locusplot.pdf", plot = conditional, device = NULL, path = dir,
       scale = 1, width = width, height=width*1/2, units = "in",
       dpi = 300, limitsize = TRUE)

####################################################################################################################
### Plot log10(BF) for AK3 FM variants
rsid <- "rs12005199"
trait <- "PLT_COUNT"

# Using revised FM allowing for 10 variants
region<-1
FM_region <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/wrapper_test/",trait,".region",region,".snp"))

locuszoom <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/",rsid,"/final_locuslist.",trait,".txt"))
FM_region$POS <- str_split_fixed(str_split_fixed(FM_region$snp,":",2)[,2],"_",2)[,1]
FM_region$POS <- as.integer(as.character(FM_region$POS))
FM_region <- merge(locuszoom,FM_region,by="POS")
FM_region <- subset(FM_region,snp_log10bf > -Inf)
FM_region$sentinel <- ifelse(FM_region$SNP =="rs409950" | 
                               FM_region$SNP =="rs12005199", "yes", "no")

fm <- ggplot(subset(FM_region,snp_log10bf > -Inf),aes(POS/(10^6),snp_log10bf)) + 
  geom_point(aes(fill=RSQR),shape=21,size=sz)  +
  pretty_plot()+
  scale_y_continuous(expand = c(0.05, 0))+
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1],name="R2") +
  guides(fill=guide_colorbar(title.vjust=0.75))+
  locustheme+
  labs(x="Position on Chromosome 9 (Mb)",y="log10(Bayes factor)") + 
  geom_point(data=subset(FM_region,sentinel=="yes"),
             aes(x=POS/(10^6),y=snp_log10bf),
             fill="yellow",shape=21,size=sz)+
  geom_label_repel(data = subset(FM_region, sentinel=="yes"),
                   aes(label = paste(SNP,RSQR,sep=" - ")),
                   size = sz,
                   force=TRUE,
                   nudge_x = 0.75,
                   nudge_y=0) 

dir="/AK3_PLT_COUNT/"
width=3
ggsave("AK3_fm_log10bf_10causalvariants.pdf", plot = fm, device = NULL, path = dir,
       scale = 1, width = width, height=width*1/2, units = "in",
       dpi = 300, limitsize = TRUE)

####  Plot log10(BF) for AK3 FM variants using original FM allowing for 5 variants
region<-130
FM_region <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/",trait,"/output/region",region,".snp"))

FM_region$POS <- str_split_fixed(str_split_fixed(FM_region$snp,":",2)[,2],"_",2)[,1]
FM_region$POS <- as.integer(as.character(FM_region$POS))
FM_region <- merge(locuszoom,FM_region,by="POS")
FM_region <- subset(FM_region,snp_log10bf > -Inf)
FM_region$sentinel <- ifelse(FM_region$SNP =="rs409950" | 
                               FM_region$SNP =="rs12005199", "yes", "no")

fm <- ggplot(subset(FM_region,snp_log10bf > -Inf),aes(POS/(10^6),snp_log10bf)) + 
  geom_point(aes(fill=RSQR),shape=21,size=sz)  +
  pretty_plot()+
  scale_y_continuous(expand = c(0.05, 0))+
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1],name="R2") +
  guides(fill=guide_colorbar(title.vjust=0.75))+
  locustheme+
  labs(x="Position on Chromosome 9 (Mb)",y="log10(Bayes factor)") + 
  geom_point(data=subset(FM_region,sentinel=="yes"),
             aes(x=POS/(10^6),y=snp_log10bf),
             fill="yellow",shape=21,size=sz)+
  geom_label_repel(data = subset(FM_region, sentinel=="yes"),
                   aes(label = paste(SNP,RSQR,sep=" - ")),
                   size = sz,
                   force=TRUE,
                   nudge_x = 0.75,
                   nudge_y=0) 

width=3
ggsave("AK3_fm_log10bf_5causalvariants.pdf", plot = fm, device = NULL, path = dir,
       scale = 1, width = width, height=width*1/2, units = "in",
       dpi = 300, limitsize = TRUE)
