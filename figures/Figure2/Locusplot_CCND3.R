library(data.table)
library(BuenColors)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Gviz)
library(ggrepel)
library(gglabeller)
library(ggplot2)
library(GGally)
library(dplyr)
library(stringr)

###################################################################################################
# CCND3
rsid <- "rs9349205"
trait <- "RBC_COUNT"
locuszoom <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/rs9349205/final_locuslist.",trait,".txt"))

# Tag the two causal variants identified from finemapping
locuszoom$sentinel <- ifelse(locuszoom$RSQR ==1 | 
                               locuszoom$SNP =="rs112233623", "yes", "no")
mycolours <- c("yes" = "yellow", "no" = "grey50")

ggplot(locuszoom,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(fill=RSQR),shape=21,size=3) + pretty_plot() +
  ggtitle(paste(rsid,trait)) + 
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1]) +
  theme(plot.title = element_text(size=14,hjust = 0.50,face="bold")) +
  labs(x="Position on Chromosome 6 (Mb)",y="-log10(P-value)") + 
  geom_point(data=subset(locuszoom,sentinel=="yes"),
             aes(x=POS/(10^6),y=-log10(PVAL)),
             fill="yellow",shape=21,size=3)+
  geom_label_repel(data = subset(locuszoom, sentinel=="yes"),
    aes(label = paste(SNP,RSQR,sep=" - ")),
    size = 3,
    force=TRUE,
    nudge_x = 0.5) 

#### Conditional Locus plot (conditioning on rs9349205)
rsid <- "rs112233623"
trait <- "RBC_COUNT"
conditional <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/rs9349205/conditional_sumstats.rs9349205.",trait,".txt")) %>% arrange(PVAL)
rsqre <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/rs112233623/locuszoom_ld.",trait,".txt"))
conditional <- merge(conditional, rsqre, by.x="SNP", by.y="snp1")
conditional$sentinel <- ifelse(conditional$rsquare ==1 | 
                                 conditional$SNP =="rs9349205", "yes", "no")

ggplot(conditional,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(fill=rsquare),shape=21,size=3) + pretty_plot() +
  ggtitle(paste(rsid,trait)) + 
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1]) +
  theme(plot.title = element_text(size=14,hjust = 0.50,face="bold")) +
  labs(x="Position on Chromosome 6 (Mb)",y="-log10(P-value)") + 
  geom_point(data=subset(conditional,sentinel=="yes"),
             aes(x=POS/(10^6),y=-log10(PVAL)),
             fill="yellow",shape=21,size=3) +
  geom_label_repel(
    data = subset(conditional, sentinel=="yes"),
    aes(label = paste(SNP,rsquare,sep=" - ")),
    size = 3,
    nudge_x = 0.5) 

# ggmatrix(list(p1,p2),nrow=1,ncol=2,
#          xAxisLabels = c("CCND3 Locus", "CCND3 Locus Conditioned on rs9349205"),
#          legend=c(1,1)) +
#   theme(plot.title = element_text(hjust = 0.5))

### Plot log10(BF) for FM variants
region<-37
FM_region <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/",trait,"/output/region",region,".snp"))
FM_region$POS <- str_split_fixed(str_split_fixed(FM_region$snp,":",2)[,2],"_",2)[,1]
FM_region$POS <- as.integer(as.character(FM_region$POS))
FM_region <- merge(locuszoom,FM_region,by="POS")
FM_region <- subset(FM_region,snp_log10bf > -Inf)
FM_region$sentinel <- ifelse(FM_region$RSQR ==1 | 
                               FM_region$SNP =="rs112233623", "yes", "no")

ggplot(subset(FM_region,snp_log10bf > -Inf),aes(POS/(10^6),snp_log10bf)) + 
  geom_point(aes(fill=RSQR),shape=21,size=3) + pretty_plot() +
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1]) +
  theme(plot.title = element_text(size=14,hjust = 0.50,face="bold")) +
  labs(x="Position on Chromosome 6 (Mb)", y="log10(Bayes factor)")+ 
  geom_point(data=subset(FM_region,sentinel=="yes"),
             aes(x=POS/(10^6),y=snp_log10bf),
             fill="yellow",shape=21,size=3) +
  geom_label_repel(data = subset(FM_region, sentinel=="yes"),
                   aes(label = paste(SNP)),
                   size = 3,
                   force=TRUE,
                   nudge_x = 0.5)

# Hard called LD finemap results
hardcallregion<-36
FM_hardcall_region <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/",trait,"/hardcalloutput/region",hardcallregion,".snp"))
FM_hardcall_region$POS <- str_split_fixed(str_split_fixed(FM_hardcall_region$snp,":",2)[,2],"_",2)[,1]
FM_hardcall_region$POS <- as.integer(as.character(FM_hardcall_region$POS))
FM_hardcall_region <- merge(locuszoom,FM_hardcall_region,by="POS")
FM_hardcall_region[FM_hardcall_region$snp_log10bf < 0,"snp_log10bf"] <- 0 
#FM_hardcall_region <- subset(FM_hardcall_region,snp_log10bf > -Inf)
FM_hardcall_region$sentinel <- ifelse(FM_hardcall_region$RSQR ==1 | 
                                        FM_hardcall_region$SNP =="rs112233623", "yes", "no")

# Plot FM log10bf with hard called variants labeled
ggplot(subset(FM_hardcall_region,snp_log10bf > -Inf),aes(POS/(10^6),snp_log10bf)) + 
  geom_point(aes(fill=RSQR),shape=21,size=3) + pretty_plot() +
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1]) +
  theme(plot.title = element_text(size=14,hjust = 0.50,face="bold")) +
  labs(x="Position on Chromosome 6 (Mb)", y="log10(Bayes factor)") + 
  geom_point(data=subset(FM_hardcall_region,sentinel=="yes"),
             aes(x=POS/(10^6),y=snp_log10bf),
             fill="yellow",shape=21,size=3) +
  geom_label_repel(data = subset(FM_hardcall_region, sentinel=="yes"),
                  aes(label = paste(SNP)),
                  size = 3,
                  force=TRUE,
                  nudge_x = 0.5,
                  nudge_y=2.5) 

# Hard called LD finemap results
locuszoom$hardcall <- ifelse( locuszoom$SNP =="rs72867133" | 
                                locuszoom$SNP =="rs11440956"| 
                                locuszoom$SNP =="rs202170071"| 
                                locuszoom$SNP =="rs141797580"| 
                                locuszoom$SNP =="rs80215763", "yes", "no")

# Replot, labeling the incorrectly predicted hard called variants
ggplot(locuszoom,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(fill=RSQR),shape=21,size=3) + pretty_plot() +
  ggtitle(paste("CCND3 hard called variants",trait)) + 
  scale_fill_gradientn(colors = jdb_palette("solar_extra")) +
  theme(plot.title = element_text(size=14,hjust = 0.5,face="bold")) +
  labs(x="Position on Chromosome 6 (Mb)") + 
  geom_label_repel(
    data = subset(locuszoom, hardcall=="yes"),
    aes(label = paste(SNP,RSQR,sep=" - ")),
    size = 3,
    force=TRUE,
    nudge_x = 0.5)
