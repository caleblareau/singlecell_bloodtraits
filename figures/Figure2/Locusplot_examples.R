library(data.table)
library(BuenColors)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Gviz)
library(ggrepel)
library(gglabeller)
library(ggplot2)
library(GGally)
library(dplyr)

###################################################################################################
# CCND3
rsid <- "rs9349205"
locuszoom <- fread("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/rs9349205/final_locuslist.RBC_COUNT.txt")
trait <- "RBC_COUNT"

# Tag the two causal variants identified from finemapping
locuszoom$sentinel <- ifelse(locuszoom$RSQR ==1 | 
                               locuszoom$SNP =="rs112233623", "yes", "no")

p1 <- ggplot(locuszoom,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(colour=RSQR),size=3) + pretty_plot() +
  ggtitle(paste(rsid,trait)) + 
  scale_color_gradientn(colors = jdb_palette("brewer_heat")[-1]) +
  theme(plot.title = element_text(size=14,hjust = 0.5,face="bold")) +
  labs(x="Position on Chromosome 4 (Mb)") + 
  geom_text_repel(
    data = subset(locuszoom, sentinel=="yes"),
    aes(label = paste(SNP,RSQR,sep=" - ")),
    size = 3,
    force=TRUE,
    nudge_x = 0.15) 

#### Conditional Locus plot
rsid <- "rs112233623"
trait <- "RBC_COUNT"
conditional <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/rs9349205/conditional_sumstats.rs9349205.",trait,".txt")) %>% arrange(PVAL)
rsqre <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/rs112233623/locuszoom_ld.RBC_COUNT.txt"))
conditional <- merge(conditional, rsqre, by.x="SNP", by.y="snp1")
conditional$sentinel <- ifelse(conditional$rsquare ==1 | 
                                 conditional$SNP =="rs9349205", "yes", "no")

p2 <- ggplot(conditional,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(colour=rsquare),size=3) + pretty_plot() +
  ggtitle(paste(rsid,trait)) + 
  scale_color_gradientn(colors = jdb_palette("brewer_heat")[-1]) +
  theme(plot.title = element_text(size=14,hjust = 0.5,face="bold")) +
  labs(x="Position on Chromosome 4 (Mb)") + 
  geom_text_repel(
    data = subset(conditional, sentinel=="yes"),
    aes(label = paste(SNP,rsquare,sep=" - ")),
    size = 3,
    force=TRUE,
    nudge_x = 0.15) 

ggmatrix(list(p1,p2),nrow=1,ncol=2,
         xAxisLabels = c("CCND3 Locus", "CCND3 Locus Conditioned on rs9349205"),
         legend=c(1,1)) +
  theme(plot.title = element_text(hjust = 0.5))

# Hard called LD finemap results
locuszoom$sentinel <- ifelse( locuszoom$SNP =="rs72867133" | 
                                locuszoom$SNP =="rs11440956"| 
                                locuszoom$SNP =="rs202170071"| 
                                locuszoom$SNP =="rs141797580"| 
                                locuszoom$SNP =="rs80215763", "yes", "no")

# Replot, labeling the incorrectly predicted hard called variants
ggplot(locuszoom,aes(POS/(10^6),-log10(PVAL))) + 
  geom_point(aes(colour=RSQR),size=3) + pretty_plot() +
  ggtitle(paste(rsid,trait)) + 
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  theme(plot.title = element_text(size=14,hjust = 0.5,face="bold")) +
  labs(x="Position on Chromosome 4 (Mb)") + 
  geom_text_repel(
    data = subset(locuszoom, sentinel=="yes"),
    aes(label = paste(SNP,RSQR,sep=" - ")),
    size = 3,
    force=TRUE,
    nudge_x = 0.15) 