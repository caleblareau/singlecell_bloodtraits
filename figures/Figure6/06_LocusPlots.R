library(dplyr)
library(magrittr)
library(data.table)
library(BuenColors)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggrepel)
library(gglabeller)
library(ggplot2)
library(GGally)
library(stringr)

# Read in all fine-mapped variants
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)

# Read in known variants
known <- read.table("../Figure6/known/knownVars.txt", sep="\t", header=F)
names(known) <- c("gene", "rsID", "trait", "chr", "pos", "PP", "region")

sz <- 1.5
width=5
height=3
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

####################################################################################################################
### SUPPLEMENTARY EXAMPLES
####################################################################################################################

makePlots <- function(rsid, trait, region_number) {
  locuszoom <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/association_plots/",rsid,"/final_locuslist.",trait,".txt"))
  
  # Tag the two causal variants identified from finemapping
  locuszoom$sentinel <- ifelse(locuszoom$SNP == rsid, "yes", "no")
  locuszoom$log10 <- -log10(locuszoom$PVAL)
  
  # Fix log10 p-values
  fix <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/MPV/nodup.snps/arghef3_update.txt"))
  vec <- fix[[2]]; names(vec) <- fix[[1]]
  locuszoom$log10 <- ifelse(locuszoom$SNP %in% names(vec), unname(vec[locuszoom$SNP]), locuszoom$log10)

  p1 <- ggplot(locuszoom,aes(POS/(10^6),log10)) + 
    geom_point(aes(fill=RSQR),shape=21,size=sz, stroke=0.5)  +
    pretty_plot()+
    scale_y_continuous(expand = c(0.05, 0))+
    scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1],name="R2") +
    guides(fill=guide_colorbar(title.vjust=0.75))+
    locustheme+
    labs(x="Position on Chromosome",y="-log10(P-value)") + 
    geom_point(data=subset(locuszoom,sentinel=="yes"),
               aes(x=POS/(10^6),y=log10),
               fill="yellow",shape=21,size=sz)
  
  # Finemap
  region<-region_number
  FM_region <- fread(paste0("/Volumes/broad_sankaranlab/ebao/FINEMAP/",trait,"/output/region",region,".snp"))
  FM_region$POS <- str_split_fixed(str_split_fixed(FM_region$snp,":",2)[,2],"_",2)[,1]
  FM_region$POS <- as.integer(as.character(FM_region$POS))
  FM_region <- merge(locuszoom,FM_region,by="POS")
  #FM_region <- subset(FM_region,snp_log10bf > -Inf)
  FM_region$snp_log10bf <- ifelse(FM_region$snp_log10bf == -Inf, -2.5, FM_region$snp_log10bf)
  FM_region$sentinel <- ifelse(FM_region$SNP == rsid, "yes", "no")
  
  fm <- ggplot(subset(FM_region,snp_log10bf >= -2.5),aes(POS/(10^6),snp_log10bf)) + 
    geom_point(aes(fill=RSQR),shape=21,size=sz, stroke=0.5)  +
    pretty_plot()+
    scale_y_continuous(expand = c(0.05, 0))+
    scale_fill_gradientn(colors = jdb_palette("solar_extra")[-1],name="R2") +
    guides(fill=guide_colorbar(title.vjust=0.75))+
    locustheme+
    labs(x="Position on Chromosome 6 (Mb)",y="log10(Bayes factor)") + 
    geom_point(data=subset(FM_region,sentinel=="yes"),
               aes(x=POS/(10^6),y=snp_log10bf),
               fill="yellow",shape=21,size=sz)
  
  return(list(p1,fm))
}

#### Run function to create locus plot and log10BF plots for each variant
# CEBPA
cebpa <- makePlots(rsid="rs78744187",trait="BASO_COUNT",region_number = 46)

# TRIM58
trim58 <- makePlots(rsid="rs3811444",trait="MCV",region_number = 2)

# ARGHEF3
arghef3 <- makePlots(rsid="rs1354034",trait="MPV",region_number = 40)

# PIK3CG
pik3g <- makePlots(rsid="rs342293",trait="MPV",region_number = 104)

# SMIM1
smim1 <- makePlots(rsid="rs1175550",trait="RETIC_COUNT",region_number = 4)

# RBM38
rbm38 <- makePlots(rsid="rs737092",trait="RBC_COUNT",region_number = 138)

# SH2B3
sh2b3 <- makePlots(rsid="rs3184504",trait="PLT_COUNT",region_number = 158)

# DNM3
dnm3 <- makePlots(rsid="rs2038479",trait="MPV",region_number = 2)


# Plot motifs for these 
BASO.motifs <- readRDS("../../data/motifbreakR/BASO_COUNT_Motifbreakr_output_PP001.rds")
MPV.motifs <- readRDS("../../data/motifbreakR/MPV_Motifbreakr_output_PP001.rds")
RETIC.motifs <- readRDS("../../data/motifbreakR/RETIC_COUNT_Motifbreakr_output_PP001.rds")
RBC.motifs <- readRDS("../../data/motifbreakR/RBC_COUNT_Motifbreakr_output_PP001.rds")
PLT.motifs <- readRDS("../../data/motifbreakR/PLT_COUNT_Motifbreakr_output_PP001.rds")

BASO.motifs[BASO.motifs$snpPos == known[1,]$pos]
MPV.motifs[MPV.motifs$snpPos == known[3,]$pos]; plotMB(MPV.motifs, "chr3:56849749:T:C")
MPV.motifs[MPV.motifs$snpPos == known[4,]$pos]
RETIC.motifs[RETIC.motifs$snpPos == known[5,]$pos]
RBC.motifs[RBC.motifs$snpPos == known[6,]$pos]
MPV.motifs[MPV.motifs$snpPos == known[8,]$pos]



