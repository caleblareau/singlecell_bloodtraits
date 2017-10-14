library(data.table)
library(dplyr)
library(BuenColors)

# Heritability Enrichment of FM variants and GW-significant variants
traits=c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

h2 <-fread("/Volumes/broad_sankaranlab/ebao/LDscore/herit/justherit/trait_h2.txt")
colnames(h2) <- c("trait","h_obs","h_obs_se")

##############################################################################################################################
# FM enrichment of heritability
# Read in enrichments for FM PP001 variants
ctHeme <- vector("list",length=16)
ct_enrichment <- vector("list",length=length(ctHeme))

celltypes <- paste0(traits,"_PP001L2")
FM_enrichment <- NULL
dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/finemap_variants/"
for (i in 1:length(ctHeme)){
  ctHeme[[i]]<-fread(paste0(dir,traits[i],".",celltypes[i],".results"))
  FM_enrichment[i] <- ctHeme[[i]]$Enrichment[1]
}
h2$PP001_enrichment <- FM_enrichment
h2$PP001_enrichment_se <- sapply(ctHeme, function(y){
  return(y[,Enrichment_std_error][1])
})

# FM PP01 Heritability Enrichments
PP01_snps <- vector("list",length=16)
PP01_enrichment <- NULL
celltypes <- paste0(traits,".",traits,"_PP0.01")

dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/FM_PP01/"
for (i in 1:length(ctHeme)){
  PP01_snps[[i]]<-fread(paste0(dir,celltypes[i],".results"))
  PP01_enrichment[i] <- PP01_snps[[i]]$Enrichment[1]
}
h2$PP01_enrichment <- PP01_enrichment
h2$PP01_enrichment_se <- sapply(PP01_snps, function(y){
  return(y[,Enrichment_std_error][1])
})

# GWAS significant SNPs heritability
gwas_snps <- vector("list",length=16)
GWAS_enrichment <- NULL
celltypes <- paste0(traits,".",traits,".sigSNPs")

dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/gwas_snps/"
for (i in 1:length(ctHeme)){
  gwas_snps[[i]]<-fread(paste0(dir,celltypes[i],".results"))
  GWAS_enrichment[i] <- gwas_snps[[i]]$Enrichment[1]
}
h2$gwas_enrichment <- GWAS_enrichment
h2$gwas_enrichment_se <- sapply(gwas_snps, function(y){
  return(y[,Enrichment_std_error][1])
})

################################################################################################
# Compare FM Enrichment to GWAS Enrichment
melted <- melt(h2[,c("trait","gwas_enrichment","PP001_enrichment","PP01_enrichment")],
               id.vars = "trait",
               measure.vars = list(c("gwas_enrichment","PP001_enrichment","PP01_enrichment")))

melted$se <- melt(h2[,c("trait","gwas_enrichment_se","PP001_enrichment_se","PP01_enrichment_se")],
               id.vars = "trait",
               measure.vars = list(c("gwas_enrichment_se","PP001_enrichment_se","PP01_enrichment_se")))$value

melted$variable <- factor(melted$variable,levels=(c("PP01_enrichment","PP001_enrichment","gwas_enrichment")))
melted$trait <- factor(melted$trait,levels=traits)
  
# Plot figure
ggplot(data=melted, aes(x=trait,y=value,fill=variable)) +
  geom_bar(stat="identity",position="dodge") +
  theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values=jdb_palette("Zissou")[c(1,3, 5)],
                    labels=c("PP01 Variants", "PP001 Variants", "GW-Significant Variants")) +
  guides(fill=guide_legend(title="",reverse=TRUE)) + 
  labs(y="Pr(h2g)/Pr(SNPs)") +
  theme(plot.title = element_text(size=10,hjust = 0.5,face="bold"), 
        panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y=element_blank())  +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2,position=position_dodge(.9),color="black")