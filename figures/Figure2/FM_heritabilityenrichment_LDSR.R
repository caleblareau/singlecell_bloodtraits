library(data.table)
library(dplyr)
library(BuenColors)

# Heritability Enrichment of FM variants and GW-significant variants
traits=c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

h2 <-fread("/Volumes/broad_sankaranlab/ebao/LDscore/herit/justherit/trait_h2.txt")
colnames(h2) <- c("trait","h_obs","h_obs_se")

##############################################################################################################################
# FM enrichment of heritability
ctHeme <- vector("list",length=16)
ct_enrichment <- vector("list",length=length(ctHeme))

# Read in enrichments for FM PP001 variants
celltypes <- paste0(traits,"_PP001L2")
FM_enrichment <- NULL
dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/finemap_variants/"
for (i in 1:length(ctHeme)){
  ctHeme[[i]]<-fread(paste0(dir,traits[i],".",celltypes[i],".results"))
  FM_enrichment[i] <- ctHeme[[i]]$Enrichment[1]
}

# GWAS significant SNPs heritability
gwas_snps <- vector("list",length=16)
celltypes <- paste0(traits,".",traits,".sigSNPs")

dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/gwas_snps/"
for (i in 1:length(ctHeme)){
  gwas_snps[[i]]<-fread(paste0(dir,celltypes[i],".results"))
  GWAS_h2[i] <- gwas_snps[[i]]$Prop._h2[1]
  GWAS_enrichment[i] <- gwas_snps[[i]]$Enrichment[1]
}

################################################################################################
# Compare FM Enrichment to GWAS Enrichment
h2$gwas_enrichment <- GWAS_enrichment
h2$FM_enrichment <- FM_enrichment
h2$gwas_enrichment_se <- sapply(gwas_snps, function(y){
  return(y[,Enrichment_std_error][1])
})
h2$FM_enrichment_se <- sapply(ctHeme, function(y){
  return(y[,Enrichment_std_error][1])
})

melted <- melt(h2[,c("trait","gwas_enrichment","FM_enrichment")],
               id.vars = "trait",
               measure.vars = list(c("gwas_enrichment","FM_enrichment")))

melted$se <- melt(h2[,c("trait","gwas_enrichment_se","FM_enrichment_se")],
               id.vars = "trait",
               measure.vars = list(c("gwas_enrichment_se","FM_enrichment_se")))$value

# Plot figure
ggplot(data=melted, aes(x=trait,y=value,fill=factor(variable))) +
  geom_bar(stat="identity",position="dodge") +
  ggtitle("Heritability Enrichments") + 
  theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values=jdb_palette("Zissou")[c(1,5)],
                    labels=c("FM Variants", "GW-Significant Variants")) +
  guides(fill=guide_legend(title="")) + 
  theme(plot.title = element_text(size=10,hjust = 0.5,face="bold"), 
        panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y=element_blank(),axis.title.x=element_blank())  +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2,position=position_dodge(.9),color="black")

