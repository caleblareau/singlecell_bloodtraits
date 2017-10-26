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

# FM PP10 Heritability Enrichments
PP10_snps <- vector("list",length=16)
PP10_enrichment <- NULL
celltypes <- paste0(traits,".",traits,".PP10")

dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/FM_PP10/"
for (i in 1:length(ctHeme)){
  PP10_snps[[i]]<-fread(paste0(dir,celltypes[i],".results"))
  PP10_enrichment[i] <- PP10_snps[[i]]$Enrichment[1]
}
h2$PP10_enrichment <- PP10_enrichment
h2$PP10_enrichment_se <- sapply(PP10_snps, function(y){
  return(y[,Enrichment_std_error][1])
})

# GWAS LD80 to Sentinel Heritability Enrichments
LD80_snps <- vector("list",length=16)
LD80_enrichment <- NULL
celltypes <- paste0(traits,".",traits,".LD80")

dir="/Volumes/broad_sankaranlab/ebao/LDscore/herit/LD80/"
for (i in 1:length(ctHeme)){
  LD80_snps[[i]]<-fread(paste0(dir,celltypes[i],".results"))
  LD80_enrichment[i] <- LD80_snps[[i]]$Enrichment[1]
}
h2$LD80_enrichment <- LD80_enrichment
h2$LD80_enrichment_se <- sapply(LD80_snps, function(y){
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
melted <- melt(h2[,c("trait","gwas_enrichment","PP001_enrichment","PP01_enrichment","PP10_enrichment","LD80_enrichment")],
               id.vars = "trait",
               measure.vars = list(c("gwas_enrichment","PP001_enrichment","PP01_enrichment","PP10_enrichment","LD80_enrichment")))

melted$se <- melt(h2[,c("trait","gwas_enrichment_se","PP001_enrichment_se","PP01_enrichment_se","PP10_enrichment_se","LD80_enrichment_se")],
               id.vars = "trait",
               measure.vars = list(c("gwas_enrichment_se","PP001_enrichment_se","PP01_enrichment_se","PP10_enrichment_se","LD80_enrichment_se")))$value

#melted$value <- factor(melted$value,levels=(c("PP10_enrichment_se","PP001_enrichment_se","gwas_enrichment_se","LD80_enrichment_se","PP01_enrichment_se")))

bedsizes <- fread("../../data/Finemap/bedsizes.txt")
melted$size <- melt(bedsizes[,c("trait","all_gwas","PP001","PP01","PP10","LD80")],id.vars="trait")$value


# Plot figure
toplot <-c("PP10_enrichment","PP01_enrichment","LD80_enrichment")
melted$variable <- factor(melted$variable,levels=(toplot))
melted$trait <- factor(melted$trait,levels=traits)

sz=4
wd = 6.5
p1 <- ggplot(data=subset(melted,variable %in% toplot),
                   aes(x=trait,y=value,fill=variable)) +
  geom_bar(stat="identity",position="dodge") +
  theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values=jdb_palette("Zissou")[c(5,3,1)],
                    labels=c("PP10 Variants","PP01 Variants","R2 > 0.80 Variants")) +
  guides(fill=guide_legend(title="",reverse=TRUE)) + 
  labs(y="Pr(h2g)/Pr(SNPs)") +
  theme(plot.title = element_text(size=10,hjust = 0.5,face="bold"), 
        panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y=element_blank(),
        axis.text = element_text(size=sz),
        axis.title.x=element_text(size=sz*1.3,margin = margin(t = 4, r = 0, b = 0, l = 0)),
        legend.position = c(.95, .98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical",
        legend.key.size = unit(wd/30, "in"),
        legend.text = element_text(size=sz*1.2))  +
  scale_y_continuous(expand = c(0.04, 0))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2,position=position_dodge(.9),color="black") + 
  geom_text(position=position_dodge(.9),aes(label=size,y=value+se), hjust=-0.5, color="black", size=sz/2.5)
  

dir="/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/LDSR/Heritabilities"
ggsave("LDSR_heritabilityenrichments_LD80_PP01_PP10.pdf", plot = p1, device = NULL, path = dir,
       scale = 1, width = wd, height=wd/1.5, units = "in",
       dpi = 300, limitsize = TRUE)


###########################################################################
# Total proportion of heritability
output <- sapply(traits, function(trait){
  PP001_herit <- sapply(ctHeme, function(y){
    return(y[,Prop._h2][1])
  })
})
all_herit <- do.call(rbind,(lapply(seq(16), function(y){
  PP001_h2<- ctHeme[[y]][,Prop._h2][1]
  PP001_h2_se <- ctHeme[[y]][,Prop._h2_std_error][1]
  
  PP01_h2<- PP01_snps[[y]][,Prop._h2][1]
  PP01_h2_se <- PP01_snps[[y]][,Prop._h2_std_error][1]
  
  gwas_h2<- gwas_snps[[y]][,Prop._h2][1]
  gwas_h2_se <- gwas_snps[[y]][,Prop._h2_std_error][1]
  
  return(cbind(PP001_h2,PP001_h2_se,
               PP01_h2, PP01_h2_se,
               gwas_h2, gwas_h2_se))
}))) %>% as.data.frame()
all_herit$trait<- traits

melted <- melt(all_herit[,c("trait","gwas_h2","PP001_h2","PP01_h2")])

melted$se <- melt(all_herit[,c("trait","gwas_h2_se","PP001_h2_se","PP01_h2_se")],
                  id.vars = "trait")$value

melted$variable <- factor(melted$variable,levels=(c("gwas_h2","PP001_h2","PP01_h2")))
melted$trait <- factor(melted$trait,levels=traits)

# Plot figure
p1 <- ggplot(data=melted, aes(x=trait,y=value,fill=variable)) +
  geom_bar(stat="identity",position="dodge") +
  theme_bw() + 
  scale_fill_manual(values=jdb_palette("Zissou")[c(1,3,5)],
                    labels=c("GW-Significant Variants", "PP001 Variants", "PP01 Variants")) +
  guides(fill=guide_legend(title="")) + 
  labs(y="Pr(h2g)/Pr(SNPs)") +
  theme(plot.title = element_text(size=10,hjust = 0.5,face="bold"), 
        panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.text = element_text(size=sz),
        axis.title.y=element_text(size=sz,margin = margin(t = 3, r = 0, b = 0, l = 0)),
        legend.position = c(.95, .98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical",
        legend.key.size = unit(width/30, "in"),
        legend.text = element_text(size=sz),
        axis.text.x = element_text(angle = 60, hjust = 1))  +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2,position=position_dodge(.9),color="black")

dir="/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/LDSR/Heritabilities"
wd = 6.5
ggsave("LDSR_prop_h2.pdf", plot = p1, device = NULL, path = dir,
       scale = 1, width = wd, height=wd, units = "in",
       dpi = 300, limitsize = TRUE)


# Make bed file for LD80 variants
# LD80 <- readRDS("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits/data/Finemap/LD80.rds")
# for (i in 1:16){
#   nonsentinel <- as.data.frame(LD80[[i]])[,c("chromosome","nonsentinel")]
#   colnames(nonsentinel) <- c("chromosome","start")
#   sentinel <- as.data.frame(LD80[[i]])[,c("chromosome","sentinel")]
#   colnames(sentinel) <- c("chromosome","start")
#   merged <- bind_rows(nonsentinel,sentinel)[!duplicated(bind_rows(nonsentinel,sentinel)),] %>% arrange(chromosome)
#   merged$chromosome <- paste0("chr",merged$chromosome)
#   merged$end <- merged$start + 1
#   write.table(merged,paste0("/Volumes/broad_sankaranlab/ebao/LDscore/beds/LD80/",traits[i],".LD80.bed"),quote = FALSE, sep = "\t", col.names = F, row.names = F)
# }
# LD80[1]
