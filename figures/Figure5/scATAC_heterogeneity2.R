library(GGally)
library(corrplot)
library(BuenColors)
library(ggrepel)
library(SummarizedExperiment)
library(ggbeeswarm)
library(chromVAR)
library(irlba)
library(data.table)
library(dplyr)
library(stringr)
library(qvalue)
library(gridExtra)

######################################################################################################
# Functions
# k-means on enrichments
ggallyplot <- function(allcells,celltype,traitstoplot,smoothed=FALSE,graph=TRUE){
  set.seed(42)
  if(smoothed) enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,29:ncol(allcells)]
  else {
    enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:28]
  }
  
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  rownames(enrichments) <- enrichments$name
  enrichments <- enrichments[,c(traitstoplot)]
  colnames(enrichments) <- abbreviate(colnames(enrichments))
  enrichments$kmeans <- as.character(kmeans(scale(as.matrix(enrichments)),2)$cluster)
  
  # test for significance of k-means clusters
  pvs <- NULL
  for (col in 1:(ncol(enrichments)-1)){
    pv <- t.test(enrichments[enrichments$kmeans %in% 1,col], enrichments[enrichments$kmeans %in% 2,col])$p.value
    pvs <- rbind(pvs,c(colnames(enrichments)[col],pv))
  }
  pvs <- cbind(pvs,p.adjust(pvs[,2],method = "fdr"))
  colnames(pvs)<-c("trait","p.value","FDR")
  print(pvs)
  
  if (graph){
    boxplots <- vector(mode="list",length=length(traitstoplot))
    for (i in 1:(ncol(enrichments)-1)){
      boxplots[[i]] <- ggplot(enrichments, aes_string(x="kmeans", y=colnames(enrichments)[i],fill="kmeans")) + 
        geom_boxplot(outlier.shape=NA) + pretty_plot() + 
        scale_fill_manual(values=as.character(jdb_color_maps2[c("Mega","Ery")])) +
        geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
        theme(axis.title.x = element_blank()) 
    }
    return(boxplots)
  } else{
    return(enrichments)
  }
  
}

# kmeans on ATAC PCs
ggallyplot_atac <- function(allcells,celltype,traitstoplot,smoothed=FALSE,graph=FALSE,numPCs=5){
  set.seed(42)
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:12]
  enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:(3+numPCs-1)]),2)$cluster)
  subenrichments<- split(enrichments, f = enrichments$kmeans )
  
  if(smoothed) {
    enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,29:ncol(allcells))])
  }
  else {
    enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
  }
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  enrichments <- enrichments[,c(traitstoplot,"kmeans")]
  colnames(enrichments) <- abbreviate(colnames(enrichments))
  colnames(enrichments)[ncol(enrichments)] <- "kmeans"
  
  # test for significance of k-means clusters
  pvs <- NULL
  for (col in 1:(ncol(enrichments)-1)){
    pv <- t.test(enrichments[enrichments$kmeans %in% 1,col], enrichments[enrichments$kmeans %in% 2,col])$p.value
    pvs <- rbind(pvs,c(colnames(enrichments)[col],pv))
  }
  pvs <- cbind(pvs,p.adjust(pvs[,2],method = "fdr"))
  colnames(pvs)<-c("trait","p.value","FDR")
  print(pvs)
  
  # Check if splitting again will create significant differences
  for (i in 1:length(subenrichments)){
    subenrichments[[i]]$kmeans<-  NULL
    subenrichments[[i]]$kmeans <- as.character(kmeans(scale(subenrichments[[i]][,3:(3+numPCs-1)]),2)$cluster)
  }
  
  if(smoothed) {
    subenrichments <- lapply(subenrichments, function(y){
      merge(y,as.data.frame(allcells)[allcells$type %in% celltype,c(1,29:ncol(allcells))])
    })
  }
  else {
    subenrichments <- lapply(subenrichments, function(y){
      merge(y,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
    })
  }
  
  subenrichments <- lapply(subenrichments, function(enrichments){
    colnames(enrichments) <- remove_smooth(colnames(enrichments))
    enrichments <- enrichments[,c(traitstoplot,"kmeans")]
    colnames(enrichments) <- abbreviate(colnames(enrichments))
    colnames(enrichments)[ncol(enrichments)] <- "kmeans"
    return(enrichments)
  })
  
  # test for significance of k-means clusters
  print("Sub-cluster p-values:")
  pvs <- NULL
  lapply(subenrichments, function(enrichments){
    for (col in 1:(ncol(enrichments)-1)){
      pv <- t.test(enrichments[enrichments$kmeans %in% 1,col], enrichments[enrichments$kmeans %in% 2,col])$p.value
      pvs <- rbind(pvs,c(colnames(enrichments)[col],pv))
    }
    pvs <- cbind(pvs,p.adjust(pvs[,2],method = "fdr"))
    colnames(pvs)<-c("trait","p.value","FDR")
    print(pvs)
  })
  
  # Plot pairwise correlations
  if (graph){
    ggpairs(enrichments,ggplot2::aes(colour=kmeans)) + 
      theme_bw() + pretty_plot()
  }
}

remove_smooth <- function(cols){
  return(gsub("_PP001","",gsub("smooth_","",cols)))
}

compare_subgroups_plot <- function(allcells,celltype,traitstoplot,numPCs=5,graph=TRUE){
  set.seed(42)
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:12]
  enrichments$kmeans <- as.character(kmeans(scale(as.matrix(enrichments[,3:(3+numPCs-1)])),2)$cluster)
  enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  
  if (graph==TRUE){
    enrichments <- enrichments[,c(traitstoplot,"kmeans")]
    boxplots <- vector(mode="list",length=length(traitstoplot))
    for (i in 1:(ncol(enrichments)-1)){
      boxplots[[i]] <- ggplot(enrichments, aes_string(x="kmeans", y=colnames(enrichments)[i],fill="kmeans")) + 
        geom_boxplot(outlier.shape=NA) + pretty_plot() + 
        scale_fill_manual(values=as.character(jdb_color_maps2[c("Mono","Ery")])) +
        geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
        theme(axis.title.x = element_blank()) 
    }
    
    # Plot PCs
    smoothed_traits <- colnames(allcells)[29:ncol(allcells)]
    smoothed_traits <- sapply(traitstoplot, function(y){
      grep(y,smoothed_traits,value=TRUE)
    })
    ct <- allcells[allcells$type %in% celltype,]
    PC_clusters <- lapply(smoothed_traits,function(y) plotPCs(ct,y,celltype=celltype))
    
    # combine box plots and PC scatter plots and plot them together
    all <- c(boxplots,PC_clusters)
    return(all)
  } else {
    enrichments <- enrichments[,c("name",traitstoplot,"kmeans")]
    return(enrichments)
  }
  
}

plotPCs <- function(graph,trait,celltype,legend=TRUE){
  if(legend){
    p<- ggplot(graph,aes(PC2,PC3)) + geom_point(aes_string(colour=trait)) +
      scale_color_gradientn(colors=jdb_palette("brewer_heat")) +
      pretty_plot() + theme_bw() + 
      ggtitle(paste(celltype,gsub("_PP001","",gsub("smooth_","",trait)))) + 
      theme(plot.title = element_text(size=14,hjust = 0.5,face="bold"),
            panel.background = element_rect(fill = "white", colour = "grey50"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text.x = element_blank(), axis.ticks.x=element_blank(),
            axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
      labs(x="",y="",fill="ChromVAR z-score") 
  } else{
    p<- ggplot(graph,aes(PC2,PC3)) + geom_point(aes_string(colour=trait)) +
      scale_color_gradientn(colors=jdb_palette("brewer_heat")) +
      pretty_plot() + theme_bw() + 
      ggtitle(paste(celltype,gsub("_PP001","",gsub("smooth_","",trait)))) + 
      theme(legend.position="none",plot.title = element_text(size=14,hjust = 0.5,face="bold"),
            panel.background = element_rect(fill = "white", colour = "grey50"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text.x = element_blank(), axis.ticks.x=element_blank(),
            axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
      labs(x="",y="")
  }
  return(p)
}

# Extract legend from a plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
######################################################################################################
#CMPs
allcells <- fread("../../data/singlecell/scATAC/sc_traitenrichments_11aug.txt")

# CMP k-means clustering
CMPtraitstoplot <- c("RBC_COUNT","PLT_COUNT","MPV","MONO_COUNT")

# Cluster and plot cluster means 
plots <- compare_subgroups_plot(allcells,celltype="CMP",traitstoplot=traitstoplot,
                                numPCs=5,graph=TRUE)
CMP_boxplots <- ggmatrix(plots[1:4],nrow=1,ncol=length(traitstoplot),
                         xAxisLabels = traitstoplot)
ggsave(CMP_boxplots, file="5B_CMP_kmeans.pdf",
       width=6,height=3)

# Extract legend and save separately
legend <- g_legend(plots[[1]])
grid.arrange(legend, ncol=1, nrow=1)

# Extract dataframe with annotation for k-means cluster by ATAC PCs
enrichments <- compare_subgroups_plot(allcells,celltype="CMP",traitstoplot=traitstoplot,
                                      numPCs=5,graph=FALSE)

# Merge with TF zscores
TF_zscores <- fread("../../data/singlecell/scATAC/tfDeviationsTable.tsv")
#CMP TFs
TFs_of_interest <- c("GATA1","KLF1","CEBPA","IRF8")

# For the TFs of interest, extract their z-scores for all single cells of a cell type and plot k-means cluster vs. z-score
TFplots <- lapply(TFs_of_interest, function(TF){
  idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
  zscores <- TF_zscores %>% dplyr::select(cellnames,idx)
  colnames(zscores) <- c("cellnames",TF)
  merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
  p <- ggplot(merged, aes_string(x="kmeans", y=colnames(merged)[ncol(merged)],fill="kmeans")) + 
    geom_boxplot(outlier.shape=NA) + pretty_plot() + 
    scale_fill_manual(values=as.character(jdb_color_maps2[c("Mono","Ery")])) +
    geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
    
    theme(axis.title.x = element_blank())
  return(p)
})

# Plot differential z-scores for TFs of interest
ggmatrix(TFplots,1,length(TFs_of_interest),
         xAxisLabels = TFs_of_interest) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(TFs, file="5B_CMP_TFzscores.pdf",width=6,height=3)

# CMP TF rank order plot
TF_differences <- readRDS("../../data/singlecell/scATAC/CMP_TFzscore_differences.rds")
TF_differences$TF_name <- str_split_fixed(TF_differences$TF, "_",n=4)[,3]
TF_differences$FDR <- qvalue(TF_differences$pval)$qvalues

# Color the GATA TFs
idx <- grep("GATA",TF_differences$TF_name)
TF_differences$highlight<- "F"
TF_differences[idx,"highlight"] <- "T"

# Take the top TFs from TFs of interest
labelidx <- sapply(TFs_of_interest, function(y) {
  grep(paste("^",y,"$", sep=""), TF_differences$TF_name)[1]
})
TF_differences$toLabel <- "F"
TF_differences[labelidx,"toLabel"] <- TF_differences[labelidx,"TF_name"]

# Rank order plot
p <- ggplot(TF_differences,aes(x=rank,y=-log10(FDR))) + 
  geom_point(shape=21,size=3.5,aes(fill=highlight)) +
  pretty_plot() +
  scale_fill_manual(values = as.character(jdb_color_maps2[c("Mono","Ery")]),
                    labels=c("Other TFs", "GATA TFs")) +
  labs(x="Rank",y="-log10(FDR)") + 
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black")) +
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text_repel(
    data = subset(TF_differences, toLabel!="F"),
    aes(label = toLabel),
    size = 2.5,
    nudge_y =2,
    nudge_x=35,
    point.padding=0.50,
    min.segment.length=0,
    segment.alpha=0.3,
    direction="both")

ggsave(p, file="5C_CMP_TFs_rankorderplot.pdf",
       width=6,height=6)

####################################################################################################  
# MEP k-means clustering
meptraits=c("RBC_COUNT","HCT","PLT_COUNT","MPV")

# K-means cluster by ATAC PCs
enrichments <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=traitstoplot,
                                      numPCs=5,graph=FALSE)
plots  <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=traitstoplot,
                                 numPCs=5,graph=TRUE)
ggmatrix(plots,1,length(meptraits),xAxisLabels = meptraits)

# Cluster MEP population by g-chromVAR
mep_chromvar <- ggallyplot(allcells,celltype="MEP",traitstoplot=meptraits,smoothed=FALSE)
mep_chromvar_plots <- ggmatrix(mep_chromvar,1,length(meptraits),
                               xAxisLabels = meptraits) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(mep_chromvar_plots, file="5D_MEP_gChromVAR_Clustering.pdf",
       width=6,height=3)

legend <- g_legend(mep_chromvar[[1]])
grid.arrange(legend, ncol=1, nrow=1)

# MEP TFs
# K-means cluster by g-ChromVAR z-scores
enrichments <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=meptraits,
                                      numPCs=5,graph=FALSE)

TFs_of_interest <- c("KLF1","GATA1","MEF2C")

# For the TFs of interest, extract their z-scores for all single cells of a cell type and plot k-means cluster vs. z-score
TFplots <- lapply(TFs_of_interest, function(TF){
  idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
  zscores <- TF_zscores %>% dplyr::select(cellnames,idx)
  colnames(zscores) <- c("cellnames",TF)
  merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
  p <- ggplot(merged, aes_string(x="kmeans", y=colnames(merged)[ncol(merged)],fill="kmeans")) + 
    geom_boxplot(outlier.shape=NA) + pretty_plot() + 
    scale_fill_manual(values=as.character(jdb_color_maps2[c("Mega","Ery")])) +
    geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
    
    theme(axis.title.x = element_blank())
  return(p)
})

# Plot differential z-scores for TFs of interest
mep_TFs <- ggmatrix(TFplots,1,length(TFs_of_interest),
         xAxisLabels = TFs_of_interest) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(mep_TFs, file="5D_MEP_TFzscores.pdf",
       width=6,height=3)

# Calculate pvalues for all TFs
TFplots <- lapply(colnames(TF_zscores)[2:ncol(TF_zscores)], function(TF){
  # idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
  zscores <- TF_zscores %>% dplyr::select(cellnames,TF)
  merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
  p <- t.test(merged[merged$kmeans %in% 1,TF], merged[merged$kmeans %in% 2,TF])$p.value
  return(p)
})

# MEP TF rank order plot
#TF_differences <- readRDS("/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/scATAC/MEP_TFzscore_differences.rds")

TF_differences <- as.data.frame(unlist(TFplots))
TF_differences$TF <- colnames(TF_zscores)[2:ncol(TF_zscores)]
colnames(TF_differences) <- c("pval","TF")
TF_differences <- arrange(TF_differences,pval)
TF_differences$logp <- -1*log10(TF_differences$pval)
TF_differences$rank <- seq(1,nrow(TF_differences),1)
TF_differences$TF_name <- str_split_fixed(TF_differences$TF, "_",n=4)[,3]
TF_differences$FDR <- qvalue(TF_differences$pval)$qvalues

# Color the GATA TFs
idx <- grep("GATA",TF_differences$TF_name)
TF_differences$highlight<- "F"
TF_differences[idx,"highlight"] <- "T"

# Take the top TFs from TFs of interest
labelidx <- sapply(TFs_of_interest, function(y) {
  grep(paste("^",y,"$", sep=""), TF_differences$TF_name)[1]
})
TF_differences$toLabel <- "F"
TF_differences[labelidx,"toLabel"] <- TF_differences[labelidx,"TF_name"]

#saveRDS(TF_differences,"/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/scATAC/MEP_gChromVAR_TFzscore_differences.rds")
# Rank order plot
ggplot(TF_differences,aes(x=rank,y=-log10(FDR))) + 
  geom_point(shape=21,size=3.5,aes(fill=highlight)) +
  pretty_plot() +
  scale_fill_manual(values = as.character(jdb_color_maps2[c("Mono","Ery")]),
                    labels=c("Other TFs", "GATA TFs")) +
  labs(x="Rank",y="-log10(FDR)") + 
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black")) +
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text_repel(
    data = subset(TF_differences, toLabel!="F"),
    aes(label = toLabel),
    size = 2.5,
    nudge_y =2,
    nudge_x=35,
    point.padding=0.50,
    min.segment.length=0,
    segment.alpha=0.3,
    direction="both")

ggsave(p, file="/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/Figure 5 - scATAC heterogeneity/5D_CMP_TFs_rankorderplot.pdf",
       width=6,height=6)


########################################################################################
# Corrplot of chromVAR enrichments

M<- cor(CMP_enrichments[,c(1:16)])
corrplot(M,method="square",order="hclust",hclust.method = "complete", tl.cex=0.5,
         cl.lim=c(-1,1),cl.cex=0.8,mar=c(0,0,2,0),tl.col="black",tl.srt=45,title="ChromVAR Correlations")

# Read unix files into R lists
setwd("/Volumes/broad_sankaranlab/ebao/LDscore/gen.correlations/constrained_UK10K/")
constrained_gencors <- vector(mode = "list", length = length(traits))

traits=c("BASO_COUNT","EO_COUNT","HCT","HGB","LYMPH_COUNT", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL","MONO_COUNT", "MPV", "NEUTRO_COUNT", "PLT_COUNT", "RBC_COUNT","RETIC_COUNT","WBC_COUNT")

constrained_gencors <- lapply(traits, function(y) read.table(paste0(y,"_UK10K_constrained.gcsummary.txt"),header=TRUE,row.names=NULL))

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

output <- make.gc.matrix(constrained_gencors,traits)[[1]]
constrained.pv.matrix <- make.gc.matrix(constrained_gencors,traits)[[2]]
p.adjust(constrained.pv.matrix[lower.tri(constrained.pv.matrix)],"fdr")

par(mfrow=c(1,2))
# Genetic correlations
reorderedtraits<- dimnames(corrplot(output,method="square",order="hclust",hclust.method = "complete", tl.cex=0.5,
                                    cl.lim=c(-1,1),cl.cex=0.6,tl.col="black",tl.srt=45,mar=c(1,0,2,0),title="Genetic Correlations"))[[1]]
celltypes <- unique(allcells$type)
# ChromVAR
M<- cor(allcells[,29:ncol(allcells)])
rownames(M)=colnames(M)=traits
corrplot(M[reorderedtraits,reorderedtraits],method="square", tl.cex=0.5,
         cl.lim=c(-1,1),cl.cex=0.6,mar=c(1,0,2,0),tl.col="black",tl.srt=45,title="ChromVAR Correlations")


M<- cor(CMP_enrichments[,1:16])
M<- cor(allcells[allcells$type %in% "MEP",29:ncol(allcells)])
rownames(M)=colnames(M)=traits
reorderedtraits<- dimnames(corrplot(M,method="square",order="hclust",hclust.method = "complete", tl.cex=0.5,
                                    cl.lim=c(-1,1),cl.cex=0.6,tl.col="black",tl.srt=45,mar=c(1,0,2,0),title="ChromVAR Correlations"))[[1]]

celltype <- "CMP"
CMP<- cor(allcells[allcells$type %in% celltype,29:ncol(allcells)])
rownames(CMP)=colnames(CMP)=traits
reorderedtraits<- dimnames(corrplot(CMP,method="square",order="hclust",hclust.method = "complete", tl.cex=0.5,
                                    cl.lim=c(-1,1),cl.cex=0.6,tl.col="black",tl.srt=45,mar=c(1,0,2,0),title=celltype))[[1]]

corrplot(M[reorderedtraits,reorderedtraits],method="square", tl.cex=0.5,
         cl.lim=c(-1,1),cl.cex=0.6,mar=c(1,0,2,0),tl.col="black",tl.srt=45,title="MEP")

par(mfrow=c(1,4))
celltype <- "HSC"
cormat<- cor(allcells[allcells$type %in% celltype,29:ncol(allcells)])
rownames(cormat)=colnames(cormat)=traits
corrplot(cormat[reorderedtraits,reorderedtraits],method="square", tl.cex=0.5,
         cl.lim=c(-1,1),cl.cex=0.6,mar=c(1,0,2,0),tl.col="black",tl.srt=45,title=celltype)

reorderedtraits<- dimnames(corrplot(cormat,method="square",order="hclust",hclust.method = "complete", tl.cex=0.5,
                                    cl.lim=c(-1,1),cl.cex=0.6,tl.col="black",tl.srt=45,mar=c(1,0,2,0),title=celltype))[[1]]

celltype <- "MPP"
cormat<- cor(allcells[allcells$type %in% celltype,29:ncol(allcells)])
rownames(cormat)=colnames(cormat)=traits
corrplot(cormat[reorderedtraits,reorderedtraits],method="square", tl.cex=0.5,
         cl.lim=c(-1,1),cl.cex=0.6,mar=c(1,0,2,0),tl.col="black",tl.srt=45,title=celltype)


par(mfrow=c(1,2))
celltype <- "CMP"
cormat<- cor(allcells[allcells$type %in% celltype,29:ncol(allcells)])
rownames(cormat)=colnames(cormat)=traits
reorderedtraits<- dimnames(corrplot(cormat,method="square",order="hclust",hclust.method = "complete", tl.cex=0.5,
                                    cl.lim=c(-1,1),cl.cex=0.6,tl.col="black",tl.srt=45,mar=c(1,0,2,0),title=celltype))[[1]]

celltype <- "MEP"
cormat<- cor(allcells[allcells$type %in% celltype,29:ncol(allcells)])
rownames(cormat)=colnames(cormat)=traits
corrplot(cormat[reorderedtraits,reorderedtraits],method="square", tl.cex=0.5,
         cl.lim=c(-1,1),cl.cex=0.6,mar=c(1,0,2,0),tl.col="black",tl.srt=45,title=celltype)

