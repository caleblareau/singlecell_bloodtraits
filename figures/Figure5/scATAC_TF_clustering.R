library(BuenColors)
jdb_color_maps2 <- c(jdb_color_maps, "Mega" = "#FF347E", "UNK" = "#8D91AD", "mDC"= "#FFD700", "MCP" = "#C390D4")
names(jdb_color_maps2)[9] <- "Mono"
library(ggrepel)
library(SummarizedExperiment)
library(ggbeeswarm)
library(irlba)
library(data.table)
library(tidyverse)
library(qvalue)
library(gridExtra)
library(fpc)
"%ni%" <- Negate("%in%")

######################################################################################################
# Functions
remove_smooth <- function(cols){
  return(gsub("_PP001","",gsub("smooth_","",cols)))
}

compare_subgroups_plot <- function(allcells,celltype,traitstoplot,numPCs=5,graph=TRUE,colors=c("Mono","Ery"),kmeans=TRUE){
  set.seed(42)
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:(numPCs+2)]
  if (kmeans==TRUE){
    enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:(numPCs+2)]),centers=2,nstart = 500, iter.max=1000)$cluster)
  } else{
    res <- pamk(scale(enrichments[,3:(numPCs+2)]),krange=2,criterion="multiasw")
    enrichments$kmeans <- res$pamobject$clustering
    enrichments$kmeans <- sapply(enrichments$kmeans, function(x) {switch(x,1,2)})
    enrichments$kmeans <- as.factor(enrichments$kmeans)
  }
  enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  
  if (graph==TRUE){
    enrichments <- enrichments[,c(traitstoplot,"kmeans")]
    boxplots <- vector(mode="list",length=length(traitstoplot))
    for (i in 1:(ncol(enrichments)-1)){
      boxplots[[i]] <- ggplot(enrichments, aes_string(x="kmeans", y=colnames(enrichments)[i],fill="kmeans")) +
        geom_boxplot(outlier.shape=NA) + pretty_plot() +
        scale_fill_manual(values=as.character(jdb_color_maps2[colors])) +
        geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
        guides(fill=guide_legend(title=paste("ATAC","clustering",sep="\n"))) +
        theme(axis.title.x = element_blank())
    }
    return(boxplots)
  } else {
    enrichments <- enrichments[,c("name",traitstoplot,"kmeans")]
    return(enrichments)
  }
  
}

# CMP TF clustering -------------------------------------------------------
#setwd("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits/figures/Figure5/")

#CMPs
allcells <- fread("../../data/singlecell/scATAC/sc_traitenrichments_11aug.txt")

# Use 100-iteration average scATAC
avg <- fread("../../data/singlecell/scATAC/scGWASenrichments_average.tsv")
avg <- dcast(avg, V1~V2)
allcells <- merge(allcells[,1:12], avg, by.x="name",by.y="V1")

# CMP k-means clustering
traitstoplot <- c("RBC_COUNT","MPV","PLT_COUNT","MONO_COUNT")

# Extract dataframe with annotation for k-medoids clustering by ATAC PCs
enrichments <- compare_subgroups_plot(allcells,celltype="CMP",traitstoplot=traitstoplot,
                                      numPCs=5,graph=FALSE,kmeans=F)

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
    scale_fill_manual(values=as.character(jdb_color_maps2[c("Ery","Mono")])) +
    geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
    
    theme(axis.title.x = element_blank())
  return(p)
})

# Plot differential z-scores for TFs of interest
CMP_TFs <- ggmatrix(TFplots,1,length(TFs_of_interest),
                    xAxisLabels = TFs_of_interest) + 
  theme(plot.title = element_text(hjust = 0.5))

#ggsave(CMP_TFs, file="5Cb_CMP_TFzscores_REVISED.pdf",width=6,height=3, useDingbats=F)


# MEP TF clustering -------------------------------------------------------
celltype <- "MEP"
meptraits=c("RBC_COUNT", "HCT","PLT_COUNT") 

# K-medoids cluster by ATAC PCs
enrichments <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=meptraits,
                                      numPCs=7,graph=FALSE,kmeans=F)

TFs_of_interest <- c("GATA1","KLF1","MEF2C")

# For the TFs of interest, extract their z-scores for all single cells of a cell type and plot k-means cluster vs. z-score
TFplots <- lapply(TFs_of_interest, function(TF){
  idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
  zscores <- TF_zscores %>% dplyr::select(cellnames,idx)
  colnames(zscores) <- c("cellnames",TF)
  merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
  p <- ggplot(merged, aes_string(x="kmeans", y=colnames(merged)[ncol(merged)],fill="kmeans")) + 
    geom_boxplot(outlier.shape=NA) + pretty_plot() + 
    scale_fill_manual(values=as.character(jdb_color_maps2[c("Ery","Mega")])) +
    geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
    theme(axis.title.x = element_blank())
  return(p)
})

# Plot differential z-scores for TFs of interest
mep_TFs <- ggmatrix(TFplots,1,length(TFs_of_interest),
                    xAxisLabels = TFs_of_interest) + 
  theme(plot.title = element_text(hjust = 0.5))

#ggsave(mep_TFs, file="5Eb_MEP_TFzscores_REVISED.pdf",width=6,height=3, useDingbats=F)
