library(GGally)
library(corrplot)
library(BuenColors)
library(ggrepel)
library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)
library(chromVAR)
library(irlba)
library(ggbeeswarm)
library(SummarizedExperiment)

##########################################################################################
# Functions
# Cluster by ChromVAR, then plot ggpairs
ggallyplot <- function(allcells,celltype,traitstoplot,smoothed=FALSE){
  if(smoothed) enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,29:ncol(allcells)]
  else {
    enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:28]
  }
  
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  enrichments <- enrichments[,c(traitstoplot)]
  colnames(enrichments) <- abbreviate(colnames(enrichments))
  enrichments$kmeans <- as.character(kmeans(scale(enrichments),2)$cluster)
  
  # test for significance of k-means clusters
  pvs <- NULL
  for (col in 1:(ncol(enrichments)-1)){
    pv <- t.test(enrichments[enrichments$kmeans %in% 1,col], enrichments[enrichments$kmeans %in% 2,col])$p.value
    pvs <- rbind(pvs,c(colnames(enrichments)[col],pv))
  }
  print(pvs)
  
  ggpairs(enrichments,ggplot2::aes(colour=kmeans)) + theme_bw() + pretty_plot()
}

# Plot ATAC PCs in scatter plot 
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
# Cluster by ATAC
ggallyplot_atac <- function(allcells,celltype,traitstoplot,smoothed=FALSE,graph=FALSE,numPCs=5){
  # Filter for the cell type of interest
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:12]
  enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:(3+numPCs-1)]),2)$cluster)
  # Split k-means into different list elements
  subenrichments<- split(enrichments, f = enrichments$kmeans )
  
  if(smoothed) {
    enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,29:ncol(allcells))])
  }
  else {
    enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
  }
  # Rename colnames so that they just have the base trait names
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  # Filter by traits of interest
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
  
  # Check if splitting by k-means again will create significant differences
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
  
  # test for significance of k-means clusters and print p-values
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
  
  # Plot pairwise correlations in a ggpairs figure
  if (graph){
    ggpairs(enrichments,ggplot2::aes(colour=kmeans)) + 
      theme_bw() + pretty_plot()
  }
}

compare_subgroups_plot <- function(allcells,celltype,traitstoplot,numPCs=5){
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:12]
  enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:(3+numPCs-1)]),2)$cluster)
  enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  enrichments <- enrichments[,c(traitstoplot,"kmeans")]
  
  boxplots <- vector(mode="list",length=length(traitstoplot))
  for (i in 1:(ncol(enrichments)-1)){
    boxplots[[i]] <- ggplot(enrichments, aes_string(x="kmeans", y=colnames(enrichments)[i],fill="kmeans")) + 
      geom_boxplot(outlier.shape=NA) + pretty_plot() + 
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
}

compare_subgroups_plot_v2 <- function(allcells,celltype,traitstoplot,numPCs=5,PCs){
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:2]
  PCs$kmeans <- as.character(kmeans(scale(PCs[,1:numPCs]),2)$cluster)
  enrichments <- merge(enrichments,PCs,by.x="name",by.y="cellnames")
  enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  enrichments <- enrichments[,c(traitstoplot,"kmeans")]
  
  boxplots <- vector(mode="list",length=length(traitstoplot))
  for (i in 1:(ncol(enrichments)-1)){
    boxplots[[i]] <- ggplot(enrichments, aes_string(x="kmeans", y=colnames(enrichments)[i],fill="kmeans")) + 
      geom_boxplot(outlier.shape=NA) + pretty_plot() + 
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
  ggmatrix(all,nrow=2,ncol=length(traitstoplot),
           xAxisLabels = traitstoplot,
           legend=c(1,1))+ 
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification = c("right", "top"))
}
remove_smooth <- function(cols){
  return(gsub("_PP001","",gsub("smooth_","",cols)))
}
##########################################################################################
# Code begins here
allcells <- fread("../../data/singlecell/scATAC/sc_traitenrichments_11aug.txt")

# CMP k-means clustering
CMPtraitstoplot <- c("RBC_COUNT","PLT_COUNT","MPV","MONO_COUNT")

# plot and get p-values for clustered differences, clustered by ATAC
ggallyplot_atac(allcells,celltype="CMP",traitstoplot=CMPtraitstoplot,smoothed=FALSE)

# Generate CMP plot 
celltype <- "CMP"
traitstoplot <- c("RBC_COUNT","PLT_COUNT","MPV","MONO_COUNT")
traitstoplot <- c("HCT","PLT_COUNT","MPV","MONO_COUNT")

#this function call does everything that the next ~50 subsequent lines do
plots <- compare_subgroups_plot(allcells,celltype="CMP",traitstoplot=traitstoplot,
                       numPCs=6)
ggmatrix(plots,nrow=2,ncol=length(traitstoplot),
         xAxisLabels = traitstoplot)

enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:7]
enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:7]),2)$cluster)
enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
colnames(enrichments) <- remove_smooth(colnames(enrichments))
enrichments <- enrichments[,c(traitstoplot,"kmeans")]

boxplots <- vector(mode="list",length=length(traitstoplot))
for (i in 1:(ncol(enrichments)-1)){
  boxplots[[i]] <- ggplot(enrichments, aes_string(x="kmeans", y=colnames(enrichments)[i],fill="kmeans")) + 
    geom_boxplot(outlier.shape=NA) + pretty_plot() + 
    geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
    theme(axis.title.x = element_blank()) 
}
ggmatrix(boxplots,1,length(traitstoplot),
         xAxisLabels = traitstoplot,
         legend=c(1,1)) + 
  ggtitle("CMP Heterogeneity")+
  theme(plot.title = element_text(hjust = 0.5))

smoothed_traits <- colnames(allcells)[29:ncol(allcells)]
smoothed_traits <- sapply(traitstoplot, function(y){
  grep(y,smoothed_traits,value=TRUE)
  })
CMPs <- allcells[allcells$type %in% celltype,]
PC_clusters <- lapply(smoothed_traits,function(y) plotPCs(CMPs,y,celltype=celltype))

# Plot PC spatial clustering
ggmatrix(PC_clusters,nrow=1,ncol=length(traitstoplot),
         xAxisLabels = traitstoplot)+ 
  ggtitle("CMP Heterogeneity")+
  theme(plot.title = element_text(hjust = 0.5))

all <- c(boxplots,PC_clusters)
ggmatrix(all,nrow=2,ncol=length(traitstoplot),
         xAxisLabels = traitstoplot,
         legend=c(1,1))+ 
  ggtitle("CMP Heterogeneity")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.justification = c("right", "top"))

# MEP k-means clustering
meptraitstoplot <- c("HCT","PLT_COUNT")
# Get p-values for ATAC clustering differences in chromVAR score
ggallyplot_atac(allcells,celltype="MEP",traitstoplot=meptraitstoplot,smoothed=FALSE,numPCs=6)
#ggallyplot(allcells,celltype="MEP",traitstoplot=meptraitstoplot,smoothed=FALSE)

# Generate PC scatter + boxplot figure for MEP
numPCs=6
celltype <- "MEP"
traitstoplot <- c("HCT","PLT_COUNT")
plots <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=traitstoplot,
                       numPCs=numPCs)

ggmatrix(plots,nrow=2,ncol=length(traitstoplot),
         xAxisLabels = traitstoplot,
         legend=c(1,1))+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.justification = c("right", "top"))

ggmatrix(plots,nrow=2,ncol=length(traitstoplot),
         xAxisLabels = traitstoplot)

library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(plots[[1]])
p <- grid.arrange(legend, ncol=1, nrow=1)

################################################################################################################################
# Check if scATAC sub-populations stratify by TF z-scores
celltype <- "CMP"
traitstoplot <- c("RBC_COUNT","PLT_COUNT","MPV","MONO_COUNT")

# K-means cluster by ATAC PCs
enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:13]
enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:7]),2)$cluster)

# Filter data set for traits of interest
enrichments <- merge(enrichments,as.data.frame(allcells)[allcells$type %in% celltype,c(1,13:28)])
colnames(enrichments) <- remove_smooth(colnames(enrichments))
enrichments <- enrichments[,c("name",traitstoplot,"kmeans")]

# Merge with TF zscores
TF_zscores <- fread("../../data/singlecell/scATAC/tfDeviationsTable.tsv")
TFs_of_interest <- c("GATA1","KLF1","CEBPA","IRF8")

TFplots <- lapply(TFs_of_interest, function(TF){
  idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
  zscores <- TF_zscores %>% dplyr::select(cellnames,idx)
  colnames(zscores) <- c("cellnames",TF)
  merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
  p <- ggplot(merged, aes_string(x="kmeans", y=colnames(merged)[ncol(merged)],fill="kmeans")) +
    geom_boxplot(outlier.shape=NA) + pretty_plot() +
    geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
    theme(axis.title.x = element_blank())
  return(p)
})

# Plot boxplots for TFs of interest
ggmatrix(TFplots,1,length(TFs_of_interest),
         xAxisLabels = TFs_of_interest,
         legend=c(1,1)) + 
  theme(plot.title = element_text(hjust = 0.5))

ggmatrix(TFplots,1,length(TFs_of_interest),
         xAxisLabels = TFs_of_interest)

TFplots <- readRDS("../../data/singlecell/scATAC/CMP_kmeans_TFdifferences.rds")
# Calculate pvalues for all TFs between CMP sub-populations
# TFplots <- lapply(colnames(TF_zscores)[2:ncol(TF_zscores)], function(TF){
#   # idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
#   zscores <- TF_zscores %>% dplyr::select(cellnames,TF)
#   merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
#   p <- t.test(merged[merged$kmeans %in% 1,TF], merged[merged$kmeans %in% 2,TF])$p.value
#   return(p)
# })
#saveRDS(TFplots,"/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/scATAC/CMP_kmeans_TFdifferences.rds")

# Calculate p-value of TF zscores difference between clusters
TF_differences <- as.data.frame(unlist(TFplots))
TF_differences$TF <- colnames(TF_zscores)[2:ncol(TF_zscores)]
colnames(TF_differences) <- c("pval","TF")
TF_differences <- arrange(TF_differences,pval)
TF_differences$FDR <- p.adjust(TF_differences$pval)
TF_differences$rank <- seq(1,nrow(TF_differences),1)
TF_differences$TF_name <- str_split_fixed(TF_differences$TF, "_",n=4)[,3]
TF_differences$toLabel <- "F"
idx <- grep("GATA",TF_differences$TF_name)
idx <- idx[idx<30]
TF_differences[idx,"toLabel"] <- TF_differences[idx,"TF_name"]
TF_differences$highlight<- "F"
TF_differences[idx,"highlight"] <- "T"

# TF Rank order plot
ggplot(TF_differences,aes(x=rank,y=-log10(FDR))) + 
  geom_point(shape=21,size=3.5,aes(fill=highlight)) +
  pretty_plot() +
  scale_fill_manual(values = c("cyan","orange"),
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
    nudge_x=30,
    point.padding=0.5,
    min.segment.length=0,
    segment.alpha=0.3,
    direction="both")
+
# theme(legend.position="none") +



ggsave(p, file="CMP_TFs_rankorderplot.pdf")
