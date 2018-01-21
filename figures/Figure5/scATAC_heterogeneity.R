library(GGally)
library(corrplot)
library(BuenColors)
jdb_color_maps2 <- c(jdb_color_maps, "Mega" = "#FF347E", "UNK" = "#8D91AD", "mDC"= "#FFD700", "MCP" = "#C390D4")
names(jdb_color_maps2)[9] <- "Mono"
library(ggrepel)
library(SummarizedExperiment)
library(ggbeeswarm)
library(irlba)
library(data.table)
library(dplyr)
library(stringr)
library(qvalue)
library(gridExtra)
library(fpc)
"%ni%" <- Negate("%in%")
######################################################################################################
# Functions
# k-means on enrichments
ggallyplot <- function(allcells,celltype,traitstoplot,smoothed=FALSE,graph=TRUE,colors=c("Mega","Ery"),kmeans=TRUE){
  set.seed(42)
  if(smoothed) enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,29:ncol(allcells)]
  else {
    enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:28]
  }
  
  colnames(enrichments) <- remove_smooth(colnames(enrichments))
  rownames(enrichments) <- enrichments$name
  enrichments <- enrichments[,c(traitstoplot)]
  colnames(enrichments) <- abbreviate(colnames(enrichments))
  if (kmeans){
    enrichments$kmeans <- as.character(kmeans(scale(as.matrix(enrichments)),centers=2,nstart = 1000,iter.max = 1000)$cluster)
  } else{
    res <- pamk(scale(enrichments),krange=2,criterion="multiasw")
    enrichments$kmeans <- res$pamobject$clustering
    enrichments$kmeans <- sapply(enrichments$kmeans, function(x) {switch(x,1,2)})
    enrichments$kmeans <- as.factor(enrichments$kmeans)
  }
  
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
        scale_fill_manual(values=as.character(jdb_color_maps2[colors])) +
        geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
        guides(fill=guide_legend(title=paste("g-ChromVAR","clustering",sep="\n"))) + 
        theme(axis.title.x = element_blank()) 
    }
    return(boxplots)
  } else{
    return(enrichments)
  }
  
}

# kmeans on ATAC PCs
ggallyplot_atac <- function(allcells,celltype,traitstoplot,smoothed=FALSE,graph=FALSE,numPCs=5,kmeans=TRUE){
  set.seed(42)
  enrichments <- as.data.frame(allcells)[allcells$type %in% celltype,1:12]
  if (kmeans==TRUE){
    enrichments$kmeans <- as.character(kmeans(scale(enrichments[,3:(3+numPCs-1)]),centers=2,nstart = 500, iter.max=1000)$cluster)
  } else{
    res <- pamk(scale(enrichments[,3:(3+numPCs-1)]),krange=2,criterion="multiasw")
    enrichments$kmeans <- res$pamobject$clustering
    enrichments$kmeans <- sapply(enrichments$kmeans, function(x) {switch(x,1,2)})
    enrichments$kmeans <- as.factor(enrichments$kmeans)
  }
  
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
    subenrichments[[i]]$kmeans <- as.character(kmeans(scale(subenrichments[[i]][,3:(3+numPCs-1)]),centers=2,nstart = 100)$cluster)
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
    # # Plot PCs
    # smoothed_traits <- colnames(allcells)[29:ncol(allcells)]
    # smoothed_traits <- sapply(traitstoplot, function(y){
    #   grep(y,smoothed_traits,value=TRUE)
    # })
    # ct <- allcells[allcells$type %in% celltype,]
    # PC_clusters <- lapply(smoothed_traits,function(y) plotPCs(ct,y,celltype=celltype))
    # 
    # # combine box plots and PC scatter plots and plot them together
    # all <- c(boxplots,PC_clusters)
    # return(all)
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
#setwd("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits/figures/Figure5/")

#CMPs
allcells <- fread("../../data/singlecell/scATAC/sc_traitenrichments_11aug.txt")
# Use 100-iteration average scATAC
avg <- fread("../../data/singlecell/scATAC/scGWASenrichments_average.tsv")
avg <- dcast(avg, V1~V2)
allcells <- merge(allcells[,1:12], avg, by.x="name",by.y="V1")

# CMP k-means clustering
traitstoplot <- c("RBC_COUNT","MPV","PLT_COUNT","MONO_COUNT")

# Cluster and plot cluster means 
numPCs=5
ggallyplot_atac(allcells,numPCs = numPCs, celltype="CMP",traitstoplot=traitstoplot,smoothed=FALSE,kmeans=FALSE)
plots <- compare_subgroups_plot(allcells,celltype="CMP",traitstoplot=traitstoplot,
                                numPCs=5,graph=TRUE,colors=c("Ery","Mono"),kmeans=F)
CMP_boxplots <- ggmatrix(plots[1:4],nrow=1,ncol=length(traitstoplot),
                         xAxisLabels = traitstoplot)
ggsave(CMP_boxplots, file="5Ca_CMP_kmedoids_REVISED.pdf",
       width=6,height=3, useDingbats=F)

# Extract legend and save separately
legend <- g_legend(plots[[1]])
grid.arrange(legend, ncol=1, nrow=1)

# Overlay PC2 and PC3 plot with k-medoids results
enrichments <- compare_subgroups_plot(allcells,celltype="CMP",traitstoplot=traitstoplot,
                                      numPCs=5,graph=FALSE,kmeans=FALSE)
others <- allcells[allcells$type %ni% "CMP",] %>% select(name,paste0(traitstoplot,"_PP001"))
colnames(others) <- remove_smooth(colnames(others))
others$kmeans <- 3
others$kmeans <- as.factor(others$kmeans)
combined <- bind_rows(enrichments,others)

pcplot <- merge(combined[,c("name","kmeans")], allcells[,1:12], by="name")

cmp_pcs <- ggplot(pcplot,aes(PC2,PC3)) + geom_point(aes(colour=kmeans)) +
  scale_color_manual(values=as.character(c(jdb_color_maps2[c("Ery","Mono")],"gray"))) +
  pretty_plot() + theme_bw() + 
  theme(plot.title = element_text(size=14,hjust = 0.5,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  labs(x="",y="")
ggsave(cmp_pcs, file="CMP_pamk_pcoverlay.pdf",
       width=5,height=4, useDingbats=F)

# CMP g-ChromVAR clustering
traitstoplot <- c("RBC_COUNT","PLT_COUNT","MPV","MONO_COUNT")
cmp_chromvar <- ggallyplot(allcells,celltype="CMP",traitstoplot=traitstoplot,smoothed=FALSE,colors=c("Mono","Ery"),kmeans=F)
cmp_chromvar_plots <- ggmatrix(cmp_chromvar,1,length(traitstoplot),
                               xAxisLabels = traitstoplot) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(cmp_chromvar_plots, file="CMP_gChromVAR_Clustering.pdf",
       width=6,height=3,useDingbats=F)

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

ggsave(CMP_TFs, file="5Cb_CMP_TFzscores_REVISED.pdf",width=6,height=3, useDingbats=F)

# CMP TF rank order plot
TF_differences <- readRDS("../../data/singlecell/scATAC/CMP_ATAC_TFzscore_differences_kmedoids.rds")
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
    nudge_y =5,
    nudge_x=100,
    point.padding=0.50,
    min.segment.length=0,
    segment.alpha=0.3)

ggsave(p, file="5D_CMP_TFs_rankorderplot.pdf",
       width=6,height=6, useDingbats=F)

####################################################################################################  
# MEP clustering
celltype <- "MEP"
meptraits=c("RBC_COUNT", "HCT","PLT_COUNT") 
numPCs=7
# K-medoids cluster by ATAC PCs
ggallyplot_atac(allcells,kmeans=F,celltype=celltype,traitstoplot=meptraits,smoothed=FALSE,numPCs = numPCs)
enrichments <- compare_subgroups_plot(allcells,celltype=celltype,traitstoplot=meptraits,
                                      numPCs=numPCs,graph=FALSE,kmeans=F)
plots  <- compare_subgroups_plot(allcells,celltype=celltype,traitstoplot=meptraits,
                                 numPCs=numPCs,graph=TRUE,colors=c("Ery","Mega"),kmeans=F)

mep_atac_plots <- ggmatrix(plots,1,length(meptraits),xAxisLabels = meptraits)
ggsave(mep_atac_plots, file="5Ea_MEP_ATAC_Clustering_4traits_REVISED.pdf",
       width=6,height=3, useDingbats=F)
legend <- g_legend(plots[[1]])
grid.arrange(legend, ncol=1, nrow=1)

# Cluster MEP population by g-chromVAR
celltype <- "MEP"
mep_chromvar <- ggallyplot(allcells,kmeans=F,celltype="MEP",traitstoplot=meptraits,smoothed=FALSE,colors=c("Ery","Mega"))
mep_chromvar_plots <- ggmatrix(mep_chromvar,1,length(meptraits),
                               xAxisLabels = meptraits) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(mep_chromvar_plots, file="5Ea_MEP_gChromVAR_Clustering_REVISED.pdf",
       width=6,height=3, useDingbats=F)

# MEP color PC2 and PC3 plot by k-medoids clustering
enrichments <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=meptraits,
                                      numPCs=numPCs,graph=FALSE,kmeans=F,colors=c("Ery","Mega"))
# enrichments <- ggallyplot(allcells,kmeans=F,celltype=celltype,traitstoplot=meptraits,smoothed=FALSE,colors = c("Mega","Ery"),graph=FALSE)
# enrichments$name <- rownames(enrichments)
others <- allcells[allcells$type %ni% "MEP",] %>% select(name,paste0(traitstoplot,"_PP001"))
colnames(others) <- remove_smooth(colnames(others))
others$kmeans <- 3
others$kmeans <- as.factor(others$kmeans)
combined <- bind_rows(enrichments,others)

pcplot <- merge(combined[,c("name","kmeans")], allcells[,1:12], by="name")

MEP_pamk_pcoverlay <- ggplot(pcplot,aes(PC2,PC3)) + geom_point(aes(colour=kmeans)) +
  scale_color_manual(values=as.character(c(jdb_color_maps2[c("Ery","Mega")],"gray"))) +
  pretty_plot() + theme_bw() + 
  theme(plot.title = element_text(size=14,hjust = 0.5,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  labs(x="",y="")

ggsave(MEP_pamk_pcoverlay, file="MEP_pamk_pcoverlay.pdf",
       width=5,height=4, useDingbats=F)

# MEP TFs
# K-medoids cluster by g-ChromVAR z-scores
enrichments <- compare_subgroups_plot(allcells,celltype="MEP",traitstoplot=meptraits,
                                      numPCs=numPCs,graph=FALSE,kmeans=F)

# enrichments <- ggallyplot(allcells,kmeans=F,celltype=celltype,traitstoplot=meptraits,
#                           smoothed=FALSE,colors = c("Ery","Mega"),graph=F)
# enrichments$name <- rownames(enrichments)

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

ggsave(mep_TFs, file="5Eb_MEP_TFzscores_REVISED.pdf",
       width=6,height=3, useDingbats=F)

# MEP TF rank order plot
TF_differences <- readRDS("../../data/singlecell/scATAC/MEP_ATAC_TFzscore_differences_kmedoids.rds")

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
  scale_fill_manual(values = as.character(jdb_color_maps2[c("Mega","Ery")]),
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

ggsave(p, file="5F_MEP_TFs_rankorderplot.pdf",width=6,height=6, useDingbats=F)

######################################################################################################
# # k-mers clustering
# allcells <- fread("../../data/singlecell/scATAC/sc_traitenrichments_11aug.txt")
# 
# kmers <- readRDS("/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/scATAC/kmer.6mer.dev.rds")
# kmer_mat <- t(assays(kmers)[["z"]])
# rownames(kmer_mat) <- colData(kmers)[,1]
# 
# kmer.pca <- prcomp(kmer_mat,center = TRUE,scale. = TRUE)
# # Take the first 10 PCs
# comp <- data.frame(kmer.pca$x[,1:10])
# comp$name <- rownames(comp)
# 
# kmers <- allcells[,1:2]
# kmers <- merge(kmers,comp,by="name")
# kmers <- merge(kmers,allcells[,-c(2:12)],by="name")
# 
# #saveRDS(kmers,"/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/scATAC/kmer_pcs.rds")
# 
# # CMP k-means clustering
# traitstoplot <- c("RBC_COUNT","PLT_COUNT","MPV","MONO_COUNT")
# 
# # Cluster by k-mers and plot cluster means 
# ggallyplot_atac(kmers,kmeans=F,celltype="CMP",traitstoplot=traitstoplot,smoothed=FALSE,numPCs = 5)
# plots <- compare_subgroups_plot(kmers,celltype="CMP",traitstoplot=traitstoplot,
#                                 numPCs=5,kmeans=F,graph=TRUE)
# ggmatrix(plots[1:4],nrow=1,ncol=length(traitstoplot),
#                          xAxisLabels = traitstoplot)
# ggsave(CMP_boxplots, file="5Ca_CMP_kmeans.pdf",
#        width=6,height=3, useDingbats=F)
# 
# # Extract dataframe with annotation for k-means cluster by ATAC PCs
# enrichments2 <- compare_subgroups_plot(kmers,celltype="CMP",traitstoplot=traitstoplot,
#                                       numPCs=5,graph=FALSE,kmeans=F)
# 
# # Merge with TF zscores
# TF_zscores <- fread("../../data/singlecell/scATAC/tfDeviationsTable.tsv")
# #CMP TFs
# TFs_of_interest <- c("GATA1","KLF1","CEBPA","IRF8")
# 
# # For the TFs of interest, extract their z-scores for all single cells of a cell type and plot k-means cluster vs. z-score
# TFplots <- lapply(TFs_of_interest, function(TF){
#   idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
#   zscores <- TF_zscores %>% dplyr::select(cellnames,idx)
#   colnames(zscores) <- c("cellnames",TF)
#   merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
#   p <- ggplot(merged, aes_string(x="kmeans", y=colnames(merged)[ncol(merged)],fill="kmeans")) + 
#     geom_boxplot(outlier.shape=NA) + pretty_plot() + 
#     scale_fill_manual(values=as.character(jdb_color_maps2[c("Mono","Ery")])) +
#     geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
#     
#     theme(axis.title.x = element_blank())
#   return(p)
# })
# 
# # Plot differential z-scores for TFs of interest
# TFs <- ggmatrix(TFplots,1,length(TFs_of_interest),
#                 xAxisLabels = TFs_of_interest) + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# # MEPs
# # K-means cluster by ATAC PCs
# meptraits=c("HGB","HCT","PLT_COUNT")
# celltype="MEP"
# numPCs=5
# ggallyplot_atac(kmers,kmeans=F,celltype=celltype,traitstoplot=meptraits,smoothed=FALSE,numPCs = numPCs)
# enrichments <- compare_subgroups_plot(kmers,celltype="MEP",traitstoplot=meptraits,
#                                       numPCs=numPCs,graph=FALSE,kmeans=F)
# plots  <- compare_subgroups_plot(kmers,celltype="MEP",traitstoplot=meptraits,
#                                  numPCs=numPCs,graph=TRUE,colors=c("Mega","Ery"),kmeans=F)
# ggmatrix(plots,1,length(meptraits),xAxisLabels = meptraits)
# ggsave(mep_atac_plots, file="MEP_ATAC_Clustering_4traits.pdf",
#        width=6,height=3, useDingbats=F)
# 
# TFs_of_interest <- c("GATA1","KLF1","MEF2C","FLI1")
# 
# enrichments <- compare_subgroups_plot(kmers,celltype="MEP",traitstoplot=meptraits,
#                                       numPCs=numPCs,graph=FALSE,kmeans=F)
# 
# # For the TFs of interest, extract their z-scores for all single cells of a cell type and plot k-means cluster vs. z-score
# TFplots <- lapply(TFs_of_interest, function(TF){
#   idx <- grep(TF,colnames(TF_zscores),value=TRUE)[1]
#   zscores <- TF_zscores %>% dplyr::select(cellnames,idx)
#   colnames(zscores) <- c("cellnames",TF)
#   merged <- merge(enrichments,zscores,by.x="name",by.y="cellnames")
#   p <- ggplot(merged, aes_string(x="kmeans", y=colnames(merged)[ncol(merged)],fill="kmeans")) + 
#     geom_boxplot(outlier.shape=NA) + pretty_plot() + 
#     scale_fill_manual(values=as.character(jdb_color_maps2[c("Mega","Ery")])) +
#     geom_quasirandom(varwidth = TRUE,alpha=0.35,size=0.5)+
#     
#     theme(axis.title.x = element_blank())
#   return(p)
# })
# 
# # Plot differential z-scores for TFs of interest
# # mep_TFs <- 
# ggmatrix(TFplots,1,length(TFs_of_interest),
#          xAxisLabels = TFs_of_interest) + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggsave(mep_TFs, file="5Eb_MEP_TFzscores.pdf",
#        width=6,height=3, useDingbats=F)
# 
# ###########################################################################################
# # RNA
# 
# genes <- fread("../../data/singlecell/imputedRNA/RNA_smoothed.genes.txt")
# RNA <- fread("../../data/singlecell/imputedRNA/RNA_smoothed.txt")
