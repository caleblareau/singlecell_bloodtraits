library(BuenColors)
library(ggplot2)
library(reshape2)
source("boostrapSplineFunctions.R")

# Import observed 
observed <- read.table("../../data/singlecell/scATAC/weightedSingleCellScores-shuffled.txt", sep = "\t", header = TRUE)
observeddf <- observed[,c(1,2,13:28)]
colnames(observeddf) <- gsub("raw_", "", colnames(observeddf))

# Load pseudotime DFs
ptDF <- lapply(list.files("../../data/pseudotime", pattern= ".txt", full.names = TRUE), function(file){
  read.table(file, header = FALSE, stringsAsFactors = FALSE)
})
names(ptDF) <- sapply(list.files("../../data/pseudotime", pattern= ".txt"), function(file) strsplit(file, split = "_")[[1]][2])
melt_ptDF <- melt(ptDF)[,c(1,3,4)]
names(melt_ptDF) <- c("name", "pseudotime", "trajectory")

# Merge data frames
mergedf <-  merge(melt_ptDF, observeddf, by = c("name"))
meltdf <- melt(mergedf, id.vars = c("name", "pseudotime", "trajectory", "type"))

# One master plot
p1 <- ggplot(meltdf, aes(pseudotime, value, color = type)) +
  stat_smooth(se = FALSE, inherit.aes = FALSE, aes(pseudotime, value)) +
  geom_point() + facet_wrap(trajectory ~ variable, scales = "free", nrow = 4, ncol = 16) +
  pretty_plot() + theme(legend.position="bottom") + 
  labs(list(title = "GWAS in scHeme", x = "Pseudotime", y = "Zscore")) +
  scale_colour_manual(values = ejc_color_maps, name = "Cell Type") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill = "white") )
ggsave(p1, file = paste0("ps_plots/", "allPseudotimeDots.pdf"), width = 30, height = 10)

# One master plot
p1 <- ggplot(meltdf, aes(pseudotime, value, color = type)) +
  stat_smooth(se = TRUE, inherit.aes = FALSE, aes(pseudotime, value)) +
  facet_wrap(trajectory ~ variable, scales = "free", nrow = 4, ncol = 16) +
  pretty_plot() + theme(legend.position="bottom") + 
  labs(list(title = "GWAS in scHeme", x = "Pseudotime", y = "Zscore")) +
  scale_colour_manual(values = jdb_color_maps2, name = "Cell Type") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill = "white") )
ggsave(p1, file = paste0("ps_plots/", "allPseudotime.pdf"), width = 30, height = 10)


# Do bootstrap estimates to get CIs and then plot per trajectory with fixed axes for proper comparison
dummyOut <- lapply(unique(meltdf$trajectory), function(traj){
  sdf <- meltdf[meltdf$trajectory == traj, ]
  ldf <- lapply(unique(sdf$variable), function(gwas){
    print(paste0(traj, "_", gwas))
    spline.cis(sdf[sdf$variable == gwas, c(2,6)], B=1000,alpha=0.05)
  })
  names(ldf) <- unique(sdf$variable)
  meltdfTrait <- melt(ldf, id.vars = c("x", "lower.ci", "upper.ci", "main.curve"))
  
  p1 <- ggplot(meltdfTrait) + pretty_plot() + 
    labs(list(title = traj, x = "Pseudotime", y = "GWAS Z Score"))  +  facet_wrap(~ L1, nrow = 4, ncol = 4) + 
    geom_linerange(data = meltdfTrait, aes(x=x, ymin=lower.ci, ymax=upper.ci), inherit.aes = FALSE, alpha=0.03) + 
    scale_y_continuous(limits = c(-1, 2.5)) +
    geom_line(data = meltdfTrait, aes(x=x, y=main.curve), inherit.aes = FALSE) + theme(legend.position="none") +
    theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill = "white") )
  ggsave(p1, file = paste0("ps_plots/", traj,"FixedAxes_n1_p25.pdf"), width = 10, height = 10)
  traj
})
