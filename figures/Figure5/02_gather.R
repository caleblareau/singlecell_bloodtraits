library(magrittr)

# Gather the results from 01 once
if(FALSE){
  files <- list.files("../../data/resample", full.names = TRUE, pattern = "*.rds$")
  gatheredList <- do.call(c, lapply(files, readRDS))
  saveRDS(gatheredList, "simResults.10k.rds")
}

# Get observed variance
observed <- read.table("../../data/singlecell/scATAC/weightedSingleCellScores-shuffled.txt", sep = "\t", header = TRUE)
observeddf <- t(observed[,c(13:28)])
colnames(observeddf) <- observed[,2]
rownames(observeddf) <- gsub("raw_", "", rownames(observeddf))
obs <- reshape2::melt(observeddf)
obs <- aggregate(. ~ Var1 + Var2, obs, sd)

# Get permuted variance mean and standard deviation
variancePermutedList <- readRDS("simResults.10k.rds")
names(variancePermutedList) <- paste0("p", as.character(1:length(variancePermutedList)))
means <- aggregate(. ~ Var1 + Var2, reshape2::melt(variancePermutedList)[,c(1,2,3)], mean)
sds   <-  aggregate(. ~ Var1 + Var2, reshape2::melt(variancePermutedList)[,c(1,2,3)], sd )

colnames(obs) <- c("TRAIT", "CELL", "OBS")
colnames(means) <- c("CELL", "TRAIT", "MEAN")
colnames(sds) <- c("CELL", "TRAIT", "SD")

mdf <- merge(means, sds) %>% merge(obs)
mdf$Z <- (mdf$OBS- mdf$MEAN)/mdf$SD
df <- mdf[order(mdf$Z, decreasing = TRUE),]
  