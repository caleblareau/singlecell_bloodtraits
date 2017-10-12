library(magrittr)
library(data.table)

makePvalueTable <- function(celltype){
  df <- read.table(paste0("data/",celltype,"_weightedchromVAReach.txt"), sep = "\t", header = TRUE)
  
  traits <- colnames(df)[13:28]
  
  df <- lapply(traits, function(trait){
    high <- df[,trait] > 0
    highidx <- which(high)
    
    distTwo <- function(df, idx){
      PCs1 <- colMeans(df[idx, paste0("PC", as.character(1:5))]) %>% unname()
      PCs2 <- colMeans(df[!(1:dim(df)[1] %in% idx), paste0("PC", as.character(1:5))]) %>% unname()
      dist <- (PCs1 - PCs2)^2 %>% sum() %>% sqrt()
      return(dist)
    }
    distObs <- distTwo(df, highidx)
    distPermutedVector <- sapply(1:100, function(i){
      distTwo(df, sample(1:dim(df)[1], length(highidx), replace = FALSE, prob = NULL))
    })
    Z <- (distObs - mean(distPermutedVector))/sd(distPermutedVector)
    p.value <- 1-pnorm(Z)
    data.frame(celltype, trait, Z, p.value, nhigh = sum(high), nlow = sum(!high))
  }) %>% rbindlist() %>% as.data.frame()
  return(df)
}

makePvalueTable("HSC")
makePvalueTable("CMP")
makePvalueTable("MEP")