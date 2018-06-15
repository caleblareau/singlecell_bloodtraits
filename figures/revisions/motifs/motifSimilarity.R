library(PWMEnrich)
library(motifbreakR)
library(data.table)
library(dplyr)
library(BuenColors)

data(hocomoco)

lapply(2:length(hocomoco@listData), function(i){
  print(i)
  lapply(1:(i-1), function(j){
    
    # Set up integer weight matrice
    mot1 <- round(hocomoco@listData[[i]]*1000); storage.mode(mot1) <- "integer"
    mot2 <- round(hocomoco@listData[[j]]*1000); storage.mode(mot2) <- "integer"
    val <- motifSimilarity(mot1, mot2, trim = 0.4, self.sim = FALSE)
    
    data.frame(motif1 = names(hocomoco@listData)[i], motif2 = names(hocomoco@listData)[j], cor = val)
  }) %>% rbindlist() %>% data.frame(stringsAsFactors = FALSE) -> df
  df
}) %>% rbindlist() %>% data.frame(stringsAsFactors = FALSE) -> df

length(unique(c(df$motif1, df$motif2)))
saveRDS(df, file = "hocomoco_similarity.rds")
