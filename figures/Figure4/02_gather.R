files <- list.files("../../data/resample", full.names = TRUE, pattern = "*.rds$")

gatheredList <- do.call(c, lapply(files, readRDS))