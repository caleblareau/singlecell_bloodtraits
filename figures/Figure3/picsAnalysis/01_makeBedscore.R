library(readxl)
p <- read_xls("masterfile9.xls")

picsnew <- data.frame(
  chr = p$chr, 
  start = p$pos,
  end = p$pos + 1,
  snp = p$SNP,
  PP = p$PICS_probability, stringsAsFactors = FALSE
)

lapply(unique(p$Disease), function(d){
  write.table(picsnew[p$Disease == d, ], file = paste0("pics_bedscore/", d, ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
})
