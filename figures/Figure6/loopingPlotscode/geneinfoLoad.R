library(BuenColors)

# Deal with exon annotation
load("hg19-geneinfo-exon.rda")

geneinfo$chrom <- paste0("chr", as.character(geneinfo$chrom))
geneinfo$strand <- ifelse(geneinfo$strand == -1, "-", "+")
colnames(geneinfo) <- c("chromosome", "start", "end", "symbol", "score", "strand")
geneinfo$gene <- geneinfo$symbol
geneinfo$transcript <- geneinfo$symbol 

map2color<-function(x, pal = as.character(jdb_palette("solar_blues")),limits=c(0,1)){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

