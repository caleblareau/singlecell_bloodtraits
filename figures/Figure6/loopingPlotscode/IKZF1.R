library(Gviz)
library(data.table)
library(GenomicRanges)
library(GenomicInteractions)
library(InteractionSet)
library(diffloop)
library(BuenColors)


source("geneinfoLoad.R") #hi-C interactions?
genome <- "hg19"

chr <- "chr7"
fromBP <- 50080000 
toBP   <- 50560000
bp <- 50343720
gene <- "IKZF1"

snps <- makeGRangesFromDataFrame(data.frame(chr = rep(chr,3),
                                            start = c(50187623, 50427982, 50497912),
                                            end = c(50187623, 50427982, 50497912)))

# Make GRange of region
g_region <- makeGRangesFromDataFrame(data.frame(chr = chr, start = fromBP, end = toBP))


# Get relevant peaks
snpsInRegion <- snps
snpsTrack <- AnnotationTrack(snpsInRegion, fill = c("black"))
geneLoci <- geneinfo[geneinfo$chromosome == chr & geneinfo$start > fromBP & geneinfo$end < toBP & geneinfo$symbol == gene,]

snp_track <- AnnotationTrack(padGRanges(snpsInRegion, pad = 1000), stacking = "dense", fill = c("#0081C9", "#8F1336", "#A65AC2"))
displayPars(snp_track) <- list( max.height = 25, stackHeight = 1, shape = "box")

# Build Interactions set
anchor.one <-  snps
anchor.two <-  makeGRangesFromDataFrame(data.frame(chr = chr, start = bp, end = bp))[rep(1,3)] 
interaction_counts<- c(5,5,5)
gi <- GenomicInteractions(anchor.one, anchor.two, counts=interaction_counts)
gi <- gi[mcols(gi)$counts > 0]
interaction_track <- InteractionTrack(gi, chromosome=chr)

displayPars(interaction_track) = list(col.interactions= c("#0081C9", "#8F1336", "#A65AC2"), 
                                      col.anchors.fill ="black",
                                      col.anchors.line = "black",
                                      interaction.dimension=100,
                                      anchor.height = 0,
                                      rotation = 0)

#availableDisplayPars(interaction_track)

itrack <- IdeogramTrack(genome = genome, chromosome = chr)

gtrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(geneLoci, genome = genome, chromosome = chr, name = "  ", transcriptAnnotation = "symbol", fill = "black")

pdf(file = paste0("../plots/",gene, ".loops.pdf"), width = 8, height = 4)

plotTracks(list(itrack, gtrack, interaction_track, snp_track, grtrack), from = fromBP, to = toBP,
           background.title = "white", sizes = c(0.05, 0.15, 0.2, 0.01, 0.05), innerMargin = 0, margin = 0)
dev.off()

