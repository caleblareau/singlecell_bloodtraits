library(Gviz)
library(data.table)
library(GenomicRanges)
library(GenomicInteractions)
library(InteractionSet)
library(diffloop)
library(BuenColors)


source("geneinfoLoad.R")
genome <- "hg19"

chr <- "chr9"
fromBP <- 34100000 
toBP   <- 34720000
bp <- 34650699
gene <- "IL11RA"

snps <- makeGRangesFromDataFrame(data.frame(chr = rep(chr,1),
                                            start = c(34195572),
                                            end = c(34195572)))

# Make GRange of region
g_region <- makeGRangesFromDataFrame(data.frame(chr = chr, start = fromBP, end = toBP))


# Get relevant peaks
snpsInRegion <- snps
snpsTrack <- AnnotationTrack(snpsInRegion, fill = c("black"))
geneLoci <- geneinfo[geneinfo$chromosome == chr & geneinfo$start > fromBP & geneinfo$end < toBP,]

snp_track <- AnnotationTrack(padGRanges(snpsInRegion, pad = 1000), stacking = "dense", fill = c("firebrick"))
displayPars(snp_track) <- list( max.height = 25, stackHeight = 1, shape = "box")

# Build Interactions set
anchor.one <-  snps
anchor.two <-  makeGRangesFromDataFrame(data.frame(chr = chr, start = bp, end = bp))
interaction_counts<- 1
gi <- GenomicInteractions(anchor.one, anchor.two, counts=interaction_counts)
gi <- gi[mcols(gi)$counts > 0]
interaction_track <- InteractionTrack(gi, chromosome=chr)

displayPars(interaction_track) = list(col.interactions= c("firebrick"), 
                                      col.anchors.fill ="black",
                                      col.anchors.line = "black",
                                      interaction.dimension=100,
                                      anchor.height = 0,
                                      rotation = 0)

#availableDisplayPars(interaction_track)

itrack <- IdeogramTrack(genome = genome, chromosome = chr)

gtrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(geneLoci, genome = genome, chromosome = chr, name = "  ", transcriptAnnotation = "symbol", fill = "black")

pdf(file = paste0("../plots/",gene, ".loops.pdf"), width = 8, height = 6)

plotTracks(list(itrack, gtrack, interaction_track, snp_track, grtrack), from = fromBP, to = toBP,
           background.title = "white", sizes = c(0.05, 0.15, 0.2, 0.01, .7), innerMargin = 0, margin = 0)
dev.off()

