library(data.table)
library(dplyr)
library(rtracklayer)
library(preprocessCore)
library(ggjoy)
library(BuenColors)
library(regioneR)
library(stringr)
"%ni%" <- Negate("%in%")

#Make fine-mapped df and gr
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)
CS.df$seqnames <- paste0("chr",CS.df$seqnames)

CS.gr <- makeGRangesFromDataFrame(CS.df, keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)

# UKBB exclusion list
exdf <- read.table("figures/revisions/exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
CS.gr <- CS.gr[CS.gr$var %ni% exdf,]

#' Read in annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
coding.gr <- import("../../data/annotations/Coding_UCSC.bed",format="bed")
intron.gr <- import("../../data/annotations/Intron_UCSC.bed",format="bed")
promoter.gr <- import("../../data/annotations/Promoter_UCSC.fixed.bed",format="bed")
UTR.gr <- import("../../data/annotations/UTR_3_UCSC.bed",format="bed")

# Do overlaps
gr_t <- unique(CS.gr)
ov_1 <- findOverlaps(gr_t, coding.gr)
ov_2 <- findOverlaps(gr_t, promoter.gr)
ov_3 <- findOverlaps(gr_t, UTR.gr)
ov_4 <- findOverlaps(gr_t, intron.gr)

# Classify each variant
class <- ifelse(1:length(gr_t) %in% queryHits(ov_1), "coding",
                ifelse(1:length(gr_t) %in% queryHits(ov_2), "promoter",
                       ifelse(1:length(gr_t) %in% queryHits(ov_3), "utr",
                              ifelse(1:length(gr_t) %in% queryHits(ov_4), "intron", "intergenic"))))

# Colormap
annotationColors <- jdb_palette("brewer_spectra")[c(1,3,4,6,7)]
names(annotationColors) <- c("coding", "promoter", "utr", "intergenic", "intron")

# Pie Chart of fine-mapped variants with Percentages
slices <- table(class) %>% as.data.frame() %>% dplyr::select(Freq) %>% unlist %>%  as.integer()
lbls <- table(class) %>% as.data.frame() %>% dplyr::select(class) %>% unlist %>%  as.character()
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # add % to labels 
pie(slices,labels = lbls, col=annotationColors)

# Read VEP table ----------------------------------------------------------
vep <- data.frame(fread("VEP/coding_variants_all_VEPoutput.txt"))
colnames(vep)[1] <- "Uploaded_variation"

vep$Uploaded_variation <- gsub("/","_",vep$Uploaded_variation) %>% str_replace(.,"_",":") 
vep_nodups <- vep[!duplicated(vep$Uploaded_variation),]
vep_nodups %>% dplyr::select(Uploaded_variation,SYMBOL,Consequence,EXON,INTRON,Existing_variation) %>% 
  right_join(.,CS.df,by=c("Uploaded_variation"="var")) -> CS.df_vepmerged
colnames(CS.df_vepmerged)[1] <- "var"
CS.df_vepmerged[is.na(CS.df_vepmerged)] <- "-"
CS.gr <- makeGRangesFromDataFrame(CS.df_vepmerged, keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)

# Coding consequences
coding_consequences <- c("missense_variant","synonymous_variant","frameshift_variant",
                         "inframe_insertion","stop_gained","stop_retained_variant",
                         "start_lost","stop_lost","coding_sequence_variant","incomplete_terminal_codon_variant")

all_consequences <- vep$Consequence %>% stringr::str_split(.,",") %>% unlist() 
coding_con <- all_consequences[all_consequences %in% coding_consequences] %>% table() %>% as.data.frame()
colnames(coding_con)[1] <- "Categorie"
coding_con$Categorie <- paste0(coding_con$Categorie,": ",round(100*coding_con$Freq/(sum(coding_con$Freq)),2),"%")

par(mar=c(5,0,4,2))
pie(coding_con$Freq, labels=NA, col=jdb_palette("brewer_spectra"),bty='L')
legend(.9, .5, as.character(coding_con$Categorie), cex=0.7, fill=jdb_palette("brewer_spectra"))

# Plotly version
# # Plot piechart and potentially export as PDF
# coding_con %>% plot_ly(labels = ~Categorie, values = ~Freq,
#                        textposition='inside',
#                        textinfo = 'percent',
#                        insidetextfont = list(color = '#FFFFFF')) %>%
#   add_pie(hole = 0) %>%
#   layout(margin=m,showlegend = T, legend=list(x=1,y=0.9),
#          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
#          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F)) #%>% export("coding_subtypes.pdf")
