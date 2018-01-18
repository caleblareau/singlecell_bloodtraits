library(data.table)
library(dplyr)
library(BuenColors)

icd <- data.frame(fread(paste0('zcat < ', "../../data/UKBB_BC_PP001/ICD10all_PP5_summaryTable.txt.gz")))

erythroid.traits <- c("HCT","HGB","MCH","MCHC","MCV","MEAN_RETIC_VOL","RBC_COUNT","RETIC_COUNT","MPV","PLT_COUNT")
myeloid.traits <- c("EO_COUNT","BASO_COUNT","NEUTRO_COUNT","MONO_COUNT")
lymph.traits <- c("LYMPH_COUNT", "WBC_COUNT")

getDF <- function(v, lin) {
  lapply(v, function(trait){
    tdf <- icd[icd$trait == trait, ]
    odf <- tdf[tdf$pval < 0.05/(dim(tdf)[1]),]
    odf
  }) %>% rbindlist() %>% as.data.frame() -> odf1
  
  df<- data.frame(
    counts = sort(table(unique(odf1[,c("chr", "pos", "V2")])$V2), decreasing = TRUE),
    lineage = lin
  )
  
  return(df)
}

ery <- getDF(erythroid.traits, "Megakaryocyte / \nErythroid")
myeloid <- getDF(myeloid.traits, "Myeloid")
lymph <- getDF(lymph.traits, "Lymphoid")
all <- getDF(c(erythroid.traits, myeloid.traits, lymph.traits), "ALL")

df <- rbind(ery,myeloid,lymph)
df$label <- "Other"
keep <- c(1:7,9,10,11,14,23) # check order of `all` for which ones to keep
df$label <- gsub("Diagnoses - main ICD10: ", "", df$counts.Var1)
df <- df[df$counts.Var1 %in% all[keep,1],]

# Manually override pleiotropic // UGLY
head(icd[icd$ICD10 == "I25",], 20) # 2 pleiotropic
head(icd[icd$ICD10 == "I21",])# 1 pleiotropic
df$counts.Freq <- ifelse(df$label == "I25 Chronic ischaemic heart disease", df$counts.Freq -2, df$counts.Freq)
df$counts.Freq <- ifelse(df$label == "I21 Acute myocardial infarction", df$counts.Freq -1, df$counts.Freq)

pdf_heart <- data.frame(
  counts.Var1 = "Diagnoses - main ICD10: I25 Chronic ischaemic heart disease",
  counts.Freq = 2, 
  lineage = "Multi-lineage",
  label = "I25 Chronic ischaemic heart disease"
)

pdf_infarc <- data.frame(
  counts.Var1 = "Diagnoses - main ICD10: I21 Acute myocardial infarction",
  counts.Freq = 1, 
  lineage = "Multi-lineage",
  label = "I21 Acute myocardial infarction"
)

finaldf <- rbind(df, pdf_heart, pdf_infarc)
sumsort <- sort(tapply(finaldf$counts.Freq, list(finaldf$label), sum))

finaldf$label <- factor(as.character(finaldf$label), (names(sumsort)))

p1 <- ggplot(finaldf, aes(x = label, y= counts.Freq, fill = lineage)) +
  geom_bar(width = 0.7, colour="black", stat = "identity", position = position_stack(reverse = FALSE)) +
  labs(x = "ICD 10 Code", y = "# PheWAS Variants Identified from PP > 0.5 Blood Variants", fill = "")+
  pretty_plot() +coord_flip() + scale_fill_manual(values = alpha(c("#8F1336","#FF5A00", "#0081C9", "#46A040"), 1)) + 
  theme(legend.position = "bottom")
ggsave(p1, file = "ICD10_phewas_hightlights.pdf",  useDingbats=FALSE)
