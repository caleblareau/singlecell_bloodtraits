library(data.table)
library(ggplot2)
library(BuenColors)
library(plotly)
library(dplyr)
library(scales)
library(webshot)

### Moderating alleles
trait <- "RBC_COUNT"
unit <- "10^6 cells/uL"
rs_alleles <- readRDS("CCND3.haplotypes.rds")

# Only look at the imputed genotypes with 0, 1, or 2 for boxplots
rs_alleles_onlyabsolutes <- as.data.frame(rs_alleles[rs_alleles$alleleA.rs112233623 == 0 |
                                                       rs_alleles$alleleA.rs112233623 == 1 | 
                                                       rs_alleles$alleleA.rs112233623 == 2,])

# Write out genotypes for rs1
allelecol1 <- "alleleA.rs2408955"
rs_alleles_onlyabsolutes$rs1_A <- "AA"
rs_alleles_onlyabsolutes[grep(0,rs_alleles_onlyabsolutes[,allelecol1]),"rs1_A"] <- "AA"
rs_alleles_onlyabsolutes[grep(1,rs_alleles_onlyabsolutes[,allelecol1]),"rs1_A"] <- "AG"
rs_alleles_onlyabsolutes[grep(2,rs_alleles_onlyabsolutes[,allelecol1]),"rs1_A"] <- "GG"
rs_alleles_onlyabsolutes$rs1_A <- factor(rs_alleles_onlyabsolutes$rs1_A, levels = c("GG","AG","AA"))

# Write out genotypes for rs2
allelecol2 <- "alleleA.rs112233623"
rs_alleles_onlyabsolutes$rs2_A <- "TT"
rs_alleles_onlyabsolutes[grep(0,rs_alleles_onlyabsolutes[,allelecol2]),"rs2_A"] <- "TT"
rs_alleles_onlyabsolutes[grep(1,rs_alleles_onlyabsolutes[,allelecol2]),"rs2_A"] <- "CT"
rs_alleles_onlyabsolutes[grep(2,rs_alleles_onlyabsolutes[,allelecol2]),"rs2_A"] <- "CC"
rs_alleles_onlyabsolutes$rs2_A <- factor(rs_alleles_onlyabsolutes$rs2_A, levels = c("CC","CT","TT"))

# Make column for all haplotype combinations
rs_alleles_onlyabsolutes$con_rs1 <- ifelse(rs_alleles_onlyabsolutes$rs1_A =="GG" & rs_alleles_onlyabsolutes$rs2_A =="CT", "GG/CT", "")
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="AA" & rs_alleles_onlyabsolutes$rs2_A=="CT","con_rs1"] <- "AA/CT"
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="AG" & rs_alleles_onlyabsolutes$rs2_A=="CT","con_rs1"] <- "AG/CT"
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="GG" & rs_alleles_onlyabsolutes$rs2_A=="CT","con_rs1"] <- "GG/CT"
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="AA" & rs_alleles_onlyabsolutes$rs2_A=="CC","con_rs1"] <- "AA/CC"
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="AG" & rs_alleles_onlyabsolutes$rs2_A=="CC","con_rs1"] <- "AG/CC"
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="GG" & rs_alleles_onlyabsolutes$rs2_A=="CC","con_rs1"] <- "GG/CC"
rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes$rs1_A=="AA" & rs_alleles_onlyabsolutes$rs2_A=="TT","con_rs1"] <- "AA/TT"
rs_alleles_onlyabsolutes$con_rs1 <- factor(rs_alleles_onlyabsolutes$con_rs1, 
                                           levels = c("GG/CC","AG/CC","AA/CC",
                                                      "GG/CT","AG/CT","AA/CT","AA/TT"))

# The baseline is the homozygous ref/ref haplotypes.
baseline <- as.data.frame(rs_alleles_onlyabsolutes[rs_alleles_onlyabsolutes[,"con_rs1"] %in% "GG/CC",trait])[,1] %>% median()
# Calculate median phenotype value for each haplotype. 
meds <- rs_alleles_onlyabsolutes %>% group_by(con_rs1) %>% summarise(median = median(RBC_COUNT, na.rm = TRUE))
meds <- meds[1:5,]

ggplot(subset(rs_alleles_onlyabsolutes,rs2_A=="CC" | rs2_A=="CT"), 
              aes_string(x="con_rs1", trait)) + 
  geom_boxplot() + pretty_plot() +
  geom_boxplot(outlier.colour=NA) + 
  labs(x="Haplotype (rs2408955/rs112233623)",
       y=paste0("RBC Count (",unit,")")) +
  geom_hline(yintercept = baseline,linetype="dashed") +
  geom_text(data = meds, aes(x = con_rs1, y = baseline, label = paste(round(median,3))), 
            size = 3, vjust = -1.8) +
  coord_cartesian(ylim = limits)

# Pie chart of allele frequencies
proportions <- data.frame(table(rs_alleles_onlyabsolutes$con_rs1))
colnames(proportions) <- c("Genotype","Freq")
proportions$Percent <- proportions$Freq / sum(proportions$Freq)

p <- plot_ly(proportions, labels = ~Genotype, values = ~Freq, type = 'pie',textposition = 'outside',textinfo = 'label+percent') %>%
  layout(title = '',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

tmpFile <- tempfile(fileext = ".pdf")
export(p, file = "/Users/erikbao/Dropbox (MIT)/HMS/Sankaran Lab/ATACSeq_GWAS/Examples/CCND3/rs9349205_rs112233623.conditional_piechart.pdf")