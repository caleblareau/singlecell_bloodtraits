library(data.table)
library(ggplot2)
library(BuenColors)
library(plotly)
library(dplyr)
library(scales)

########################################################################################################################
# CCND3 Bar plots
trait <- "RBC_COUNT"
rs_alleles <- fread("../../../data/examples/rawphenos_mixed.rs9349205_rs112233623.txt")
colnames(rs_alleles)[2] <- "alleleA.rs2408955"

# Round imputed genotypes to nearest integer
rs_alleles_onlyabsolutes <- as.data.frame(rs_alleles)
rs_alleles_onlyabsolutes$alleleA.rs112233623 <- round(rs_alleles_onlyabsolutes$alleleA.rs112233623,0)
rs_alleles_onlyabsolutes$alleleB.rs112233623 <- 2-rs_alleles_onlyabsolutes$alleleA.rs112233623
rs_alleles_onlyabsolutes$alleleB.rs2408955 <- 2-rs_alleles_onlyabsolutes$alleleA.rs2408955

rs_alleles_onlyabsolutes <- rs_alleles_onlyabsolutes[,c("alleleB.rs2408955","alleleB.rs112233623","RBC_COUNT")]

# Make bar plot with mean and se of RBC_COUNT
eval_CCND3 <- function(rs112233623, rs112233623_l, rs9349205, rs9349205_l){
  v <- subset(rs_alleles_onlyabsolutes,alleleB.rs112233623==rs112233623 & alleleB.rs2408955==rs9349205)
  if (nrow(v)<10){
    v <- NULL
  }
  data.frame(name = paste0("rs112233623: ", rs112233623_l, "\nrs9349205: ", rs9349205_l), estimate = mean(v$RBC_COUNT), se = sd(v$RBC_COUNT)/sqrt(length(v$RBC_COUNT)))
}

CCND3plotdf <- rbind(
  eval_CCND3(rs112233623 = 0, rs112233623_l = "CC", rs9349205 = 0, rs9349205_l = "GG"),
  eval_CCND3(rs112233623 = 0, rs112233623_l = "CC", rs9349205 = 1, rs9349205_l = "AG"),
  eval_CCND3(rs112233623 = 0, rs112233623_l = "CC", rs9349205 = 2, rs9349205_l = "AA"),
  eval_CCND3(rs112233623 = 1, rs112233623_l = "TC", rs9349205 = 0, rs9349205_l = "GG"),
  eval_CCND3(rs112233623 = 1, rs112233623_l = "TC", rs9349205 = 1, rs9349205_l = "AG"),
  eval_CCND3(rs112233623 = 1, rs112233623_l = "TC", rs9349205 = 2, rs9349205_l = "AA"),
  eval_CCND3(rs112233623 = 2, rs112233623_l = "TT", rs9349205 = 0, rs9349205_l = "GG"),
  eval_CCND3(rs112233623 = 2, rs112233623_l = "TT", rs9349205 = 1, rs9349205_l = "AG"),
  eval_CCND3(rs112233623 = 2, rs112233623_l = "TT", rs9349205 = 2, rs9349205_l = "AA")
)

CCND3plotdf <- CCND3plotdf[complete.cases(CCND3plotdf),]

limits <- c(as.numeric(quantile(rs_alleles_onlyabsolutes[,trait],0.3)),
            as.numeric(quantile(rs_alleles_onlyabsolutes[,trait],0.7)))

ccnd3_phenos<- ggplot(CCND3plotdf, aes(x = name, y = estimate)) + 
  geom_bar(stat = "identity", color = "black", fill = "dodgerblue4") + pretty_plot() +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.1) +
  labs(x = "", y = "RBC_COUNT (10^6 cells/uL)") + coord_cartesian(ylim = limits)

cowplot::ggsave(ccnd3_phenos, file = "CCND3_gwas_phenotypes.pdf", height = 5, width = 8)

########################################################################################################################
# AK3
# Bar plots
rs_alleles <- fread("../../../data/examples/rawphenos_mixed.rs12005199_rs409950.txt")
trait <- "PLT_COUNT"
# Round imputed genotypes to nearest integer
rs_alleles_onlyabsolutes <- as.data.frame(rs_alleles)
rs_alleles_onlyabsolutes$alleleA.rs12005199 <- round(rs_alleles_onlyabsolutes$alleleA.rs12005199,0)
rs_alleles_onlyabsolutes$alleleB.rs12005199 <- 2-rs_alleles_onlyabsolutes$alleleA.rs12005199
rs_alleles_onlyabsolutes$alleleA.rs409950 <- round(rs_alleles_onlyabsolutes$alleleA.rs409950,0)
rs_alleles_onlyabsolutes$alleleB.rs409950 <- 2-rs_alleles_onlyabsolutes$alleleA.rs409950

rs_alleles_onlyabsolutes <- rs_alleles_onlyabsolutes[,c("alleleB.rs12005199","alleleB.rs409950",trait)]

# Make bar plot with mean and se of RBC_COUNT
eval_AK3 <- function(rs409950, rs409950_l, rs12005199, rs12005199_l,trait="PLT_COUNT"){
  v <- subset(rs_alleles_onlyabsolutes,alleleB.rs409950==rs409950 & alleleB.rs12005199==rs12005199)
  if (nrow(v)<10){
    v <- NULL
  }
  data.frame(name = paste0("rs409950: ", rs409950_l, "\nrs12005199: ", rs12005199_l), estimate = mean(v[,trait]), se = sd(v[,trait])/sqrt(length(v[,trait])))
}

AK3plotdf <- rbind(
  eval_AK3(rs409950 = 2, rs409950_l = "AA", rs12005199 = 2, rs12005199_l = "AA"),
  eval_AK3(rs409950 = 1, rs409950_l = "AC", rs12005199 = 2, rs12005199_l = "AA"),
  eval_AK3(rs409950 = 0, rs409950_l = "CC", rs12005199 = 2, rs12005199_l = "AA"),
  eval_AK3(rs409950 = 2, rs409950_l = "AA", rs12005199 = 1, rs12005199_l = "AG"),
  eval_AK3(rs409950 = 1, rs409950_l = "AC", rs12005199 = 1, rs12005199_l = "AG"),
  eval_AK3(rs409950 = 0, rs409950_l = "CC", rs12005199 = 1, rs12005199_l = "AG"),
  eval_AK3(rs409950 = 2, rs409950_l = "AA", rs12005199 = 0, rs12005199_l = "GG"),
  eval_AK3(rs409950 = 1, rs409950_l = "AC", rs12005199 = 0, rs12005199_l = "GG"),
  eval_AK3(rs409950 = 0, rs409950_l = "CC", rs12005199 = 0, rs12005199_l = "GG")
)

AK3plotdf <- AK3plotdf[complete.cases(AK3plotdf),]

limits <- c(as.numeric(quantile(rs_alleles_onlyabsolutes[,trait],0.3)),
            as.numeric(quantile(rs_alleles_onlyabsolutes[,trait],0.7)))

ak3_phenos <- ggplot(AK3plotdf, aes(x = name, y = estimate)) + 
  geom_bar(stat = "identity", color = "black", fill = "firebrick") + pretty_plot() +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.1) +
  labs(x = "", y = "Platelet Count (10^6 cells/uL)") + coord_cartesian(ylim = limits)
cowplot::ggsave(ak3_phenos, file = "AK3_gwas_phenotypes.pdf", height = 5, width = 10)

########################################################################################################################
# CDK6 - EO_COUNT
# Bar plots
rs_alleles <- fread("../../../data/examples/rawphenos_mixed.rs8_rs445.txt")
trait <- "EO_COUNT"

colnames(rs_alleles)[1:2] <- c("alleleA.rs1","alleleA.rs2")
# rs1 = rs445, rs2 = rs8

# Round imputed genotypes to nearest integer
rs_alleles_onlyabsolutes <- as.data.frame(rs_alleles)
rs_alleles_onlyabsolutes$alleleA.rs1 <- round(rs_alleles_onlyabsolutes$alleleA.rs1,0)
rs_alleles_onlyabsolutes$alleleB.rs1 <- 2-rs_alleles_onlyabsolutes$alleleA.rs1
rs_alleles_onlyabsolutes$alleleA.rs2 <- round(rs_alleles_onlyabsolutes$alleleA.rs2,0)
rs_alleles_onlyabsolutes$alleleB.rs2 <- 2-rs_alleles_onlyabsolutes$alleleA.rs2

rs_alleles_onlyabsolutes <- rs_alleles_onlyabsolutes[,c("alleleB.rs1","alleleB.rs2",trait)]

# Make bar plot with mean and se of trait
eval_gene <- function(rs2, rs2_l, rs1, rs1_l,trait="EO_COUNT"){
  v <- subset(rs_alleles_onlyabsolutes,alleleB.rs2==rs2 & alleleB.rs1==rs1)
  if (nrow(v)<10){
    v <- NULL
  }
  data.frame(name = paste0("rs8: ", rs2_l, "\nrs445: ", rs1_l), estimate = mean(v[,trait]), se = sd(v[,trait])/sqrt(length(v[,trait])))
}
CDK6plotdf <- rbind(
  eval_gene(rs2 = 2, rs2_l = "TT", rs1 = 2, rs1_l = "TT"),
  eval_gene(rs2 = 1, rs2_l = "CT", rs1 = 2, rs1_l = "TT"),
  eval_gene(rs2 = 0, rs2_l = "CC", rs1 = 2, rs1_l = "TT"),
  eval_gene(rs2 = 2, rs2_l = "TT", rs1 = 1, rs1_l = "CT"),
  eval_gene(rs2 = 1, rs2_l = "CT", rs1 = 1, rs1_l = "CT"),
  eval_gene(rs2 = 0, rs2_l = "CC", rs1 = 1, rs1_l = "CT"),
  eval_gene(rs2 = 2, rs2_l = "TT", rs1 = 0, rs1_l = "CC"),
  eval_gene(rs2 = 1, rs2_l = "CT", rs1 = 0, rs1_l = "CC"),
  eval_gene(rs2 = 0, rs2_l = "CC", rs1 = 0, rs1_l = "CC")
)

CDK6plotdf <- CDK6plotdf[complete.cases(CDK6plotdf),]

limits <- c(as.numeric(quantile(rs_alleles_onlyabsolutes[,trait],0.5)),
            as.numeric(quantile(rs_alleles_onlyabsolutes[,trait],0.65)))

cdk6_phenos <- ggplot(CDK6plotdf, aes(x = name, y = estimate)) + 
  geom_bar(stat = "identity", color = "black", fill = "firebrick") + pretty_plot() +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.1) +
  labs(x = "", y = "Eosinophil Count (10^6 cells/uL)") + coord_cartesian(ylim = limits)
cowplot::ggsave(cdk6_phenos, file = "CDK6_eo_count_phenotypes.pdf", height = 5, width = 10)


#########
# Piechart of haplotype frequencies for the two AK3 SNPs
rs_alleles_onlyabsolutes$con_rs1 <- paste(rs_alleles_onlyabsolutes$alleleB.rs409950,rs_alleles_onlyabsolutes$alleleB.rs12005199,sep="/")
map <- setNames(c("AA/AA","AC/AA","CC/AA",
                  "AA/AG","AC/AG","CC/AG",
                  "AA/GG", "AC/GG", "CC/GG"),
                c("2/2", "1/2", "0/2",
                  "2/1","1/1","0/1",
                  "2/0","1/0","0/0"))
rs_alleles_onlyabsolutes$haps <- map[rs_alleles_onlyabsolutes$con_rs1]
proportions <- data.frame(table(rs_alleles_onlyabsolutes$haps))
colnames(proportions) <- c("Genotype","Freq")

proportions$Percent <- proportions$Freq / sum(proportions$Freq)

p <- plot_ly(proportions, labels = ~Genotype, values = ~Freq, type = 'pie',textposition = 'outside',textinfo = 'label+percent') %>%
  layout(title = '',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

tmpFile <- tempfile(fileext = ".pdf")
export(p, file = "AK3.haplotype_piechart.pdf")
