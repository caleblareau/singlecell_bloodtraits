library(BuenColors)
library(dplyr)
source("summarySE.R")

df <- read.table("reporter_experiments_both.txt", header = TRUE)

## A/A is reference
# rs409950 alt is C; rs12005199 alt is G

AK3 <- df %>% filter(type %in% c("pGL","AK3_AA", "AK3_AG", "AK3_CA", "AK3_CG"))

AK3r <- df %>% filter(type %in% c("AK3_AA", "AK3_AG", "AK3_CA", "AK3_CG"))
AK3r$XG <- ifelse(AK3r$type %in% c("AK3_AG", "AK3_CG"),1,0)
AK3r$CX <- ifelse(AK3r$type %in% c("AK3_CA", "AK3_CG"),1,0)
AK3mod <- lm(value ~ Experiment + XG + CX + XG*CX, data = AK3r)

# Get empty promoter
pGL <- df %>% filter(type %in% c("pGL")) %>% summarySE(measurevar="value")

# Evaluate at the different haplotyeps
eval_AK3 <- function(rs409950, rs409950_l, rs12005199, rs12005199_l){
  v <- predict(AK3mod, data.frame(Experiment = "Exp1", XG = rs12005199, CX = rs409950), se.fit = TRUE)
  data.frame(name = paste0("rs409950: ", rs409950_l, "\nrs12005199: ", rs12005199_l), estimate = unname(v$fit), se = v$se.fit)
}

AK3plotdf <- rbind(
  data.frame(name = "pGL", estimate = pGL$value, se = pGL$se),
  eval_AK3(rs409950 = 0, rs409950_l = "A", rs12005199 = 0, rs12005199_l = "A"),
  eval_AK3(rs409950 = 1, rs409950_l = "C", rs12005199 = 0, rs12005199_l = "A"),
  eval_AK3(rs409950 = 0, rs409950_l = "A", rs12005199 = 1, rs12005199_l = "G"),
  eval_AK3(rs409950 = 1, rs409950_l = "C", rs12005199 = 1, rs12005199_l = "G")
)

ak3gg <- ggplot(AK3plotdf, aes(x = name, y = estimate)) + 
  geom_bar(stat = "identity", color = "black", fill = "firebrick") + pretty_plot() +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.1) +
  labs (x = "", y = "Luciferase Activity (A.U.)")
cowplot::ggsave(ak3gg, file = "plots/AK3_luciferase_regression.pdf", height = 5, width = 8)

##################################################################################################################

# C/G is reference
# rs112233623 alt is T; rs9349205 alt is A

CCND3 <- df %>% filter(type %in% c("pGL","CCND3_CG", "CCND3_TG", "CCND3_CA", "CCND3_TA"))
ggplot(CCND3, aes(x = type, y = value)) + 
  geom_bar(aes(fill = Experiment), position = "dodge", stat="identity") + pretty_plot()

CCND3r <- df %>% filter(type %in% c("CCND3_CG", "CCND3_TG", "CCND3_CA", "CCND3_TA"))
CCND3r$XA <- ifelse(CCND3r$type %in% c("CCND3_CA", "CCND3_TA"),1,0)
CCND3r$TX <- ifelse(CCND3r$type %in% c("CCND3_TG", "CCND3_TA"),1,0)
CCND3mod  <- lm(value ~ Experiment + XA + TX + XA*TX, data = CCND3r)

# Evaluate at the different haplotyeps
eval_CCND3 <- function(rs112233623, rs112233623_l, rs9349205, rs9349205_l){
  v <- predict(CCND3mod, data.frame(Experiment = "Exp1", TX = rs112233623, XA = rs9349205), se.fit = TRUE)
  data.frame(name = paste0("rs112233623: ", rs112233623_l, "\nrs9349205: ", rs9349205_l), estimate = unname(v$fit), se = v$se.fit)
}

CCND3plotdf <- rbind(
  data.frame(name = "pGL", estimate = pGL$value, se = pGL$se),
  eval_CCND3(rs112233623 = 0, rs112233623_l = "C", rs9349205 = 0, rs9349205_l = "G"),
  eval_CCND3(rs112233623 = 0, rs112233623_l = "C", rs9349205 = 1, rs9349205_l = "A"),
  eval_CCND3(rs112233623 = 1, rs112233623_l = "T", rs9349205 = 0, rs9349205_l = "G"),
  eval_CCND3(rs112233623 = 1, rs112233623_l = "T", rs9349205 = 1, rs9349205_l = "A")
)

ccnd3gg <- ggplot(CCND3plotdf, aes(x = name, y = estimate)) + 
  geom_bar(stat = "identity", color = "black", fill = "dodgerblue4") + pretty_plot() +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.1) +
  labs (x = "", y = "Luciferase Activity (A.U.)")
cowplot::ggsave(ccnd3gg, file = "plots/CCND3_luciferase_regression.pdf", height = 5, width = 8)

