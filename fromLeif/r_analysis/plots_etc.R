library(BuenColors)
library(dplyr)

df <- read.table("reporter_experiments_both.txt", header = TRUE)

## A/A is reference

AK3 <- df %>% filter(type %in% c("pGL","AK3_AA", "AK3_AG", "AK3_CA", "AK3_CG"))
ggplot(AK3, aes(x = type, y = value)) + 
  geom_bar(aes(fill = Experiment), position = "dodge", stat="identity") + pretty_plot()

AK3r <- df %>% filter(type %in% c("AK3_AA", "AK3_AG", "AK3_CA", "AK3_CG"))
AK3r$XG <- ifelse(AK3r$type %in% c("AK3_AG", "AK3_CG"),1,0)
AK3r$CX <- ifelse(AK3r$type %in% c("AK3_CA", "AK3_CG"),1,0)
summary(lm(value ~ Experiment + XG + CX + XG*CX, data = AK3r))

# C/G is reference

CCND3 <- df %>% filter(type %in% c("pGL","CCND3_CG", "CCND3_TG", "CCND3_CA", "CCND3_TA"))
ggplot(CCND3, aes(x = type, y = value)) + 
  geom_bar(aes(fill = Experiment), position = "dodge", stat="identity") + pretty_plot()

CCND3r <- df %>% filter(type %in% c("CCND3_CG", "CCND3_TG", "CCND3_CA", "CCND3_TA"))
CCND3r$XA <- ifelse(CCND3r$type %in% c("CCND3_CA", "CCND3_TA"),1,0)
CCND3r$TX <- ifelse(CCND3r$type %in% c("CCND3_TG", "CCND3_TA"),1,0)
summary(lm(value ~ Experiment + XA + TX + XA*TX, data = CCND3r))
