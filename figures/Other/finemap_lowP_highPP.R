# Read in UKBB CS
CS.df <- read.table("../../data/Finemap/UKBB_BC_v3.bed")
names(CS.df) <- c("seqnames","end","start","annotation","PP")
CS.df[,c("trait","var","region")] <- str_split_fixed(CS.df$annotation, "-", 3)
CS.df <- CS.df %>% dplyr::select(-annotation)

# Read in FM PP01 SNPs with GWAS p-value information
fm_snps <- readRDS("../../data/Finemap/all.finemap.snps.0.01.rds")

# Merge PPs with p-values
merged <- merge(CS.df, fm_snps,by.x=c("var","trait"),by.y=c("SNP","trait"))
# Extract variants with high PP, and low p-values
lowP_highPP <- merged[merged$PP.x > 0.75 & merged$p.value > 5e-8,c("var","PP.x","trait","p.value")]
lowP_highPP$var %>% unique %>% length

# See how many total FM variants have high PP (with no restriction on p-value)
merged[merged$PP.x > 0.75,"var"] %>% unique %>% length
CS.df[CS.df$PP>0.75,"var"] %>% unique %>% length
