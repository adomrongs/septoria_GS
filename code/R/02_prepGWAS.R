library(tidyverse)
library(lme4)
library(rrBLUP)

load("data/septoria_genotype.Rdata")

# prepare phenotypic data
raw_pheno <- read_csv("data/raw_phenotypes.csv")
raw_pheno <- raw_pheno %>% 
  mutate(Year = unname(sapply(strsplit(Picture, "_"), function(x) x[[2]])),
         Isolate = unname(sapply(strsplit(Picture, "_"), function(x) paste0(x[[2]], "_", x[[3]], "_", x[[4]]))))

phenotype <- raw_pheno %>% 
  dplyr::select(Isolate, Year, Leaf = leave_id, Rep = REP, Wheat_line = varieties, 
                PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) %>% 
  mutate(across(c(PLACL, pycnidiaPerCm2Lesion, pycnidiaPerCm2Leaf))) %>% 
  filter(Leaf == 2)

phenotype$Isolate[phenotype$Isolate=="22_Conil_Fer"]<-"22_ConilFer_L1"
phenotype$Isolate[phenotype$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"
phenotype$Isolate[phenotype$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2"

intersect_isolates <- intersect(phenotype$Isolate, genotype[,"Isolate"])

# adjust both genotype and phenotype
genotype <- genotype[genotype[,"Isolate"] %in% intersect_isolates,]
phenotype_clean <- phenotype %>% 
  filter(Isolate %in% intersect_isolates) %>% 
  mutate(across(c(Isolate, Year, Leaf, Rep, Wheat_line), as.factor))

# extract BLUPs
extract_BLUPs <- function(phenotype, trait){
  model <- lmer(formula(paste(trait, "~ Rep + Wheat_line + Year + (1|Isolate)")), data = phenotype)
  BLUPs <- ranef(model)$Isolate[, 1] 
  names <- rownames(ranef(model)$Isolate)
  BLUPs <- cbind(names, BLUPs)
  return(BLUPs)
}

traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
BLUPs_list <- list()
for(i in seq_along(traits)){
  trait <- traits[[i]]
  blup <- extract_BLUPs(phenotype_clean, trait)
  BLUPs_list[[i]] <- blup
}

BLUPs_df <- data.frame(do.call(cbind, BLUPs_list))
BLUPs_df <- BLUPs_df %>% 
  mutate(across(c(2,4,6), as.numeric)) %>% 
  dplyr::select(c(1,2,4,6)) 
colnames(BLUPs_df) <- c("Isolate", traits)

# impute data and create kinship
KI <- A.mat(genotype[,-1], return.imputed = T, impute.method = "EM")

geno_imputed <- data.frame(KI$imputed)
geno_imputed <- cbind(genotype[,"Isolate"], geno_imputed)
colnames(geno_imputed)[1] <- "Isolate"

K <- data.frame(KI$A)
K <- K %>% 
  mutate(Isolate = genotype[,"Isolate"]) %>% 
  dplyr::select(Isolate, everything())
K_sep <- K
colnames(K) <- rownames(K) <- NULL

# create map file
map <- data.frame(
  ID = colnames(genotype)[-1],
  CHROM = gsub("X", "", sapply(strsplit(colnames(genotype)[-1], "_"), function(x) x[1])),
  POS = sapply(strsplit(colnames(genotype)[-1], "_"), function(x) x[2])
)

dim(K); dim(geno_imputed); dim(BLUPs_df); dim(map)

save(geno_imputed, BLUPs_df, map, K_sep, file = "data/septoria_data.Rdata")

