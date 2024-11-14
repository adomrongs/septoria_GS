library(tidyverse)
library(vcfR)
library(snpReady)
library(lme4)
library(rrBLUP)

vcf <- read.vcfR("data/septoria_ale_clean.vcf.gz")
numeric_gt <- t(extract.gt(vcf, as.numeric = T)) # septoria is an haploid, thus we can just convert the genotypes to numeric
maf <- colMeans(numeric_gt, na.rm = T) # calculate minor allele frequency
range(maf) # it is correctly filtered because it does not have any vakllues below 0.05 (or above 0.95). But we have top turn around the values with maf > 0.5 (thats not minor, its major alelle!!)
numeric_gt[,which(maf > 0.5)] <- (numeric_gt[,which(maf > 0.5)] -1) * -1 # turn 0 into 1 and 1 into 0
maf_clean <- colMeans(numeric_gt, na.rm = T) 
range(maf_clean) # now maf ranges from 0.05 to 0.5

# Now I would like to change the chr names from enssembl to ncbi notation
chr_names <- read_tsv("data/chr_names.tsv")
chr_names <- chr_names %>% 
  dplyr::select(`RefSeq seq accession`, `Chromosome name`)

marker_names <- colnames(numeric_gt) 
chr <- sapply(strsplit(marker_names, "_"), function(x) paste0(x[[1]], "_", x[[2]])) #split marker names in chr and position
pos <- sapply(strsplit(marker_names, "_"), function(x) x[[3]])

new_markers <- data.frame(Old = chr) %>% 
  left_join(chr_names, by = c("Old" = "RefSeq seq accession"))
new_markers_names <- paste0(new_markers$`Chromosome name`, "_", pos) # make the change and assign the new names
colnames(numeric_gt) <- new_markers_names

# Create the map file ans dubset the markers in chr 1-13 and then readjust geno
map <- data.frame(CHROM = new_markers$`Chromosome name`,
                  POS = pos,
                  ID = colnames(numeric_gt))
map <- map %>% 
  filter(CHROM <= 13)

genotype_septoria <- numeric_gt[, colnames(numeric_gt) %in% map$ID]

# adjust rownames to convert from fasta to isolates
samples <- substr(rownames(genotype_septoria), 1, nchar(rownames(genotype_septoria)) - 
                    ifelse(grepl("merge", rownames(genotype_septoria)), 27, 13))
samples[samples == "S25_combined_EKDN23H000001-1A"] <- "S25_combined"
info_strains <- read_csv("data/INFO_STRAINS.csv")
info_strains <- info_strains %>% 
  mutate(Fasta = sapply(strsplit(Fasta, "_"), function(x) paste(x[1:2], collapse = "_")))

length(intersect(samples, info_strains$Fasta))
rownames_df <- data.frame(Fasta = samples) %>% 
  left_join(info_strains, by = "Fasta")
rownames(genotype_septoria) <- rownames_df$Isolate
  
# Create the Isolate column and create map file
genotype_septoria_clean <- data.frame(genotype_septoria) %>% 
  mutate(Isolate = rownames(genotype_septoria)) %>% 
  dplyr::select(Isolate, everything())
rownames(genotype_septoria_clean) <- NULL

# Check map and geno dimensions
dim(genotype_septoria_clean); dim(map)

# Now lets start with the phenotype
raw_pheno <- read_csv("data/raw_phenotypes.csv")
# extract year, and isolate information from picture column
raw_pheno <- raw_pheno %>% 
  mutate(Year = unname(sapply(strsplit(Picture, "_"), function(x) x[[2]])),
         Isolate = unname(sapply(strsplit(Picture, "_"), function(x) paste0(x[[2]], "_", x[[3]], "_", x[[4]]))))
# select the relevant columns and transform them to numeric and factors
phenotype <- raw_pheno %>% 
  dplyr::select(Isolate, Year, Leaf = leave_id, Rep = REP, Wheat_line = varieties, 
                PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) %>% 
  mutate(across(c(PLACL, pycnidiaPerCm2Lesion, pycnidiaPerCm2Leaf), as.numeric)) %>% 
  filter(Leaf ==2)
# adjust a few lines needed to martch with the genotype
phenotype$Isolate[phenotype$Isolate=="22_Conil_Fer"]<-"22_ConilFer_L1"
phenotype$Isolate[phenotype$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"
phenotype$Isolate[phenotype$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2"
# final adjustment
phenotype <- phenotype %>% 
  mutate(across(c(Isolate, Year, Leaf, Rep, Wheat_line), as.factor)) %>% 
  filter(Isolate %in% genotype_septoria_clean[,"Isolate"])
str(phenotype)
# check the intersect
intersect_isolates <- intersect(phenotype$Isolate, genotype_septoria_clean[,"Isolate"])
# now calculate blups
extract_BLUPs <- function(phenotype, trait) {
  # Ajustar el modelo lineal mixto
  model <- lmer(formula(paste(trait, "~ Year + Rep + Wheat_line + (1|Isolate)")), data = phenotype)
  blups <- data.frame(ranef(model)$Isolate) 
  return(blups)
}

traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
blups_list <- list()
for(i in seq_along(traits)){
  trait <- traits[[i]]
  blups_list[[i]] <- extract_BLUPs(phenotype, trait)
}
names(blups_list) <- traits
blups_df <- do.call(cbind, blups_list)
blups_df <- data.frame(blups_df) %>% 
  rownames_to_column("Isolate") %>% 
  dplyr::select(Isolate, everything())
colnames(blups_df) <- c("Isolate", traits)

# Impute genotype and calculate kinship
KI <- A.mat(genotype_septoria_clean[,-1], return.imputed = T, impute.method = "EM")
imputed_genotype_septoria <- data.frame(KI$imputed)
imputed_genotype_septoria <- imputed_genotype_septoria %>% 
  mutate(Isolate = genotype_septoria_clean[,"Isolate"]) %>% 
  dplyr::select(Isolate, everything())

K <- data.frame(KI$A)
K <- K %>% 
  mutate(Isolate = genotype_septoria_clean[,"Isolate"]) %>% 
  dplyr::select(Isolate, everything())

# Renaming and checking dimensions
septoria_geno <- imputed_genotype_septoria
septoria_pheno <- blups_df
septoria_kinship <- K
septoria_map <- map

dim(septoria_geno); dim(septoria_pheno); dim(septoria_kinship); dim(septoria_map)
save(septoria_geno, septoria_pheno, septoria_kinship, septoria_map,
     file = "data/1_septoria_data.Rdata")







