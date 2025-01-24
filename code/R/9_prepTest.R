library(tidyverse)
library(vcfR)
library(rrBLUP)
library(factoextra)
library(ggiraph)
library(ggrepel)
source("code/R/function_septoria_GS.R")

# we have to subset the snps selected for the train in the test so that we can intergate fairly 
# the snps form the test set and also plot the PCA joined
load("data/modified_data/1_septoria.Rdata")
snps <- colnames(genotype_septoria)[-1]

chr <- map_chr(snps, ~strsplit(.x, "_")[[1]][1])
chr_df <- data.frame(v1 = chr, 
                     v2 = as.numeric(gsub("X", "", chr)))
chr_names <- read_tsv("data/raw_data/chr_names.tsv")
chr_names <- chr_names[, c("Chromosome name", "RefSeq seq accession")]
chr_df <- chr_df %>% 
  left_join(chr_names, by = c("v2" = "Chromosome name"))

pos <- map_chr(snps, ~strsplit(.x, "_")[[1]][2])
snp_list <- cbind(chr_df, pos) %>% 
  dplyr::select(3,4)

write_tsv(snp_list, file = "data/modified_data/snp_subset.txt", col_names = F)

# subset both samples and snps from the raw_snps.vcf.gz containing 119 samples x 3 M snps

#==============================================================================
# Prepare genotype
#==============================================================================

# load the vcf file correspondent to the test samples
raw_test_vcf <- read.vcfR("data/raw_data/clean_test.vcf.gz")
raw_allele_mat <- t(extract.gt(raw_test_vcf, return.alleles = T)) # extract allele matrix
# asjust marker names
chr_names <- read_tsv("data/raw_data/chr_names.tsv")
chr_names <- chr_names %>% 
  dplyr::select(`Chromosome name`, `RefSeq seq accession`)

marker_names <- colnames(raw_allele_mat) 
chr <- sapply(strsplit(marker_names, "_"), function(x) paste0(x[[1]], "_", x[[2]])) #split marker names in chr and position
pos <- sapply(strsplit(marker_names, "_"), function(x) x[[3]])

marker_names <- data.frame(old_chr = chr, pos = pos)
marker_names <- marker_names %>% 
  left_join(chr_names, by = c("old_chr" = "RefSeq seq accession"))
marker_names <- marker_names %>% 
  mutate(new_marker_names = paste0("X", `Chromosome name`, "_", pos))

colnames(raw_allele_mat) <- marker_names$new_marker_names

# subset the markers included in the train matrix
clean_allele_mat <- raw_allele_mat[, colnames(raw_allele_mat) %in% colnames(genotype_septoria)]

# change the rownames to isolate format
samples <- substr(rownames(clean_allele_mat), 1, nchar(rownames(clean_allele_mat)) - 
                    ifelse(grepl("merge", rownames(clean_allele_mat)), 27, 13))
info_strains <- read_csv("data/raw_data/INFO_STRAINS.csv")
info_strains <- info_strains %>% 
  mutate(Fasta = sapply(strsplit(Fasta, "_"), function(x) paste(x[1:2], collapse = "_")))
length(intersect(samples, info_strains$Fasta))

rownames_df <- data.frame(Fasta = samples) %>% 
  left_join(info_strains, by = "Fasta")
fasta_to_isolate <- setNames(rownames_df$Isolate, rownames_df$Fasta)
rownames(clean_allele_mat) <- unname(fasta_to_isolate[samples])
rownames(clean_allele_mat) <- gsub("\\.", "_", rownames(clean_allele_mat))

clean_allele_mat[clean_allele_mat == "."] <- NA
M <- allele2numeric(clean_allele_mat)
# add Issolate column

genotype_test <- data.frame(M) %>% 
  mutate(Isolate = rownames(clean_allele_mat)) %>% 
  dplyr::select(Isolate, everything())
ncol(genotype_test) == ncol(genotype_septoria)

genotype_combined <- rbind(genotype_septoria, genotype_test)
M_imputed <- A.mat(genotype_combined[,-1], return.imputed = T, impute.method = "EM")
M_impute <- M_imputed$imputed
M_impute <- data.frame(M_impute) %>% 
  mutate(Isolate = genotype_combined[,1]) %>% 
  dplyr::select(Isolate, everything())

genotype_all <- M_impute
k_all <- A.mat(as.matrix(genotype_all[,-1]))
rownames(k_all) <- colnames(k_all) <- genotype_all[,1]

#==============================================================================
# Prepare phenotype
#==============================================================================

raw_septoria_phenotype <- read_csv("data/raw_data/raw_phenotypes.csv")

factor_cols <- c("REP", "varieties")
septoria_phenotype <- raw_septoria_phenotype %>%
  mutate(
    Isolate = if_else(
      varieties == "Don Ricardo",
      gsub("\\.", "_", sapply(str_split(Picture, "_"), 
                              function(x) paste0(x[3], "_", x[4], "_", x[5]))),
      gsub("\\.", "_", sapply(str_split(Picture, "_"), 
                              function(x) paste0(x[2], "_", x[3], "_", x[4])))
    ),
    Year = if_else(
      varieties == "Don Ricardo",
      sapply(str_split(Picture, "_"), function(x) x[3]),
      sapply(str_split(Picture, "_"), function(x) x[2])
    ),
    BRep = if_else(
      varieties == "Don Ricardo",
      sapply(str_split(Picture, "_"), function(x) x[6]),
      sapply(str_split(Picture, "_"), function(x) x[5])
    ), 
    across(all_of(factor_cols), as.factor)
  ) %>% 
  dplyr::select(Isolate, Line = varieties, Trial = REP, Leaf = leave_id, Year, BRep,
                PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) |> 
  mutate(across(Isolate:BRep, as.factor))

septoria_phenotype$Isolate[septoria_phenotype$Isolate=="22_Conil_Fer"]<-"22_ConilFer_L1"
septoria_phenotype$Isolate[septoria_phenotype$Isolate=="22_EcijaSecSha_L1"]<-"22_EcijaSecSah_L1"
septoria_phenotype$Isolate[septoria_phenotype$Isolate=="22_EcijaSecSim_L1"]<-"22_EcijaSecSim_L2"

cleaned_septoria_phenotype_1 <- septoria_phenotype %>% 
  filter(Isolate %in% genotype_all[,1]) |> 
  filter(PLACL < 101) |> 
  droplevels()

cleaned_septoria_phenotype_2 <- cleaned_septoria_phenotype_1 %>% 
  filter(Leaf == 2) %>% 
  dplyr::select(-Leaf) |> 
  droplevels()

cleaned_septoria_phenotype_2 <- remove_outliers(cleaned_septoria_phenotype_2,
                                               cols = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"))

write_csv(cleaned_septoria_phenotype_1, file = 'data/modified_data/gh1_clean_pheno.csv')


#==============================================================================
# Prepare test phenotype 
#==============================================================================

raw_test_pheno <- readxl::excel_sheets("data/raw_data/Output_Exp5greenhouse2023_19strains_Rep1and2.xlsx") %>%
  map_df(~ readxl::read_xlsx("data/raw_data/Output_Exp5greenhouse2023_19strains_Rep1and2.xlsx", 
                     col_names = TRUE, sheet = .x) %>% dplyr::select(-1))

test_pheno <- raw_test_pheno %>% 
  dplyr::select(strain, dpi, cultivar, replicate, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) %>% 
  mutate(
    across(c(strain, dpi, cultivar, replicate), as.factor),
    across(c(PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion), as.numeric)
  ) %>% 
  dplyr::select(Isolate = strain, Line = cultivar, TRep = dpi, BRep = replicate,
                PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion)

clean_test_pheno <- remove_outliers(test_pheno, c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"))

save(genotype_all, k_all, cleaned_septoria_phenotype_1, cleaned_septoria_phenotype_2,
     clean_test_pheno, file = "data/modified_data/5_predictions.Rdata")

#==============================================================================
# Plot combined PCA
#==============================================================================

info_strains$Isolate <- gsub("\\.", "_", info_strains$Isolate )
pca_df <- data.frame(Isolate = genotype_all[,1]) %>% 
  left_join(info_strains, by = "Isolate") %>% 
  dplyr::select(Isolate, Region = Partition)

colors <- c("Train" = "#DD5129FF", 
            "Test" = "#0F7BA2FF")

shapes <- c("Train" = 16,
            "Test" = 16)

# calculate PCA plot for the combined dataset
PCA_all <- plotPCA(genotype = genotype_all,
        regions = pca_df$Region,
        colors = colors,
        names = pca_df$Isolate, 
        interactive = T,
        shapes = shapes)

png(paste0("outputs/plots/PCA_all_sep.png"), width = 3000, height = 3000, res = 400)
PCA_all$pca_plot
dev.off()






























