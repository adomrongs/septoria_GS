library(tidyverse)

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


