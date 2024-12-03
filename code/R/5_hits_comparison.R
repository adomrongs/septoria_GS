library(tidyverse)
library(gt)

#lets group the results from the GWAS
directories <- c("outputs/GWAS_wheat/All_leaves", "outputs/GWAS_wheat/Leaf2/")
results_list <- list()
for(i in seq_along(directories)){
  dir <- directories[[i]]
  results_list[[i]] <- grepGAPITres(dir)
}

all_leaves <- do.call(rbind, results_list[[1]])[,-1]
leaf2 <- do.call(rbind, results_list[[2]])[,-1]
dim(all_leaves); dim(leaf2)

write_csv(all_leaves, file ="outputs/GWAS_wheat/All_leaves/results_all_leaves.csv")
write_csv(leaf2, file ="outputs/GWAS_wheat/Leaf2/results_leaf2.csv")

# lets compare how many hits are shared bwetween all leaves and just leaf2
columns_to_compare <- c("SNP", "Chr", "Pos", "traits")
all_leaves_subset <- all_leaves[, columns_to_compare]
leaf2_subset <- leaf2[, columns_to_compare]

common_rows <- dplyr::intersect(all_leaves_subset, leaf2_subset)
nrow(common_rows) # there is a total of 4 hits shared by the 2 

# which model do the hits belong to? In both scenarios. Also, how many times does the snp appear per scenario?
shared_snps <- map(common_rows$SNP, function(snp) {
  tmp_df <- all_leaves %>% 
    filter(SNP == snp)
  return(tmp_df)
})
shared_snps_all_df <- do.call(rbind, shared_snps)

shared_snps_leaf2 <- map(common_rows$SNP, function(snp) {
  tmp_df <- leaf2 %>% 
    filter(SNP == snp)
  return(tmp_df)
})
shared_snps_leaf2_df <- do.call(rbind, shared_snps_leaf2)