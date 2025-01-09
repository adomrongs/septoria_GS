library(tidyverse)
library(biomaRt)
library(seqinr)
source('code/R/function_septoria_GS.R')

hits <- list.files("outputs/GWAS_sep/", full.names = T)
cultivar <-  c(rep("Athoris", 2), rep("Don", 2), rep("Sculptur", 1), rep("Svevo", 2), rep("All", 2))
pcs <- c(1, 2, 2, 3, 3, 2, 1, 2, 2)
leaf <- c(2, 3, 2, 3, 2, 2, 3, 2, 3)

ref_df <- data.frame(Cultivar = cultivar, PCs = pcs, Leaf = leaf)

df_list <- list()
for(i in seq_along(hits)){
  file <- hits[[i]]
  tmp_file <- read_csv(file)
  df_list[[i]] <- tmp_file |> 
    mutate(PCs = pcs[[i]], Leaf = leaf[[i]], Cultivar = cultivar[[i]])
}


df_hits <- do.call(rbind, df_list)
snps <- df_hits |> 
  mutate(Trait = gsub("^BLINK.", "", traits),
         Chrom = as.numeric(str_extract(SNP, "(?<=X)(\\d+)(?=_)"))) |> 
  arrange(Trait, Chrom) |> 
  dplyr::select(Trait, SNP, Chrom, PCs, Leaf, Cultivar)

sort(table(snps$Chrom), decreasing = T) # chromosomes 1,9 and 8 are the chromocmes with the largest amount of hits
final_hits <- snps |> 
  group_by(SNP) |> 
  filter(n() == 1) |> 
  ungroup() 

dim(final_hits) # we will be looking for candidate genes affected by these 28 snps
# it would be correct to select the hits that are being detected twice beacuse this would mean that 
# we would be paying attention just to those ones that affect multiple traits

mart <- useEnsemblGenomes(biomart = "fungi_mart", dataset = "ztritici_eg_gene")
attributes <- c("ensembl_gene_id", "start_position", "end_position", "strand", "description", "peptide")
filters <- c("chromosome_name", "start", "end")
distances <- list(500, 1000, 2000) #base pairs (bp)
chr <- final_hits$Chrom
traits <- final_hits$Trait
markers <- as.numeric(map(strsplit(final_hits$SNP, "_"), \(x) x[[2]]))
out_dir <- 'outputs/FASTA_sep'

candidate_genes <- find_genes(mart, attributes, filters, distances, chr, traits, out_dir)

# summary abolut genes and hits
# Number of hits and genes per trait
table1 <- candidate_genes |> 
  group_by(trait) |> 
  summarise(
    Hits = n_distinct(marker), 
    Genes = n_distinct(gene)
  ) |> 
  arrange(desc(Hits), desc(Genes))

# Number of Genes and hit per hit 
table2 <- candidate_genes |> 
  group_by(marker) |>
  summarise(Genes = n_distinct(gene), Trait = n_distinct(trait)) 
mean(table2$Genes); mean(table2$Trait)

# Number of genes and hit per chromosome 
table3 <- candidate_genes |> 
  mutate(chr = as.numeric(str_extract(marker, "(?<=X)(\\d+)(?=_)"))) |> 
  group_by(chr) |>
  summarise(Hits = n_distinct(marker), Genes = n_distinct(gene)) |> 
  arrange(desc(Hits), desc(Genes))



