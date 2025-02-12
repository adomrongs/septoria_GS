library(tidyverse)
library(biomaRt)
library(seqinr)
library(CMplot)
library(here)
library(LDheatmap)
library(janitor)
library(GO.db)
library(paletteer)
library(hrbrthemes)
library(ggdist)
library(gridExtra)
library(scales)
source('code/R/function_septoria_GS.R')

load("data/modified_data/1_septoria.Rdata")

hits <- grep("results", list.files("outputs/GWAS_sep/", full.names = T), value = T)
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
df_hits <- do.call(rbind, df_list) |> 
  mutate(traits = gsub("BLINK.", "", traits)) |> 
  dplyr::select(-1) |> 
  relocate(Cultivar, Leaf, PCs) |> 
  arrange(Cultivar, PCs, Chr, Pos)

for (leaf in unique(df_hits$Leaf)) {
  df_hits_leaf <- df_hits |> filter(Leaf == leaf)
  output_file <- paste0("outputs/postGWAS_sep/cultivars_l", leaf, ".csv")
  fileConn <- file(output_file, "w")
  cultivars <- unique(df_hits_leaf$Cultivar)
  for (cultivar in cultivars) {
    subset_df <- df_hits_leaf[df_hits_leaf$Cultivar == cultivar, ]
    write.table(subset_df, file = fileConn, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)
    cat("\n", file = fileConn)
  }
  close(fileConn)
}

df_hits_l2 <- df_hits |> 
  filter(Leaf ==2)

snps <- df_hits_l2 |> 
  mutate(Trait = gsub("^BLINK.", "", traits),
         Chrom = as.numeric(str_extract(SNP, "(?<=X)(\\d+)(?=_)"))) |> 
  arrange(Trait, Chrom) |> 
  dplyr::select(Trait, SNP, Chrom, PCs, Leaf, Cultivar)

sort(table(snps$Chrom), decreasing = T) # chromosomes 1,9 and 8 are the chromocmes with the largest amount of hits
final_hits <- snps |> 
  group_by(SNP) |> 
  filter(n() == 1) |> 
  ungroup() 

dim(final_hits) # we will be looking for candidate genes affected by these 18 snps
# it would be correct to select the hits that are being detected twice beacuse this would mean that 
# we would be paying attention just to those ones that affect multiple traits

mart <- useEnsemblGenomes(biomart = "fungi_mart", dataset = "ztritici_eg_gene")
attributes <- c("ensembl_gene_id", "start_position", "end_position", "strand", "description", "peptide")
filters <- c("chromosome_name", "start", "end")
distances <- list(500, 1000, 2000) #base pairs (bp)
chr <- final_hits$Chrom
traits <- final_hits$Trait
markers <- as.numeric(map(strsplit(final_hits$SNP, "_"), \(x) x[[2]]))
out_dir <- 'outputs/postGWAS_sep/FASTA_sep/'

candidate_genes <- find_genes(mart, attributes, filters, distances, chr, traits, out_dir, markers)
write_csv(candidate_genes, 'outputs/postGWAS_sep/genes_l2.csv')

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

#===============================================================================
# Post GWAS for the hits identified using all cultivars and leaf 2
#===============================================================================

# generate both gwas and genes tables
all_l2 <- snps |> 
  filter(Leaf == 2, 
         Cultivar == "All") 

chr_all <- all_l2$Chrom
traits_all <- all_l2$Trait
markers_all <- as.numeric(map(strsplit(all_l2$SNP, "_"), \(x) x[[2]]))
out_dir_all <- 'outputs/postGWAS_sep/FASTA_all'
candidate_genes_all <- find_genes(mart = mart,
                                  attributes = attributes,
                                  filters = filters,
                                  distances = distances,
                                  chr = chr_all,
                                  traits = traits_all,
                                  out_dir = out_dir_all,
                                  markers = markers_all)

gwas_sep_table_all <- df_hits |> 
  filter(Cultivar == "All", 
         Leaf == 2) |> 
  dplyr::select(-c(Cultivar, Leaf, PCs)) |> 
  relocate(Trait = traits) |> 
  mutate(MAF = round(2* MAF, 2)) |> 
  left_join(candidate_genes_all, by = c('SNP' = 'marker')) |> 
  dplyr::select(Trait, SNP, Chr, Pos, P.value, MAF, Distance = distance, Gene = gene) |> 
  distinct()

genes_sep_table_all <- candidate_genes_all |> 
  dplyr::select(Trait = trait, Marker = marker, Distance = distance, Gene = gene, 
                Protein = accesion, GO_id, GO_name = Go_name) |> 
  group_by(Marker, Distance, Gene, Protein, GO_id, GO_name) |> 
  summarise(Traits = paste(unique(Trait), collapse = " and "), .groups = "drop") |> 
  relocate(Traits)

genes_sep <- genes_sep_table_all |>
  mutate(
    Traits = ifelse(duplicated(Traits), NA, Traits),
    Marker = ifelse(duplicated(Marker), NA, Marker),
    Distance = ifelse(duplicated(Distance), NA, Distance),
    Gene = ifelse(duplicated(Gene), NA, Gene),
    Protein = ifelse(duplicated(Protein), NA, Protein)
  ) |> 
  remove_empty(which = "rows")

go_ids_sep <- genes_sep_table_all |> 
  dplyr::select(c(Gene, Protein, GO_id, GO_name)) |> 
  mutate(GO_class = Ontology(GO_id),
         Gene = ifelse(duplicated(Gene), NA, Gene),
         Protein = ifelse(duplicated(Protein), NA, Protein)
  ) 

write_csv(gwas_sep_table_all, file = 'outputs/postGWAS_sep/gwas_sep_table_all.csv')
write_csv(genes_sep, file = 'outputs/postGWAS_sep/genes_sep_all.csv')
write_csv(go_ids_sep, file = 'outputs/postGWAS_sep/go_ids_sep_all.csv')


# plot manhattan and qqplot
leaf <- read_csv('outputs/GWAS_sep/GAPIT.Association.GWAS_Results.BLINK.pycnidiaPerCm2Leaf.csv') # file corresponding to all cultivars leaf2
lesion <- read_csv('outputs/GWAS_sep/GAPIT.Association.GWAS_Results.BLINK.pycnidiaPerCm2Lesion.csv')

pvalues <- leaf |> 
  dplyr::select(1:4) |> 
  bind_cols(lesion[,4]) |> 
  dplyr::rename(
    pycnidiaPerCm2Leaf = 4,
    pycnidiaPerCm2Lesion = 5
  )

setwd('outputs/postGWAS_sep/')
CMplot(pvalues,
       col = c('#0F7BA2FF', '#43B284FF'),
       plot.type="m",
       multraits=TRUE,
       threshold= 0.05/nrow(pvalues),
       threshold.lty =1,
       threshold.lwd = 1,
       threshold.col=c("black"),
       amplify=TRUE,
       bin.size=1e6,
       signal.cex=1,
       file="jpg",
       file.name="",
       dpi=300,
       file.output= TRUE,
       verbose=TRUE,
       points.alpha=250,
       legend.ncol=1,
       legend.pos="left")
setwd(here())

plotCMqq(leaf[,1:4], "Leaf", '#0F7BA2FF', 'outputs/postGWAS_sep/')
plotCMqq(lesion[,1:4], "Lesion", '#43B284FF', 'outputs/postGWAS_sep/')

# plot LDheatmap
snps <- colnames(genotype_septoria)[-1]
snps_df <- data.frame(SNPs = snps, 
                      Chromosome = as.numeric(str_extract(snps, 'X([0-9])_', group = T)),
                      Position = as.numeric(str_extract(snps, '_(.*)', group = T)))

hits <- plotLD(gwas_table = gwas_sep_table_all, genotype = genotype_septoria[, -1], df_info = snps_df)


# plot boxplot of allelic diffferences
load("data/modified_data/GWASdataBefGapit.RData")

blups_sep <- BLUPS |> dplyr::select(1:4)
names_conversion <- data.frame(geno = genotype_septoria[,1], 
                               lower_case = tolower(genotype_septoria[,1])) |> 
  left_join(blups_sep, by = c("lower_case" = "Isolate")) |> 
  na.omit() |> 
  dplyr::select(-c(PLACL:pycnidiaPerCm2Lesion))

blups_sep_ready <- blups_sep |> 
  left_join(names_conversion, by = c("Isolate" = "lower_case")) |> 
  dplyr::select(Isolate = geno, PLACL:pycnidiaPerCm2Lesion)
genotype_sep_ready <- genotype_septoria[genotype_septoria[,1] %in% blups_sep_ready$Isolate, ]


boxplot_list <- list()
for(i in 1:nrow(gwas_sep_table_all)){
  marker <- gwas_sep_table_all$SNP[i]
  trait <- gwas_sep_table_all$Trait[i]
  print(paste0("working on marker: ", marker))
  
  phenotype <- blups_sep_ready %>%  dplyr::select(Isolate, trait)
  genotype <- genotype_sep_ready[genotype_sep_ready[, "Isolate"] %in% phenotype[,1], ]
  
  boxplot_list[[i]] <- plot_allelic_diff_sep(phenotype = phenotype,
                                       genotype =  genotype,
                                       marker = marker,
                                       trait = trait)
}

png(paste0("outputs/postGWAS_sep/boxplot_sep.png"), width = 3000, height = 4000, res = 400)
grid.arrange(grobs = boxplot_list)
dev.off()

# marker density plot
setwd("outputs/postGWAS_sep/")
CMplot(pvalues[,1:4],
       bin.breaks= seq(0, 650, 100),
       plot.type="d",
       bin.size=1e4,
       chr.den.col=c("#43B284FF", "#FAB255FF", "#DD5129FF"),
       file="jpg",
       file.name="",
       dpi=300,
       main="Marker_density",
       file.output=T,
       verbose=F,
       width=9,
       height=6)
setwd(here())


# blast results

genes <- c('Mycgr3G107386',
           'Mycgr3G67942',
           'Mycgr3G47177',
           'Mycgr3G37599',
           'Mycgr3G30712',
           'Mycgr3G30347',
           'Mycgr3G107385',
           'Mycgr3G101090')
protein <- c('hypoothetical protein', 
             'related to sialidase',
             'voltage gated chloride channel like protein',
             'short-chain dehydrogenase like protein',
             'DNA polymerase epsilon subunit C',
             'putative 6-phosphogluconate dehydrogenase, NADP-binding, 6-phosphogluconate dehydrogenase',
             'cell wall biogenesis protein phosphatase Ssd1',
             'alpha/beta-Hydrolase like protein'
)
species <- c('Zymoseptoria brevis',
             'Ramularia collo-cygni',
             'Zymoseptoria brevis', 
             'Zymoseptoria brevis',
             'Pseudocercospora fuligena',
             'Septoria linicola',
             'Zymoseptoria brevis',
             'Zymoseptoria brevis'
)
percentage_identity <- c('95.85',
                         '71.23',
                         '97.66',
                         '52.84',
                         '84.21',
                         '65.62',
                         '99.06',
                         '95.50')
accesion <- c('KJY01215.1',
              'XP_023622909.1',
              'KJY02505.1',
              'KJX93753.1',
              'KAF7196700.1',
              'KAI5364816.1',
              'KJY01214.1',
              'KJY02094.1')

blast_df <- data.frame(zt_genes = genes,
                       protein = protein,
                       species = species,
                       percentage_identity = percentage_identity,
                       accesion = accesion
)

new_names <- c('Z.tritici genes', 'Protein', 'Species', '%Identity', 'Accesion')
names(blast_df) <- new_names
write_csv(blast_df, file = 'outputs/postGWAS_sep/blast_table.csv')


#===============================================================================
# Cultivar-sepecific hits 
#===============================================================================

hits <- read_csv('outputs/postGWAS_sep/gwas_sep_table_all.csv') |> 
  distinct(SNP)

dirs <- as.list(list.dirs('outputs/GWAS_sep', full.names = T, recursive = F))

tmp_dfs <- map(dirs, \(x) {
  tmp <- read_csv(paste0(x, '/GAPIT.Association.Filter_GWAS_results.csv')) |> 
    mutate(Cultivar = basename(x)) |> 
    dplyr::select(-1)
})

dfs_cultivars <- bind_rows(tmp_dfs) |> 
  mutate(traits = gsub('BLINK.', '', traits)) |> 
  filter(!grepl('log', traits))

pvalues_dfs <- map(dirs, \(x) {
  cultivar <- basename(x)
  tmp_files <- list.files(x, full.names = TRUE)
  tmp_files <- tmp_files[!grepl('Filter', tmp_files)]
  tmp_dfs <- map(tmp_files, read_csv) |> 
    bind_cols()
  
  # Seleccionar columnas según la cantidad de archivos
  cols_to_select <- if (length(tmp_files) == 1) {
    1:5
  } else if (length(tmp_files) == 2) {
    c(1:5, 12)
  } else if (length(tmp_files) == 3) {
    c(1:5, 12, 20)
  } else {
    1:5  
  }
  
  traits <- map(tmp_files, \(x) gsub(paste0('outputs/GWAS_sep/', cultivar,'/GAPIT.Association.GWAS_Results.BLINK.'), '', x))
  traits <- unlist(map(traits, \(x) gsub('.csv', '', x)))
  def_traits <- paste0(cultivar, '_', traits)
  tmp_dfs <- tmp_dfs |> 
    dplyr::select(all_of(cols_to_select)) |> 
    dplyr::select(-5)
  
  colnames(tmp_dfs) <- c('SNP', 'Chr', 'Pos', def_traits)
  
  return(tmp_dfs)
})

final_df <- pvalues_dfs |> 
  purrr::reduce(left_join, by = c('SNP', 'Chr', 'Pos'))
names(final_df) <- gsub('pycnidiaPerCm2Leaf', 'PCm2Leaf', names(final_df))
names(final_df) <- gsub('pycnidiaPerCm2Lesion', 'PCm2Lesion', names(final_df))

colors_trigo <- c("#8E6B3D", "#8E6B3D", "#D8B06A", "#D8B06A", "#D8B06A", "#5B4C44", "#5B4C44", "#9E5B40", "#9E5B40")
shapes <- c(15,16,17,15,16,15,16,15,16)

setwd("outputs/postGWAS_sep/")
CMplot(final_df,
       col = colors_trigo,
       plot.type="m",
       multraits=TRUE,
       threshold= 0.05/nrow(final_df),
       threshold.lty =1,
       threshold.lwd = 1,
       threshold.col=c("black"),
       amplify=F,
       pch = shapes, 
       bin.size=1e6,
       signal.cex=1,
       file="jpg",
       file.name="cultivars",
       dpi=300,
       file.output= T,
       verbose=TRUE,
       points.alpha=250,
       legend.ncol=5,
       legend.pos="left")
setwd(here())

gwas_cultivars_table <- dfs_cultivars |> 
  mutate(traits = case_when(
    traits == "pycnidiaPerCm2Leaf" ~ "PCm2Leaf",
    traits == "pycnidiaPerCm2Lesion" ~ "PCm2Lesion",
    TRUE ~ traits  # Mantener el valor original si no coincide con ninguna condición
  )) |> 
  dplyr::relocate(Cultivar, Traits = traits)

write_csv(gwas_cultivars_table, file = 'outputs/postGWAS_sep/cultivars_hit.csv')

# See which strains has the significant snps 
# identify all markers 
snps <- unique(c(all_l2$SNP, dfs_cultivars$SNP))
search_info <- function(df, marker) {
  new_df <- df |> 
    dplyr::select(marker, Isolate, Region) |> 
    filter(!!sym(marker) == 1) |> 
    mutate(Marker = marker) |> 
    dplyr::select(-!!sym(marker))
  return(new_df)
}

snps_df <- genotype_septoria[, colnames(genotype_septoria) %in% snps] |> 
  mutate(Isolate = genotype_septoria[, 1]) |> 
  left_join(info) |> 
  dplyr::select(Isolate, Region, everything()) 

search_results <- purrr::map_dfr(colnames(snps_df)[grepl("^X", colnames(snps_df))], 
                                 ~ search_info(snps_df, .))

isolates_per_variant <- search_results |> 
  group_by(Isolate)  |> 
  summarize(
    Count = n(),                            # Contar cuántas veces aparece cada Isolate
    Markers = toString(unique(Marker))       # Concatenar los nombres de los Markers para cada Isolate
  ) |> 
  arrange(desc(Count))

region_per_variant <- search_results |> 
  distinct(Isolate, Region) |>
  group_by(Region) |> 
  summarize(n = n())


