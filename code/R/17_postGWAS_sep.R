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
library(patchwork)
library(ape)
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
#===============================================================================
# Overall Hits on Leaf2
#===============================================================================
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


#===============================================================================
# Plot Manhattan and QQplot
#===============================================================================
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

#===============================================================================
# Plot LDHeatmap
#===============================================================================
snps <- colnames(genotype_septoria)[-1]
snps_df <- data.frame(SNPs = snps, 
                      Chromosome = as.numeric(str_extract(snps, 'X([0-9])_', group = T)),
                      Position = as.numeric(str_extract(snps, '_(.*)', group = T)))

hits <- plotLD(gwas_table = gwas_sep_table_all, genotype = genotype_septoria[, -1], df_info = snps_df)

#===============================================================================
# Plot Boxplot allelic differences
#===============================================================================
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

#===============================================================================
# Marker Density plot
#===============================================================================
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


#===============================================================================
# Blast Results
#===============================================================================

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
#===============================================================================
# Cultivar-sepecific hits 
#===============================================================================
#===============================================================================

# obtain all the hits for each cultivar specific analysis and arrange them in a df
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

# transform the dataframes in order to join the differente pvalues for the traits identififed in each 
# cultivar, and them arranged them to be able to plot a multi trait manhattan plot
pvalues_dfs <- map(dirs, process_pvalues)

# bind all the lists from the previous steps
final_df <- pvalues_dfs |> 
  purrr::reduce(left_join, by = c('SNP', 'Chr', 'Pos'))

# Renaming to match the paper traits names
names(final_df) <- gsub('pycnidiaPerCm2Leaf', 'PCm2Leaf', names(final_df))
names(final_df) <- gsub('pycnidiaPerCm2Lesion', 'PCm2Lesion', names(final_df))

final_df_plot <- final_df |> 
  dplyr::select(-matches('Effect')) |> 
  dplyr::select(-matches('MAF')) 

traits <- names(final_df_plot)[-c(1:3)]  # Exclude the first 3 columns
traits <- gsub('.*_', '', traits)
dir <- 'outputs/postGWAS_sep/'

map2(4:ncol(final_df_plot), traits, \(x, y) {
  fastqqplot(x, final_df_plot, y, dir)
})


#===============================================================================
# Post Gwas on the cultivar-sepcific traits
#===============================================================================

gwas_cultivars_table <- dfs_cultivars |> 
  mutate(traits = case_when(
    traits == "pycnidiaPerCm2Leaf" ~ "PCm2Leaf",
    traits == "pycnidiaPerCm2Lesion" ~ "PCm2Lesion",
    TRUE ~ traits  # Mantener el valor original si no coincide con ninguna condición
  )) |> 
  dplyr::relocate(Cultivar, Traits = traits)
  
# add also effects
effects <- final_df |> 
  dplyr::select(c(SNP, matches('_Effect'))) |> 
  filter(SNP %in% gwas_cultivars_table$SNP) |> 
  pivot_longer(cols = -c(SNP),   # Keep the first 3 columns
               names_to = "Traits",         # Name of the new column for traits
               values_to = "Effect") |> 
  mutate(Traits = gsub('_Effect$', '', Traits))  |> 
  separate(Traits, into = c("Cultivar", "Traits"), sep = "_") |> 
  mutate(Traits = case_when(
    Traits == "pycnidiaPerCm2Lesion" ~ "PCm2Lesion", 
    Traits == "pycnidiaPerCm2Leaf" ~ "PCm2Leaf",
    TRUE ~ Traits  # Keeps other Trait values unchanged if needed
  )) |> 
  left_join(gwas_cultivars_table, by = c('Cultivar', 'Traits', 'SNP')) |> 
  na.omit()
  
# define everything you need to obtain genes and GOs associated to those markers
attributes <- c("ensembl_gene_id", "start_position", "end_position", "strand", "description", "peptide")
filters <- c("chromosome_name", "start", "end")
distances <- list(500, 1000, 2000) #base pairs (bp)
chr_all <- gwas_cultivars_table$Chr
traits_all <- gwas_cultivars_table$Traits
markers_all <- as.numeric(map(strsplit(gwas_cultivars_table$SNP, "_"), \(x) x[[2]]))
out_dir_all <- 'outputs/postGWAS_sep/FASTA_cultivars'
cultivars_genes <- find_genes(mart = mart,
                                  attributes = attributes,
                                  filters = filters,
                                  distances = distances,
                                  chr = chr_all,
                                  traits = traits_all,
                                  out_dir = out_dir_all,
                                  markers = markers_all)

inter_table <- effects |> 
  left_join(cultivars_genes, by = c('SNP' = "marker")) |> 
  mutate(GO_class = Ontology(GO_id))

# Create GWas table
gwas_table <- inter_table |> 
  dplyr::select(-start, end_position, GO_id:GO_class) |> 
  dplyr::select(Cultivar, Traits, SNP, P.value, Effect, MAF, Distance = distance, Gene = gene, Protein = accesion) |> 
  distinct()
# Create genes table
genes_table <- inter_table |> 
  dplyr::select(-c(trait, SNP, distance, Cultivar, Traits, Pos, P.value, MAF, pos, GO_id, GO_class, Go_name)) |> 
  distinct() |> 
  dplyr::relocate(Gene = gene, Chr) |> 
  arrange((as.numeric(Chr))) |> 
  na.omit()

# save tables
write_csv(gwas_table, file = 'outputs/postGWAS_sep/cultivars_hit.csv')

# Manhattan plot 

mahattan_plot <- final_df_plot |> 
  filter(SNP %in% unique(gwas_table$SNP) | row_number() %% 10 == 0) |> 
  mutate(across(-c(SNP, Chr, Pos), as.numeric)) |>
  rowwise() |> 
  mutate(
    Athoris = min(c_across(starts_with("Athoris")), na.rm = TRUE),
    `Don Ricardo` = min(c_across(starts_with("Don Ricardo")), na.rm = TRUE),
    Sculptur = min(c_across(starts_with("Sculptur")), na.rm = TRUE),
    Svevo = min(c_across(starts_with("Svevo")), na.rm = TRUE)
  ) |> 
  ungroup() |> 
  dplyr::select(SNP, Chr, Pos, Athoris, `Don Ricardo`, Sculptur, Svevo)

add_snp <- final_df_plot |> 
  filter(SNP == 'X11_849459') |> 
  dplyr::select(SNP, Chr, Pos, Athoris = Athoris_PCm2Leaf) |> 
  mutate(SNP = 'X11_849460')
mahattan_plot <- bind_rows(mahattan_plot, add_snp)

colors_cultivars <- c('Athoris' = "#8E6B3D",
                      'Don Ricardo' = "#D8B06A",
                      'Sculptur' = "#5B4C44",
                      'Svevo' = "#9E5B40")
df_snps <- gwas_table |> 
  distinct(Cultivar, SNP, Traits) |> 
  mutate(color = colors_cultivars[Cultivar]) |> 
  arrange(Cultivar)
df_snps[1, 'SNP'] <- 'X11_849460'

snp_list <- df_snps %>%
  group_by(Cultivar) %>%
  summarise(SNPs = list(SNP), .groups = 'drop') %>%
  deframe()
color_list <- df_snps |> 
  group_by(Cultivar) |> 
  summarise(color = list(unique(color)), .groups = 'drop') %>%
  deframe()


shape_vector <- df_snps %>%
  mutate(Shape = case_when(
    Traits == "PCm2Leaf" ~ 16,
    Traits == "PCm2Lesion" ~ 17,
    Traits == "PLACL" ~ 18
  )) |> 
  group_by(Cultivar) |> 
  summarise(Shape = list((Shape)), .groups = 'drop') %>%
  deframe()


setwd("outputs/postGWAS_sep/")
CMplot(mahattan_plot,
       col = c('#98FB98', '#3CB371'),
       type = 'p',
       plot.type="c",
       outward = T, 
       multraits=TRUE,
       H = 5, 
       LOG10 = T, 
       chr.labels=paste("Chr",c(1:13),sep=""), 
       threshold= 0.05/nrow(final_df_plot),
       threshold.lty =1,
       threshold.lwd = 1,
       threshold.col=c("black"),
       amplify=T,
       highlight = snp_list, 
       highlight.col = color_list, 
       highlight.pch = 16,
       highlight.cex = 1,
       signal.cex = 0,
       file="jpg",
       file.name="cultivars",
       dpi=300,
       file.output= T,
       points.alpha=250,
       legend.ncol=5,
       legend.pos="left", 
       cir.axis.col = 'black')
setwd(here())


#===============================================================================
# Identify which strains have each significant snps
#===============================================================================
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

#===============================================================================
# merge table with new JGI annotation
#===============================================================================

genes_table <- genes_table |> 
  mutate(Gene = gsub('G', '_', Gene))

translate_genes <- read_csv('data/raw_data/new_annotations.csv')
translate_genes <- translate_genes |>
  separate(JGI_genes, into = c('delete', 'Gene'), sep = '\\|') |> 
  filter(Gene %in% genes_table$Gene) |> 
  janitor::clean_names() |> 
  dplyr::select(gene:seq_length, gene_2) |> 
  dplyr::relocate(Ensembl = gene_2)

fasta_file <- "data/raw_data/z.tritici.IP0323.reannot.proteins.fasta"
fasta <- read.fasta(file = fasta_file, seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
names(fasta) <- gsub('\\.1$', '', names(fasta))
subset <- fasta[names(fasta) %in% translate_genes$gene]
write.fasta(subset, names(subset), file.out="outputs/postGWAS_sep/FASTA_cultivars/new_proteins.fasta")

effectorp <- readxl::read_excel('outputs/postGWAS_sep/EffectorP.xlsx') |> 
  janitor::remove_empty() |> 
  janitor::clean_names() |> 
  mutate(number_identifier = gsub('\\.1$', '', number_identifier)) |> 
  dplyr::rename(gene = number_identifier) |> 
  dplyr::select(-c(non_effector, prediction))
signalp <- read_table('outputs/postGWAS_sep/SignalP.txt', col_names = T) |> 
  janitor::clean_names() |> 
  mutate(id = gsub('\\.1$', '', id),
         sp_sec_spi = round(sp_sec_spi, 2)) |> 
  dplyr::select(gene = id, SignalP = sp_sec_spi)

lines <- readLines('outputs/postGWAS_sep/TMRs.gff3')
tmrs_df <- tibble(
    line = lines
  ) |> 
  filter(str_detect(line, "Number of predicted TMRs")) |> 
  mutate(
    gene = str_extract(line, "#\\s*(\\S+)"),  # Extrae el ID de la proteína
    tmrs = as.numeric(str_extract(line, "(?<=Number of predicted TMRs:\\s)\\d+"))  # Extrae el número de TMRs
  ) |> 
  dplyr::select(gene, tmrs) |> 
  mutate(gene = gsub('# ', '', gene))

fasta_files <- as.list(grep('.faa', list.files('outputs/postGWAS_sep/FASTA_cultivars/', full.names = T), value = T))
aa_length <- map(fasta_files, count_characters)
length_df <- data.frame(Ensembl = unlist(map(fasta_files, \(x) gsub('.faa', '', basename(x)))),
                        Ensembl_length = unlist(aa_length)) |> 
  mutate(Ensembl = gsub('G', '_', Ensembl))

IPR <- readxl::read_xlsx('outputs/postGWAS_sep/IPR_cultivars.xlsx')


final_gene_table <- translate_genes |> 
  left_join(effectorp) |> 
  left_join(signalp) |> 
  left_join(length_df) |> 
  left_join(tmrs_df) |> 
  left_join(IPR) |> 
  dplyr::relocate(Ensembl, Ensembl_length)
  
write_csv(final_gene_table, file = 'outputs/postGWAS_sep/cultivars_genes.csv')

#===============================================================================
# BBoxplot per region + pca
#===============================================================================
load('data/modified_data/8_pca.Rdata')
gh1_pheno <- read_csv('data/modified_data/gh1_clean_pheno.csv') |> 
  left_join(info) |> 
  mutate(across(c(Line, Region), as.factor)) |> 
  filter(Leaf == 2)
info <- read_csv('data/modified_data/info_strain_complete.csv') |> 
  dplyr::select(Isolate, Region)
pca <- pca_data |> 
  left_join(info, by = c('GenoID' = 'Isolate'))

colors <- c("Córdoba" = "#DD5129FF",
            "Cádiz" = "#0F7BA2FF",
            "Sevilla" = "#43B284FF",
            "Huelva" = "#FAB255FF", 
            'Non_selected' = 'grey')

load('data/modified_data/blues_cultivars/athoris.RData')
athoris_blues <- BLUPS
load('data/modified_data/blues_cultivars/don_Ricardo.RData')
donR_blues <- BLUPS
load('data/modified_data/blues_cultivars/Sculptur.RData')
sculptur_blues <- BLUPS
load('data/modified_data/blues_cultivars/Svevo.RData')
svevo_blues <- BLUPS

translate_names <- data.frame(Isolate = genotype_septoria[,1], 
                              lower_names = tolower(genotype_septoria[,1])) |> 
  left_join(info) |> 
  dplyr::rename(true_names = 'Isolate')

blues_df <- list(athoris_blues, donR_blues, sculptur_blues, svevo_blues) |> 
  map(\(x) x |> dplyr::select(1:4)) |> 
  map(\(x) x |> 
        left_join(translate_names, by = c('Isolate' = 'lower_names')) |> 
        dplyr::select(Isolate, PLACL:pycnidiaPerCm2Lesion, Region)
        ) 

blues_list <- list(blues_df[[1]], blues_df[[1]],
                   blues_df[[2]], blues_df[[2]], blues_df[[2]], blues_df[[2]], blues_df[[2]], blues_df[[2]], blues_df[[2]], blues_df[[2]],blues_df[[2]],
                   blues_df[[3]], blues_df[[3]], 
                   blues_df[[4]], blues_df[[4]])
# 
# plot <- pmap(list(snps, traits, cultivars, blues_list), \(x,y,z,k) pca_boxplot(pca, genotype, colors, x, y, z, k))
# out_dir <- 'outputs/postGWAS_sep/pca_boxplot'
# dir.create(out_dir)
# 
# map2(plot, snps, \(x, y) {
#   png(paste0(out_dir, '/', y, ".png"), width = 6000, height = 4000, res = 400)
#   print(x)  # Imprimir el gráfico para que se guarde
#   dev.off()  # Cerrar el dispositivo gráfico
# })

#===============================================================================
# Boxplot per cultivar
#===============================================================================

snps <- as.list((dfs_cultivars$SNP))
cultivars <- as.list((dfs_cultivars$Cultivar))
traits <- as.list((dfs_cultivars$traits))
genotype <- genotype_septoria[, colnames(genotype_septoria) %in% c('Isolate', snps)]
genotype[, -1] <- round(genotype[,-1], 0)
genotype$Isolate <- tolower(genotype$Isolate)

cultivars_boxplot <- pmap(list(snps, traits, cultivars, blues_list), \(x,y,z,k) boxplot_cultivars(genotype, colors, x, y, z, k))

png(paste0("outputs/postGWAS_sep/cultivar_pt1.png"), width = 5000, height = 10000, res = 400)
grid.arrange(grobs = cultivars_boxplot[1:15], ncol = 3)
dev.off()






