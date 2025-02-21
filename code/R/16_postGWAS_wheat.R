library(tidyverse)
library(ape)
library(BiocManager)
library(GO.db)
library(paletteer)
library(hrbrthemes)
library(gridExtra)
library(biomaRt)
library(ggdist)
library(CMplot)
library(here)
library(LDcorSV)
source("code/R/function_septoria_GS.R")

#==============================================================================
# Find GENES related to HITS from the GWAS
#==============================================================================

gff <- ape::read.gff('data/raw_data/Triticum_aestivum.IWGSC.59.gff3') %>%
  filter(type == 'gene') %>% # select only genes
  dplyr::select(seqid, start, end, attributes) %>% 
  mutate(Gene = map_chr(strsplit(attributes, ";"), ~ { # subset gene ID from attribute list
    strsplit(.x[grep("^gene_id=", .x)], "=")[[1]][2]}),
    Chr = paste0("Chr", seqid), # create a character chromosome name
    across(c(start, end), as.numeric)
  ) %>% 
  dplyr::select(Gene, Chr, start, end) %>%
  arrange(Chr, start)

results_gwas_4PC <- read_csv("outputs/GWAS_wheat/PC4/GAPIT.Association.Filter_GWAS_results.csv")
hits <- results_gwas_4PC %>% 
  mutate(traits_clean = gsub("BLINK.", "", traits)) %>% 
  dplyr::select(-c(...1, traits))

max_distance <- 5e4
results <- data.frame(SNP = character(), Chr = integer(), Pos = integer(),
                      Gene = character(), Distance = numeric(), traits = character(),
                      stringsAsFactors = FALSE)
for (i in 1:nrow(hits)) {
  # Subset the genes on the same chromosome as the current GWAS hit
  chr_genes <- gff %>% filter(Chr == hits$Chr[i])
  # Loop through each gene on the chromosome
  for (j in 1:nrow(chr_genes)) {
    gene <- chr_genes[j, ]
    # Calculate the distance from the SNP to the gene
    if (hits$Pos[i] >= gene$start & hits$Pos[i] <= gene$end) {
      # If the SNP is within the gene, distance is 0
      distance <- 0
    } else if (hits$Pos[i] < gene$start) {
      # SNP is upstream of the gene, distance is the difference from start of gene
      distance <- gene$start - hits$Pos[i]
    } else if (hits$Pos[i] > gene$end) {
      # SNP is downstream of the gene, distance is the difference from end of gene
      distance <- hits$Pos[i] - gene$end
    }
    
    # Keep genes with distance 0 (SNP inside the gene) or within the max distance
    if (distance == 0 || distance <= max_distance) {
      results <- rbind(results, data.frame(SNP = hits$SNP[i],
                                           Chr = hits$Chr[i],
                                           Pos = hits$Pos[i],
                                           Gene = gene$Gene,
                                           Start = gene$start,
                                           End = gene$end,
                                           Distance = distance, 
                                           traits = hits$traits_clean[i]))
    }
  }
}

#==============================================================================
# Find GOs related to GENES
#==============================================================================

unique_genes <- unique(results$Gene)
GOs <- list()
remove_by_exp <- c()

for (i in seq_along(unique_genes)) {
  url_gene <- unique_genes[i]
  if (!i %% 25) { print(i) }
  
  if (!url_gene %in% remove_by_exp) {
    url_init = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgo%2Cprotein_name%2Cprotein_existence&format=tsv&query=%28"
    url_fin = "%29"
    url <- paste0(url_init, url_gene, url_fin)
    
    r <- httr::GET(url)
    
    # Usar tryCatch para manejar errores en GOwide2long
    GOs[[i]] <- tryCatch(
      {
        data.frame(Gene = url_gene, GOwide2long(content2table(r)))
      },
      error = function(e) {
        # Agregar el gene a remove_by_exp en caso de error
        remove_by_exp <<- c(remove_by_exp, url_gene)
        return(NULL)  # Retorna NULL si hay un error
      }
    )
  }
}
GO_final <- do.call('rbind', GOs)
rownames(GO_final) = NULL

#==============================================================================
# Find classification of the GOs found
#==============================================================================

GO_final$GO_code <- sapply(GO_final$GO_code, function(x) gsub("^\\[(.*)\\]$", "\\1", x))
GO_final$GO_class <- sapply(GO_final$GO_code, function(x) Ontology(x))

results_df <- hits %>% 
  left_join(results, by = "SNP") %>% 
  left_join(GO_final, by = "Gene") %>% 
  dplyr::select(traits_clean, SNP, Chr.x, Pos.x, P.value,
                Gene, Entry, Prot_name, GO_code, GO_name, GO_class) %>% 
  group_by(SNP, traits_clean) %>% 
  distinct()

plotGO_class <- results_df %>% 
  na.omit() %>% 
  mutate(GO_class = factor(GO_class)) %>% 
  group_by(traits_clean, GO_class) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  group_by(traits_clean) %>%
  mutate(percentage = count / sum(count) * 100)

colors <- c(
  "BP" = "#795548FF", 
  "CC" = "#F4A261FF", 
  "MF" = "#2A9D8FFF"  
)

ggplot(plotGO_class, aes(x = traits_clean, y = percentage, fill = GO_class)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +  # Adjust bars to sum to 100%
  scale_y_continuous(labels = scales::percent_format()) +  # Show percentages on the Y-axis
  geom_text(
    aes(label = count),  # Column that contains the counts
    position = position_fill(vjust = 0.5),  # Position the labels in the center of each segment
    color = "black",  # Color of the text
    size = 7  # Text size
  ) +
  labs(
    title = "Gene Ontology (GO) categories",
    x = NULL,
    y = "Percentage",
    fill = "GO_class"  # This will be used for the fill legend title
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # Elimina las líneas mayores del eje X
    panel.grid.minor.x = element_blank(),  # Elimina las líneas menores del eje X
    panel.grid.major.y = element_line(color = "gray", size = 0.5),  # Mantiene las líneas menores del eje Y
    plot.title = element_blank(),  # Elimina el título del gráfico
    axis.title.x = element_text(size = 18, face = "bold", vjust = -2, margin = margin(b = 40), hjust = 0.5),  # Centra el texto del título X y ajusta el margen
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20), hjust = 0.5),  # Centra el texto del título Y
    axis.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),  # Centra el texto de los ticks del eje X
    axis.text.y = element_text(size = 18, hjust = 0.5),  # Centra el texto de los ticks del eje Y
    legend.position = "top",  # Coloca la leyenda en la parte superior
    legend.title = element_blank(),  # Elimina el título de la leyenda
    legend.text = element_text(size = 22)  # Tamaño del texto de la leyenda
  ) +
  scale_fill_manual(
    values = colors) +
  guides(fill = guide_legend(title = NULL)) +
  coord_flip()

#==============================================================================
# Final Tables
#==============================================================================

gwas_table <- hits %>% 
  left_join(results, by = "SNP") %>% 
  dplyr::select(Trait = traits_clean, SNP, 
                Chr = Chr.x, Pos = Pos.x, P.value, MAF, Gene, Distance) %>%
  distinct()

effect_list <- list()
for(i in 1:nrow(gwas_table)){
  row <- gwas_table[i, ]
  marker <- row$SNP
  trait <- row$Trait
  model <- "BLINK"
  parent_directory <- "outputs/GWAS_wheat/PC4/"
  effect_list[[i]] <- grepEffect(marker = marker,
                                 trait = trait,
                                 model = model,
                                 parent_directory = parent_directory)
}
effect_df <- do.call(rbind, effect_list)
effect_df <- effect_df %>% distinct()

gwas_table <- gwas_table %>% 
  left_join(effect_df, by = c("SNP", "Trait"))

# Change the name of the markers from IBW... to its assigned names
markers_names <- read_delim("data/raw_data/90K.CSv2.Ann.info.txt", delim = "\t")
markers_names <- markers_names %>% dplyr::select(ID, Name)
gwas_table2 <- gwas_table %>% 
  left_join(markers_names, by = c("SNP" = "ID")) %>% 
  dplyr::select(Trait, SNP = Name, Chr, Pos, P.value, MAF, Effect, Gene, Distance)

write_csv(gwas_table2, file = "outputs/postGWAS_wheat/gwas_wheat_table.csv")

attr <- c("ensembl_gene_id", "description")
triticum <- useEnsemblGenomes(biomart = "plants_mart", dataset = "taestivum_eg_gene")
description <- getBM(attributes = attr,
                     filters = "ensembl_gene_id",
                     values = GO_final$Gene,
                     mart = triticum)

genes_GO_table <- results %>% 
  left_join(GO_final, by = "Gene") %>% 
  left_join(description, by = c("Gene" = "ensembl_gene_id")) %>% 
  dplyr::select(Trait = traits, SNP, Gene, Entry, Prot_name, Start, End, GO_code, GO_name, GO_class) %>% 
  distinct()

genes_GO_table <- genes_GO_table %>% 
  left_join(markers_names, by = c("SNP" = "ID")) %>% 
  dplyr::select(Trait, SNP = Name, Gene, Protein = Entry, Name = Prot_name, Start, End, GO_code, GO_name, GO_class)

write_csv(genes_GO_table, file = "outputs/postGWAS_wheat/genes_GO.csv")

#==============================================================================
# Plot Allelic Diff
#==============================================================================

load("data/modified_data/3_wheat_GWAS.Rdata")

boxplot_list <- list()
for(i in 1:nrow(gwas_table)){
  marker <- gwas_table$SNP[i]
  trait <- gwas_table$Trait[i]
  print(paste0("working on marker: ", marker))
  
  phenotype <- blues_wheat %>%  dplyr::select(GenoID, trait)
  genotype <- genotype_wheat[genotype_wheat[, "GenoID"] %in% phenotype[,1], ]
  
  boxplot_list[[i]] <- plotAllelicdiff(phenotype = phenotype,
                                       genotype =  genotype,
                                       marker = marker,
                                       trait = trait)
}

png(paste0("outputs/postGWAS_wheat/boxplot_wheat.png"), width = 5000, height = 2000, res = 400)
grid.arrange(grobs = boxplot_list, ncol = 2)
dev.off()

#==============================================================================
# Plot Manhattan and QQplot
#==============================================================================

pvalues <- read_csv("outputs/GWAS_wheat/PC4/GAPIT.Association.GWAS_Results.BLINK.PLACL.csv")
pvalues <- pvalues[,1:4]
png(paste0("outputs/plots/boxplot_wheat"), width = 5000, height = 2000, res = 400)

setwd('outputs/postGWAS_wheat/')
CMplot(pvalues,
       plot.type = "m",
       multraits = FALSE,
       col = c('grey', 'darkgrey'),
       threshold = 0.05/nrow(pvalues),
       threshold.lty = c(1, 2), 
       threshold.lwd = c(1, 1),
       threshold.col = c("black", "grey"),
       amplify = TRUE,
       bin.size = 1e6,
       chr.den.col = NULL,
       highlight = unique(gwas_table$SNP),
       highlight.text = unique(gwas_table2$SNP), 
       highlight.col = '#F4A261FF', 
       signal.cex = 1,
       file = "jpg",
       file.name = 'PLACL',
       dpi = 300,
       file.output = T,
       verbose = FALSE,
       points.alpha = 250,
       legend.ncol = 2,
       legend.pos = "middle",
       main = 'PLACL')
setwd(here())

# QQ plot
plotCMqq(
  df = pvalues,
  name = "PLACL",
  color = '#F4A261FF',
  dir = "outputs/postGWAS_wheat/")

#==============================================================================
# Histogram for BLUEs
#==============================================================================

traits <- c("PLACL" = "#F4A261FF",
            "pycnidiaPerCm2Leaf" = "#795548FF",
            "pycnidiaPerCm2Lesion" = "#2A9D8FFF")

hist <- plotHist(blues_wheat,
                 columns = colnames(blues_wheat)[-1],
                 color = traits)

png(paste0("outputs/plots/hist_blues.png"), width = 5000, height = 2000, res = 400)
grid.arrange(grobs = hist, ncol = 3)
dev.off()

#==============================================================================
# Which individuals have the snps
#==============================================================================

IWB11991_indv <- genotype_wheat[, c('GenoID', colnames(genotype_wheat)[colnames(genotype_wheat) %in% gwas_table$SNP])] |> 
  dplyr::select(1,2) |> 
  filter(IWB11991 == 1 | IWB11991 == 2) |> 
  pull(GenoID)

IWB74799_indv <-  genotype_wheat[, c('GenoID', colnames(genotype_wheat)[colnames(genotype_wheat) %in% gwas_table$SNP])] |> 
  dplyr::select(1,3) |> 
  filter(IWB74799 == 0) |> 
  pull(GenoID)

cultivars <- c('PyrSep_105', 'PyrSep_107', 'PyrSep_115')
cross1 <- IWB11991_indv[IWB11991_indv %in% cultivars]
cross2 <- IWB74799_indv[IWB74799_indv %in% cultivars]

 