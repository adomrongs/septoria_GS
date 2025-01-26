library(sf)
library(ggmap)
library(tidyverse)
library(factoextra)
library(ggrepel)
library(RColorBrewer)
library(ggiraph)
library(vcfR)
library(ggside)
library(scales)
source("code/R/function_septoria_GS.R")

load('data/modified_data/1_septoria.Rdata')

#===============================================================================
# Create data to plot map
#===============================================================================

# extract sf
my_sf <- read_sf("data/raw_data/spain-provinces.geojson")
andalucia <- my_sf |> 
  filter(cod_ccaa == "01") # keep just provinces in andalucia

# prepare isolate to match info
places <- genotype_septoria[,1] |> 
  str_split("_") |> 
  map_chr(\(x) x[[2]]) |>  # Extraemos la segunda parte después de "_"
  map(\(name) {
    name <- str_remove_all(name, "[0-9]")  # Elimina cualquier número del nombre
    uppercase_positions <- str_locate_all(name, "[A-Z]")[[1]]  # Encuentra todas las mayúsculas
    if (nrow(uppercase_positions) < 2) return(c(name, ""))  # Si no hay 2 mayúsculas, dejar como está
    split_point <- uppercase_positions[2, 1]  # Toma la posición de la segunda mayúscula
    c(str_sub(name, 1, split_point - 1), str_sub(name, split_point, nchar(name)))  # Divide en dos partes
  }) |> 
  map_chr(\(x) x[[1]])

# prepare regions
regions <- readxl::read_xlsx('data/raw_data/regions_strain.xlsx') |> 
  dplyr::select(1, 3, 4) |> 
  distinct() |> 
  mutate(Strain = map_chr(Strain, \(name) {
    name <- str_remove_all(name, "[0-9]")  # Elimina números
    uppercase_positions <- str_locate_all(name, "[A-Z]")[[1]]  # Encuentra posiciones de mayúsculas
    if (nrow(uppercase_positions) < 2) return(name)  # Si no hay 2 mayúsculas, dejar igual
    split_point <- uppercase_positions[2, 1]  # Segunda mayúscula
    str_sub(name, 1, split_point - 1)  # Extraer solo la primera parte
  })) |> 
  distinct()
# match isolates, egions and localidades
names_correct <- data.frame(table(places)) |> 
  dplyr::rename(Places = places, n = Freq) |> 
  left_join(regions, by = c("Places" = "Strain")) |> 
  distinct() |> 
  mutate(Place = dplyr::recode(Place, 
                               "Finca Montera (Utrera)" = "Utrera",
                               "Finca Vega del Vivero (Jerez)" = "Jerez",
                               "Finca Montana (SANLÚC. B.)" = "Sanlucar de Barrameda")) |> 
  dplyr::rename(Localidad = Place) |> 
  dplyr::slice(-c(6,7))

# extract coordinates for Locations
geo_data <- geocode(names_correct$Localidad, output = "more")
geo_data[9, ] <- geocode("21890 Manzanilla, Huelva", output = "more")
geo_data[4, ] <- geocode("Cordoba, Spain", output = "more")
final_df <- cbind(names_correct, geo_data) |> 
  mutate(Region = as.factor(Region))

colors <- c("Córdoba" = "#DD5129FF",
            "Cádiz" = "#0F7BA2FF",
            "Sevilla" = "#43B284FF",
            "Huelva" = "#FAB255FF")

#===============================================================================
# plot boxplot per regions
#===============================================================================
# load info
info <- read_csv('data/modified_data/info_strain_complete.csv')

# associate each isolate to a region
gh1_pheno <- read_csv('data/modified_data/gh1_clean_pheno.csv') |> 
  left_join(info) |> 
  mutate(across(c(Line, Region), as.factor))

gh1_pheno |>
  dplyr::select(Isolate, Region) |> 
  distinct() |> 
  count(Region)
# count of how many isolates per region
counts <- gh1_pheno |>
  dplyr::select(Isolate, Region) |> 
  distinct() |> 
  count(Region) |> 
  collapse(Isolate, n)
# count for map too 
count2 <- final_df |> 
  group_by(Localidad, lon, lat, Region) |> 
  summarise(n = n()) |> 
  mutate(colors = colors[Region])

#plot map ( we needed the count)
map <- ggplot(andalucia) +
  geom_sf(fill = "lightgrey", color = "white") +
  theme_void() + 
  geom_point(data = final_df, 
             aes(x = lon, y = lat, size = n, fill = Region),  # Color de relleno según Region
             shape = 21,      # Círculo con borde
             stroke = 0.3,    # Grosor del borde negro
             color = "black") +  # Color del borde de los puntos
  geom_text_repel(data = count2,
                  aes(x = lon, y = lat, label = Localidad), 
                  size = 4,           # Tamaño del texto
                  box.padding = 1.5,  # Espaciado alrededor de la etiqueta
                  point.padding = 0.4,  # Espaciado entre el punto y la etiqueta
                  segment.size = 0.5,   # Grosor del segmento
                  segment.curvature = 0.3,  # Curvatura del segmento (puedes ajustarlo)
                  min.segment.length = 0.1,  # Longitud mínima del segmento
                  force = 5, show.legend = F)  +  
  scale_fill_manual(values = colors,  # Usa 'fill' en vez de 'color'
                    labels = paste0(counts$Region," (n = ", counts$n, ")"), 
                      guide = guide_legend(ncol = 4, 
                                           title = NULL,
                                           position = "top", 
                                           color = NULL, 
                                           override.aes = list(size = 6, color = NULL), 
                                           theme = theme(
                                             legend.text = element_text(size = 17), 
                                             legend.key.size = unit(4, "cm")
                                           )
                      )
  ) + 
  scale_size_continuous(name = "Number of Isolates", 
                        range = c(2, 10),  
                        guide = guide_legend(ncol = 3,
                                             position = 'bottom', 
                                             title.position = "top")) +
  theme(
    legend.box = "vertical",
    legend.box.just = "center"
  )

png(paste0("outputs/plots/map.png"), width = 4000, height = 3000, res = 400)
map
dev.off()


# plot boxplot
plot_regions <- ggplot(gh1_pheno) +
  geom_boxplot(aes(x = Region, y = PLACL, fill = Region), width = 0.2, size = 0.5) +
  scale_fill_manual(values = colors) +  
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.4),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16, margin = margin(t = 10)),  # Espacio para etiquetas
    axis.text.y = element_text(size = 13), 
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_discrete(labels = paste0(counts$Region," (n = ", counts$n, ")"))

png(paste0("outputs/plots/boxplot_regions.png"), width = 5000, height = 1500, res = 400)
plot_regions
dev.off()

#===============================================================================
# plot PCA with side density plot
#===============================================================================
load("data/modified_data/1_septoria.Rdata")
shapes2 <- c("21" = 17,
             "22" = 16)
pca_df <- info |> 
  filter(Isolate %in% genotype_septoria[,1]) # just for the first 100 isolates

# perform PCA
PCA <- computePCA(genotype = genotype_septoria,
                          regions = pca_df$Region,
                          colors = colors,
                          shapes = shapes2,
                          names = pca_df$Isolate,
                          interactive = TRUE, 
                          shape_col = pca_df$year
)

PCA <- tibble(PCA) |> 
  dplyr::select(PC1, PC2, GenoID, regions, shape, Names, other) |> 
  mutate(across(regions:other, as.factor))
  
# plot PCA with density sides
pca_plot <- ggplot(data = PCA) +
  geom_point(aes(x = PC1, y = PC2, color = regions, shape = other), size = 3) +
  geom_xsidedensity(aes(x = PC1, y = after_stat(density), color = regions), 
                    alpha = 1, size = 0.9, show.legend = F, outline.type = "upper") +
  geom_ysidedensity(aes(y = PC2, x = after_stat(density), color = regions), 
                    alpha = 1, size = 0.9, show.legend = F, outline.type = "upper") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "PC1 (3.5%)",
    y = "PC2 (2.5%)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 18),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 17),
    strip.text = element_text(size = 10, face = "plain", color = "black", hjust = 0.5),
    strip.background = element_rect(fill = "lightgray"),
    panel.border = element_rect(colour = "darkgrey", fill = NA, size = 0.5),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10),
    ggside.panel.border = element_blank(),
    ggside.axis.text.y = element_blank(),
    ggside.axis.text.x = element_blank(), 
    ggside.panel.scale = 0.2
    ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +  # Para que la densidad tenga los mismos colores
  scale_shape_manual(values = shapes2, 
                     labels = c("2021", "2022"))


png(paste0("outputs/plots/PCA_sep_complete2.png"), width = 3000, height = 3000, res = 400)
pca_plot
dev.off()

#===============================================================================
# Compute LD decay
#===============================================================================
dir.create('outputs/lists')

# We want to run code/bash/compute_ld.sh. There, will compute there the LD between all  
# markers pairs in 20kb windows separately for each Region. But first we have to 
# create lists containing the individuals contained in each Region but in fasta format
samples_fasta <- read_delim('outputs/lists/samples.txt', delim = '\t', col_names = F)
samples_fasta <- samples_fasta |> 
  mutate(Fasta = gsub("_L.*", "", X1),
         Fasta = gsub("_m.*", "", Fasta)) |> 
  left_join(info) |> 
  dplyr::rename(full = X1)

regions <- unique(samples_fasta$Region)
for(i in seq_along(regions)){
  region <- regions[[i]]
  tmp_df <- samples_fasta |> 
    filter(Region == region) |> 
    dplyr::select(full)
  
  write.table(tmp_df, file = paste0('outputs/lists/', region,'.list'), col.names = F, row.names = F, quote = F)
}

# Now, usinfg plink we would have generated .ld files for each Region
tables <- list.files('outputs/ld', full.names = T)
tables_data <- map(tables, \(x) read.table(x, header = T, stringsAsFactors = F))
names <- list("Córdoba", "Cádiz", "Sevilla", "Huelva")
# here we will apply the HW function to calculate LD decay for Region and be able to plot it
ld_dfs <- map2(tables_data, names, computeLD)
ld_combined <- bind_rows(ld_dfs) |> 
  mutate(Region = as.factor(Region))

# where (in distance) does LD decay reach 0.2 R2
ld_threshold <- ld_combined %>%
  filter(LD <= 0.2) %>%            # Filtrar valores donde LD <= 0.2
  group_by(Region) %>%             # Agrupar por región
  slice_min(distance) %>%          # Extraer la primera distancia mínima
  ungroup()

#plot LD decay per Region
ld_decay_plot <- ggplot(ld_combined) +
  geom_line(aes(x = distance, y = LD, color = Region), linewidth = 1, show.legend = FALSE) +  
  geom_point(data = ld_threshold, aes(x = distance, y = LD, color = Region), size = 3) +  
  scale_color_manual(values = colors,
                     name = "Distance to R² < 0.2",  # Cambia el título de la leyenda
                     labels = paste0(ld_threshold$Region, ': ', ld_threshold$distance)) +  # Personaliza las etiquetas de la leyenda
  labs(title = "LD Decay by Region", x = "Distance (bp)", y = "Predicted R²") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "darkgrey"),  
    axis.ticks = element_line(color = "darkgrey"),  
    axis.ticks.length = unit(6, "pt"),  
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = c(0.75, 0.75),  # Cambia esto para mover la leyenda dentro del gráfico (ajusta las coordenadas)
    legend.background = element_rect(color = NA),  # Fondo semi-transparente
    legend.key = element_blank(),  # Quita el fondo de los ítems de la leyenda
    legend.text = element_text(size = 12),
    panel.border = element_blank()  # Elimina el borde del panel entero
  )


png(paste0("outputs/plots/LD_decay.png"), width = 3000, height = 3000, res = 400)
ld_decay_plot
dev.off()

#===============================================================================
# plot SNP eff distribution 
#===============================================================================
# First, we have to run code/bash/snpeff.slurm code to annotate snps

vcf <- read.vcfR('data/modified_data/septoria_ale_clean_ann.vcf.gz')
fixed <- data.frame(vcf@fix)
fix_df <- fixed |> 
  mutate(snp = paste0('X', CHROM, '_', POS),
         putative_impact = map_chr(strsplit(INFO, "\\|"), \(x) x[[3]])) |> 
  filter(snp %in% colnames(genotype_septoria)) |> 
  dplyr::select(snp, CHROM, putative_impact)

snps_df <- fix_df |> 
  group_by(CHROM, putative_impact) |> 
  summarise(n = n()) |> 
  mutate(CHROM = as.numeric(CHROM)) |> 
  arrange(CHROM)

colors2 <- c("HIGH" = "#9A348EFF", 
             "LOW" = "#1E88E5FF",  
             "MODERATE" = "#FFC107FF",
             "MODIFIER" = "#2E7D32FF")

snps_plot <- ggplot(snps_df) +
  geom_col(aes(x = factor(CHROM), y = n, fill = putative_impact), position = "stack", alpha = 0.8) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +  # Expresa en miles con "K"
  labs(x = "Chromosome", y = "Count (Thousands)", fill = "Putative Impact") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.major.y = element_line(color = "darkgrey"), 
    axis.ticks.x.bottom = element_blank(),
    axis.text = element_text(size = 14), 
    axis.text.x = element_text(margin = margin(t = -10)),  # Reduce el margen superior
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),  # Adjust top margin for x-axis title
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    legend.position = 'top',
    legend.title.position = "top",  # Move legend title to the top
    legend.title.align = 0.5,
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  scale_fill_manual(values = colors2) +
  guides(fill = guide_legend(ncol = 4)) 
snps_plot

png(paste0("outputs/plots/snps_plot.png"), width = 3000, height = 2000, res = 400)
snps_plot
dev.off()

# add marker putative impact to gwas table
table <- read_csv('outputs/postGWAS_sep/gwas_sep_table_all.csv') |> 
  left_join(fix_df, by = c("SNP" = 'snp')) |> 
  dplyr::rename('Putative impact' = 'putative_impact')

write_csv(table, file = 'outputs/postGWAS_sep/gwas_sep_table_all.csv')

#===============================================================================
# Additiional genomic info
#===============================================================================
info_extra <- read_tsv('~/Downloads/sequence_report.tsv') 

table_impact <- snps_df |> 
  group_by(putative_impact) |> 
  summarise(n_new = sum(n)) |> 
  mutate(percentage = (n_new/sum(n_new))*100)

map_info <- map_septoria |> 
  group_by(Chromosome) |> 
  summarize(Markers = n())
impact_df <- snps_df |> 
  pivot_wider(names_from = putative_impact, values_from = n) |> 
  dplyr::rename(Chromosome = CHROM)
genomic_info  <- info_extra |> 
  dplyr::select(Chromosome = `Chromosome name`, Length = 'Seq length') |> 
  left_join(map_info) |> 
  mutate(Density = Length/Markers) |> 
  left_join(impact_df) |> 
  mutate(
    HIGH = paste0(HIGH, " (", round((HIGH / Markers) * 100, 2), "%)"),
    LOW = paste0(LOW, " (", round((LOW / Markers) * 100, 2), "%)"),
    MODERATE = paste0(MODERATE, " (", round((MODERATE / Markers) * 100, 2), "%)"),
    MODIFIER = paste0(MODIFIER, " (", round((MODIFIER / Markers) * 100, 2), "%)"),
    Density = round(Density, 2)
  ) |> 
  dplyr::slice(1:13) |> 
  dplyr::select(Chromosome, Length, Markers, Density, High = HIGH, Low = LOW, Moderate = MODERATE, Modifier = MODIFIER)
  
write_csv(genomic_info, file = 'outputs/postGWAS_sep/genomic_info_table.csv')

# overall density 
bases <- info_extra |> 
  dplyr::select(length = 'Seq length') |> 
  mutate(length = as.numeric(length)) |>
  dplyr::slice(1:13) |> 
  pull(length) |> 
  sum()
density <- bases/(dim(genotype_septoria)[2]-1)



