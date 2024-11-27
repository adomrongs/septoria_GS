library(tidyverse)
library(factoextra)
library(ggrepel)
library(RColorBrewer)
source("code/R/function_septoria_GS.R")

#-------------------------------------------------------------------------------
# PCA and GRM Zymoseptoria
#-------------------------------------------------------------------------------

# load and read data
load("data/modified_data/1_septoria.Rdata")
info_strains <- read_csv("data/raw_data/INFO_STRAINS.csv")
regions <- readxl::read_xlsx("data/raw_data/regions_strain.xlsx")
regions <- regions %>% dplyr::select(Region, Strain)

#change the name of the strains to match with the ones from the marker matrix
new_names <- sapply(regions$Strain, function(nombre) {
  # Verifica si el nombre comienza con un número usando una expresión regular
  if (grepl("^[0-9]+_", nombre)) {
    return(nombre)  # Si comienza con un número, lo deja igual
  } else {
    partes <- strsplit(nombre, "_")[[1]]  # Divide el nombre por "_"
    nuevo_nombre <- paste(partes[2], partes[1], partes[3], sep = "_")  # Reordena y une
    return(nuevo_nombre)
  }
})
new_names <- gsub("\\.", "_", new_names)
regions$Strain <- unname(new_names)

# Assign its region to each isolate
pca_df <- data.frame(Isolate = genotype_septoria[,"Isolate"]) %>% 
  left_join(regions, by = c("Isolate" = "Strain"))
# Fix manually the ones that did not match the information from the file
missing_regions <- pca_df[is.na(pca_df$Region), ]
missing_regions <- missing_regions %>% 
  mutate(Region = case_when(
    grepl("^\\d+_Cor", Isolate) ~ "Córdoba",
    grepl("^\\d+_Jer", Isolate) ~ "Cádiz",
    grepl("^\\d+_Con", Isolate) ~ "Cádiz",
    grepl("^\\d+_Eci", Isolate) ~ "Sevilla",
    grepl("^\\d+_Esc", Isolate) ~ "Huelva",
    TRUE ~ NA_character_  # Para cualquier otro caso, deja <NA>
  ))
pca_df[is.na(pca_df$Region), "Region"] <-  missing_regions$Region
sum(is.na(pca_df$Region))

# assign colors
colors <- c("Córdoba" = "#DD5129FF",
            "Cádiz" = "#0F7BA2FF",
            "Sevilla" = "#43B284FF",
            "Huelva" = "#FAB255FF")
# assign shapes
shapes <- c("Córdoba" = 16,
            "Cádiz" = 16,
            "Sevilla" = 16,
            "Huelva" = 16)

# adjust pca_df to match the isolates in genotype_septoria
pca_df <- pca_df %>% 
  mutate(color = colors[Region],
                 shape = shapes[Region]) %>% 
  filter(Isolate %in% genotype_septoria[, "Isolate"]) %>% 
  distinct()

# create pca and screeplot 
PCA_plot_sep <- plotPCA(genotype = genotype_septoria,
                          regions = pca_df$Region,
                          colors = colors,
                          shapes = shapes)

# save results
png(paste0("outputs/plots/PCA_sep.png"), width = 3000, height = 3000, res = 400)
PCA_plot_sep$pca_plot
dev.off()

png(paste0("outputs/plots/screeplot_sep.png"), width = 3000, height = 3000, res = 400)
PCA_plot_sep$scree_plot
dev.off()

# create heatmap for the GRM
k_sep_mat <- as.matrix(k_septoria[,-1])
png(paste0("outputs/plots/GRM_sep.png"), width = 5000, height = 5000, res = 400)
heatmap(k_sep_mat, col = colorRampPalette(brewer.pal(8, "Oranges"))(25), 
        RowSideColors = pca_df$color,
        labRow = k_septoria[,1],
        labCol = k_septoria[,1])
dev.off()
#-------------------------------------------------------------------------------
# PCA and GRM WHEAT
#-------------------------------------------------------------------------------

#load information
wheat_information <- readxl::read_xlsx("data/raw_data/wheat_families.xlsx")
wheat_information <- wheat_information %>% 
  dplyr::select(Sample_Name, INIA.name, groups, genotypes) %>% 
  mutate(Sample_Name = gsub("\\.", "_", Sample_Name))

# assign checks
wheat_information[which(grepl("Ckeck", wheat_information$groups)), "genotypes"] <- "CHECK"
colnames(wheat_information) <- c("GenoID", "Name", "nothing", "region")

# assign colors
colors <- c("CHECK" = "#DD5129FF",
            "Landraces" = "#0F7BA2FF",
            "Lines" = "#43B284FF",
            "Cultivars" = "#FAB255FF")
#asign shapes
shapes <- c("CHECK" = 17,
            "Landraces" = 16,
            "Lines" = 16,
            "Cultivars" = 16)

# create data.frame for the information requirted for the plotting
load(file = "data/modified_data/2_wheat.Rdata")
regions_df <- genotype_wheat %>% 
  left_join(wheat_information) %>% 
  dplyr::select(GenoID, Name, region) %>% 
  mutate(color = colors[region],
         shape = shapes[region])

# create pca and screeplot
PCA_plot_wheat <- plotPCA(genotype = genotype_wheat,
               regions = regions_df$region,
               colors = colors,
               names = regions_df$Name,
               shapes = shapes)
# save rsults
png(paste0("outputs/plots/PCA_wheat.png"), width = 3000, height = 3000, res = 400)
PCA_plot_wheat$pca_plot
dev.off()

png(paste0("outputs/plots/screeplot_wheat.png"), width = 3000, height = 3000, res = 400)
PCA_plot_wheat$scree_plot
dev.off()

# create heatmap for the GRM
k_wheat_mat <- as.matrix(k_wheat[,-1])
png(paste0("outputs/plots/GRM_wheat.png"), width = 5000, height = 5000, res = 400)
heatmap(k_wheat_mat, col = colorRampPalette(brewer.pal(8, "Oranges"))(25), 
        RowSideColors = regions_df$color,
        labRow = k_wheat[,1],
        labCol = k_wheat[,1])
dev.off()

