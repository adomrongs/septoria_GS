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
shared_snps_all_df %>%
  gt() %>%
  # Poner los nombres de las columnas en negrita
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels()
  ) %>%
  # Añadir un color de fondo alternado a las filas 1 y 2 (rojo claro) y filas 4, 5, y 6 (verde claro)
  tab_style(
    style = list(cell_fill(color = "lightcoral")),
    locations = cells_body(rows = c(1, 2)) # Filas 1 y 2
  ) %>%
  tab_style(
    style = list(cell_fill(color = "lightgreen")),
    locations = cells_body(rows = c(4, 5, 6)) # Filas 4, 5 y 6
  ) %>%
  # Cambiar el color del texto de toda la tabla a negro
  tab_style(
    style = list(cell_text(color = "black")),
    locations = cells_body()
  ) %>%
  # Añadir un título
  tab_header(
    title = md("**HITS ALL_LEAVES**") # Título en negrita
  ) %>%
  # Ajustar espaciado y estilos generales
  tab_options(
  )
table(shared_snps_all_df$SNP)


shared_snps_leaf2 <- map(common_rows$SNP, function(snp) {
  tmp_df <- leaf2 %>% 
    filter(SNP == snp)
  return(tmp_df)
})
shared_snps_leaf2_df <- do.call(rbind, shared_snps_leaf2)
shared_snps_leaf2_df %>%
  gt() %>%
  # Poner los nombres de las columnas en negrita
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels()
  ) %>%
  # Añadir un color de fondo alternado a las filas 1 y 2 (rojo claro) y filas 4, 5, y 6 (verde claro)
  tab_style(
    style = list(cell_fill(color = "lightcoral")),
    locations = cells_body(rows = c(1, 2, 5, 6)) # Filas 1 y 2
  ) %>%
  # Cambiar el color del texto de toda la tabla a negro
  tab_style(
    style = list(cell_text(color = "black")),
    locations = cells_body()
  ) %>%
  # Añadir un título
  tab_header(
    title = md("**HITS LEAF2**") # Título en negrita
  )

