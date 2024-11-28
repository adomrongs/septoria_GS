plotPCA <- function(genotype, regions = NULL, names = NULL, colors = NULL, shapes = NULL){
  if (!is.matrix(genotype) && !is.data.frame(genotype)) {
    stop("Input must be a matrix or data frame.")
  }
  
  # Check if regions is provided and has the correct length
  if (!is.null(regions)) {
    if (length(regions) != nrow(genotype)) {
      stop("The length of 'regions' must match the number of rows in 'genotype'.")
    }
  }
  
  # double center the genotypic matrix for the PCA
  row_means <- rowMeans(genotype[, -1])
  col_means <- colMeans(genotype[, -1])
  overall_mean <- mean(as.matrix(genotype[, -1]))
  genotype[, -1] <- genotype[, -1] - row_means + overall_mean
  genotype[, -1] <- t(t(genotype[, -1]) - col_means)
  
  # compute the PC analysis
  PCA <- prcomp(genotype[, -1], center = F)
  PCs <- PCA$x
  scree_plot <- fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 100)) # extract scree plot
  
  # Prepare the data for PCA plot
  pca_data <- as.data.frame(PCs) %>% 
    mutate(GenoID = genotype[,1])
  if (!is.null(regions)) {
    pca_data$regions <- as.factor(regions)  # Add regions as a factor
    pca_data$shape <- as.factor(ifelse(pca_data$regions == "CHECK", 17, 16))
  } else {
    pca_data$regions <- as.factor("All")  # Default value if no regions are provided
  }
  
  if(!is.null(names)){
    pca_data$Names <- names
  } else {
    pca_data$Names <- NA
  }
  
    pca_plot <- ggplot(data = pca_data) +
      geom_point(aes(x = PC1, y = PC2, color = regions, shape = regions), size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_text_repel(
        data = pca_data %>% filter(regions == "CHECK"),  # Filter data where regions == "CHECK"
        aes(x = PC1, y = PC2, label = Names),  # Ensure x and y are mapped to PC1 and PC2
        size = 4,  # Text size
        box.padding = 1.5, max.overlaps = Inf # Push text farther vertically
      ) +
      labs(
        x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""),
        y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = "")
      ) +  # Axis labels
      theme_minimal() +  # Clean theme
      theme(
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",  # Move legend to the top
        legend.title = element_blank(),  # Remove legend title
        legend.key = element_blank(),  # Optional: remove key background
        legend.text = element_text(size = 18),  # Increase text size of legend
        plot.title = element_blank(),  # Remove title from the plot
        axis.title.x = element_text(size = 18),  # Axis label sizes
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 17),  # Axis text sizes
        strip.text = element_text(size = 10, face = "plain", color = "black", hjust = 0.5),
        strip.background = element_rect(fill = "lightgray"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        plot.margin = margin(t = 10, r = 40, b = 10, l = 10) # Panel border
      ) +
      scale_color_manual(values = colors) +
      scale_shape_manual(values = shapes)
  
  return(list(scree_plot = scree_plot, pca_plot = pca_plot))
}


plotGrid <- function(phenotype_long, trait, colors) {
  plot <- ggplot(data = phenotype_long) +
    stat_halfeye(aes_string(y = "Value", fill = trait), 
                 adjust = 20,
                 width = 5,
                 justification = -0.7,
                 .width = 0,
                 point_colour = NA) +
    # Dot plot
    stat_dots(aes_string(y = "Value"),
              side = "left",
              justification = 10,
              binwidth = NA,
              dotsize = 3,
              overflow = "compress",
              scale = 0.4)  +
    # Boxplot
    geom_boxplot(aes_string(y = "Value", fill = trait), 
                 width = 4,
                 outlier.shape = NA,
                 color = "black") + 
    scale_fill_manual(values = colors) +
    theme(
      panel.background = element_blank(),  # Background color
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 18),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 17),
      strip.text = element_text(size = 10, face = "plain", color = "black", hjust = 0.5),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 10),
      legend.position = "none",  # Hide legend
      axis.title.x = element_text(size = 15),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.ticks = element_blank(),
      text = element_text(family = "Arial") # Optional: Set a custom font family
    ) + 
    # Faceting: 3 rows for Trait, 7 columns for trait
    facet_grid(Trait ~ .data[[trait]], scales = "free_y") +  # Updated for dynamic faceting
    xlab(label = paste0(trait))  # Labeling x-axis with trait name
  return(plot)  # Return the plot
}

plotHist <- function(df, columns, colors) {
  # Select specified columns from the dataframe
  df_selected <- df %>% dplyr::select(all_of(columns))
  
  # Convert to long format
  df_long <- df_selected %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

  plots <- list()  # Initialize an empty list to store plots
  for(i in seq_along(columns)){
    col <- columns[i]  # Obtener el nombre de la columna
    color <- colors[[i]]
    
    # Crear un histograma para la variable actual
    range_col <- range(df[[col]], na.rm = TRUE)  # Obtener el rango de la columna (min y max)
    binwidth <- diff(range_col) / 10  # Definir el ancho del bin
    
    p <- ggplot(data = df_long %>% filter(variable == col), aes(x = value)) +
      geom_histogram(binwidth = binwidth, fill = color, color = "black", alpha = 0.7) +  # Usar el color directamente
      labs(title = paste(col), x = "Value", y = "Count") +  # Agregar etiquetas
      theme_ipsum() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(colour = "black", size = 3),
        legend.position = "top",
        plot.title = element_blank(),  # Elimina el título del gráfico
        strip.text = element_text(size = 10, face = "plain", color = "black", hjust = 0.5),
        strip.background = element_rect(fill = "lightgray"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text.x = element_text(size = 11, margin = margin(t = 10)),  # Tamaño del texto del eje x aumentado a 18
        axis.text.y = element_text(size = 11, margin = margin(r = 10)),  # Tamaño del texto del eje y aumentado a 18
        axis.title.x = element_text(size = 15),  # Tamaño del texto del título del eje x a 15
        axis.title.y = element_text(size = 15)   # Tamaño del texto del título del eje y a 15
      ) + 
      facet_grid(. ~ variable)  # Crear facetas para cada variable
    
    plots[[i]] <- p  # Almacenar el gráfico en la lista
  }
  
  return(plots)  # Return the list of plots
}

remove_outliers <- function(df, cols) {
  df_clean <- df
  for (col in cols) {
    # Get the outliers for the column
    outliers <- boxplot.stats(df[[col]])$out
    # Replace outliers with NA
    df_clean[[col]][df[[col]] %in% outliers] <- NA
    print(paste0(length(outliers), " outliers have been removed for ", col))
  }
  df_clean <- na.omit(df_clean)
  return(df_clean)
}

extract_blups <- function(phenotype, trait, formula) {
  model <- lmer(as.formula(paste(trait, formula)), data = phenotype)
  blups <- data.frame(ranef(model)$Plant) %>% 
    rownames_to_column("GenoID") %>% 
    dplyr::select(GenoID, everything())
  
  colnames(blups) <- c("GenoID", trait)
  rownames(blups) <- NULL
  return(blups)
}

extract_blups_df <- function(phenotype, traits, formula) {
  blups_list <- list()
  
  # Populate blups_list and capture dimensions
  dims_list <- vector("list", length(traits))
  for (i in seq_along(traits)) {
    trait <- traits[i]
    blups <- extract_blups(phenotype, trait, formula)
    blups_list[[i]] <- blups
    dims_list[[i]] <- dim(blups)
  }
  
  # Check if dimensions match
  if (length(unique(sapply(dims_list, paste, collapse = "x"))) == 1) {
    message("Dimensions match. Returning a combined blups data frame.")
    blups_df <- purrr::reduce(blups_list, ~ dplyr::left_join(.x, .y, by = "GenoID"))
    return(blups_df)
  } else {
    message("Dimensions do not match. Returning blups as a list.")
    return(blups_list)
  }
}

grepGAPITres <- function(parent_directory){
  subdirs <- list.dirs(parent_directory, recursive = F)
  results <- list()
  
  for(i in seq_along(subdirs)){
    message("Working on file ", dir)
    dir <- subdirs[[i]]
    file <- paste0(dir, "/GAPIT.Association.Filter_GWAS_results.csv")
    results_df <- read_csv(file) %>% 
      mutate(File = dir)
    results[[i]] <- results_df
  }
  
  return(results)
}