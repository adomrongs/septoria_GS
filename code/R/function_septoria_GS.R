allele2numeric <- function(matrix, maf.filter = NULL, heter.filter = NULL) {
  # Internal function to convert genotypes to numeric values
  convert <- function(col) {
    # Split the alleles
    alleles <- strsplit(col, "/")
    alleles <- do.call(rbind, alleles) # Convert the list into a two-column matrix
    
    # Handle NA: if any entry is NA, keep it as NA
    alleles[is.na(alleles)] <- NA
    
    # Frequency table of the alleles
    allele_counts <- sort(table(alleles), decreasing = TRUE)
    
    # Identify the second most frequent allele (minor allele)
    if (length(allele_counts) < 2) return(rep(NA, length(col))) # If there's only one allele, return NA
    
    minor_allele <- names(allele_counts)[2]
    
    # Convert alleles into 0/1 while preserving NA for missing values
    numeric_values <- ifelse(alleles == minor_allele, 1, 
                             ifelse(alleles != minor_allele, 0, NA))
    
    # Sum by row while preserving NA values
    numeric_values_sum <- rowSums(numeric_values, na.rm = F)
    
    return(numeric_values_sum)
  }
  
  # Internal function to calculate the Minor Allele Frequency (MAF)
  extract.maf <- function(col) {
    # Exclude NA values for frequency calculation
    col_no_na <- col[!is.na(col)]
    
    # Frequency table of numeric values
    allele_counts <- sort(table(col_no_na), decreasing = TRUE)
    n <- length(col_no_na) # Total number of non-NA observations
    
    if (length(allele_counts) < 2) return(0) # If no second allele exists, MAF is zero
    
    maf_value <- allele_counts[2] / n
    return(maf_value)
  }
  
  # Apply the conversion to each column
  M <- apply(matrix, 2, convert)
  
  # Calculate MAF for each column
  maf <- apply(M, 2, extract.maf)
  
  # Filter columns based on MAF if maf.filter is defined
  if (!is.null(maf.filter)) {
    valid_markers <- which(maf >= maf.filter) # Indices of valid columns
    M <- M[, valid_markers, drop = FALSE]
    cat("Markers removed with Minor Allele Frequency below", maf.filter, ":",
        ncol(matrix) - length(valid_markers), "\n")
  }
  
  if(!is.null(heter.filter)){
    het_snps <- apply(M, 2, function(column) sum(column == 1, na.rm = T))
    het_snps_per_sample <- het_snps / nrow(M)
    het_snps <- het_snps_per_sample[het_snps_per_sample > heter.filter]
    
    M <- M[, !colnames(M) %in% names(het_snps)]
    cat("Markers removed with Heterozygosity over", heter.filter, ":",
        length(het_snps), "\n")
  }
  
  return(M)
}

plotPCA <- function(genotype, regions = NULL, names = NULL, colors = NULL, shapes = NULL, interactive = NULL){
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
      geom_point_interactive(aes(x = PC1, y = PC2, color = regions, shape = regions, tooltip = Names, data_id = Names),
                             size = 3) +
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
    
    if(!is.null(interactive)){
      interactive_plot <- girafe(ggobj = pca_plot, 
                                 options = list(opts_hover(css = "stroke:#000; stroke-width: 0.5px; transition: all 0.3s ease;"),
                                                opts_hover_inv("opacity:0.5;filter:saturate(10%);")))
      return(list(scree_plot = scree_plot, pca_plot = pca_plot, interactive = interactive_plot))
    } else {
      return(list(scree_plot = scree_plot, pca_plot = pca_plot))
    }
}

plotPCA2 <- function(genotype, regions = NULL, names = NULL, colors = NULL, shape_col = NULL,  shapes = NULL, interactive = NULL){
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
  
  if(!is.null(shape_col)){
    pca_data$other <- shape_col
  } else {
    pca_data$other <- NA
  }
  
  pca_plot <- ggplot(data = pca_data) +
    geom_point_interactive(aes(x = PC1, y = PC2, color = regions, shape = other, tooltip = Names, data_id = Names),
                           size = 3) +
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
  
  if(!is.null(interactive)){
    interactive_plot <- girafe(ggobj = pca_plot, 
                               options = list(opts_hover(css = "stroke:#000; stroke-width: 0.5px; transition: all 0.3s ease;"),
                                              opts_hover_inv("opacity:0.5;filter:saturate(10%);")))
    return(list(scree_plot = scree_plot, pca_plot = pca_plot, interactive = interactive_plot))
  } else {
    return(list(scree_plot = scree_plot, pca_plot = pca_plot))
  }
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

extract_blues <- function(phenotype, trait, formula) {
  model <- lm(as.formula(paste(trait, formula)), data = phenotype)
  blues <- coef(model)[grep("Plant", names(coef(model)), value = T)]
  names(blues) <- gsub("Plant", "", names(blues))
  blues_df <- data.frame(blues) %>% 
    rownames_to_column("GenoID") %>% 
    dplyr::select(GenoID, everything())
  
  colnames(blues_df) <- c("GenoID", trait)
  rownames(blues_df) <- NULL
  return(blues_df)
}

extract_blues_df <- function(phenotype, traits, formula) {
  blues_list <- list()
  
  # Populate blues_list and capture dimensions
  dims_list <- vector("list", length(traits))
  for (i in seq_along(traits)) {
    trait <- traits[i]
    blups <- extract_blues(phenotype, trait, formula)
    blues_list[[i]] <- blups
    dims_list[[i]] <- dim(blups)
  }
  
  # Check if dimensions match
  if (length(unique(sapply(dims_list, paste, collapse = "x"))) == 1) {
    message("Dimensions match. Returning a combined blues data frame.")
    blues_df <- purrr::reduce(blues_list, ~ dplyr::left_join(.x, .y, by = "GenoID"))
    return(blues_df)
  } else {
    message("Dimensions do not match. Returning blues as a list.")
    return(blues_list)
  }
}

extract_blues_adapted <- function(phenotype, trait, formula, colname) {
  model <- lm(as.formula(paste(trait, formula)), data = phenotype)
  blues <- coef(model)[grep(colname, names(coef(model)), value = T)]
  names(blues) <- gsub(colname, "", names(blues))
  blues_df <- data.frame(blues) %>% 
    rownames_to_column(colname) %>% 
    dplyr::select(colname, everything())
  
  colnames(blues_df) <- c(colname, trait)
  rownames(blues_df) <- NULL
  return(blues_df)
}

extract_blues_df_adapted <- function(phenotype, traits, formula, colname) {
  blues_list <- list()
  
  # Populate blues_list and capture dimensions
  dims_list <- vector("list", length(traits))
  for (i in seq_along(traits)) {
    trait <- traits[i]
    blups <- extract_blues_adapted(phenotype, trait, formula, colname)
    blues_list[[i]] <- blups
    dims_list[[i]] <- dim(blups)
  }
  
  # Check if dimensions match
  if (length(unique(sapply(dims_list, paste, collapse = "x"))) == 1) {
    message("Dimensions match. Returning a combined blues data frame.")
    blues_df <- purrr::reduce(blues_list, ~ dplyr::left_join(.x, .y, by = colname))
    return(blues_df)
  } else {
    message("Dimensions do not match. Returning blues as a list.")
    return(blues_list)
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

runS1 <- function(trait, Kw, Kmix, pheno, genoW, map, wtest, formula, wModel = FALSE){
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  
  ptrain <- pheno
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  if (!wModel) {
    Xwtrain <- model.matrix(~ 1, data = ptrain)
    results <- data.frame()
  } 
  if (wModel) {
    bonferroni <- 0.05/nrow(map)
    # prepare train data
    gwas_geno <- data.frame(genoW[genoW[,1] %in% wtrain, ]) # subset genotype
    # gwas_k <- Kw %>% 
    #   filter(rownames(.) %in% gwas_geno[,1]) %>% 
    #   dplyr::select(all_of(gwas_geno[,1])) %>% 
    #   rownames_to_column("GenoID") %>% 
    #   dplyr::select(GenoID, everything())
    gwas_pheno <- extract_blues_df_adapted(ptrain, trait, formula, "Plant")
    colnames(gwas_pheno) <- c("GenoID", trait)
    
    dim(gwas_geno); dim(gwas_pheno); dim(map)#; dim(gwas_k)
    # run GWAS
    tmp <- capture.output({
      scores <- GAPIT(Y = gwas_pheno,
                      GD = gwas_geno,
                      GM = map,
                      KI = NULL,
                      CV = NULL,
                      PCA.total = 2,
                      model = "Blink",
                      file.output = F)
    })
    rm(tmp)
    gc()
    setwd(here())
    
    results <- scores[["GWAS"]] %>% 
      arrange(P.value)
    
    hits_bonferroni <- results %>% filter(P.value <= bonferroni)
    
    if (nrow(hits_bonferroni) == 0) {
      sSNPs <- results$SNP[1]
    }
    if (nrow(hits_bonferroni) >= 3) {
      sSNPs <- results$SNP[1:3]
    }
    if (nrow(hits_bonferroni) == 2) {
      sSNPs <- results$SNP[2]
    }
    if (nrow(hits_bonferroni) == 1) {
      sSNPs <- results$SNP[2]
    }
    sSNPs_data <- genoW[, sSNPs, drop = FALSE] %>%
      data.frame() %>% 
      mutate(ID = genoW[,1])
    
    Xwtrain <- model.matrix(~ 1, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  K12_wheat <- Zwtrain %*% as.matrix(Kw) %*% t(Zwtrain)
  K12_mix <- Zmixtrain %*% Kmix %*% t(Zmixtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model1 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2 <- computeH2(model1)
  
  model1_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2_I <- computeH2(model1_I, interaction = T)
  
  return(list(
    S1 = model1, # results model without interaction term
    S1_I = model1_I, # results with interaction term
    wtest = wtest, # partition test 
    wtrain = wtrain,
    hits = results,
    H_list = list(H2, H2_I)
  ))
}

eval_S1 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$S1)
  predictions_I <- predict(strategy$S1_I)
  
  total_results <- data.frame(Plant = phenotype$Plant, 
                              Strain = phenotype$Strain,
                              Observed = phenotype[[trait]], 
                              Predicted = predictions,
                              Predicted_I = predictions_I)
  test_results <- total_results %>% 
    filter(Plant %in% strategy$wtest)
  # overall correlation
  lines_test <- phenotype$Plant %in% strategy$wtest
  pheno_test <- phenotype[lines_test,]
  cor <- cor(predictions[lines_test],
             pheno_test[[trait]],
             use = "complete.obs")
  cor_I <- cor(predictions_I[lines_test],
               pheno_test[[trait]],
               use = "complete.obs")
  
  # within strain correlation for both model and model_I
  cor_withinStrain <- test_results %>%
    split(.$Strain) %>%  # Divide el dataframe por cada valor único de 'Strain'
    map(~ cor(.x$Observed, .x$Predicted, use = "complete.obs"))
  cor_withinStrain_I <- test_results %>%
    split(.$Strain) %>%  # Divide el dataframe por cada valor único de 'Strain'
    map(~ cor(.x$Observed, .x$Predicted_I, use = "complete.obs"))
  # prediction accuracy calculation as the prediction ability divided by srt of h2
  accuracy_withinStrain <- map(cor_withinStrain, ~ .x/sqrt(strategy$H_list[[1]]))
  accuracy_withinStrain_I <- map(cor_withinStrain, ~ .x/sqrt(strategy$H_list[[2]]))
  # save resulst
  cor_results <- list(
    cor = cor,
    cor_I = cor_I,
    cor_withinStrain = cor_withinStrain,
    cor_withinStrain_I = cor_withinStrain_I,
    hits <- strategy$hits,
    accuracy = accuracy_withinStrain,
    accuracy_I = accuracy_withinStrain_I
  )
  
  return(cor_results)
}

runS2 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, formula, wModel = FALSE) {
  
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~ 1, data = ptrain)
    results <- data.frame()
  } 
  if (wModel) {
    bonferroni <- 0.05/nrow(map)
    gwas_geno <- genoW
    gwas_pheno <- extract_blues_df_adapted(ptrain, trait, formula, "Plant")
    colnames(gwas_pheno) <- c("GenoID", trait)
    # gwas_k <- data.frame(Kw) %>% 
    #   rownames_to_column("GenoID") %>% 
    #   dplyr::select(GenoID, everything())
    
    dim(gwas_pheno); dim(gwas_geno); dim(map)#; dim(gwas_k)
    tmp <- capture.output({
      scores <- GAPIT(Y = gwas_pheno,
                      GD = gwas_geno,
                      GM = map,
                      KI = NULL,
                      CV = NULL,
                      PCA.total = 2,
                      model = "Blink",
                      file.output = F)
    })
    rm(tmp)
    gc()
    setwd(here())
    
    results <- scores[["GWAS"]] %>% 
      arrange(P.value)
    hits_bonferroni <- results %>% filter(P.value <= bonferroni)
    
    if (nrow(hits_bonferroni) == 0) {
      sSNPs <- results$SNP[1]
    }
    if (nrow(hits_bonferroni) >= 3) {
      sSNPs <- results$SNP[1:3]
    }
    if (nrow(hits_bonferroni) == 2) {
      sSNPs <- results$SNP[2]
    }
    if (nrow(hits_bonferroni) == 1) {
      sSNPs <- results$SNP[2]
    }
    sSNPs_data <- genoW[, sSNPs, drop = FALSE] %>%
      data.frame() %>% 
      mutate(ID = genoW[,1])
    
    Xwtrain <- model.matrix(~ 1, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  
  K12_mix <- Zmixtrain %*% as.matrix(Kmix) %*% t(Zmixtrain)
  K12_wheat <- Zwtrain %*% as.matrix(Kw) %*% t(Zwtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model2 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2 <- computeH2(model2)
  
  model2_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2_I <- computeH2(model2_I, interaction = T)
  
  return(list(model2 = model2,
              model2_I = model2_I,
              sMix = sMix,
              H_list = list(H2, H2_I),
              hits = results))
}

eval_S2 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model2)
  predictions_I <- predict(strategy$model2_I)
  
  total_results <- data.frame(Plant = phenotype$Plant, 
                              Strain = phenotype$Strain, 
                              Observed = phenotype[[trait]], 
                              Predicted = predictions, 
                              Predicted_I = predictions_I)
  test_results <- total_results %>% 
    filter(Strain %in% strategy$sMix)
  
  cor <- cor(test_results$Observed, test_results$Predicted)
  cor_I <- cor(test_results$Observed, test_results$Predicted_I)
  
  accuracy <- cor/sqrt(strategy$H_list[[1]])
  accuracy_I <- cor_I/sqrt(strategy$H_list[[2]])
  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I,
    accuracy = accuracy,
    accuracy_I = accuracy_I,
    hits = strategy$hits
  )
  
  return(list(CorrelationResults = correlationResults))
}

run_S3 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, formula, wtest, wModel = FALSE) {
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~ 1, data = ptrain)
    results <- data.frame()
  } 
  if (wModel) {
    bonferroni <- 0.05/nrow(map)
    gwas_pheno <- extract_blues_df_adapted(ptrain, trait, formula, "Plant")
    colnames(gwas_pheno) <- c("GenoID", trait)
    gwas_geno <- data.frame(genoW[genoW[,1] %in% gwas_pheno$GenoID, ])
    # gwas_k <- Kw %>% 
    #   filter(rownames(.) %in% gwas_geno[,1]) %>% 
    #   dplyr::select(all_of(gwas_geno[,1])) %>% 
    #   rownames_to_column("GenoID") %>% 
    #   dplyr::select(GenoID, everything())
    
    dim(gwas_geno); dim(gwas_pheno); dim(map)#; dim(gwas_k)
    tmp <- capture.output({
      scores <- GAPIT(Y = gwas_pheno,
                      GD = gwas_geno,
                      GM = map,
                      KI = NULL,
                      CV = NULL,
                      PCA.total = 2,
                      model = "Blink",
                      file.output = F)
    })
    rm(tmp)
    gc()
    setwd(here())
    
    results <- scores[["GWAS"]] %>% 
      arrange(P.value)
    hits_bonferroni <- results %>% filter(P.value <= bonferroni)
    
    if (nrow(hits_bonferroni) == 0) {
      sSNPs <- results$SNP[1]
    }
    if (nrow(hits_bonferroni) >= 3) {
      sSNPs <- results$SNP[1:3]
    }
    if (nrow(hits_bonferroni) == 2) {
      sSNPs <- results$SNP[2]
    }
    if (nrow(hits_bonferroni) == 1) {
      sSNPs <- results$SNP[2]
    }
    sSNPs_data <- genoW[, sSNPs, drop = FALSE] %>%
      data.frame() %>% 
      mutate(ID = genoW[,1])
    
    Xwtrain <- model.matrix(~ 1, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  
  K12_mix <- Zmixtrain %*% as.matrix(Kmix) %*% t(Zmixtrain)
  K12_wheat <- Zwtrain %*% as.matrix(Kw) %*% t(Zwtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model3 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2 <- computeH2(model3)
  
  model3_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2_I <- computeH2(model3_I, interaction = T)
  
  return(list(model3 = model3,
              model3_I = model3_I,
              sMix = sMix,
              wtest = wtest,
              H_list = list(H2, H2_I),
              hits = results))
}

eval_S3 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model3)
  predictions_I <- predict(strategy$model3_I)
}

scenario1 <- function(allResults) {
  # Initialize lists to store data
  cor_strain_list <- cor_strain_I_list <- acc_strain_list <- acc_strain_I_list <- list()
  
  for (result_name in names(allResults)) {
    s1_models <- grep("Scenario1", names(allResults[[result_name]]), value = TRUE)
    
    for (model in s1_models) {
      # Process cor_strain
      tmp_df <- data.frame(allResults[[result_name]][[model]][[3]]) %>% 
        mutate(Matrix = result_name,
               Model = ifelse(grepl("1w", model), "Weighted", "Normal"),
               Interaction = "")
      cor_strain_list[[paste(result_name, model, "noI", sep = "_")]] <- tmp_df
      
      tmp_acc <- data.frame(allResults[[result_name]][[model]][[6]]) %>% 
        mutate(Matrix = result_name,
               Model = ifelse(grepl("1w", model), "Weighted", "Normal"),
               Interaction = "")
      acc_strain_list[[paste(result_name, model, "noI", sep = "_")]] <- tmp_acc
      
      # Process cor_strain_I
      tmp_df_I <- data.frame(allResults[[result_name]][[model]][[4]]) %>% 
        mutate(Matrix = result_name,
               Model = ifelse(grepl("1w", model), "Weighted", "Normal"),
               Interaction = "+I")
      cor_strain_I_list[[paste(result_name, model, "I", sep = "_")]] <- tmp_df_I
      
      tmp_acc_I <- data.frame(allResults[[result_name]][[model]][[7]]) %>% 
        mutate(Matrix = result_name,
               Model = ifelse(grepl("1w", model), "Weighted", "Normal"),
               Interaction = "+I")
      acc_strain_I_list[[paste(result_name, model, "I", sep = "_")]] <- tmp_acc_I
    }
  }
  
  # Combine results from lists
  cor_strain <- bind_rows(cor_strain_list)
  cor_strain_I <- bind_rows(cor_strain_I_list)
  acc_strain <- bind_rows(acc_strain_list)
  acc_strain_I <- bind_rows(acc_strain_I_list)
  
  # Combine both datasets and add Info column
  s1_df <- bind_rows(cor_strain, cor_strain_I) %>% 
    mutate(Info = paste0(Model, Interaction)) %>% 
    dplyr::select(-Model, -Interaction)
  
  acc_df <- bind_rows(acc_strain, acc_strain_I) %>% 
    mutate(Info = paste0(Model, Interaction)) %>% 
    dplyr::select(-Model, -Interaction)
  
  # Pivot longer for final format
  s1_df <- s1_df %>% 
    pivot_longer(cols = -c(Info, Matrix), names_to = "Mix", values_to = "Cor") %>% 
    mutate(across(c(Info, Mix, Matrix), as.factor)) %>% 
    dplyr::select(Info, Matrix, Mix, Cor)
  
  acc_df <- acc_df %>% 
    pivot_longer(cols = -c(Info, Matrix), names_to = "Mix", values_to = "Cor") %>% 
    mutate(across(c(Info, Mix, Matrix), as.factor)) %>% 
    dplyr::select(Info, Matrix, Mix, Cor)
  
  
  return(list(ability = s1_df,  accuracy = acc_df))
}

scenario2 <- function(allResults){
  cor_strain_list <- cor_strain_I_list <- acc_strain_list <- acc_strain_I_list <- list()
  
  for (result_name in names(allResults)) {
    s2_models <- grep("Scenario2", names(allResults[[result_name]]), value = TRUE)
    for (model in s2_models) {
      tmp_df <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[1]]) %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("2w", model), "Weighted", "Normal"),
               Interaction = "")
      colnames(tmp_df)[1] <- "Cor"
      cor_strain_list[[paste(result_name, model, "noI", sep = "_")]] <- tmp_df
      
      tmp_acc <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[3]]) %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("2w", model), "Weighted", "Normal"),
               Interaction = "")
      colnames(tmp_acc)[1] <- "Cor"
      acc_strain_list[[paste(result_name, model, "noI", sep = "_")]] <- tmp_acc
      
      tmp_df_I <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[2]]) %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("2w", model), "Weighted", "Normal"),
               Interaction = "+I")
      colnames(tmp_df_I)[1] <- "Cor"
      cor_strain_I_list[[paste(result_name, model, "I", sep = "_")]] <- tmp_df_I
      
      tmp_acc_I <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[4]])  %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("2w", model), "Weighted", "Normal"),
               Interaction = "+I")
      colnames(tmp_acc_I)[1] <- "Cor"
      acc_strain_I_list[[paste(result_name, model, "I", sep = "_")]] <- tmp_acc_I
    }
  }
  
  cor_strain <- bind_rows(cor_strain_list)
  cor_strain_I <- bind_rows(cor_strain_I_list)
  acc_strain <- bind_rows(acc_strain_list)
  acc_strain_I <- bind_rows(acc_strain_I_list)
  
  s2_df <- bind_rows(cor_strain, cor_strain_I) %>% 
    mutate(Info = paste0(Model, Interaction)) %>% 
    dplyr::select(Info, Matrix, Mix, Cor) %>% 
    mutate(across(c(Info, Mix, Matrix), as.factor))
  
  acc_df <- bind_rows(acc_strain, acc_strain_I) %>% 
    mutate(Info = paste0(Model, Interaction)) %>% 
    dplyr::select(Info, Matrix, Mix, Cor) %>% 
    mutate(across(c(Info, Mix, Matrix), as.factor))
  
  return(list(ability = s2_df,  accuracy = acc_df))
}

scenario3 <- function(allResults){
  cor_strain_list <- cor_strain_I_list <- acc_strain_list <- acc_strain_I_list <- list()
  
  for (result_name in names(allResults)) {
    s3_models <- grep("Scenario3", names(allResults[[result_name]]), value = TRUE)
    for (model in s3_models) {
      tmp_df <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[1]]) %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("3w", model), "Weighted", "Normal"),
               Interaction = "")
      colnames(tmp_df)[1] <- "Cor"
      cor_strain_list[[paste(result_name, model, "noI", sep = "_")]] <- tmp_df
      
      tmp_acc <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[3]])%>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("3w", model), "Weighted", "Normal"),
               Interaction = "")
      colnames(tmp_acc)[1] <- "Cor"
      acc_strain_list[[paste(result_name, model, "noI", sep = "_")]] <- tmp_acc
      
      tmp_df_I <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[2]]) %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("3w", model), "Weighted", "Normal"),
               Interaction = "+I")
      colnames(tmp_df_I)[1] <- "Cor"
      cor_strain_I_list[[paste(result_name, model, "I", sep = "_")]] <- tmp_df_I
      
      tmp_acc_I <- data.frame(allResults[[result_name]][[model]][["CorrelationResults"]][[4]])  %>% 
        mutate(Mix = sub(".* ", "", model), 
               Matrix = result_name,
               Model = ifelse(grepl("3w", model), "Weighted", "Normal"),
               Interaction = "+I")
      colnames(tmp_acc_I)[1] <- "Cor"
      acc_strain_I_list[[paste(result_name, model, "I", sep = "_")]] <- tmp_acc_I
    }
  }
  
  cor_strain <- bind_rows(cor_strain_list)
  cor_strain_I <- bind_rows(cor_strain_I_list)
  acc_strain <- bind_rows(acc_strain_list)
  acc_strain_I <- bind_rows(acc_strain_I_list)
  
  s3_df <- bind_rows(cor_strain, cor_strain_I) %>% 
    mutate(Info = paste0(Model, Interaction)) %>% 
    dplyr::select(Info, Matrix, Mix, Cor) %>% 
    mutate(across(c(Info, Mix, Matrix), as.factor))
  
  acc_df <- bind_rows(acc_strain, acc_strain_I) %>% 
    mutate(Info = paste0(Model, Interaction)) %>% 
    dplyr::select(Info, Matrix, Mix, Cor) %>% 
    mutate(across(c(Info, Mix, Matrix), as.factor))
  
  return(list(ability = s3_df,  accuracy = acc_df))
}

combine_elements <- function(data_list, element_name) {
  # Check if the element exists in the first list item for safety
  if (!element_name %in% names(data_list[[1]])) {
    stop(paste("Element", element_name, "not found in the list items."))
  }
  
  # Extract and combine the specified element from all sublists
  combined_df <- data_list %>%
    purrr::map(element_name) %>%
    dplyr::bind_rows()
  
  return(combined_df)
}

plotCV <- function(results, colors, stat, strategy, subheader){
  subheader_main <- "The boxplots are colored based on the different Genomic Relationship Matrices employed in the calculations. These are described in the legend following the structure Wheat_GRM/Mixes_GRM.
Additionally, 'G' denotes the GRM calculated from the marker matrix, while 'I' represent the identity."
  
  
  p <- ggplot(results, aes(x = Mix, y = Cor, fill = Matrix)) +
    geom_boxplot(width = 0.6) +
    scale_fill_manual(values = colors) + 
    scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 0.1)) +
    labs(
      title = paste0("Strategy ", strategy),
      subtitle = subheader_main,
      caption = subheader, 
      x = NULL,
      y = paste("prediction", stat)
    ) +
    theme(
      plot.subtitle = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 10, b = 10)),
      legend.title = element_blank(), 
      legend.position = "top",
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "lightgray", linewidth = 0.3),
      plot.title = element_text(hjust = 0, size = 18, face = "bold", family = "Arial"),
      strip.text = element_text(size = 10, color = "black", family = "Arial"),
      strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.title.y = element_text(size = 12, family = "Arial"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, family = "Arial"),
      plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 20, b = 20))
    ) + 
    facet_grid(. ~ Info)
  
  return(p)
}

results2plot <- function(df, name, colors, stat, strategy, subheader){
  plot <- plotCV(df, colors, stat, strategy, subheader)
  # Guardar el gráfico
  png(paste0("outputs/plots/", name, ".png"), width = 6000, height = 3000, res = 400)
  print(plot)
  dev.off()
}

computeH2 <- function(model, interaction = NULL){
  varGw <- model$ETA[[2]]$varU
  varGm <- model$ETA[[3]]$varU
  varE <- model$varE
  
  if(is.null(interaction)){
    Hw <- varGw/(varGw + (varE/16))
    Hm <- varGm/(varGm + (varE/800))
  }else{
    varGE <- model$ETA[[4]]$varU
    Hw <- varGw/(varGw + (varE/16) + (varGE/800))
    Hm <- varGm/(varGm + (varE/800) + (varGE/800))
  }
  
  H2 <- (Hw + Hm)/2
  return(H2)
}

cv_septoria <- function(genotype, phenotype, kinship, map, test, trait, blues_all,  wModel = FALSE){
  #===============================================
  # Create data for train and test
  # ==============================================
  train <- setdiff(rownames(kinship), test)
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  ptrain <- phenotype %>% filter(Isolate %in% train)
  blues <- blues_all %>% filter(Isolate %in% train)
  gtrain <- genotype %>% filter(genotype[,1] %in% blues$Isolate)
  ktrain <- kinship %>%
    filter(rownames(.) %in% blues$Isolate) %>%
    dplyr::select(all_of(blues$Isolate)) %>%
    rownames_to_column("Isolate") %>%
    dplyr::select(Isolate, everything())
  
  dim(blues); dim(gtrain); dim(map); dim(ktrain)
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  if (!wModel) {
    formula_blups <- as.formula(paste(trait, "~ Line + Year + Trial + Leaf + BRep"))
  } 
  if (wModel) {
    tmp <- capture.output({
      scores <- GAPIT(Y = blues,
                      GD = gtrain,
                      GM = map,
                      KI = ktrain,
                      CV = NULL,
                      PCA.total = 3,
                      model = "Blink",
                      file.output = F)
    })
    rm(tmp)
    gc()
    setwd(here())
    
    results <- scores[["GWAS"]] %>% 
      arrange(P.value)
    hits_bonferroni <- results %>% filter(P.value <= bonferroni)
    
    if (nrow(hits_bonferroni) == 0) {
      sSNPs <- results$SNP[1]
    }
    if (nrow(hits_bonferroni) >= 3) {
      sSNPs <- results$SNP[1:3]
    }
    if (nrow(hits_bonferroni) == 2) {
      sSNPs <- results$SNP[2]
    }
    if (nrow(hits_bonferroni) == 1) {
      sSNPs <- results$SNP[2]
    }
    sSNPs_data <- gtrain[, sSNPs, drop = FALSE] %>%
      mutate(Isolate = gtrain[,1])
    
    ptrain <- ptrain %>% 
      left_join(sSNPs_data)
    
    formula_blups <- as.formula(
      paste0(trait, "~ Line + Year + Trial + Leaf + BRep ", paste(sSNPs, collapse = " + "))
    )
  }
  
  model <- mmer(formula_blups,
                random = ~ vsr(Isolate, Gu = as.matrix(kinship)),
                rcov = ~ units,
                data = ptrain,
                verbose = TRUE)
  H2 <- h2_sommer(model, n = 12)
  
  blups_test <- data.frame(Isolate = names(model$U[[1]][[1]])) %>%
    mutate(!!trait := model$U[[1]][[1]]) %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  
  blues_test <- blues_all %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  
  ability <- cor(blups_test[,trait], blues_test[,trait])
  accuracy <- ability/sqrt(H2)
  
  results <- list(ability = ability, accuracy = accuracy)
  return(results)
  
}

h2_sommer <- function(model, n){
  varG <- model$sigma[[1]]
  varE <- model$sigma[[2]]
  h2 <- varG/(varG+(varE/n))
  
  return(h2)
}

combine_trait <- function(trait_name) {
  # Crear una lista vacía para almacenar resultados
  results <- list()
  
  # Recorrer cada pheno y modelo, extrayendo el trait correspondiente
  for (pheno_name in names(blups_list)) {
    for (model_name in names(blups_list[[pheno_name]])) {
      # Extraer la data frame del trait
      trait_data <- blups_list[[pheno_name]][[model_name]][[trait_name]]
      
      # Verificar si el trait_data existe
      if (!is.null(trait_data)) {
        # Extraer la columna BLUP y añadirla a los resultados
        results[[paste0(pheno_name, "_", model_name)]] <- trait_data$BLUP
        # Conservar los rownames
        names(results[[paste0(pheno_name, "_", model_name)]]) <- rownames(trait_data)
      }
    }
  }
  
  # Combinar las columnas en un data frame
  results_df <- as.data.frame(results)
  
  # Mantener los rownames
  rownames(results_df) <- rownames(trait_data)
  
  return(results_df)
}

process_and_filter <- function(trait_name, test) {
  combine_trait(trait_name) %>%
    filter(rownames(.) %in% test) %>%
    arrange(rownames(.))
}

correlation_dfs <- function(ref_df, ref_col, df) {
  ref_values <- ref_df[[ref_col]]
  # Calculate the correlation with each column in df
  correlations <- sapply(colnames(df), function(col) {
    cor(ref_values, df[[col]], use = "complete.obs", method = "pearson")
  })
  # Return a data frame with column names and corresponding correlation values
  result <- data.frame(Model = names(correlations), Correlation = correlations)
  rownames(result) <- NULL
  return(result)
}

content2table = function(r){
  a = strsplit(httr::content(r, encoding = 'UTF-8'), '\n')[[1]]
  b = data.frame(do.call('rbind',sapply(a, strsplit, split = '\t')))
  rownames(b) = NULL
  colnames(b) = b[1,]
  b = b[-1,]
  for (i in 1:nrow(b)){
    if (b[i,1] == b[i,2]){b[i,2] = NA}
  }
  return(b)
}

GOwide2long <- function(b){
  output = list()
  for (i in 1:nrow(b)){
    c = str_match(strsplit(b[i,2], '; ')[[1]], '(.+) (.+)')
    ########
    if(length(c[,2]) > 0){
      GO_name = c[,2]
    } else {
      GO_name = NA
    }
    
    if(length(c[,3]) > 0){
      GO_code = c[,3]
    } else {
      GO_code = NA
    }
    ########
    # output[[i]] = data.frame(Entry = unname(b[i,1]),
    #                          Prot_name = b[i,3], 
    #                          GO_name = c[,2], GO_code = c[,3])
    output[[i]] = data.frame(Entry = unname(b[i,1]),
                             Prot_name = b[i,3], 
                             Prot_existence = b[i,4],
                             GO_name = GO_name, GO_code = GO_code)
  }
  output = do.call('rbind', output)
  return(output)
}

grepEffect <- function(marker, trait, model, parent_directory){
  file <- paste0(parent_directory, "/GAPIT.Association.GWAS_Results.", model,".", trait,".csv")
  pvalues <- read_csv(file)
  effect <- pvalues %>% 
    filter(SNP == marker) %>% 
    pull(Effect)
  
  tmp_df <- data.frame(SNP = marker, Effect = effect, Trait = trait)
  return(tmp_df)
}

plotAllelicdiff <- function(phenotype, genotype, marker, trait){
  # Seleccionar la columna del marcador de genotipo
  tmp_geno <- genotype %>% 
    dplyr::select(GenoID, all_of(marker)) %>%
    mutate(across(all_of(marker), round))
  
  tmp_pheno <- phenotype %>%
    left_join(tmp_geno, by = "GenoID")
  
  clean_trait <- gsub("^BLINK.", "", trait)
  
  tmp_pheno <- tmp_pheno %>%
    filter(.data[[marker]] %in% c(0, 1, 2)) %>%
    mutate(SNP = marker) %>%
    dplyr::select(GenoID, all_of(marker), SNP, clean_trait)
  
  
  if(length(unique(tmp_pheno[[marker]])) == 2) {
    # Realiza un t-test para los dos valores únicos
    test_results <- t.test(
      tmp_pheno[[trait]] ~ as.factor(tmp_pheno[[marker]]), 
      data = tmp_pheno
    )
    p_value <- test_results$p.value
    pvalue_df <- data.frame(level_1= names(table(tmp_geno[[marker]])[1]), 
                            level_2 = names(table(tmp_geno[[marker]])[2]),
                            pvalue = p_value) %>% 
      mutate(significance = ifelse(pvalue < 0.05, "*", "NS")) %>% 
      distinct()
  } else if(length(unique(tmp_pheno[[marker]])) == 3) {
    # Realiza una ANOVA para los tres valores
    tmp_pheno[[marker]] <- as.factor(tmp_pheno[[marker]])
    anova_formula <- as.formula(paste0("`", trait, "` ~ `", marker, "`"))
    Anova <- aov(anova_formula, data = tmp_pheno)
    tukey_results <- TukeyHSD(x = Anova)
    p_value <- tukey_results[[marker]][, "p adj"]
    pvalue_df <- data.frame(pvalue = p_value) %>% 
      rownames_to_column("levels") %>% 
      mutate(level_1 = sapply(str_split(levels, "-"), function(x) x[[1]]),
             level_2 = sapply(str_split(levels, "-"), function(x) x[[2]])) %>% 
      dplyr::select(level_1, level_2, pvalue) %>% 
      mutate(significance = ifelse(pvalue < 0.05, "*", "NS"))
  } else {
    message("El marcador no tiene exactamente 2 o 3 valores únicos. No se realizó ningún análisis.")
  }
  
  # Calcular los conteos y eliminar clases con menos de 2 valores
  # valid_classes <- tmp_pheno %>%
  #   group_by(!!sym(marker)) %>%
  #   summarise(n = n(), .groups = 'drop') %>%
  #   filter(n >= 2) %>%  # Mantener solo clases con al menos 2 valores
  #   pull(!!sym(marker)) # Extraer los valores válidos del marcador
  # 
  # # Filtrar tmp_pheno para mantener solo clases válidas
  # tmp_pheno <- tmp_pheno %>%
  #   filter(.data[[marker]] %in% valid_classes)
  
  # Recalcular los conteos para las clases válidas
  counts <- tmp_pheno %>%
    group_by(!!sym(marker)) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(label = paste0(as.factor(!!sym(marker)), "\n", "n = ", n))
  
  
  
  # Add a dataframe for significance annotations with larger spacing
  annotation_df <- pvalue_df %>%
    mutate(
      x = as.numeric(factor(level_1, levels = levels(as.factor(tmp_pheno[[marker]])))),  # Ajusta las coordenadas x
      xend = as.numeric(factor(level_2, levels = levels(as.factor(tmp_pheno[[marker]])))), # Ajusta las coordenadas xend
      y_base = max(tmp_pheno[[clean_trait]], na.rm = TRUE) + 0.1 * diff(range(tmp_pheno[[clean_trait]])), # Base y-value
      y_offset = seq(0, by = 0.1 * diff(range(tmp_pheno[[clean_trait]])), length.out = nrow(pvalue_df)), # Más espacio entre líneas
      y = y_base + y_offset,  # Calcula el y final
      label_y = y + 0.03 * diff(range(tmp_pheno[[clean_trait]])), # Etiqueta ligeramente por encima del segmento
      label = significance
    )
  
  markers_names <- read_delim("data/raw_data/90K.CSv2.Ann.info.txt", delim = "\t")
  markers_names <- markers_names %>% dplyr::select(ID, Name) %>% 
    filter(ID %in% marker)
  colnames(tmp_pheno)[2] <- markers_names$Name
  tmp_pheno$SNP <- markers_names$Name
  marker <- markers_names$Name
  # Create the plot
  p <- ggplot(tmp_pheno, aes(x = as.factor(.data[[marker]]),    
                             y = .data[[clean_trait]],
                             fill = as.factor(.data[[marker]]))) +
    
    # Half violin plot
    stat_halfeye(adjust = .5, width = .3, justification = -.7, .width = 0, point_colour = NA) +
    
    # Dot plot
    stat_dots(side = "left",
              justification = 1.5,
              binwidth = NA,
              dotsize = 0.8,
              overflow = "compress",
              scale = 0.4)  +
    
    # Boxplot
    geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", aes(fill = as.factor(.data[[marker]]))) +
    
    # Escala de colores
    scale_fill_manual(values = c("#DD5129FF", "#0F7BA2FF", "#43B284FF")) +
    scale_x_discrete(
      labels = counts$label,
      expand = c(0,0)) +# Ajusta los márgenes izquierdo (0.05) y derecho (0.1)
    # Estilo de tema
    theme_ipsum() + 
    theme(
      plot.margin = unit(c(1, 1, 1, 0), "lines"),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_text(size = 13, face = "plain", color = "black", hjust = 0.5),
      strip.background = element_rect(fill = "lightgray"),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.title.x = element_text(size = 18, face = "plain"),
      axis.title.y = element_text(size = 18, face = "plain"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )  +
    labs(x = NULL, y = trait) +
    facet_grid(. ~ SNP) +
    
    # Add significance annotations
    geom_segment(data = annotation_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
    geom_text(data = annotation_df, aes(x = (x + xend) / 2, y = label_y, label = label), inherit.aes = FALSE, size = 5) +
    
    # Adjust the plot range to include annotations
    coord_cartesian(ylim = c(min(tmp_pheno[[clean_trait]], na.rm = TRUE), 
                             max(annotation_df$label_y + 0.1 * diff(range(tmp_pheno[[clean_trait]])))),
                    xlim = c(min = 0.3, max = ifelse(length(unique(tmp_pheno[[marker]])) == 3, 3.5, 2.5)))
  
  return(p)
}

plotCMManhattan <- function(df, name, color, dir) {
  setwd(dir)
  CMplot(df,
         plot.type = "m",
         multraits = FALSE,
         col = color,
         threshold = 0.05/nrow(df),
         threshold.lty = c(1, 2), 
         threshold.lwd = c(1, 1),
         threshold.col = c("black", "grey"),
         amplify = TRUE,
         bin.size = 1e6,
         chr.den.col = NULL,
         signal.col = NULL,
         signal.cex = 1,
         file = "jpg",
         file.name = name,
         dpi = 300,
         file.output = TRUE,
         verbose = FALSE,
         points.alpha = 250,
         legend.ncol = 2,
         legend.pos = "middle",
         main = name)
  setwd(here())
}

plotCMqq <- function(df, name, color, dir){
  setwd(dir)
  CMplot(df,
         main= name,
         plot.type="q",
         col = color,
         box= T,
         file="jpg",
         file.name= name,
         dpi=300,
         conf.int=TRUE,
         conf.int.col=NULL,
         threshold.col="red",
         threshold.lty=2,
         amplify = TRUE,
         signal.pch=13, 
         file.output=T,
         verbose=TRUE,
         width=5,
         height=5)
  setwd(here())
}
