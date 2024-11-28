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

runS1 <- function(trait, Kw, Kmix, pheno, genoW, map, wtest, wModel = NULL){
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  
  ptrain <- pheno
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  # prepare train data
  gtrain <- data.frame(genoW[genoW[,1] %in% wtrain, ]) # subset genotype
  ptrain_GWAS <- pheno[pheno$Plant %in% gtrain[,1], ] # subset pheno
  Ktrain <- data.frame(A.mat(gtrain[,-1])) 
  colnames(Ktrain) <- gtrain[,1]
  Ktrain <- Ktrain %>% #subset K
    mutate(GenoID = gtrain[,1]) %>% 
    dplyr::select(GenoID, everything())
  
  formula <- "~ Strain + Rep + Leaf + (1|Plant)"
  BLUPs <- extract_blups_df(ptrain_GWAS, "PLACL", formula) # obtain blups
  gtrain_GWAS <- gtrain %>% # adjust genotype to blups lines
    filter(GenoID %in% BLUPs$GenoID)
  Ktrain_GWAS <- Ktrain %>% # adjust kinship to blups lines
    filter(GenoID %in% BLUPs$GenoID) %>% 
    dplyr::select(c(GenoID, which(colnames(Ktrain) %in% BLUPs$GenoID)))
  
  dim(BLUPs); dim(gtrain_GWAS); dim(map); dim(Ktrain_GWAS)
  # run GWAS
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUPs,
                    GD = gtrain_GWAS,
                    GM = map,
                    KI = Ktrain_GWAS,
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
  
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  if (!wModel) {
    Xwtrain <- model.matrix(~1, data = ptrain) 
  } 
  if (wModel) {
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
    
    Xwtrain <- model.matrix(~1, data = ptrain) %>%
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
  
  model1_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  return(list(
    S1 = model1, # results model without interaction term
    S1_I = model1_I, # results with interaction term
    wtest = wtest, # partition test 
    wtrain = wtrain,
    hits = results
  ))
}

eval_S1 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$S1)
  predictions_I <- predict(strategy$S1_I)
  
  lines_test <- phenotype$Plant %in% strategy$wtest
  
  pheno_test <- phenotype[lines_test,]
  
  cor <- cor(predictions[lines_test],
             pheno_test[[trait]],
             use = "complete.obs")
  cor_I <- cor(predictions_I[lines_test],
               pheno_test[[trait]],
               use = "complete.obs")
  
  withinStrainCorrelations <- function(predictions, phenotypeData, testLocations) {
    data <- data.frame(predictions = predictions[testLocations],
                       traits = phenotypeData[[trait]][testLocations],
                       Strain = phenotypeData$Strain[testLocations])
    
    listData <- split(data, data$Strain)
    lapply(listData, function(df) cor(df$predictions, df$traits, use = "complete.obs"))
  }
  
  cor_withinStrain <- withinStrainCorrelations(predictions, phenotype, lines_test)
  cor_withinStrain_I <- withinStrainCorrelations(predictions_I, phenotype, lines_test)
  
  cor_results <- list(
    cor = cor,
    cor_I = cor_I,
    cor_withinStrain = cor_withinStrain,
    cor_withinStrain_I = cor_withinStrain_I,
    hits <- strategy$hits
  )
  
  return(cor_results)
}

runS2 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, wModel = NULL) {
  #===============================================
  # Run GWAS 
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  Kw <- Kw %>% 
    rownames_to_column("GenoID") %>% 
    dplyr::select(GenoID, everything())
  
  
  formula <- "~ Strain + Rep + Leaf + (1|Plant)"
  BLUPs <- extract_blups_df(phenotype, "PLACL", formula) 
  
  dim(BLUPs); dim(genoW); dim(map); dim(Kw)
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUPs,
                    GD = genoW,
                    GM = map,
                    KI = Kw,
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
  
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~1, data = ptrain) 
  } 
  if (wModel) {
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
    
    Xwtrain <- model.matrix(~1, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Xwtrain <- model.matrix(~1, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  
  K12_mix <- Zmixtrain %*% as.matrix(Kmix) %*% t(Zmixtrain)
  K12_wheat <- Zwtrain %*% as.matrix(Kw[,-1]) %*% t(Zwtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model2 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  model2_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  return(list(model2 = model2,
              model2_I = model2_I,
              sMix = sMix))
}

eval_S2 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model2)
  predictions_I <- predict(strategy$model2_I)
  
  mix_test <- phenotype$Strain %in% strategy$sMix
  
  ptest <- phenotype[mix_test,]
  ptestTrait <- ptest[[trait]]
  
  cor <- cor(predictions[mix_test], ptestTrait)
  cor_I <- cor(predictions_I[mix_test], ptestTrait)
  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I 
  )
  
  return(list(CorrelationResults = correlationResults))
}

run_S3 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, wtest, wModel = NULL) {
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  gtrain <- data.frame(genoW[genoW[,1] %in% wtrain, ])
  
  ptrain_GWAS <- phenotype[phenotype$Plant %in% gtrain[,1], ]
  Ktrain <- data.frame(A.mat(gtrain[,-1]))
  colnames(Ktrain) <- rownames(Ktrain) <- gtrain[,1]
  Ktrain <- Ktrain %>% 
    rownames_to_column("GenoID") %>% 
    dplyr::select(GenoID, everything())
  rownames(Ktrain) <- NULL
  
  formula <- "~ Strain + Rep + Leaf + (1|Plant)"
  BLUPs <- extract_blups_df(ptrain_GWAS, "PLACL", formula) 
  gtrain_GWAS <- gtrain %>% # adjust genotype to blups lines
    filter(GenoID %in% BLUPs$GenoID)
  Ktrain_GWAS <- Ktrain %>% # adjust kinship to blups lines
    filter(GenoID %in% BLUPs$GenoID) %>% 
    dplyr::select(c(GenoID, which(colnames(Ktrain) %in% BLUPs$GenoID)))
  
  dim(BLUPs); dim(gtrain_GWAS); dim(map); dim(Ktrain_GWAS)
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUPs,
                    GD = gtrain_GWAS,
                    GM = map,
                    KI = Ktrain_GWAS,
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
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~1, data = ptrain) 
  } 
  if (wModel) {
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
    
    Xwtrain <- model.matrix(~1, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Xwtrain <- model.matrix(~1, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  
  K12_mix <- Zmixtrain %*% as.matrix(Kmix) %*% t(Zmixtrain)
  K12_wheat <- Zwtrain %*% as.matrix(Kw) %*% t(Zwtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model3 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  model3_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  return(list(model3 = model3,
              model3_I = model3_I,
              sMix = sMix,
              wtest = wtest))
}

eval_S3 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model3)
  predictions_I <- predict(strategy$model3_I)
  
  mix_test <- phenotype$Strain %in% strategy$sMix
  wheat_test <- phenotype$Plant %in% strategy$wtest
  
  all_test <- mix_test & wheat_test
  
  ptest <- phenotype[all_test,]
  ptestTrait <- ptest[[trait]]
  
  cor <- cor(predictions[all_test], ptestTrait)
  cor_I <- cor(predictions_I[all_test], ptestTrait)
  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I
  )
  
  return(list(CorrelationResults = correlationResults))
}