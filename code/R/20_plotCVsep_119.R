library(tidyverse)
library(ggdist)
library(gt)
library(factoextra)
library(ggiraph)
library(ggrepel)
library(patchwork)
library(ggbeeswarm)
source("code/R/function_septoria_GS.R")

load("data/modified_data/1_septoria.Rdata")
load("data/modified_data/cv_septoria/iter_4.Rdata")
load("data/modified_data/6_CV_septoria.Rdata")
load("data/modified_data/5_predictions.Rdata")

models <- c("normal", "weighted")
traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")

final_list <- list()
partitions <- list()
for(i in 1:30){
  file_path <- paste0("data/modified_data/cv_septoria_119/iter_", i, ".Rdata")
  load(file_path)
  names(results_lits) <- c("normal", "weighted")
  for(model in models){
    for(trait in traits){
      final_list[[paste0("iter", i)]][[model]][[trait]] <- extract_predictions(results_lits, model, trait, i)
    }
  } 
  partitions[[i]] <- results_lits |> pluck(3)
}

PLACL <- prepare_traits(final_list, "PLACL")
pycnidiaLeaf <- prepare_traits(final_list, "pycnidiaPerCm2Leaf")
pycnidiaLesion <- prepare_traits(final_list, "pycnidiaPerCm2Lesion")
plot_df <- rbind(PLACL, pycnidiaLeaf, pycnidiaLesion)

colors <- c("PLACL.normal" = "#DD5129FF",
            "PLACL.weighted" = "#F59B6AFF", 
            "pycnidiaPerCm2Leaf.normal" = "#0F7BA2FF",
            "pycnidiaPerCm2Leaf.weighted" = "#58A9D4FF", 
            "pycnidiaPerCm2Lesion.normal" = "#43B284FF",
            "pycnidiaPerCm2Lesion.weighted" = "#75D1A0FF")

plot_cv_sep <- ggplot(plot_df) +
  # Boxplot with color by combination of Trait and Model
  geom_boxplot(aes(x = Trait, y = Accuracy, fill = interaction(Trait, Model)),  
               width = 0.5) +
  geom_jitter(aes(x = Trait, y = Accuracy, fill = interaction(Trait, Model)),
              size = 2, alpha = 0.9, 
              color = "black",   # Black border
              stroke = 0.5,      # Border width
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5),
              shape = 21)  +   # Jitter y dodge juntos + 
  ggdist::stat_slab(data = plot_df[plot_df$Model == "weighted",] ,
                    aes(x = Trait, y = Accuracy, fill = interaction(Trait, Model),
                        side = "right"),
                    alpha = 0.6,
                    width = 0.4,
                    position = position_nudge(x = 0.1)) +
  ggdist::stat_slab(data = plot_df[plot_df$Model == "normal",] ,
                    aes(x = Trait, y = Accuracy, fill = interaction(Trait, Model),
                        side = "left"),
                    alpha = 0.6,
                    width = 0.4,
                    position = position_nudge(x = -0.1)) +
  scale_fill_manual(values = colors, 
                    breaks = c("PLACL.normal", "PLACL.weighted", "pycnidiaPerCm2Leaf.normal",
                               "pycnidiaPerCm2Leaf.weighted", "pycnidiaPerCm2Lesion.normal", "pycnidiaPerCm2Lesion.weighted")) +
  scale_x_discrete(expand = c(0.2, 0.2)) + 
  labs(
    x = NULL, 
    y = "Accuracy"
  ) +
  scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 0.2)) +
  # Estilo y temas
  theme(
    plot.subtitle = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 10, b = 10)),
    legend.title = element_blank(), 
    legend.position = 'top',
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "lightgray", linewidth = 0.3),
    plot.title = element_text(hjust = 0, size = 18, face = "bold", family = "Arial"),
    strip.text = element_text(size = 10, color = "black", family = "Arial"),
    strip.background = element_rect(fill = "lightgray", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 12, family = "Arial"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, family = "Arial"),
    plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 20, b = 20)),
    plot.margin = margin(2, 30, 2, 30)
  )

png(paste0("outputs/plots/cv_sep_119_accurac.png"), width = 3000, height = 3000, res = 400)
plot_cv_sep
dev.off()

plot_df_summary <- plot_df %>%
  group_by(Trait, Model) %>%
  summarise(Mean = mean(Accuracy, na.rm = TRUE), 
            SD = sd(Accuracy, na.rm = TRUE)) %>%
  mutate(Final = paste0(round(Mean, 3), " ± ", round(SD, 2))) %>%
  dplyr::select(Trait, Model, Final)


load("data/modified_data/1_septoria.Rdata")
bonferroni <- 0.05/nrow(map_septoria)

hits_list <- list()  # Inicializas la lista vacía
for(i in 1:30) {
  if(i != 20) {
    file_path <- paste0("data/modified_data/cv_septoria_119//iter_", i, ".Rdata")
    load(file_path)
    for(trait in traits) {
      # Llenar la lista con el valor de hits para cada trait e iteración
      hits_list[[paste0("iter", i)]][[trait]] <- extract_hits_cv_sep(results_lits, trait, bonferroni)
    }
  }
}

hits_list_all <- map(traits, ~ hits_df(hits_list, .x))
hits_df <- bind_rows(hits_list_all)

hits_cv_sep <- ggplot(hits_df) +
  geom_histogram(aes(x = Hits)) +
  facet_grid(.~ Trait)

png(paste0("outputs/plots/hist_cv_sep_119.png"), width = 3000, height = 1500, res = 400)
hits_cv_sep
dev.off()

genotype <- genotype_all
row_means <- rowMeans(genotype[, -1])
col_means <- colMeans(genotype[, -1])
overall_mean <- mean(as.matrix(genotype[, -1]))
genotype[, -1] <- genotype[, -1] - row_means + overall_mean
genotype[, -1] <- t(t(genotype[, -1]) - col_means)

# compute the PC analysis
PCA <- prcomp(genotype[, -1], center = F)
PCs <- PCA$x
# Prepare the data for PCA plot
pca_data <- as.data.frame(PCs) %>% 
  mutate(GenoID = genotype[,1], 
         shape = str_extract(GenoID, "([0-9]{2})_", group = T))

right_group <- pca_data |> 
  filter(PC1 > 100) |> 
  pull(GenoID)

extra <- list()
for(i in 1:30){
  test <- partitions[[i]]
  extra[[i]] <- sum(test %in% right_group)
  new_df <- pca_data |> 
    mutate(regions = ifelse(GenoID %in% test, "Test", "Train"))
  
  pca_plot <- ggplot(data = new_df) +
    geom_point_interactive(aes(x = PC1, y = PC2, color = regions, shape = shape, tooltip = GenoID, data_id = GenoID),
                           size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""),
      y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = "")
    ) +  # Axis labels
    theme_minimal() +  # Clean theme
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",  # Move legend to the top
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
    scale_color_manual(values = c("Train" = "grey", 
                                  "Test" = "red")) +
    scale_shape_manual(values = c("21" = 17,
                                  "22" = 16)) +
    scale_x_continuous(
      limits = c(-50, 200),  # Set x-axis limits
      breaks = seq(-50, 200, by = 50)  # Set breaks at intervals of 50
    )
  
  new_df2 <- plot_df |> 
    mutate(color = ifelse(Iter == i, "color", "nothing"),
           interaction_ordered = factor(interaction(Trait, Model), 
                                        levels = c("PLACL.normal",
                                                   "PLACL.weighted", 
                                                   "pycnidiaPerCm2Leaf.normal",
                                                   "pycnidiaPerCm2Leaf.weighted", 
                                                   "pycnidiaPerCm2Lesion.normal",
                                                   "pycnidiaPerCm2Lesion.weighted")))
  
  accuracy_plot <- ggplot(new_df2) +
    # Boxplot with color by combination of Trait and Model
    geom_boxplot(aes(x = interaction_ordered, y = Accuracy),  
                 width = 0.5, fill = "darkgrey") +
    geom_point(aes(x = interaction_ordered, y = Accuracy, fill = color),
               size = 3, alpha = 0.9,
               stroke = 0.3,  # Border width
               shape = 21, color = "black") +  # Replace jitter with fixed points
    scale_fill_manual(values = c("color" = "red", 
                                 "nothing" = "grey")) +
    scale_x_discrete(expand = c(0.2, 0.2)) + 
    labs(
      x = NULL, 
      y = "Accuracy"
    ) +
    scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 0.2)) +
    # Estilo y temas
    theme(
      plot.subtitle = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 10, b = 10)),
      legend.title = element_blank(), 
      legend.position = 'none',
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "lightgray", linewidth = 0.3),
      plot.title = element_text(hjust = 0, size = 18, face = "bold", family = "Arial"),
      strip.text = element_text(size = 10, color = "black", family = "Arial"),
      strip.background = element_rect(fill = "lightgray", size = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.title.y = element_text(size = 12, family = "Arial"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10, family = "Arial"),
      plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 20, b = 20)),
      plot.margin = margin(2, 30, 2, 30)
    ) + 
    coord_flip()
  
  
  patchwork <- pca_plot + accuracy_plot +
    plot_annotation(
      title = paste("Iteration", i),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  # Save Output
  png(paste0("outputs/pca_accuracy_119/iter", i, ".png"), width = 7000, height = 4000, res = 400)
  print(patchwork) # Ensure patchwork is printed
  dev.off()
  
}

df <- data.frame(Iter = 1:30,
                 n_right_group = unlist(extra)) 

extra_df <- df |> 
  left_join(plot_df, by = "Iter") |>  # Join the dataframes by Iter
  mutate(Trait_Model = paste0(Trait, "_", Model)) |>  # Combine Trait and Model into a new column
  dplyr::select(-c(Model, Ability, Trait)) |>  # Remove Model, Ability, and Trait columns
  pivot_wider(
    names_from = Trait_Model,  # Column names will come from the Trait_Model
    values_from = Accuracy  # Values will come from the Accuracy column
  ) |> 
  group_by(n_right_group) |>  # Group by n_right_group
  summarise(
    n = n(),  # Count the number of observations for each group
    across(-c(Iter), list(mean = ~mean(.x, na.rm = TRUE), 
                          sd = ~sd(.x, na.rm = TRUE)), 
           .names = "{.col}_{.fn}")  # Create new columns for mean and SD
  ) |> 
  ungroup() |>  
  arrange(desc(PLACL_normal_mean))

ggplot(extra_df) + 
  geom_point(aes(x = n_right_group, y = PLACL_normal_mean, size = n), color = "black") + 
  geom_line(aes(x = n_right_group, y = PLACL_normal_mean), color = "black") +
  geom_errorbar(aes(x = n_right_group, 
                    ymin = PLACL_normal_mean - PLACL_normal_sd, 
                    ymax = PLACL_normal_mean + PLACL_normal_sd), 
                width = 0.2, color = "blue") +
  theme_minimal() +
  labs(
    x = "n_right_group", 
    y = "PLACL Normal (Mean ± SD)", 
    title = "Plot of PLACL Normal with Standard Deviation",
    size = "Sample Size (n)" # Etiqueta para la leyenda del tamaño
  ) +
  scale_size_continuous(
    breaks = seq(min(extra_df$n), max(extra_df$n), by = 1) # Forzar números enteros en la leyenda
  )


