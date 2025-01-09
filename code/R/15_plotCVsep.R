library(tidyverse)
library(ggdist)
library(gt)
source("code/R/function_septoria_GS.R")

load("data/modified_data/cv_septoria/iter_4.Rdata")

models <- c("normal", "weighted")
traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")

final_list <- list()
for(i in 1:30){
  if(i != 20){
    file_path <- paste0("data/modified_data/cv_septoria/iter_", i, ".Rdata")
    load(file_path)
    names(results_lits) <- c("normal", "weighted")
    for(model in models){
      for(trait in traits){
        final_list[[paste0("iter", i)]][[model]][[trait]] <- extract_predictions(results_lits, model, trait)
      }
    } 
  }
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

png(paste0("outputs/plots/cv_sep_accurac.png"), width = 3000, height = 3000, res = 400)
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
    file_path <- paste0("data/modified_data/cv_septoria/iter_", i, ".Rdata")
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

png(paste0("outputs/plots/hist_cv_sep.png"), width = 3000, height = 1500, res = 400)
hits_cv_sep
dev.off()





