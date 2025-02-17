library(tidyverse)
library(scales)
library(ggh4x)

files <- list.files('data/modified_data/cv_cultivars/', full.names = T)
dfs <- map(files, function(x) {
  env <- new.env()
  load(x, envir = env)
  get(ls(env), envir = env)
})

dfs <- map(dfs, as.data.frame)

accuracy_df <-  bind_rows(dfs) |> 
  rownames_to_column('Cultivar') |> 
  mutate(Cultivar = factor(str_remove(Cultivar, "\\.\\.\\..*"))) |> 
  pivot_longer(cols = -Cultivar, names_to = 'Trait', values_to = 'Accuracy') |> 
  mutate(Trait = factor(Trait))

# Definir los colores
colors_wheat <- c("Athoris" = alpha("#8E6B3D", 0.7), 
                  "Don Ricardo" = alpha("#D8B06A", 0.7), 
                  "Sculptur" = alpha("#5B4C44", 0.7),
                  "Svevo" = alpha("#9E5B40", 0.7))
colors_traits <- c("PLACL" = "#8A0000", 
                   "PCm2Leaf" = "#B5A300", 
                   "PCm2Lesion" = "#A8C300")

# Definir límites para el eje y
max_y <- round(min(max(accuracy_df$Accuracy, na.rm = TRUE) + 0.1, 1), 1)
min_val <- round(min(accuracy_df$Accuracy, na.rm = TRUE), 1)
min_y <- ifelse(min_val < 0, min_val - 0.1, 0)

strip <- strip_themed(background_x = elem_list_rect(fill = colors_wheat))
# Crear el gráfico
cv_cultivars <- ggplot(accuracy_df, aes(x = Trait, y = Accuracy, fill = Trait)) + 
  stat_summary(
    fun.data = mean_sdl, 
    fun.args = list(mult = 1), 
    geom = "errorbar", 
    width = 0.15, 
    color = 'black', 
    position = position_dodge(width = 0.7)
  ) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 23, 
    size = 3, 
    color = 'black', 
    position = position_dodge(width = 0.7)
  ) +
  scale_y_continuous(
    breaks = seq(min_y, max_y, by = 0.1),
    labels = label_number(accuracy = 0.1)
  ) + 
  scale_fill_manual(values = colors_traits) +
  theme(
    plot.subtitle = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 10, b = 10)),
    legend.title = element_blank(), 
    legend.position = "top",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "lightgray", linewidth = 0.3),
    plot.title = element_text(hjust = 0, size = 18, face = "bold", family = "Arial"),
    strip.text = element_text(size = 10, color = "black", family = "Arial"),
    strip.background = element_rect(fill = 'white', colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 12, family = "Arial"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.2, family = "Arial", margin = margin(t = 20, b = 20)),
    axis.ticks.x = element_blank()
  ) + 
  facet_wrap2(~ Cultivar, strip = strip, ncol = 4) 

png(paste0("outputs/plots/cv_cultivars.png"), width = 4000, height = 1500, res = 400)
cv_cultivars
dev.off()
