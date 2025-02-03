library(tidyverse)
library(hrbrthemes)
library(ggplot2)
library(extrafont)
source("code/R/function_septoria_GS.R")

# Results including both prediction abiliyt and accuracy

s1 <- s2 <- s3 <- list()

for( i in 1:30){
  load(paste0("data/modified_data/cv/iter_", i, ".Rdata"))
  names(allResults) <- c("G/G", "G/I", "I/G", "I/I")
  s1[[i]] <- scenario1(allResults = allResults)
  s2[[i]] <- scenario2(allResults = allResults)
  s3[[i]] <- scenario3(allResults = allResults)
}

s1_ability <- combine_elements(s1, "ability")
s1_accuracy <- combine_elements(s1, "accuracy")

s2_ability <- combine_elements(s2, "ability")
s2_accuracy <- combine_elements(s2, "accuracy")

s3_ability <- combine_elements(s3, "ability")
s3_accuracy <- combine_elements(s3, "accuracy")

colors <- c("G/G" = "#DD5129FF",
            "G/I" = "#0F7BA2FF",
            "I/G" = "#43B284FF",
            "I/I" = "#FAB255FF")



ability_list <- list(s1_ability, s2_ability, s3_ability)
accuracy_list <- list(s1_accuracy, s2_accuracy, s3_accuracy)
results <- list(ability_list, accuracy_list)

names <- c("Strategy_1", "Strategy_2", "Strategy_3")
stats <- c("ability", "accuracy")
subheaders <- c("This strategy consist of a 5-Fold Cross Validation where only 80% of the lines were observed while the 20% left is predicted. All mixes (environment) were observed",
               "This strategy consist of a Leave One Out CV in which a single mix was employed for model validation. All lines were observed",
               "This strategy is a combination between a 5-Fold and LOO CV. The validation set was composed by a mix and 20% of the lines. Thus, not all lines were observed. So, we will predict unobserved phenotypes in unobserved environments")

for(i in seq_along(results)){
  list <- results[[i]]
  stat_name <- stats[[i]]
  
  for(j in seq_along(list)){
    strategy <- list[[j]]
    subheader <- subheaders[[j]]
    name <- paste0(stat_name, "_", names[[j]])
    results2plot(df = strategy,
                 name = name,
                 colors = colors,
                 stat = stat_name,
                 strategy = j,
                 subheader = subheader)
  }
}

hits_s1 <- hits_w(s1, "Strategy 1")
hits_s2 <- hits_w(s2, "Strategy 2")
hits_s3 <- hits_w(s3, "Strategy 3")

hits_df <- bind_rows(hits_s1, hits_s2, hits_s3)
hist_cv_wheat <- ggplot(hits_df) +
  geom_histogram(aes(x = Hits)) +
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
  )+
  facet_grid(.~ Strategy) 

png(paste0("outputs/plots/hist_cv_wheat.png"), width = 3000, height = 1500, res = 400)
hist_cv_wheat
dev.off()


accuracy_table <- do.call(rbind, accuracy_list) |> 
  mutate(Strategy = c(rep("S1", 1920), rep("S2", 1920), rep("S3", 1920))) |> 
  group_by(Strategy, Info, Matrix, Mix) |> 
  summarize(
    Mean = mean(Cor, na.rm = TRUE),  # Calcular la media de 'Accuracy'
    sd = sd(Cor, na.rm = TRUE)       # Calcular la desviación estándar de 'Accuracy'
  ) |> 
  mutate(
    Mean_sd = paste0(round(Mean, 2), " ± ", round(sd, 2))  # Crear la columna 'Mean_sd' con el formato deseado
  ) |> 
  dplyr::select(-c(Mean, sd)) |> 
  pivot_wider(names_from = Mix, values_from = Mean_sd) |> 
  filter(Strategy == "S1")

accuracy_table2 <- do.call(rbind, accuracy_list) |> 
  mutate(Strategy = c(rep("S1", 1920), rep("S2", 1920), rep("S3", 1920))) |>
  filter(Strategy == "S1") |> 
  group_by(Strategy, Info, Matrix, Mix) |> 
  summarize(
    Max = max(Cor, na.rm = TRUE), 
    Min = min(Cor, na.rm = TRUE)
  )


h2_strategy <- function(ability, accuracy){
  tmp_df <- cbind(ability, accuracy)[, c(1:4, 8)]  # Evitas dplyr::select()
  colnames(tmp_df) <- c('Info', 'Matrix', 'Mix', 'Ability', 'Accuracy')
  
  df <- tmp_df |> 
    dplyr::mutate(h2 = sqrt(Ability / Accuracy))
  
  return(df)
}

h2_plots <- map2(ability_list, accuracy_list, \(x,y)  h2_strategy(x,y))

map(h2_plots, \(df) {
  ggplot(df) +
    geom_boxplot(aes(x = Mix, y = h2, fill = Matrix)) +
    facet_grid(. ~ Info)
})


pheno <- readxl::read_xlsx('data/raw_data/Output_reps_unidas_12cm_mock_checks.xlsx')  

pheno |> 
  mutate(across(Set:Rep, as.factor)) |> 
  group_by(Plant) |> 
  summarize(
    nSet = n_distinct(Set),      # Cuenta valores únicos en la columna Set
    nStrain = n_distinct(Strain),# Cuenta valores únicos en la columna Strain
    nBrep = n_distinct(Leaf),    # Cuenta valores únicos en la columna Leaf
    nRep = n_distinct(Rep),
    n = n()                      # Cuenta el número total de observaciones por grupo
  )

pheno |> 
  mutate(across(Set:Rep, as.factor)) |> 
  group_by(Strain) |> 
  summarize(
    nSet = n_distinct(Set),      # Cuenta valores únicos en la columna Set
    nPlant = n_distinct(Plant),# Cuenta valores únicos en la columna Strain
    nBrep = n_distinct(Leaf),    # Cuenta valores únicos en la columna Leaf
    nRep = n_distinct(Rep),
    n = n()                      # Cuenta el número total de observaciones por grupo
  )
  








