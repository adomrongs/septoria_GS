library(tidyverse)
library(hrbrthemes)
library(ggplot2)
library(extrafont)
source("code/R/function_septoria_GS.R")

# Results including both prediction abiliyt and accuracy

s1 <- s2 <- s3 <- list()

for( i in 1:30){
  if(i != 4){
    load(paste0("data/modified_data/cv/iter_", i, ".Rdata"))
    names(allResults) <- c("G/G", "G/I", "I/G", "I/I")
    s1[[i]] <- scenario1(allResults = allResults)
    s2[[i]] <- scenario2(allResults = allResults)
    s3[[i]] <- scenario3(allResults = allResults)
  }
}


s1_ability <- combine_elements(s1, "ability")
s1_accuracy <- combine_elements(s1, "accuracy")

s2_ability <- combine_elements(s2, "ability")
s2_accuracy <- combine_elements(s2, "accuracy")

s3_ability <- combine_elements(s3, "ability")
s3_accuracy <- combine_elements(s3, "accuracy")

colors <- c("G/G" = "#72A8D4",  # Azul claro
            "G/I" = "#F7D77A",  # Amarillo claro
            "I/G" = "#A67C52",  # Marrón
            "I/I" = "#7D3243") 


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

# Genetic parameters

files <- list.files('data/modified_data/cv', full.names = T)

extract_parameters <- function(file){
  load(file)
  table <- map(allResults, \(x) x |> pluck('Scenario1', 'parameters')) |> bind_rows()
  table_I <-  map(allResults, \(x) x |> pluck('Scenario1', 'parameters_I')) |> bind_rows()
  table_w <- map(allResults, \(x) x |> pluck('Scenario1w', 'parameters')) |> bind_rows()
  table_wI <- map(allResults, \(x) x |> pluck('Scenario1w', 'parameters_I')) |> bind_rows()
  
  no_I_table <- bind_rows(table, table_w)
  I_table <- bind_rows(table_I, table_wI) |> 
    mutate(approach = case_when(
      approach == 'normal' ~ 'normal_I',
      approach == 'weighted' ~ 'weighted_I',
      TRUE ~ approach  # Esta línea cubre cualquier otro valor que no sea 'normal' ni 'weighted'
    ))
  
  final_table <- bind_rows(no_I_table, I_table) |> 
    dplyr::relocate(approach, kw, km)
  return(final_table)
}

paramteres <- map(files, \(x) extract_parameters(x))

parameters_df <- bind_rows(paramteres) |> 
  as.tibble() |> 
  mutate(across(c(kw, km), ~ ifelse(grepl('Diagonal', .x), 'I', .x)),
         Matrix = paste0(kw, '/', km),
         vari = as.numeric(ifelse(is.na(vari), '0', vari)),
         H2_wheat = varw/(varw + varm + vare + vari),
         H2_mix = varm/(varw + varm + vare + vari)) |> 
  dplyr::select(-c(kw, km)) |> 
  group_by(approach, Matrix) |> 
  summarize(across(contains('var'), ~ round(mean(.x, na.rm = TRUE), 2))) |> 
  mutate(total = rowSums(across(c(varw, varm, vari, vare)), na.rm = TRUE),
         across(contains('var'), ~ paste0(.x, ' (', round((.x / total) * 100, 2), '%)'))) |> 
  mutate(vari = ifelse(grepl('0%', vari), '0', vari)) |> 
  dplyr::select(Approach = approach, Matrix, Varw = varw, Varm = varm, VarI = vari, VarE = vare)

write_csv(parameters_df, file = 'outputs/plots/genetic_parameters.csv')



#-------------------------------------------------------------------------------
# Leaf predictions
#-------------------------------------------------------------------------------

s1 <- s2 <- s3 <- list()

for( i in 1:30){
  if(i != 4){
    load(paste0("data/modified_data/cv_leaf/iter_", i, ".Rdata"))
    names(allResults) <- c("G/G", "G/I", "I/G", "I/I")
    s1[[i]] <- scenario1(allResults = allResults)
    s2[[i]] <- scenario2(allResults = allResults)
    s3[[i]] <- scenario3(allResults = allResults)
  }
}


s1_ability <- combine_elements(s1, "ability")
s1_accuracy <- combine_elements(s1, "accuracy")

s2_ability <- combine_elements(s2, "ability")
s2_accuracy <- combine_elements(s2, "accuracy")

s3_ability <- combine_elements(s3, "ability")
s3_accuracy <- combine_elements(s3, "accuracy")

colors <- c("G/G" = "#72A8D4",  # Azul claro
            "G/I" = "#F7D77A",  # Amarillo claro
            "I/G" = "#A67C52",  # Marrón
            "I/I" = "#7D3243") 


ability_list <- list(s1_ability, s2_ability, s3_ability)
accuracy_list <- list(s1_accuracy, s2_accuracy, s3_accuracy)
results <- list(ability_list, accuracy_list)

names <- c("Strategy_1_leaf", "Strategy_2_leaf", "Strategy_3_leaf")
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

#-------------------------------------------------------------------------------
# Lesion predictions
#-------------------------------------------------------------------------------

s1 <- s2 <- s3 <- list()

for( i in 1:30){
  if(i != 4){
    load(paste0("data/modified_data/cv_lesion//iter_", i, ".Rdata"))
    names(allResults) <- c("G/G", "G/I", "I/G", "I/I")
    s1[[i]] <- scenario1(allResults = allResults)
    s2[[i]] <- scenario2(allResults = allResults)
    s3[[i]] <- scenario3(allResults = allResults)
  }
}


s1_ability <- combine_elements(s1, "ability")
s1_accuracy <- combine_elements(s1, "accuracy")

s2_ability <- combine_elements(s2, "ability")
s2_accuracy <- combine_elements(s2, "accuracy")

s3_ability <- combine_elements(s3, "ability")
s3_accuracy <- combine_elements(s3, "accuracy")

colors <- c("G/G" = "#72A8D4",  # Azul claro
            "G/I" = "#F7D77A",  # Amarillo claro
            "I/G" = "#A67C52",  # Marrón
            "I/I" = "#7D3243") 


ability_list <- list(s1_ability, s2_ability, s3_ability)
accuracy_list <- list(s1_accuracy, s2_accuracy, s3_accuracy)
results <- list(ability_list, accuracy_list)

names <- c("Strategy_1_lesion", "Strategy_2_lesion", "Strategy_3_lesion")
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














