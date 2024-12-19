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
    results2plot(strategy, name, colors, stat_name, j, subheader)
  }
}

