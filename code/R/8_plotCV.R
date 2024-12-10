library(tidyverse)
library(hrbrthemes)
library(ggplot2)
source("code/R/function_septoria_GS.R")

s1 <- s2 <- s3 <- list()

for( i in 1:10){
  load(paste0("data/modified_data/cv/all/iter_", i, ".Rdata"))
  names(allResults) <- c("G/G", "G/I", "I/G", "I/I")
  s1[[i]] <- scenario1(allResults = allResults)
  s2[[i]] <- scenario2(allResults = allResults)
  s3[[i]] <- scenario3(allResults = allResults)
}

colors <- c("G/G" = "#DD5129FF",
            "G/I" = "#0F7BA2FF",
            "I/G" = "#43B284FF",
            "I/I" = "#FAB255FF")

results_list <- list(s1, s2, s3)
names <- c("S1", "S2", "S3")

for(i in seq_along(results_list)){
  list <- results_list[[i]]
  name <- names[[i]]
  results2plot(list, name, colors)
}


# New results introducing what I think is correctly the fixed effects

s1 <- s2 <- s3 <- list()

for( i in 1:10){
  load(paste0("data/modified_data/cv/all2/iter_", i, ".Rdata"))
  names(allResults) <- c("G/G", "G/I", "I/G", "I/I")
  s1[[i]] <- scenario1(allResults = allResults)
  s2[[i]] <- scenario2(allResults = allResults)
  s3[[i]] <- scenario3(allResults = allResults)
}

colors <- c("G/G" = "#DD5129FF",
            "G/I" = "#0F7BA2FF",
            "I/G" = "#43B284FF",
            "I/I" = "#FAB255FF")

results_list <- list(s1, s2, s3)
names <- c("S1_2", "S2_2", "S3_2")

for(i in seq_along(results_list)){
  list <- results_list[[i]]
  name <- names[[i]]
  results2plot(list, name, colors)
}

