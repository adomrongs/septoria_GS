library(sommer)
library(tidyverse)
library(here)
source('code/R/function_septoria_GS.R')

load('data/modified_data/9_cvseptoria2.Rdata')

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
traits <- list('PLACL', 'pycnidiaPerCm2Leaf', 'pycnidiaPerCm2Lesion')
test <- sample(rownames(kinship), ceiling(0.2 * nrow(kinship)))

dir <- "data/modified_data/cv_septoria2/"
dir.create(dir, recursive = T, showWarnings = FALSE)

result <- map(traits, function(trait) {
  map2(cultivars_phenotype, cultivars_blues, ~ cv_cultivar(.x, kinship, test, trait, .y))
})

df <- bind_rows(result) |> 
  t()
colnames(df) <- c('PLACL', 'PCm2Leaf', 'PCm2Lesion')

save(df, paste0(dir, "iter_", i, ".Rdata"))
