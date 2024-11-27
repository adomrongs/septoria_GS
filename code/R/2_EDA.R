library(tidyverse)
library(gt)
library(ggdist)
library(gridExtra)
library(hrbrthemes)
source("code/R/function_septoria_GS.R")

#load phenotype
raw_phenotype <- readxl::read_xlsx("data/raw_data/Output_reps_unidas_12cm_mock_checks.xlsx")
# define columns to mutate
factor_cols <- c("Set", "Plant", "Strain", "Leaf", "Rep")
numeric_cols <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
# adjust phenotype
phenotype <- raw_phenotype %>% 
  mutate(Plant = paste0("PyrSep_", sprintf("%03d", Plant)),
         across(all_of(factor_cols), as.factor),
         across(all_of(numeric_cols), as.numeric)) %>% 
  dplyr::select(all_of(c(factor_cols, numeric_cols)))

write_csv(phenotype, file = "data/modified_data/clean_phenotype.csv")

# Exploratory data analysis
summary(phenotype)
n_counts <- phenotype %>%
  group_by(Strain) %>%
  summarise(
    PLACL = sum(!is.na(PLACL)),
    pycnidiaPerCm2Leaf = sum(!is.na(pycnidiaPerCm2Leaf)),
    pycnidiaPerCm2Lesion = sum(!is.na(pycnidiaPerCm2Lesion))
  )
# lets generate a few table to visualize the missing values per Strain
na_counts <- phenotype %>% 
  group_by(Strain) %>% 
  summarise(
    NA_PLACL = sum(is.na(PLACL)),
    NA_pycnidiacm2Lesion = sum(is.na(pycnidiaPerCm2Lesion)),
    NA_pycnidiacm2Leaf = sum(is.na(pycnidiaPerCm2Lesion))) %>% 
  gt()

# Create long phenotype to make it easier for plotting
phenotype_long <- phenotype %>%
  pivot_longer(cols = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
               names_to = "Trait",
               values_to = "Value")

# plot each trait based on the strain (which is the treatment)
colors <- c("check1" = "#DD5129FF",
            "check2" = "#0F7BA2FF",
            "mix1" = "#43B284FF",
            "mix2" = "#FAB255FF",
            "mix3" = "#A64C75FF",
            "mix4" = "#6A4F97FF",
            "MOCK" = "#FF5733FF")
# Create the plot for each trait per strain 
trait_strain <- plotGrid(phenotype_long = phenotype_long,
                         trait = "Strain",
                         colors = colors)

colors2 <- c("1" = "#DD5129FF",
            "2" = "#0F7BA2FF",
            "3" = "#43B284FF")
trait_leaf <- plotGrid(phenotype_long = phenotype_long,
                         trait = "Leaf",
                         colors = colors2)

colors3 <- c("set1" = "#DD5129FF",
            "set2" = "#0F7BA2FF",
            "set3" = "#43B284FF",
            "set4" = "#FAB255FF",
            "set5" = "#A64C75FF",
            "set6" = "#6A4F97FF",
            "set7" = "#FF5733FF", 
            "set8" = "#FFD700")
trait_set <- plotGrid(phenotype_long = phenotype_long,
                       trait = "Set",
                       colors = colors3)

colors4 <- c("R1" = "#DD5129FF",
             "R2" = "#0F7BA2FF",
             "R3" = "#43B284FF",
             "R4" = "#FAB255FF")
trait_rep <- plotGrid(phenotype_long = phenotype_long,
                      trait = "Rep",
                      colors = colors4)

# save results
png(paste0("outputs/plots/trait_strain.png"), width = 4000, height = 3000, res = 400)
trait_strain
dev.off()
png(paste0("outputs/plots/trait_leaf.png"), width = 3000, height = 3000, res = 400)
trait_leaf
dev.off()
png(paste0("outputs/plots/trait_set.png"), width = 5000, height = 3000, res = 400)
trait_set
dev.off()
png(paste0("outputs/plots/trait_rep.png"), width = 4000, height = 3000, res = 400)
trait_rep
dev.off()


# finally lets plot an histogram of the total values
traits <- c("PLACL" = "#0F7BA2FF",
             "pycnidiaPerCm2Leaf" = "#43B284FF",
             "pycnidiaPerCm2Lesion" = "#FAB255FF")

hist <- plotHist(phenotype, columns = numeric_cols, color = traits)
png(paste0("outputs/plots/hist_traits.png"), width = 3000, height = 3000, res = 400)
grid.arrange(grobs = hist, ncol = 3)
dev.off()


