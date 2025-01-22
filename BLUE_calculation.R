# Load necessary libraries
library(tidyverse)
library(lme4)
library(emmeans)

# Define a function to calculate average values by group
calculate_average <- function(data, group_col, value_col) {
  data %>%
    group_by({{ group_col }}) %>%
    summarise(avgValue = mean({{ value_col }}, na.rm = TRUE), .groups = "drop")
}

# Define a function to calculate BLUEs
calculate_blue <- function(data, trait, fixed_effects, random_effects) {
  formula <- as.formula(paste(trait, "~", paste(fixed_effects, collapse = " + "), "+ (1|", random_effects, ")"))
  model <- lmer(formula, data = data)
  blue <- emmeans(model, ~ !!sym(fixed_effects[1]))
  as.data.frame(summary(blue))
}

# Define a function to calculate and save BLUEs for multiple traits
calculate_and_save_blue <- function(data, traits, fixed_effects, random_effects, output_prefix) {
  blue_results <- list()
  for (trait in traits) {
    cat("Calculating BLUE for trait:", trait, "\n")
    blue_results[[trait]] <- calculate_blue(data, trait, fixed_effects, random_effects)
    write_csv(blue_results[[trait]], paste0(output_prefix, "_", trait, ".csv"))
  }
  blue_results
}

# Define a function to prepare data subsets by year
subset_data_by_year <- function(data, year_col, years) {
  if (!is.null(years)) {
    data <- data[data[[year_col]] %in% years, ]
  }
  data
}

# Main script workflow
process_phenotype_data <- function() {
  # Read phenotype data
  df <- read.csv("phenotype_long.csv", header = TRUE)
  df[, 1:3] <- lapply(df[, 1:3], as.factor)

  # Subset data by year
  df_2019_2022 <- subset_data_by_year(df, "year", c(2019, 2020, 2022))
  
  # Calculate and save BLUEs for PSN, FSN, and PTD
  traits <- c("PSN", "FSN", "PTD")
  calculate_and_save_blue(df_2019_2022, traits, fixed_effects = c("name"), random_effects = "year", output_prefix = "blue_2019_2020")
  calculate_and_save_blue(df, traits, fixed_effects = c("name"), random_effects = "year", output_prefix = "blue_all")
}

# Extended functionality: Calculate BLUEs for all traits in a dataset
process_gwas_traits <- function(data_path, output_dir) {
  data <- read.csv(data_path)
  traits <- names(data)[4:ncol(data)]
  
  # Calculate BLUEs for each trait and save results
  blue_results <- calculate_and_save_blue(data, traits, fixed_effects = c("genotype"), random_effects = "rep", output_prefix = paste0(output_dir, "/blue"))
  
  # Combine and save all BLUE results
  combined_blue_results <- do.call(rbind, lapply(names(blue_results), function(trait) {
    trait_results <- blue_results[[trait]]
    trait_results$trait <- trait
    trait_results
  }))
  
  write.csv(combined_blue_results, file = paste0(output_dir, "/combined_blue_results.csv"), row.names = FALSE)
}

# Example usage
# Step 1: Process phenotype data
process_phenotype_data()

# Step 2: Process GWAS traits data
process_gwas_traits("../../phenotype/GWAS_traits_Yield_roop.csv", "results")
