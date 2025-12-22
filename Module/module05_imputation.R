# Module 05: Missing Value Imputation
# Purpose: Impute missing values for standardized data (mainly for NoCat group)
# Author: CodeNorm Pipeline
# Date: 2024

#' Module 05: Missing Value Imputation
#' 
#' @description
#' Impute missing values on standardized data using the Perseus approach.
#' Primarily targets groups with CatalyticGroup = NoCat; Cat groups can optionally be imputed.
#' 
#' @param dir_config Directory configuration list (from Module 1)
#' @param standardized_data_list Standardized data list (from Module 4)
#' @param sampleGroup Sample grouping information (from Module 2)
#' @param impute_cat_mean Whether to impute Cat groups with mean when n_valid = 2 (default FALSE)
#' @param random_seed Random seed for reproducibility (default 123)
#' 
#' @return List containing:
#'   - imputed_data_list: all imputed versions
#'   - imputation_params: imputation parameter info
#' 
#' @details
#' Imputation strategy:
#' 1. Perseus parameters: imputation_mean = col_mean - 1.8 * col_sd; imputation_sd = 0.3 * col_sd
#' 2. NoCat rules:
#'    - n_valid == 2: impute with group mean
#'    - n_valid < 2: sample from Perseus parameters
#' 3. Cat group: no imputation by default; optional mean imputation when n_valid == 2
#' 4. Output before/after comparison boxplot
#' 
#' @export
module05_imputation <- function(dir_config, 
                                standardized_data_list,
                                sampleGroup,
                                impute_cat_mean = FALSE,
                                random_seed = 123) {
  
  # 0. Load required packages ####
  require(tidyverse)
  require(dplyr)
  require(tidyr)
  
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("Module 05: Missing Value Imputation\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # 1. Validate input ####
  cat("\n[1] Validate input data...\n")
  
  if (!all(c("reference", "output") %in% names(dir_config))) {
    stop("❌ dir_config must include reference and output paths")
  }
  
  if (!is.list(standardized_data_list) || length(standardized_data_list) == 0) {
    stop("❌ standardized_data_list must be a non-empty list")
  }
  
  required_cols <- c("FinalName", "bioGroup", "CatalyticGroup")
  if (!all(required_cols %in% names(sampleGroup))) {
    stop("❌ sampleGroup must include: ", paste(required_cols, collapse = ", "))
  }
  
  cat("✓ Input validation passed\n")
  cat(sprintf("  - Datasets to process: %d\n", length(standardized_data_list)))
  cat(sprintf("  - Cat-group strategy: %s\n", ifelse(impute_cat_mean, "impute when n_valid=2", "no imputation")))
  cat(sprintf("  - Random seed: %d\n", random_seed))
  
  # 2. Prepare grouping info ####
  cat("\n[2] Prepare grouping info...\n")
  
  # Identify data and annotation columns from first dataset
  first_data <- standardized_data_list[[1]]
  data_cols <- sampleGroup$FinalName
  anno_cols <- setdiff(names(first_data), c("Gene", data_cols))
  
  cat(sprintf("✓ Identified %d data columns\n", length(data_cols)))
  cat(sprintf("✓ Identified %d annotation columns: %s\n", 
              length(anno_cols), 
              paste(anno_cols, collapse = ", ")))
  
  # Identify NoCat and Cat groups
  cat_groups <- sampleGroup %>%
    filter(CatalyticGroup == "Cat") %>%
    pull(bioGroup) %>%
    unique()
  
  nocat_groups <- sampleGroup %>%
    filter(CatalyticGroup == "NoCat") %>%
    pull(bioGroup) %>%
    unique()
  
  cat(sprintf("✓ Cat groups (%d): %s\n", 
              length(cat_groups), 
              paste(cat_groups, collapse = ", ")))
  cat(sprintf("✓ NoCat groups (%d): %s\n", 
              length(nocat_groups), 
              paste(nocat_groups, collapse = ", ")))
  
  # 3. Define imputation function ####
  impute_missing_values <- function(data, data_name, nocat_groups, cat_groups, impute_cat) {
    
    cat(sprintf("\n  Processing: %s\n", data_name))
    
    # Calculate imputation parameters (Perseus)
    cat("    - Calculating Perseus imputation parameters...\n")
    
    perseus_params <- data %>%
      select(all_of(data_cols)) %>%
      pivot_longer(everything(), names_to = "sample", values_to = "intensity") %>%
      filter(!is.na(intensity)) %>%
      group_by(sample) %>%
      summarise(
        col_mean = mean(intensity),
        col_sd = sd(intensity),
        .groups = "drop"
      ) %>%
      mutate(
        imputation_mean = col_mean - 1.8 * col_sd,
        imputation_sd = 0.3 * col_sd
      )
    
    # Map samples to bioGroup
    sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
    
    # Impute values
    cat("    - Performing imputation...\n")
    
    # Set seed for reproducibility
    set.seed(random_seed)
    
    data_imputed <- data %>%
      pivot_longer(
        cols = all_of(data_cols),
        names_to = "sample",
        values_to = "intensity",
        values_drop_na = FALSE
      ) %>%
      mutate(group = sample_to_group[sample]) %>%
      group_by(Gene, group) %>%
      mutate(
        n_valid = sum(!is.na(intensity)),
        mean_valid_within_group = mean(intensity, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      left_join(perseus_params, by = "sample") %>%
      rowwise() %>%
      mutate(
        intensity_imputed = case_when(
          # NoCat: n_valid=2 use mean; <2 use Perseus sampling
          group %in% nocat_groups & n_valid == 2 & is.na(intensity) ~ mean_valid_within_group,
          group %in% nocat_groups & n_valid < 2 & is.na(intensity) ~ rnorm(1, mean = imputation_mean, sd = imputation_sd),
          # Cat: optional mean imputation when n_valid=2
          impute_cat & group %in% cat_groups & n_valid == 2 & is.na(intensity) ~ mean_valid_within_group,
          # Otherwise keep original
          TRUE ~ intensity
        )
      ) %>%
      ungroup() %>%
      select(Gene, sample, intensity_imputed, all_of(anno_cols)) %>%
      distinct() %>%
      pivot_wider(
        names_from = sample,
        values_from = intensity_imputed
      ) %>%
      select(Gene, all_of(data_cols), all_of(anno_cols))
    
    # Imputation stats
    n_imputed <- sum(is.na(data[, data_cols])) - sum(is.na(data_imputed[, data_cols]))
    cat(sprintf("    ✓ Imputed %d missing values\n", n_imputed))
    
    return(data_imputed)
  }
  
  # 4. Impute across all datasets ####
  cat("\n[3] Run imputation...\n")
  
  imputed_data_list <- list()
  
  for (i in seq_along(standardized_data_list)) {
    data_name <- names(standardized_data_list)[i]
    data <- standardized_data_list[[i]]
    
    # Perform imputation
    data_imputed <- impute_missing_values(data, data_name, nocat_groups, cat_groups, impute_cat_mean)
    
    # Save to list
    imputed_name <- paste0(data_name, "_Imputed")
    imputed_data_list[[imputed_name]] <- data_imputed
    
    # Save CSV
    csv_file <- file.path(dir_config$output, paste0("Module05_", imputed_name, ".csv"))
    write.csv(data_imputed, csv_file, row.names = FALSE)
    cat(sprintf("    ✓ Saved: %s\n", basename(csv_file)))
  }
  
  cat(sprintf("\n✓ Processed %d datasets\n", length(imputed_data_list)))
  
  # 5. Generate comparison boxplots ####
  cat("\n[4] Generate before/after boxplots...\n")
  
  # Assign colors per bioGroup (aligned with Module 4)
  sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
  col_groups <- sample_to_group[data_cols]
  unique_groups <- unique(col_groups)
  n_groups <- length(unique_groups)
  
  if (n_groups <= 12) {
    group_colors <- scales::hue_pal()(n_groups)
  } else {
    group_colors <- rainbow(n_groups)
  }
  
  names(group_colors) <- unique_groups
  box_colors <- group_colors[col_groups]
  
  pdf_file <- file.path(dir_config$output, "Module05_Imputation_comparison_boxplot.pdf")
  pdf(pdf_file, width = 12, height = 6)
  
  for (i in seq_along(standardized_data_list)) {
    original_name <- names(standardized_data_list)[i]
    imputed_name <- names(imputed_data_list)[i]
    
    original_data <- standardized_data_list[[i]]
    imputed_data <- imputed_data_list[[i]]
    
    # Two-panel layout
    par(mfrow = c(1, 2))
    
    # Before imputation
    boxplot(original_data[, data_cols], 
            cex.axis = 0.4, 
            las = 2, 
            main = paste0(original_name, "\n(Before Imputation)"),
            col = box_colors,
            border = "black")
    
    # After imputation
    boxplot(imputed_data[, data_cols], 
            cex.axis = 0.4, 
            las = 2, 
            main = paste0(imputed_name, "\n(After Imputation)"),
            col = box_colors,
            border = "black")
    
    par(mfrow = c(1, 1))
  }
  
  dev.off()
  cat(sprintf("✓ Saved: %s\n", basename(pdf_file)))
  
  # 6. Validate annotation columns ####
  cat("\n[5] Validate annotation columns...\n")
  
  for (i in seq_along(imputed_data_list)) {
    original_data <- standardized_data_list[[i]]
    imputed_data <- imputed_data_list[[i]]
    
    for (anno_col in anno_cols) {
      if (!all(original_data[[anno_col]] == imputed_data[[anno_col]], na.rm = TRUE)) {
        warning(sprintf("⚠ Column %s may have changed in %s", 
                       names(imputed_data_list)[i], anno_col))
      }
    }
  }
  
  cat("✓ Annotation columns validated\n")
  
  # 7. Summarize imputation parameters ####
  cat("\n[6] Summarize imputation parameters...\n")
  
  imputation_params <- list(
    nocat_groups = nocat_groups,
    cat_groups = cat_groups,
    impute_cat_mean = impute_cat_mean,
    random_seed = random_seed,
    data_cols = data_cols,
    anno_cols = anno_cols
  )
  
  # 8. Return result ####
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("Module 05 complete\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("\nGenerated datasets:\n")
  for (name in names(imputed_data_list)) {
    cat(sprintf("  - %s\n", name))
  }
  cat("\n")
  
  return(list(
    imputed_data_list = imputed_data_list,
    imputation_params = imputation_params
  ))
}

