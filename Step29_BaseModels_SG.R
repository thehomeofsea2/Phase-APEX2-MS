# Academy of Military Sciences SGs Sample Analysis
#1. Load Packages #### R4.1
library(readr)
library(tidyverse)
library(readxl)
library(purrr)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(gridExtra)
library(pheatmap)
library(pROC)
library(openxlsx)
library(preprocessCore)
library(limma)
library(colourpicker)
library(RColorBrewer)

#2. Working Directory and Data Loading ####
setwd("D:/Bioinfomatics/MS/CodeNorm_Nuc_Batch1")
#source("../202501017 3DAUC models base.R")
## Replace with any environment that can complete the main pipeline to FC/Scatter plot
load("D:/Bioinfomatics/MS/CodeNorm_Nuc_Batch1/Module14_workspace.RData") 

str(ForStep16)
names(ForStep16)

## Group Structure
group_info <- list(
  E7A2B4 = list(
    samples = c("E7A2B4_Light_LFQ_1", "E7A2B4_Light_LFQ_2", "E7A2B4_Light_LFQ_3"),
    logFC  = "E7A2B4_vs_K19_logFC",
    logFC_FDR="E7A2B4_vs_K19_adj.P.Val"
  ),
  E7C2 = list(
    samples = c("E7C2_Light_LFQ_1", "E7C2_Light_LFQ_2", "E7C2_Light_LFQ_3"),
    logFC  = "E7C2_vs_K19_logFC",
    logFC_FDR="E7C2_vs_K19_adj.P.Val"
  ),
  C2 = list(
    samples = c("C2_H2O2_LFQ_1", "C2_H2O2_LFQ_2", "C2_H2O2_LFQ_3"),
    logFC  = "C2_vs_J37_logFC",
    logFC_FDR="C2_vs_J37_adj.P.Val"
  )
)

for (i in seq_along(ForStep16)) {
  ForStep16[[i]]=ForStep16[[i]] %>% mutate(MultiBait_Localization=Sub_HPA_Localization)
}

str(ForStep16)

unique(ForStep16[[1]]$MultiBait_Localization)
myTP_vector=c("Nucleolus","Other&Nucleolus","Nuclear&Nucleolus","Nuclear&Cytosol&Nucleolus","Cytosol&Nucleolus")
myTP_vector %in% unique(ForStep16[[1]]$MultiBait_Localization)
ForStep19=ForStep16

# Define myTP_vector2 for Fig9 TPR curve - includes all categories containing "TP"
myTP_vector2 <- unique(ForStep19[[1]]$MultiBait_Localization)[
  grepl("Nucleolus", unique(ForStep19[[1]]$MultiBait_Localization), ignore.case = FALSE)
]
cat("myTP_vector2 (for TPR calculation):", paste(myTP_vector2, collapse = ", "), "\n")
myTP_vector2==myTP_vector
str(ForStep19)

# Spatial Control Group Interface: Configure each control group with column name matching rules and output names
spatial_control_groups <- list(
  J37 = list(
    pattern = "J37_H2O2_LFQ",          # Column name matching pattern (regex) for finding corresponding LFQ columns
    display_name = "J37 (H2O2-NLS)"    # Name in charts and outputs
  ),
  K19 = list(
    pattern = "K19_Light_LFQ",
    display_name = "K19 (Light-NLS)"
  )
)

unique(ForStep19[[1]]$MultiBait_Localization)
MULTIBAIT_COLORS <- c(
  "Nucleolus" = "#DB6968",           # Red - bottom layer
  "Nuclear&Nucleolus" = "#B379B4",             # Purple - second layer
  "Nuclear&Cytosol&Nucleolus" = "#B379B4",             # Purple - same as SGs&Other
  "Other&Nucleolus" = "#4D97CD",           # Blue - third layer
  "Cytosol&Nucleolus" = "#EFCA72",   # Yellow - top layer
  "Other" = "grey89",
  "Mitochondrion" = "grey89",
  "Nuclear" = "grey89",
  "Nuclear_Cytosol" = "grey89",
  "Cytosol" = "grey89"
)

# MULTIBAIT_LEVEL_ORDER: Controls ggplot stacking order (from first element in vector stacked at bottom), facilitates alignment with schematic diagram
MULTIBAIT_LEVEL_ORDER <- c(
  "Nucleolus",
  "Nuclear&Nucleolus",
  "Nuclear&Cytosol&Nucleolus",
  "Other&Nucleolus",
  "Cytosol&Nucleolus",
  "Other",
  "Mitochondrion",
  "Nuclear",
  "Nuclear_Cytosol",
  "Cytosol"
)

# MULTIBAIT_SG_CATEGORIES: Specifies which MultiBait categories to highlight (label) in plots, label text filtered by this set
MULTIBAIT_SG_CATEGORIES <- c(
  "Nucleolus",
  "Nuclear&Nucleolus",
  "Nuclear&Cytosol&Nucleolus",
  "Other&Nucleolus",
  "Cytosol&Nucleolus"
)

# Quantitative Column Non-NA Count Threshold Interface
#    - PARAM_MIN_LFQ_NONNA: Phase3 main pipeline and MultiBait version require at least how many non-NA values
#      in the corresponding experimental group's LFQ columns before entering process_replicates() (default 1 is sufficient).
#    - PARAM_MIN_LFQ_NONNA_TPR: Fig9 TPR curve candidate sorting specific threshold, controls minimum non-NA count
#      before TP/FP calculation (default 2, ensures curve based on stable replicates).
PARAM_MIN_LFQ_NONNA <- 2
PARAM_MIN_LFQ_NONNA_TPR <- 2



#Step29 Multi-Model Complex Requirements #####

# ===================================================================
#
###    Phase 1: Data Calculation and Result Export #####
#
# ===================================================================

# -------------------------------------------------------------------
# 1. Setup and Configuration
# -------------------------------------------------------------------

# Load necessary libraries
# Ensure packages are installed: install.packages(c("dplyr", "tidyr", "pROC", "STRINGdb", "igraph", "future", "furrr", "progressr"))
library(dplyr)
library(tidyr)
library(pROC)
library(STRINGdb)
library(igraph)
library(future)
library(furrr)
library(progressr)

# Define a functional error handler as per your preference
my_error_handler <- function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
}

# Setup parallel processing
# Changed to sequential processing for stability due to repeated parallel processing failures
plan(sequential)  # Use sequential processing to avoid parallel environment issues
cat("Switched to sequential processing mode to ensure analysis stability\n")
cat("Processing time will increase, but success rate will significantly improve\n")

# --- Key Parameters ---
# This parameter can be adjusted for the local PPI calculation
local_ppi_window_size <- 10

# Create output directories
output_dir <- "Step29_Analysis_Phase1_DataExport"
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "Cytoscape_Files"), showWarnings = FALSE)

# -------------------------------------------------------------------
# 2. Input Data Preparation
# -------------------------------------------------------------------
# This section assumes 'ForStep19' and 'group_info' are loaded in your environment.

# -------------------------------------------------------------------
# 3. PPI Data Pre-processing (Local STRING Database Version) ####
# -------------------------------------------------------------------
cat("Starting PPI data preprocessing (using local STRING files)...\n")

# === Configuration: Local STRING File Paths ===
resolve_stringdb_dir <- function(base_dir = getwd()) {
  candidate <- file.path(base_dir, "Reference", "StringDb")
  if (!dir.exists(candidate)) {
    stop(paste0("Error: Reference/StringDb directory not found in working directory (", base_dir, ")."))
  }
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

LOCAL_STRING_DIR <- resolve_stringdb_dir()
LOCAL_STRING_ALIAS_FILE <- file.path(LOCAL_STRING_DIR, "9606.protein.aliases.v12.0.txt")
LOCAL_STRING_LINKS_FILE <- file.path(LOCAL_STRING_DIR, "9606.protein.links.detailed.v12.0.txt")
SCORE_THRESHOLD <- 400  # Consistent with original online version

# Verify file paths
cat(paste("STRING file directory:", LOCAL_STRING_DIR, "\n"))
cat(paste("Alias file:", LOCAL_STRING_ALIAS_FILE, "\n"))
cat(paste("Interaction file:", LOCAL_STRING_LINKS_FILE, "\n"))

# === Step 1: Load STRING Alias Mapping Table (load only once) ====
cat("Loading STRING alias mapping table...\n")
if (!file.exists(LOCAL_STRING_ALIAS_FILE)) {
  stop(paste("Error: Mapping file not found", LOCAL_STRING_ALIAS_FILE))
}

string_aliases <- read_tsv(LOCAL_STRING_ALIAS_FILE, 
                           col_names = c("STRING_id", "alias", "source"),
                           col_types = "ccc",
                           skip = 1,
                           show_col_types = FALSE)

# Filter aliases from gene name sources (keep multiple sources to improve mapping rate)
string_gene_map <- string_aliases %>%
  filter(source %in% c("BioMart_HUGO", "Ensembl_gene_name", "BLAST_UniProt_GN", 
                       "Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)")) %>%
  select(STRING_id, gene = alias) %>%
  distinct()

cat(paste("  ✅ Loaded", nrow(string_gene_map), "gene mapping records\n"))

# === Step 2: Load STRING Interaction Network (load only once) ====
cat("Loading STRING interaction network...\n")
if (!file.exists(LOCAL_STRING_LINKS_FILE)) {
  stop(paste("Error: Interaction file not found", LOCAL_STRING_LINKS_FILE))
}

# STRING files use space separation, need to specify explicitly
string_links <- read.table(LOCAL_STRING_LINKS_FILE, 
                           header = TRUE,           # First row is column names
                           sep = " ",               # Space separated
                           stringsAsFactors = FALSE,
                           comment.char = "",       # Don't process comments
                           quote = "")              # Don't process quotes

# Display actual column names for debugging
cat(paste("  File actual column names:", paste(colnames(string_links)[1:5], "..."), "\n"))
cat(paste("  Total", ncol(string_links), "columns,", format(nrow(string_links), big.mark=","), "rows\n"))

# Automatically identify column names (STRING database may use different column name formats)
col_names <- colnames(string_links)

# Identify first protein column (usually column 1)
protein1_col <- col_names[1]
# Identify second protein column (usually column 2)  
protein2_col <- col_names[2]
# Identify combined_score column (may be called combined_score or other names)
score_cols <- col_names[grepl("combined", col_names, ignore.case = TRUE)]
if (length(score_cols) == 0) {
  # If no combined column, try to find score column
  score_cols <- col_names[grepl("score", col_names, ignore.case = TRUE)]
  if (length(score_cols) == 0) {
    stop("Error: Score column not found in file")
  }
}
score_col <- score_cols[1]

cat(paste("  Using columns: [", protein1_col, "] - [", protein2_col, "] - [", score_col, "]\n"))

# Extract and rename columns
string_links <- string_links %>%
  select(protein1 = all_of(protein1_col), 
         protein2 = all_of(protein2_col), 
         combined_score = all_of(score_col)) %>%
  mutate(combined_score = as.numeric(combined_score))

# Apply score threshold filtering
string_links_filtered <- string_links %>%
  filter(combined_score >= SCORE_THRESHOLD)

cat(paste("  ✅ Loaded", format(nrow(string_links_filtered), big.mark=","), 
          "interaction records (score ≥", SCORE_THRESHOLD, ")\n"))

# === Step 3: Process PPI for each data source ====
ppi_data <- list(dfs = ForStep19, interactions = list(), ppi_global = list())

for(df_name in names(ForStep19)) {
  cat(paste("Processing PPI information for data source:", df_name, "\n"))
  
  tryCatch({
    # Get protein list
    protein_list <- unique(ForStep19[[df_name]]$Gene)
    cat(paste("  - Data source contains", length(protein_list), "proteins\n"))
    
    # Map gene names to STRING ID
    mapped_proteins <- string_gene_map %>%
      filter(gene %in% protein_list) %>%
      group_by(gene) %>%
      slice(1) %>%  # Keep only the first STRING ID for each gene
      ungroup()
    
    mapping_rate <- round(nrow(mapped_proteins) / length(protein_list) * 100, 1)
    cat(paste("  - Successfully mapped", nrow(mapped_proteins), "/", length(protein_list), 
              "proteins (", mapping_rate, "%)\n"))
    
    if (nrow(mapped_proteins) == 0) {
      cat("  ⚠️  Warning: No proteins successfully mapped to STRING ID, skipping this data source\n")
      ppi_data$interactions[[df_name]] <- data.frame(
        protein1 = character(0), 
        protein2 = character(0), 
        combined_score = numeric(0)
      )
      ppi_data$ppi_global[[df_name]] <- data.frame(
        Gene = ForStep19[[df_name]]$Gene, 
        ppi_score_global = 0
      )
      next
    }
    
    # Filter interactions (keep only mapped proteins)
    valid_ids <- mapped_proteins$STRING_id
    interactions_filtered <- string_links_filtered %>%
      filter(protein1 %in% valid_ids & protein2 %in% valid_ids)
    
    # Convert STRING ID back to gene names and scale scores to 0-1
    ppi_data$interactions[[df_name]] <- interactions_filtered %>%
      mutate(combined_score = combined_score / 1000) %>%  # Scale to 0-1
      left_join(mapped_proteins %>% select(STRING_id, gene), 
                by = c("protein1" = "STRING_id"), 
                relationship = "many-to-many") %>%
      rename(protein1_gene = gene) %>%
      left_join(mapped_proteins %>% select(STRING_id, gene), 
                by = c("protein2" = "STRING_id"), 
                relationship = "many-to-many") %>%
      rename(protein2_gene = gene) %>%
      filter(!is.na(protein1_gene) & !is.na(protein2_gene)) %>%
      select(protein1 = protein1_gene, protein2 = protein2_gene, combined_score) %>%
      distinct()
    
    cat(paste("  - Obtained", format(nrow(ppi_data$interactions[[df_name]]), big.mark=","), 
              "interaction relationships\n"))
    
    # Calculate global PPI score (total interaction strength for each protein)
    scores1 <- ppi_data$interactions[[df_name]] %>% 
      group_by(protein1) %>% 
      summarise(total_score = sum(combined_score)) %>% 
      rename(Gene = protein1)
    
    scores2 <- ppi_data$interactions[[df_name]] %>% 
      group_by(protein2) %>% 
      summarise(total_score = sum(combined_score)) %>% 
      rename(Gene = protein2)
    
    ppi_data$ppi_global[[df_name]] <- bind_rows(scores1, scores2) %>%
      group_by(Gene) %>%
      summarise(ppi_score_global = sum(total_score))
    
    cat(paste("  - Calculated global PPI scores for", nrow(ppi_data$ppi_global[[df_name]]), 
              "proteins\n"))
    
  }, error = my_error_handler)
}

cat("✅ PPI data preprocessing completed (local version).\n\n")


# -------------------------------------------------------------------
# 4. Core Worker Function for Parallel Execution
# -------------------------------------------------------------------

process_group <- function(df_name, group, ppi_data) {
  
  # Load libraries within each worker for safety
  library(dplyr)
  library(pROC)
  
  cat(paste("Starting processing", df_name, "-", group, "...\n"))
  
  # --- Unpack data for the current task ---
  df <- ppi_data$dfs[[df_name]]
  interactions_df <- ppi_data$interactions[[df_name]]
  ppi_feature_global <- ppi_data$ppi_global[[df_name]]
  
  # Safety check
  if(is.null(df)) {
    cat("❌ Data frame is empty\n")
    return(NULL)
  }
  
  if(is.null(interactions_df) || nrow(interactions_df) == 0) {
    cat("⚠️  PPI data is empty, creating empty PPI features\n")
    interactions_df <- data.frame(protein1 = character(0), protein2 = character(0), combined_score = numeric(0))
    ppi_feature_global <- data.frame(Gene = df$Gene, ppi_score_global = 0)
  }
  
  # --- Data Preparation ---
  samples <- group_info[[group]]$samples
  logFC_col <- group_info[[group]]$logFC
  
  df_processed <- df %>%
    mutate(mean_group = rowMeans(select(., all_of(samples)), na.rm = TRUE)) %>%
    mutate(mean_group_medianNorm = mean_group - median(mean_group, na.rm = TRUE)) %>%
    mutate(
      # Define ground truth based on MultiBait_Localization
      plot_group = if_else(
        MultiBait_Localization %in% myTP_vector, 
        "TP", "Other"
      ),
      # Set factor levels to ensure consistent modeling direction
      plot_group = factor(plot_group, levels = c("Other", "TP"))
    )
  
  # Filter for proteins with valid data points for modeling
  # Corrected the !is.na call to be more robust
  plot_data <- df_processed %>% 
    filter(!is.na(mean_group_medianNorm) & !is.na(.data[[logFC_col]]) & !is.na(plot_group))
  
  # If no data remains after filtering, exit early
  if(nrow(plot_data) == 0) return(NULL)
  
  # --- Local PPI Feature Calculation (Simplified) ---
  sorted_genes <- plot_data %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  
  if(nrow(interactions_df) > 0) {
    local_ppi_scores <- sapply(seq_along(sorted_genes), function(i) {
      current_gene <- sorted_genes[i]
      start_index <- max(1, i - local_ppi_window_size)
      end_index <- min(length(sorted_genes), i + local_ppi_window_size)
      neighbors <- sorted_genes[start_index:end_index]
      
      # Sum scores of interactions with neighbors in the window
      local_score <- interactions_df %>% 
        filter((protein1 == current_gene & protein2 %in% neighbors) | (protein2 == current_gene & protein1 %in% neighbors)) %>% 
        pull(combined_score) %>% 
        sum()
      return(local_score)
    })
  } else {
    # If no PPI data, set to 0
    local_ppi_scores <- rep(0, length(sorted_genes))
  }
  
  ppi_feature_local <- data.frame(Gene = sorted_genes, ppi_score_local = local_ppi_scores)
  
  # --- Merge all features and perform log-transformation ---
  plot_data_full <- plot_data %>%
    left_join(ppi_feature_global, by = "Gene") %>% 
    left_join(ppi_feature_local, by = "Gene") %>%
    mutate(
      # Convert NA PPI scores to 0 (for proteins with no known interactions)
      ppi_score_global = ifelse(is.na(ppi_score_global), 0, ppi_score_global),
      ppi_score_local = ifelse(is.na(ppi_score_local), 0, ppi_score_local),
      # Log1p transform to handle skewness, as discussed
      log_ppi_score_global = log1p(ppi_score_global),
      log_ppi_score_local = log1p(ppi_score_local)
    )
  
  # --- Model Building, Evaluation, and Protein Scoring ---
  all_aucs <- list()
  model_predictions <- data.frame(protein_id = plot_data_full$Gene)
  all_proteins_scored_list <- list()
  
  # Ensure there are at least two groups to compare for modeling
  if (n_distinct(plot_data_full$plot_group) < 2 || any(table(plot_data_full$plot_group) < 5)) {
    warning(paste("Skipping modeling for", df_name, "-", group, "due to insufficient data in one or both classes."))
    return(NULL)
  }
  
  tryCatch({
    # Reload necessary libraries in parallel environment
    library(dplyr)
    library(pROC)
    
    # Define basic models to be tested (simplified for debugging)
    models_to_run <- list(
      "1D_Abundance" = list(type="1D", predictor="mean_group_medianNorm"),
      "1D_logFC" = list(type="1D", predictor=logFC_col),
      "1D_PPI_Global" = list(type="1D", predictor="log_ppi_score_global"),
      #"1D_PPI_Local" = list(type="1D", predictor="log_ppi_score_local"),
      "2D_Abundance_logFC" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col))
      #"2D_Abundance_PPI_Global" = list(type="GLM", formula="plot_group ~ mean_group_medianNorm + log_ppi_score_global"),
      #"2D_Abundance_PPI_Local" = list(type="GLM", formula="plot_group ~ mean_group_medianNorm + log_ppi_score_local"),
      #"2D_logFC_PPI_Global" = list(type="GLM", formula=paste("plot_group ~", logFC_col, "+ log_ppi_score_global")),
      #"2D_logFC_PPI_Local" = list(type="GLM", formula=paste("plot_group ~", logFC_col, "+ log_ppi_score_local")),
      #"3D_with_Global_PPI" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col, "+ log_ppi_score_global")),
      #"3D_with_Local_PPI" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col, "+ log_ppi_score_local"))
    )
    
    cat(paste("Will test", length(models_to_run), "base models\n"))
    
    for(model_name in names(models_to_run)) {
      cat(paste("  Training model:", model_name, "\n"))
      model_info <- models_to_run[[model_name]]
      predictor_values <- NULL
      
      tryCatch({
        # Train model and get predictor values
        if(model_info$type == "1D") {
          model <- glm(as.formula(paste("plot_group ~", model_info$predictor)), data = plot_data_full, family = "binomial")
          predictor_values <- plot_data_full[[model_info$predictor]]
          model_predictions[[paste0("prob_", model_name)]] <- predict(model, type = "response")
          cat(paste("    ✅ 1D model training successful\n"))
        } else { # GLM for 2D/3D
          model <- glm(as.formula(model_info$formula), data = plot_data_full, family = "binomial")
          predictor_values <- predict(model, type = "response")
          model_predictions[[paste0("prob_", model_name)]] <- predictor_values
          cat(paste("    ✅ GLM model training successful\n"))
        }
      }, error = function(e) {
        cat(paste("    ❌ Model training failed:", e$message, "\n"))
        next
      })
      
      # Perform ROC analysis with specified direction
      roc_obj <- roc(response = plot_data_full$plot_group, predictor = predictor_values, 
                     levels = c("Other", "TP"), quiet = TRUE, direction = "<")
      all_aucs[[model_name]] <- as.numeric(pROC::auc(roc_obj))
      
      # Get best threshold using Youden index
      best_coords <- coords(roc_obj, "best", best.method = "youden", ret = "threshold")
      threshold_value <- best_coords$threshold
      
      # Score all proteins based on the model and create the detailed output table
      # Preserving all requested columns
      scored_proteins_df <- plot_data_full %>%
        mutate(
          model_score = if(model_info$type == "1D") .data[[model_info$predictor]] else predict(model, type = "response"),
          Model_best_threshold = threshold_value,
          PassThreshold = if_else(model_score > threshold_value, "Pass", "Fail"),
          isCandidate = if_else(plot_group == "Other" & PassThreshold == "Pass", "TP_candidate", NA_character_)
        ) %>%
        select(
          protein_id = Gene, 
          group = plot_group,
          contains("Localization"),
          abundance = mean_group_medianNorm, 
          logFC = all_of(logFC_col),
          ppi_score_global, ppi_score_local,
          log_ppi_score_global, log_ppi_score_local,
          model_score, 
          Model_best_threshold,
          PassThreshold, 
          isCandidate
        )
      all_proteins_scored_list[[model_name]] <- scored_proteins_df
    }
  }, error = function(e) { 
    warning(paste("Model evaluation failed for", df_name, "-", group, ":", e$message))
    cat("Detailed error information:\n")
    cat("Error type:", class(e)[1], "\n")
    cat("Error message:", e$message, "\n")
    if(exists("plot_data_full")) {
      cat("Valid data rows:", nrow(plot_data_full), "\n")
      if("plot_group" %in% names(plot_data_full)) {
        cat("Group distribution:\n")
        print(table(plot_data_full$plot_group, useNA = "ifany"))
      }
    }
  })
  
  # --- Package all results for return ---
  return(list(
    df_name = df_name, 
    group_name = group,
    auc_scores = all_aucs, 
    # This node_attributes table contains all calculated features for plotting
    node_attributes = plot_data_full,
    # This table contains the final probability scores from each model
    model_predictions = model_predictions, 
    # This list contains the fully annotated dataframes for each model
    candidates_full_list = all_proteins_scored_list,
    # Pass along interactions for this group
    interactions = interactions_df
  ))
}


# -------------------------------------------------------------------
# 5. Execute Parallel Analysis
# -------------------------------------------------------------------
cat("Starting parallel analysis for all combinations...\n")

# Enable progress bar tracking
handlers(global = TRUE)
handlers("progress")

# Create a grid of all tasks to run
tasks <- expand_grid(df_name = names(ForStep19), group = names(group_info))

with_progress({
  # Setup the progressor
  p <- progressr::progressor(steps = nrow(tasks))
  
  # Use sequential processing instead of parallel processing
  all_results <- vector("list", nrow(tasks))
  for(i in 1:nrow(tasks)) {
    df_name <- tasks$df_name[i]
    group <- tasks$group[i]
    cat(paste("Processing", df_name, "-", group, "...\n"))
    
    # Run the main worker function
    result <- process_group(df_name, group, ppi_data)
    all_results[[i]] <- result
    
    # Update progress
    p(message = paste("Completed", df_name, "-", group))
    cat(paste("✅ Completed", df_name, "-", group, "\n"))
  }
})

cat("\nAll sequential computation tasks completed.\n")


# -------------------------------------------------------------------
# 6. Consolidate Results and Export Files
# -------------------------------------------------------------------
cat("Consolidating results and exporting files...\n")

final_summary_list <- list()
final_candidates_list <- list()

# Get the task list again for matching failed tasks
tasks_check <- expand_grid(df_name = names(ForStep19), group = names(group_info))

for(i in seq_along(all_results)) {
  result <- all_results[[i]]
  
  # Check for failed tasks which return NULL
  if (is.null(result) || !is.list(result)) {
    failed_task <- tasks_check[i, ]
    warning(paste("Task", failed_task$df_name, "-", failed_task$group, "failed or did not return valid results, skipped."))
    next
  }
  
  df_name <- result$df_name
  group <- result$group_name
  
  # --- Consolidate AUC summary ---
  temp_row <- data.frame(DataAnalysisMethod = df_name, ExperimentGroup = group)
  aucs <- result$auc_scores
  if(length(aucs) > 0) { 
    temp_row <- cbind(temp_row, as.data.frame(aucs)) 
  }
  group_counts <- table(result$node_attributes$plot_group)
  temp_row$N_TP <- ifelse("TP" %in% names(group_counts), group_counts["TP"], 0)
  temp_row$N_Other <- ifelse("Other" %in% names(group_counts), group_counts["Other"], 0)
  final_summary_list[[paste(df_name, group)]] <- temp_row
  
  # --- Consolidate scored protein list ---
  candidates_by_model <- result$candidates_full_list
  if (length(candidates_by_model) > 0) {
    for (model_name in names(candidates_by_model)) {
      candidates_df <- candidates_by_model[[model_name]]
      if (nrow(candidates_df) > 0) {
        candidates_df$DataAnalysisMethod <- df_name
        candidates_df$ExperimentGroup <- group
        candidates_df$Model <- model_name
        final_candidates_list[[paste(df_name, group, model_name)]] <- candidates_df
      }
    }
  }
  
  # --- Export files for Cytoscape ---
  cytoscape_path_prefix <- file.path(output_dir, "Cytoscape_Files", paste0(df_name, "_", group))
  logFC_col <- group_info[[group]]$logFC
  
  # Node attributes file
  node_attributes_for_export <- result$node_attributes %>%
    left_join(result$model_predictions, by = c("Gene" = "protein_id")) %>%
    select(
      protein_id = Gene, 
      group = plot_group,
      contains("Localization"),
      abundance = mean_group_medianNorm, 
      logFC = all_of(logFC_col),
      log_ppi_score_global, 
      log_ppi_score_local,
      starts_with("prob_")
    )
  write.csv(node_attributes_for_export, 
            paste0(cytoscape_path_prefix, "_nodes.csv"), 
            row.names = FALSE,quote=FALSE)
  
  # Edge attributes files
  group_genes <- result$node_attributes$Gene
  group_interactions <- result$interactions %>% 
    filter(protein1 %in% group_genes & protein2 %in% group_genes)
  
  if(nrow(group_interactions) > 0){
    # 1. Edges with STRING Score
    edges_string <- group_interactions %>% 
      select(source = protein1, target = protein2, weight = combined_score)
    write.csv(edges_string, paste0(cytoscape_path_prefix, "_edges_STRING.csv"), 
              row.names = FALSE,quote=FALSE)
    
    # 2. Edges with Single Factor weights (Abundance, logFC)
    node_data_for_edges <- node_attributes_for_export %>% select(protein_id, abundance, logFC)
    edges_single_factor <- group_interactions %>% 
      select(source = protein1, target = protein2) %>%
      left_join(node_data_for_edges, by = c("source" = "protein_id")) %>% rename(source_ab = abundance, source_lf = logFC) %>%
      left_join(node_data_for_edges, by = c("target" = "protein_id")) %>% rename(target_ab = abundance, target_lf = logFC) %>%
      mutate(
        weight_Abundance = (source_ab + target_ab) / 2, 
        weight_logFC = (source_lf + target_lf) / 2
      ) %>%
      select(source, target, weight_Abundance, weight_logFC)
    write.csv(edges_single_factor, 
              paste0(cytoscape_path_prefix, "_edges_SingleFactor.csv"), 
              row.names = FALSE,quote=FALSE)
    
    # 3. Edges with Model Score weights
    node_probs_for_edges <- node_attributes_for_export %>% select(protein_id, starts_with("prob_"))
    edges_model_scores <- group_interactions %>% 
      select(source = protein1, target = protein2) %>% 
      left_join(node_probs_for_edges, by = c("source" = "protein_id"))
    
    # Prepare target df with renamed columns for joining
    target_df <- node_probs_for_edges
    colnames(target_df) <- c("protein_id", paste0(colnames(target_df)[-1], "_target"))
    
    edges_model_scores <- edges_model_scores %>% 
      left_join(target_df, by = c("target" = "protein_id"))
    
    # Calculate edge weights by multiplying node probabilities
    for (model_prob_col in colnames(node_probs_for_edges)[-1]) {
      target_prob_col <- paste0(model_prob_col, "_target")
      new_weight_col_name <- paste0("weight_", sub("prob_", "", model_prob_col))
      edges_model_scores <- edges_model_scores %>%
        mutate(!!new_weight_col_name := .data[[model_prob_col]] * .data[[target_prob_col]])
    }
    
    edges_model_scores <- edges_model_scores %>% select(source, target, starts_with("weight_"))
    write.csv(edges_model_scores, 
              paste0(cytoscape_path_prefix, "_edges_ModelScores.csv"), 
              row.names = FALSE,quote=FALSE)
  }
}

# --- Save final summary files and RData object ---

# Save model comparison summary
final_summary_df <- bind_rows(final_summary_list)
write.csv(final_summary_df, file.path(output_dir, "full_model_comparison_summary.csv"), row.names = FALSE)
cat(paste("\nModel performance summary table saved to:", file.path(output_dir, "full_model_comparison_summary.csv"), "\n"))

# Save detailed scored protein summary
if(length(final_candidates_list) > 0) {
  final_candidates_summary_df <- bind_rows(final_candidates_list) %>%
    select(DataAnalysisMethod, ExperimentGroup, Model, everything())
  write.csv(final_candidates_summary_df, file.path(output_dir, "all_proteins_scored_summary.csv"), row.names = FALSE)
  cat(paste("Protein scoring details summary table saved to:", file.path(output_dir, "all_proteins_scored_summary.csv"), "\n"))
}

# Save the essential R objects for Phase 2 (Plotting)
# Corrected 'analysis_results' to 'all_results'
save(all_results, group_info, ppi_data,
     file = file.path(output_dir, "analysis_phase1_output.RData"))
cat(paste("R objects for Phase 2 plotting saved to:", file.path(output_dir, "analysis_phase1_output.RData"), "\n"))


# -------------------------------------------------------------------
# 7. Shutdown
# -------------------------------------------------------------------

# Gracefully terminate the parallel backend
plan(sequential)
cat("\nParallel session terminated. Script execution completed.\n")



# ===================================================================
#
### Phase 2: Result Visualization (v1.1 - Corrected) ####
#
# ===================================================================

# -------------------------------------------------------------------
# 1. Setup and Configuration
# -------------------------------------------------------------------

# Load necessary libraries
# Ensure packages are installed: install.packages(c("ggplot2", "dplyr", "tidyr", "pROC", "ggrepel", "patchwork", "RColorBrewer", "viridis", "MASS", "FNN"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(MASS)
library(FNN)


# Load the data object from Phase 1
cat("Loading Phase 1 calculation results 'analysis_phase1_output.RData'...\n")
load("Step29_Analysis_Phase1_DataExport/analysis_phase1_output.RData")

# --- Create Output Directories ---
output_dir_plots <- "Step29_Analysis_Phase2_Plots_v1.2"
dir.create(output_dir_plots, showWarnings = FALSE)
plot_subdirs <- c(
  "1_ROC_Curves",
  "2_Youden_Index",
  "3_Decision_Boundaries",
  "4_3D_Like_Scatters",
  "5_PPI_Heatmaps",
  "7_Candidate_PPI_Heatmaps",
  "13_JointProb_Heatmaps",
  "15_JointProb_Heatmaps_Weighted"
)
for (subdir in plot_subdirs) {
  dir.create(file.path(output_dir_plots, subdir), showWarnings = FALSE)
}

# -------------------------------------------------------------------
# 2.Global Plotting Configuration (Modify all plot appearances here) ####
# -------------------------------------------------------------------
cat("Applying global plotting configuration...\n")

# --- General Theme ---
theme_set(theme_bw() + theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size=14),
  plot.subtitle = element_text(hjust = 0.5, size=11),
  legend.position = "bottom"
))

# --- Color Palettes ---
TP_vs_other_colors <- c("TP" = "#e41a1c", "Other" = "#bdbdbd")
# 1. ROC curve color
roc_line_color <- "#e41a1c" 
# 2. Youden Index curve color
youden_line_color <- "#1f78b4"
# 4. 3D-like scatter plot color scheme
scatter_3d_color_palette <- c("#bdbdbd","yellow","#e41a1c")
# 5. PPI heatmap color scheme
heatmap_color_palette <- c("#bdbdbd","#e41a1c")
# 6. Density plot color scheme
density_color_palette <- c("blue", "yellow", "red")

# --- Plot-specific Parameters ---
# 2. Youden Index Y-axis range
youden_y_axis_range <- c(-0.2, 1)
# 5. PPI heatmap interaction score range
heatmap_score_range <- c(0.4, 1.0)
# 6. Density plot clustering number (set to 0 to disable clustering)
density_kmeans_k <- 0


# --- NEW: Data-driven Axis Range Interfaces ---
# Note: Set as vector like c(-5, 5) to manually override range; keep as NULL to use auto-calculated global range.
axis_range_abundance    <- c(-5,10) # e.g. c(-4, 8)
axis_range_logfc          <- c(-3,8) # e.g. c(-3, 6)
axis_range_ppi_global   <- c(0,6.5) # e.g. c(0, 10)
axis_range_ppi_local      <- c(0,4) # e.g. c(0, 5)

# --- Density normalization & scaling interfaces (Scheme A) ---
# Coordinate mode: "rank01" uses [0,1] normalized rank coordinates with fixed bandwidth; "rank" keeps original 1..n rank domain with default bandwidth
density_coord_mode <- "rank01"



# Legacy code interfaces, do not delete to avoid dependency errors, not used in calculations
#-----
#----Protected area, do not modify
density_bandwidth <- c(0.1, 0.1)
density_grid_size <- 100
use_normalized_density <- TRUE
enable_weighted_density <- TRUE
weighted_density_sample_size <- 10000
pm_enable <- TRUE
pm_grid_mode <- "fixed"      
pm_grid_bins_fixed <- 50L    
pm_lower_left_percents <- c(0.25, 0.5, 0.75)
pm_weight_modes <- c("none", "combined_score_tp", "combined_score", "tp")
pm_min_quantile <- 0.01
pm_candidates_min_quantile <- 0.01
pm_weight_power <- 2
pm_max_quantile <- 0.7
pm_candidates_max_quantile <- 0.7
weighted_density_max_quantile <- 0.7
weighted_density_candidates_max_quantile <- 0.7
weighted_density_mode <- "combined_score_tp"
weighted_density_suffix <- "_Weighted"
density_max_quantile <- 0.7
density_candidates_max_quantile <- 0.7
#-----
#-----
#-----



#-------------------------------------------------------------------------------------
#----Configurable area---------------------------------------------------------Configurable area
tp_weight_both <- 2
tp_weight_single <- 1.0
tp_weight_none <- 0.25

# Override coordinate mode based on master switch
density_coord_mode <- if (use_normalized_density) "rank01" else "rank"


# JP: Joint rank probability heatmap (unweighted)
jp_enable <- TRUE
jp_pair_mode <- "directed"   # "directed"|"undirected" (currently using directed)
jp_grid_mode <- "fixed"
jp_grid_bins_fixed <- 50L
# Optional: Lower quantile (recommended 0~0.1), 0 means no background enhancement
jp_min_quantile <- 0.01
jp_candidates_min_quantile <- 0.01
jp_max_quantile <- 0.7
jp_candidates_max_quantile <- 0.7

# JP (weighted)
jp_weight_enable <- TRUE
jp_weight_modes <- c("combined_score_tp", "combined_score", "tp")
jp_weight_power <- 2
jpw_min_quantile <- 0.01
jpw_candidates_min_quantile <- 0.01
jpw_max_quantile <- 0.7
jpw_candidates_max_quantile <- 0.7


# -------------------------------------------------------------------
# 3. Pre-computation for Unified Plot Scales
# -------------------------------------------------------------------
# (This section remains the same as v1.1)
cat("Pre-computing unified axis and color ranges for all plots...\n")
plot_scales_config <- list()
safe_range <- function(x) { range(x[is.finite(x)], na.rm = TRUE) }

all_node_data <- lapply(all_results, function(res) {
  if (is.null(res)) return(NULL)
  logfc_col_name <- group_info[[res$group_name]]$logFC
  res$node_attributes %>%
    dplyr::select(
      mean_group_medianNorm, logFC = all_of(logfc_col_name),
      log_ppi_score_global, log_ppi_score_local
    )
}) %>% bind_rows()

plot_scales_config$abundance_range <- safe_range(all_node_data$mean_group_medianNorm)
plot_scales_config$logfc_range <- safe_range(all_node_data$logFC)
plot_scales_config$ppi_global_range <- safe_range(all_node_data$log_ppi_score_global)
plot_scales_config$ppi_local_range <- safe_range(all_node_data$log_ppi_score_local)

cat("Pre-computing global maximum density for density plots...\n")
global_max_density <- 0
for (res in all_results) {
  if (is.null(res)) next
  node_attr <- res$node_attributes
  group_col <- if ("plot_group" %in% colnames(node_attr)) "plot_group" else "group_Localization"
  group_map <- setNames(node_attr[[group_col]], node_attr$Gene)
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  interactions_df <- res$interactions
  n_genes <- length(sorted_genes)
  gene_map <- setNames(1:n_genes, sorted_genes)
  
  # Filter and symmetrize interaction relationships (consistent with plotting function)
  filtered_edges <- interactions_df %>%
    dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  
  # Symmetrization
  symmetric_edges <- bind_rows(
    filtered_edges,
    filtered_edges %>% 
      rename(protein1 = protein2, protein2 = protein1)
  ) %>%
  distinct()
  
  # Use same column names as plotting function for consistency
  coords_df_precalc <- symmetric_edges %>%
    mutate(row = gene_map[protein1], col = gene_map[protein2])
  coords <- coords_df_precalc %>% dplyr::select(row, col) %>% as.matrix()
  
  if (nrow(coords) > 10) {
    tryCatch({
      # Note: Coordinate order must match plotting function (coords[,2], coords[,1])
      # Grid size #####
      density_data <- MASS::kde2d(coords[,2], coords[,1], n = density_grid_size, lims = c(1, n_genes, 1, n_genes))
      global_max_density <- max(global_max_density, max(density_data$z, na.rm = TRUE))
    }, error = function(e) { warning(paste("KDE2D failed for", res$df_name, res$group_name, ":", e$message)) })
  }
}
plot_scales_config$max_density <- global_max_density
cat(paste("All Fig6 density plots will use unified maximum density value:", round(global_max_density, 8), "\n"))

# --- Override with normalized/quantile-based max if requested (Scheme A) ---
if (identical(tolower(density_coord_mode), "rank01") || density_max_quantile < 1 || density_max_quantile > 1) {
  density_max_vals <- c()
  for (res in all_results) {
    if (is.null(res)) next
    node_attr <- res$node_attributes
    sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    node_attr <- res$node_attributes
    group_col <- if ("plot_group" %in% colnames(node_attr)) "plot_group" else "group_Localization"
    group_map <- setNames(node_attr[[group_col]], node_attr$Gene)
    interactions_df <- res$interactions
    n_genes <- length(sorted_genes)
    if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered_edges <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    symmetric_edges <- bind_rows(filtered_edges, filtered_edges %>% rename(protein1 = protein2, protein2 = protein1)) %>% distinct()
    coords_df_precalc <- symmetric_edges %>% mutate(row = gene_map[protein1], col = gene_map[protein2])
    coords <- coords_df_precalc %>% dplyr::select(row, col) %>% as.matrix()
    if (nrow(coords) > 10) {
      col_norm <- coords[,2] / n_genes
      row_norm <- coords[,1] / n_genes
      density_data <- MASS::kde2d(col_norm, row_norm, n = density_grid_size, lims = c(0, 1, 0, 1), h = density_bandwidth)
      zmax <- suppressWarnings(max(density_data$z, na.rm = TRUE))
      if (is.finite(zmax)) density_max_vals <- c(density_max_vals, zmax)
    }
  }
  if (length(density_max_vals) > 0) {
    if (is.finite(density_max_quantile) && density_max_quantile > 1) {
      plot_scales_config$max_density <- as.numeric(max(density_max_vals, na.rm = TRUE) * density_max_quantile)
    } else {
      plot_scales_config$max_density <- as.numeric(stats::quantile(density_max_vals, probs = max(0, min(1, density_max_quantile)), na.rm = TRUE))
    }
    cat(paste("Fig6 global maximum density changed to normalized/quantile:", round(plot_scales_config$max_density, 8),
              paste0(" [q=", density_max_quantile, "]\n"), sep = ""))
  }
}

# Pre-compute global maximum density for candidate protein density plots (Fig8)
cat("Pre-computing global maximum density for candidate protein density plots (Fig8)...\n")
global_max_density_candidates <- 0

for (res in all_results) {
  if (is.null(res)) next
  
  interactions_df <- res$interactions
  
  # Iterate through all models for this result
  for (model_name in names(res$candidates_full_list)) {
    tryCatch({
      # Get candidate proteins that passed threshold
      candidates_df <- res$candidates_full_list[[model_name]]
      pass_proteins <- candidates_df %>% 
        filter(PassThreshold == "Pass") %>% 
        pull(protein_id)
      
      if (length(pass_proteins) < 10) next
      
      # Sort by abundance
      sorted_candidates <- res$node_attributes %>%
        filter(Gene %in% pass_proteins) %>%
        arrange(desc(mean_group_medianNorm))
      
      sorted_genes <- sorted_candidates$Gene
      n_genes <- length(sorted_genes)
      gene_map <- setNames(1:n_genes, sorted_genes)
      
      # Filter interactions between candidate proteins
      filtered_edges <- interactions_df %>%
        filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
      
      if (nrow(filtered_edges) < 20) next
      
      # Symmetrization
      symmetric_edges <- bind_rows(
        filtered_edges,
        filtered_edges %>% rename(protein1 = protein2, protein2 = protein1)
      ) %>% distinct()
      
      coords_df_precalc <- symmetric_edges %>%
        mutate(row = gene_map[protein1], col = gene_map[protein2])
      coords <- coords_df_precalc %>% dplyr::select(row, col) %>% as.matrix()
      
      # Coordinate validity check
      coords_valid <- coords[coords[,1] >= 1 & coords[,1] <= n_genes & 
                             coords[,2] >= 1 & coords[,2] <= n_genes, , drop = FALSE]
      
      if (nrow(coords_valid) > 10) {
      density_data <- MASS::kde2d(coords_valid[,2], coords_valid[,1], n = density_grid_size, 
                                   lims = c(1, n_genes, 1, n_genes))
        global_max_density_candidates <- max(global_max_density_candidates, 
                                              max(density_data$z, na.rm = TRUE))
      }
    }, error = function(e) { })
  }
}

plot_scales_config$max_density_candidates <- global_max_density_candidates
cat(paste("All Fig8 candidate protein density plots will use unified maximum density value:", round(global_max_density_candidates, 8), "\n\n"))

# ================== Weighted density global maxima (separate) ==================
cat("Pre-computing global maximum density/quantile for weighted density plots (Fig6)...\n")
if (isTRUE(enable_weighted_density)) {
weighted_density_max_vals <- c()
for (res in all_results) {
  if (is.null(res)) next
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  interactions_df <- res$interactions
  n_genes <- length(sorted_genes)
  if (n_genes < 2) next
  gene_map <- setNames(1:n_genes, sorted_genes)

  filtered_interactions <- interactions_df %>%
    dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  symmetric_interactions <- bind_rows(
    filtered_interactions,
    filtered_interactions %>% rename(protein1 = protein2, protein2 = protein1)
  ) %>% distinct()
  if (nrow(symmetric_interactions) < 20) next

  coords_df_precalc <- symmetric_interactions %>%
    mutate(row = gene_map[protein1], col = gene_map[protein2])
  coords <- coords_df_precalc %>% dplyr::select(row, col) %>% as.matrix()
  if (nrow(coords) <= 10) next

  # Get weights (combined_score), ensure non-negative and non-NA
  tp1 <- group_map[symmetric_interactions$protein1] == "TP"
  tp2 <- group_map[symmetric_interactions$protein2] == "TP"
  tp1[is.na(tp1)] <- FALSE
  tp2[is.na(tp2)] <- FALSE
  tp_w <- ifelse(tp1 & tp2, tp_weight_both, ifelse(xor(tp1, tp2), tp_weight_single, tp_weight_none))
  base_w <- symmetric_interactions$combined_score
  base_w[!is.finite(base_w) | is.na(base_w) | base_w < 0] <- 0
  w <- tp_w * base_w
  if (sum(w) <= 0) next

  # Weighted sampling
  set.seed(123)
  s <- min(weighted_density_sample_size, nrow(coords))
  idx <- sample.int(n = nrow(coords), size = s, replace = TRUE, prob = w)
  coords_s <- coords[idx, , drop = FALSE]

  # KDE calculation (by mode)
  if (identical(tolower(density_coord_mode), "rank01")) {
    col_norm <- coords_s[,2] / n_genes
    row_norm <- coords_s[,1] / n_genes
    density_data <- MASS::kde2d(col_norm, row_norm, n = density_grid_size, lims = c(0, 1, 0, 1), h = density_bandwidth)
  } else {
    density_data <- MASS::kde2d(coords_s[,2], coords_s[,1], n = density_grid_size, lims = c(1, n_genes, 1, n_genes))
  }
  zmax <- suppressWarnings(max(density_data$z, na.rm = TRUE))
  if (is.finite(zmax)) weighted_density_max_vals <- c(weighted_density_max_vals, zmax)
}
plot_scales_config$max_density_weighted <- if (length(weighted_density_max_vals) > 0) {
  if (is.finite(weighted_density_max_quantile) && weighted_density_max_quantile > 1) {
    as.numeric(max(weighted_density_max_vals, na.rm = TRUE) * weighted_density_max_quantile)
  } else {
    as.numeric(stats::quantile(weighted_density_max_vals, probs = max(0, min(1, weighted_density_max_quantile)), na.rm = TRUE))
  }
} else { 0 }
cat(paste("Fig6 weighted density plots will use unified maximum density (quantile)",
          paste0("[q=", weighted_density_max_quantile, "]:"),
          round(plot_scales_config$max_density_weighted, 8), "\n"))

cat("Pre-computing global maximum density/quantile for weighted candidate density plots (Fig8)...\n")
weighted_density_candidates_max_vals <- c()
for (res in all_results) {
  if (is.null(res)) next
  interactions_df <- res$interactions
  for (model_name in names(res$candidates_full_list)) {
    candidates_df <- res$candidates_full_list[[model_name]]
    pass_proteins <- tryCatch({
      candidates_df %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id)
    }, error = function(e) character(0))
    if (length(pass_proteins) < 10) next

    sorted_candidates <- res$node_attributes %>%
      dplyr::filter(Gene %in% pass_proteins) %>%
      dplyr::arrange(desc(mean_group_medianNorm))
    sorted_genes <- sorted_candidates$Gene
    n_genes <- length(sorted_genes)
    if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)

    filtered_edges <- interactions_df %>%
      dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    if (nrow(filtered_edges) < 20) next
    symmetric_interactions <- bind_rows(
      filtered_edges,
      filtered_edges %>% rename(protein1 = protein2, protein2 = protein1)
    ) %>% dplyr::distinct()

    coords_df_precalc <- symmetric_interactions %>%
      dplyr::mutate(row = gene_map[protein1], col = gene_map[protein2])
    coords <- coords_df_precalc %>% dplyr::select(row, col) %>% as.matrix()
    valid_idx <- which(coords[,1] >= 1 & coords[,1] <= n_genes & coords[,2] >= 1 & coords[,2] <= n_genes)
    coords_valid <- coords[valid_idx, , drop = FALSE]
    if (nrow(coords_valid) <= 10) next

    tp1 <- group_map[symmetric_interactions$protein1] == "TP"
    tp2 <- group_map[symmetric_interactions$protein2] == "TP"
    tp1[is.na(tp1)] <- FALSE
    tp2[is.na(tp2)] <- FALSE
    tp_w <- ifelse(tp1 & tp2, tp_weight_both, ifelse(xor(tp1, tp2), tp_weight_single, tp_weight_none))
    base_w <- symmetric_interactions$combined_score
    if (length(base_w) != nrow(coords)) base_w <- rep(1, nrow(coords))
    base_w[!is.finite(base_w) | is.na(base_w) | base_w < 0] <- 0
    w <- tp_w * base_w
    w_valid <- w[valid_idx]
    if (sum(w_valid) <= 0) next

    set.seed(123)
    s <- min(weighted_density_sample_size, nrow(coords_valid))
    idx <- sample.int(n = nrow(coords_valid), size = s, replace = TRUE, prob = w_valid)
    coords_s <- coords_valid[idx, , drop = FALSE]

    if (identical(tolower(density_coord_mode), "rank01")) {
      col_norm <- coords_s[,2] / n_genes
      row_norm <- coords_s[,1] / n_genes
      density_data <- MASS::kde2d(col_norm, row_norm, n = density_grid_size, lims = c(0, 1, 0, 1), h = density_bandwidth)
    } else {
      density_data <- MASS::kde2d(coords_s[,2], coords_s[,1], n = density_grid_size, lims = c(1, n_genes, 1, n_genes))
    }
    zmax <- suppressWarnings(max(density_data$z, na.rm = TRUE))
    if (is.finite(zmax)) weighted_density_candidates_max_vals <- c(weighted_density_candidates_max_vals, zmax)
  }
}
plot_scales_config$max_density_candidates_weighted <- if (length(weighted_density_candidates_max_vals) > 0) {
  if (is.finite(weighted_density_candidates_max_quantile) && weighted_density_candidates_max_quantile > 1) {
    as.numeric(max(weighted_density_candidates_max_vals, na.rm = TRUE) * weighted_density_candidates_max_quantile)
  } else {
    as.numeric(stats::quantile(weighted_density_candidates_max_vals, probs = max(0, min(1, weighted_density_candidates_max_quantile)), na.rm = TRUE))
  }
} else { 0 }
cat(paste("Fig8 weighted candidate density plots will use unified maximum density (quantile)",
          paste0("[q=", weighted_density_candidates_max_quantile, "]:"),
          round(plot_scales_config$max_density_candidates_weighted, 8), "\n\n"))
}

# --- Override candidate max with normalized/quantile-based (Scheme A) ---
if (identical(tolower(density_coord_mode), "rank01") || density_candidates_max_quantile < 1 || density_candidates_max_quantile > 1) {
  density_max_vals_candidates <- c()
  for (res in all_results) {
    if (is.null(res)) next
    interactions_df <- res$interactions
    for (model_name in names(res$candidates_full_list)) {
      candidates_df <- res$candidates_full_list[[model_name]]
      pass_proteins <- tryCatch({
        candidates_df %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id)
      }, error = function(e) character(0))
      if (length(pass_proteins) < 10) next
      sorted_candidates <- res$node_attributes %>% dplyr::filter(Gene %in% pass_proteins) %>% dplyr::arrange(desc(mean_group_medianNorm))
      sorted_genes <- sorted_candidates$Gene
      n_genes <- length(sorted_genes)
      if (n_genes < 2) next
      gene_map <- setNames(1:n_genes, sorted_genes)
      filtered_edges <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
      if (nrow(filtered_edges) < 20) next
      symmetric_edges <- bind_rows(filtered_edges, filtered_edges %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
      coords_df_precalc <- symmetric_edges %>% dplyr::mutate(row = gene_map[protein1], col = gene_map[protein2])
      coords <- coords_df_precalc %>% dplyr::select(row, col) %>% as.matrix()
      coords_valid <- coords[coords[,1] >= 1 & coords[,1] <= n_genes & coords[,2] >= 1 & coords[,2] <= n_genes, , drop = FALSE]
      if (nrow(coords_valid) > 10) {
        col_norm <- coords_valid[,2] / n_genes
        row_norm <- coords_valid[,1] / n_genes
        density_data <- MASS::kde2d(col_norm, row_norm, n = density_grid_size, lims = c(0, 1, 0, 1), h = density_bandwidth)
        zmax <- suppressWarnings(max(density_data$z, na.rm = TRUE))
        if (is.finite(zmax)) density_max_vals_candidates <- c(density_max_vals_candidates, zmax)
      }
    }
  }
  if (length(density_max_vals_candidates) > 0) {
    if (is.finite(density_candidates_max_quantile) && density_candidates_max_quantile > 1) {
      plot_scales_config$max_density_candidates <- as.numeric(max(density_max_vals_candidates, na.rm = TRUE) * density_candidates_max_quantile)
    } else {
      plot_scales_config$max_density_candidates <- as.numeric(stats::quantile(density_max_vals_candidates, probs = max(0, min(1, density_candidates_max_quantile)), na.rm = TRUE))
    }
    cat(paste("Fig8 candidate density upper limit changed to normalized/quantile:", round(plot_scales_config$max_density_candidates, 8),
              paste0(" [q=", density_candidates_max_quantile, "]\n\n"), sep = ""))
  }
}


# ================== PM/JP Global Color Upper Limit Pre-computation ==================
cat("Pre-computing PM/JP global color upper limits...\n")

# Utility function: Get TP/Other grouping column name
getGroupColName <- function(df) {
  if ("plot_group" %in% colnames(df)) {
    "plot_group"
  } else if ("group_Localization" %in% colnames(df)) {
    "group_Localization"
  } else {
    NULL
  }
}

# PM weight function (reads global tp_weight_* and combined_score)
getPmWeightsFromMode <- function(symmetric_interactions, node_attr, weight_mode) {
  n <- nrow(symmetric_interactions)
  if (n == 0) return(numeric(0))
  base_w <- symmetric_interactions$combined_score
  base_w[!is.finite(base_w) | is.na(base_w) | base_w < 0] <- 0
  gcol <- getGroupColName(node_attr)
  if (!is.null(gcol)) {
    group_map <- setNames(node_attr[[gcol]], node_attr$Gene)
  } else {
    group_map <- setNames(rep(NA_character_, nrow(node_attr)), node_attr$Gene)
  }
  tp1 <- group_map[symmetric_interactions$protein1] == "TP"
  tp2 <- group_map[symmetric_interactions$protein2] == "TP"
  tp1[is.na(tp1)] <- FALSE; tp2[is.na(tp2)] <- FALSE
  tp_w <- ifelse(tp1 & tp2, tp_weight_both,
                 ifelse(xor(tp1, tp2), tp_weight_single, tp_weight_none))
  if (identical(weight_mode, "none")) {
    w <- rep(1, n)
  } else if (identical(weight_mode, "combined_score")) {
    w <- base_w
  } else if (identical(weight_mode, "tp")) {
    w <- tp_w
  } else { # combined_score_tp
    w <- tp_w * base_w
  }
  w[!is.finite(w) | is.na(w) | w < 0] <- 0
  # Power enhancement
  if (exists("pm_weight_power", inherits = TRUE) && is.finite(pm_weight_power) && pm_weight_power != 1) {
    w <- w ^ pm_weight_power
  }
  w
}

# Map [0,1] coordinates to [1..bins]
.binIndex <- function(vals01, bins) {
  v <- pmin(1, pmax(0, vals01))
  idx <- ceiling(v * bins)
  idx[idx < 1] <- 1L
  idx[idx > bins] <- bins
  as.integer(idx)
}

# Calculate PM maximum single-cell probability (all)
.pmMaxForResult <- function(res, bins, weight_mode) {
  if (is.null(res)) return(NA_real_)
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes); if (n_genes < 2) return(NA_real_)
  gene_map <- setNames(1:n_genes, sorted_genes)
  interactions_df <- res$interactions
  filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) return(NA_real_)
  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
  w <- getPmWeightsFromMode(sym, node_attr, weight_mode)
  if (length(w) == 0 || sum(w) <= 0) return(NA_real_)
  df <- data.frame(xb = xb, yb = yb, w = w)
  agg <- df %>% group_by(yb, xb) %>% summarise(v = sum(w), .groups="drop")
  total <- sum(agg$v)
  if (total <= 0) return(NA_real_)
  max(agg$v / total, na.rm = TRUE)
}

# Calculate PM maximum single-cell probability (candidates, take maximum across models)
.pmMaxForCandidates <- function(res, bins, weight_mode) {
  if (is.null(res)) return(NA_real_)
  node_attr <- res$node_attributes
  interactions_df <- res$interactions
  bests <- c()
  for (model_name in names(res$candidates_full_list)) {
    cand <- res$candidates_full_list[[model_name]]
    pass <- tryCatch({ cand %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id) }, error = function(e) character(0))
    if (length(pass) < 2) next
    sorted_genes <- node_attr %>% dplyr::filter(Gene %in% pass) %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    n_genes <- length(sorted_genes); if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
    if (nrow(sym) < 1) next
    coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
    row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
    xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
    w <- getPmWeightsFromMode(sym, node_attr, weight_mode)
    if (length(w) == 0 || sum(w) <= 0) next
  df <- data.frame(xb = xb, yb = yb, w = w)
  agg <- df %>% group_by(yb, xb) %>% summarise(v = sum(w), .groups="drop")
  total <- sum(agg$v); if (total <= 0) next
  bests <- c(bests, max(agg$v / total, na.rm = TRUE))
}
  if (length(bests) == 0) NA_real_ else max(bests, na.rm = TRUE)
}

# Calculate JP maximum single-cell probability (all, unweighted, directed)
.jpMaxForResult <- function(res, bins) {
  if (is.null(res)) return(NA_real_)
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes); if (n_genes < 2) return(NA_real_)
  gene_map <- setNames(1:n_genes, sorted_genes)
  interactions_df <- res$interactions
  filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) return(NA_real_)
  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
  # Observation count
  obs <- data.frame(xb=xb, yb=yb) %>% group_by(yb, xb) %>% summarise(cnt=n(), .groups="drop")
  # Bin protein counts on axes (directed)
  genes_idx <- 1:n_genes
  cbin <- .binIndex(genes_idx / n_genes, bins)
  rbin <- .binIndex(genes_idx / n_genes, bins)
  nx <- as.integer(tabulate(cbin, nbins=bins))
  ny <- as.integer(tabulate(rbin, nbins=bins))
  if (nrow(obs) == 0) return(NA_real_)
  obs <- obs %>% mutate(possible = nx[yb] * ny[xb],
                        prob = ifelse(possible > 0, cnt / possible, 0))
  max(obs$prob, na.rm = TRUE)
}

# Calculate JP maximum single-cell probability (candidates, unweighted)
.jpMaxForCandidates <- function(res, bins) {
  if (is.null(res)) return(NA_real_)
  node_attr <- res$node_attributes
  interactions_df <- res$interactions
  bests <- c()
  for (model_name in names(res$candidates_full_list)) {
    cand <- res$candidates_full_list[[model_name]]
    pass <- tryCatch({ cand %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id) }, error = function(e) character(0))
    if (length(pass) < 2) next
    sorted_genes <- node_attr %>% dplyr::filter(Gene %in% pass) %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    n_genes <- length(sorted_genes); if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
    if (nrow(sym) < 1) next
    coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
    row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
    xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
    obs <- data.frame(xb=xb, yb=yb) %>% group_by(yb, xb) %>% summarise(cnt=n(), .groups="drop")
    genes_idx <- 1:n_genes
    cbin <- .binIndex(genes_idx / n_genes, bins)
    rbin <- .binIndex(genes_idx / n_genes, bins)
    nx <- as.integer(tabulate(cbin, nbins=bins))
    ny <- as.integer(tabulate(rbin, nbins=bins))
    if (nrow(obs) == 0) next
    obs <- obs %>% mutate(possible = nx[yb] * ny[xb], prob = ifelse(possible > 0, cnt / possible, 0))
    bests <- c(bests, max(obs$prob, na.rm = TRUE))
  }
  if (length(bests) == 0) NA_real_ else max(bests, na.rm = TRUE)
}

# Unify bins
pm_bins <- pm_grid_bins_fixed
jp_bins <- jp_grid_bins_fixed

# PM: Calculate upper quantile for all/candidates by different weight modes
pm_max <- setNames(rep(0, length(pm_weight_modes)), pm_weight_modes)
pm_cand_max <- setNames(rep(0, length(pm_weight_modes)), pm_weight_modes)

for (mode in pm_weight_modes) {
  vals <- c()
  for (res in all_results) {
    if (is.null(res)) next
    z <- suppressWarnings(.pmMaxForResult(res, pm_bins, mode))
    if (is.finite(z)) vals <- c(vals, z)
  }
  if (length(vals) > 0) {
    pm_max[[mode]] <- as.numeric(stats::quantile(vals, probs = max(0, min(1, pm_max_quantile)), na.rm = TRUE))
  }

  vals_c <- c()
  for (res in all_results) {
    if (is.null(res)) next
    zc <- suppressWarnings(.pmMaxForCandidates(res, pm_bins, mode))
    if (is.finite(zc)) vals_c <- c(vals_c, zc)
  }
  if (length(vals_c) > 0) {
    pm_cand_max[[mode]] <- as.numeric(stats::quantile(vals_c, probs = max(0, min(1, pm_candidates_max_quantile)), na.rm = TRUE))
  }
}
plot_scales_config$pm_max_prob_by_mode <- pm_max
plot_scales_config$pm_candidates_max_prob_by_mode <- pm_cand_max
cat("PM global upper limit (quantile) calculated:\n"); print(pm_max)
cat("PM candidates upper limit (quantile) calculated:\n"); print(pm_cand_max)

# PM: Minimum single-cell probability quantile (for color lower limit), based on non-zero cells
pm_min <- setNames(rep(0, length(pm_weight_modes)), pm_weight_modes)
pm_cand_min <- setNames(rep(0, length(pm_weight_modes)), pm_weight_modes)

.pmProbValsForResult <- function(res, bins, weight_mode) {
  if (is.null(res)) return(numeric(0))
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes); if (n_genes < 2) return(numeric(0))
  gene_map <- setNames(1:n_genes, sorted_genes)
  interactions_df <- res$interactions
  filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) return(numeric(0))
  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
  w <- getPmWeightsFromMode(sym, node_attr, weight_mode)
  if (length(w) == 0 || sum(w) <= 0) return(numeric(0))
  df <- data.frame(x_bin = xb, y_bin = yb, w = w)
  agg <- df %>% group_by(y_bin, x_bin) %>% summarise(v = sum(w), .groups = "drop")
  total <- sum(agg$v); if (total <= 0) return(numeric(0))
  prob <- agg$v / total
  prob[is.finite(prob) & prob > 0]
}

.pmProbValsForCandidates <- function(res, bins, weight_mode) {
  if (is.null(res)) return(numeric(0))
  node_attr <- res$node_attributes
  interactions_df <- res$interactions
  vals <- c()
  for (model_name in names(res$candidates_full_list)) {
    cand <- res$candidates_full_list[[model_name]]
    pass <- tryCatch({ cand %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id) }, error = function(e) character(0))
    if (length(pass) < 2) next
    sorted_genes <- node_attr %>% dplyr::filter(Gene %in% pass) %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    n_genes <- length(sorted_genes); if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
    if (nrow(sym) < 1) next
    coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
    row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
    xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
    w <- getPmWeightsFromMode(sym, node_attr, weight_mode)
    if (length(w) == 0 || sum(w) <= 0) next
    df <- data.frame(x_bin = xb, y_bin = yb, w = w)
    agg <- df %>% group_by(y_bin, x_bin) %>% summarise(v = sum(w), .groups = "drop")
    total <- sum(agg$v); if (total <= 0) next
    prob <- agg$v / total
    vals <- c(vals, prob[is.finite(prob) & prob > 0])
  }
  vals
}

for (mode in pm_weight_modes) {
  all_vals <- c(); cand_vals <- c()
  for (res in all_results) {
    if (is.null(res)) next
    all_vals <- c(all_vals, suppressWarnings(.pmProbValsForResult(res, pm_bins, mode)))
    cand_vals <- c(cand_vals, suppressWarnings(.pmProbValsForCandidates(res, pm_bins, mode)))
  }
  if (length(all_vals) > 0 && is.finite(pm_min_quantile) && pm_min_quantile > 0) {
    pm_min[[mode]] <- as.numeric(stats::quantile(all_vals, probs = max(0, min(1, pm_min_quantile)), na.rm = TRUE))
  }
  if (length(cand_vals) > 0 && is.finite(pm_candidates_min_quantile) && pm_candidates_min_quantile > 0) {
    pm_cand_min[[mode]] <- as.numeric(stats::quantile(cand_vals, probs = max(0, min(1, pm_candidates_min_quantile)), na.rm = TRUE))
  }
}
plot_scales_config$pm_min_prob_by_mode <- pm_min
plot_scales_config$pm_candidates_min_prob_by_mode <- pm_cand_min


# JP: Unweighted (all/candidates)
jp_vals <- c()
for (res in all_results) {
  if (is.null(res)) next
  z <- suppressWarnings(.jpMaxForResult(res, jp_bins))
  if (is.finite(z)) jp_vals <- c(jp_vals, z)
}
plot_scales_config$jp_max_prob <- if (length(jp_vals) > 0) {
  as.numeric(stats::quantile(jp_vals, probs = max(0, min(1, jp_max_quantile)), na.rm = TRUE))
} else { 0 }

jp_cvals <- c()
for (res in all_results) {
  if (is.null(res)) next
  zc <- suppressWarnings(.jpMaxForCandidates(res, jp_bins))
  if (is.finite(zc)) jp_cvals <- c(jp_cvals, zc)
}
plot_scales_config$jp_candidates_max_prob <- if (length(jp_cvals) > 0) {
  as.numeric(stats::quantile(jp_cvals, probs = max(0, min(1, jp_candidates_max_quantile)), na.rm = TRUE))
} else { 0 }

cat(paste("JP all maximum quantile:", round(plot_scales_config$jp_max_prob, 8), "\n"))
cat(paste("JP candidates maximum quantile:", round(plot_scales_config$jp_candidates_max_prob, 8), "\n"))

# JP: Minimum single-cell probability quantile (based on non-zero cells)
.jpProbValsForResult <- function(res, bins) {
  if (is.null(res)) return(numeric(0))
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes); if (n_genes < 2) return(numeric(0))
  gene_map <- setNames(1:n_genes, sorted_genes)
  interactions_df <- res$interactions
  filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) return(numeric(0))
  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
  obs <- data.frame(x_bin = xb, y_bin = yb) %>% group_by(y_bin, x_bin) %>% summarise(cnt = n(), .groups = "drop")
  genes_idx <- 1:n_genes
  cbin <- .binIndex(genes_idx / n_genes, bins)
  rbin <- .binIndex(genes_idx / n_genes, bins)
  nx <- as.integer(tabulate(cbin, nbins = bins))
  ny <- as.integer(tabulate(rbin, nbins = bins))
  if (nrow(obs) == 0) return(numeric(0))
  prob <- with(obs, ifelse(nx[y_bin] * ny[x_bin] > 0, cnt / (nx[y_bin] * ny[x_bin]), 0))
  prob[is.finite(prob) & prob > 0]
}

.jpProbValsForCandidates <- function(res, bins) {
  if (is.null(res)) return(numeric(0))
  node_attr <- res$node_attributes
  interactions_df <- res$interactions
  vals <- c()
  for (model_name in names(res$candidates_full_list)) {
    cand <- res$candidates_full_list[[model_name]]
    pass <- tryCatch({ cand %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id) }, error = function(e) character(0))
    if (length(pass) < 2) next
    sorted_genes <- node_attr %>% dplyr::filter(Gene %in% pass) %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    n_genes <- length(sorted_genes); if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
    if (nrow(sym) < 1) next
    coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
    row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
    xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
    obs <- data.frame(x_bin = xb, y_bin = yb) %>% group_by(y_bin, x_bin) %>% summarise(cnt = n(), .groups = "drop")
    genes_idx <- 1:n_genes
    cbin <- .binIndex(genes_idx / n_genes, bins)
    rbin <- .binIndex(genes_idx / n_genes, bins)
    nx <- as.integer(tabulate(cbin, nbins = bins))
    ny <- as.integer(tabulate(rbin, nbins = bins))
    if (nrow(obs) == 0) next
    prob <- with(obs, ifelse(nx[y_bin] * ny[x_bin] > 0, cnt / (nx[y_bin] * ny[x_bin]), 0))
    vals <- c(vals, prob[is.finite(prob) & prob > 0])
  }
  vals
}

jp_all_vals <- c(); jp_cand_vals <- c()
for (res in all_results) {
  if (is.null(res)) next
  jp_all_vals <- c(jp_all_vals, suppressWarnings(.jpProbValsForResult(res, jp_bins)))
  jp_cand_vals <- c(jp_cand_vals, suppressWarnings(.jpProbValsForCandidates(res, jp_bins)))
}
plot_scales_config$jp_min_prob <- if (length(jp_all_vals) > 0 && is.finite(jp_min_quantile) && jp_min_quantile > 0) {
  as.numeric(stats::quantile(jp_all_vals, probs = max(0, min(1, jp_min_quantile)), na.rm = TRUE))
} else { 0 }
plot_scales_config$jp_candidates_min_prob <- if (length(jp_cand_vals) > 0 && is.finite(jp_candidates_min_quantile) && jp_candidates_min_quantile > 0) {
  as.numeric(stats::quantile(jp_cand_vals, probs = max(0, min(1, jp_candidates_min_quantile)), na.rm = TRUE))
} else { 0 }


# ================== JP Weighted Global Color Upper/Lower Limit Pre-computation ==================
cat("Pre-computing JP weighted global color upper/lower limits...\n")

getJpNormalizedWeights <- function(symmetric_interactions, node_attr, weight_mode) {
  n <- nrow(symmetric_interactions)
  if (n == 0) return(numeric(0))
  # Normalize combined_score to [0,1]
  base_w <- symmetric_interactions$combined_score
  base_w[!is.finite(base_w) | is.na(base_w)] <- 0
  base_w[base_w < 0] <- 0
  base_w[base_w > 1] <- 1
  # TP weights
  gcol <- getGroupColName(node_attr)
  if (!is.null(gcol)) {
    group_map <- setNames(node_attr[[gcol]], node_attr$Gene)
  } else {
    group_map <- setNames(rep(NA_character_, nrow(node_attr)), node_attr$Gene)
  }
  tp1 <- group_map[symmetric_interactions$protein1] == "TP"
  tp2 <- group_map[symmetric_interactions$protein2] == "TP"
  tp1[is.na(tp1)] <- FALSE; tp2[is.na(tp2)] <- FALSE
  tp_w_raw <- ifelse(tp1 & tp2, tp_weight_both,
                 ifelse(xor(tp1, tp2), tp_weight_single, tp_weight_none))
  max_tp_w <- max(tp_weight_both, tp_weight_single, tp_weight_none, na.rm = TRUE)
  if (!is.finite(max_tp_w) || max_tp_w <= 0) max_tp_w <- 1
  tp_w <- pmax(0, pmin(1, tp_w_raw / max_tp_w))

  if (identical(weight_mode, "combined_score")) {
    w <- base_w
  } else if (identical(weight_mode, "tp")) {
    w <- tp_w
  } else { # combined_score_tp
    # Normalize to <=1: multiply then divide by max tp weight (equivalent to base_w * (tp_w_raw/max_tp_w))
    w <- base_w * tp_w
  }
  # Power enhancement
  if (exists("jp_weight_power", inherits = TRUE) && is.finite(jp_weight_power) && jp_weight_power != 1) {
    w <- w ^ jp_weight_power
  }
  w[!is.finite(w) | is.na(w) | w < 0] <- 0
  w[w > 1] <- 1
  w
}

.jpwMaxForResult <- function(res, bins, weight_mode) {
  if (is.null(res)) return(NA_real_)
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes); if (n_genes < 2) return(NA_real_)
  gene_map <- setNames(1:n_genes, sorted_genes)
  interactions_df <- res$interactions
  filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) return(NA_real_)
  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
  w <- getJpNormalizedWeights(sym, node_attr, weight_mode)
  if (length(w) == 0) return(NA_real_)
  obs <- data.frame(x_bin = xb, y_bin = yb, w = w) %>% group_by(y_bin, x_bin) %>% summarise(ws = sum(w), .groups = "drop")
  genes_idx <- 1:n_genes
  cbin <- .binIndex(genes_idx / n_genes, bins)
  rbin <- .binIndex(genes_idx / n_genes, bins)
  nx <- as.integer(tabulate(cbin, nbins = bins))
  ny <- as.integer(tabulate(rbin, nbins = bins))
  if (nrow(obs) == 0) return(NA_real_)
  obs <- obs %>% mutate(possible = nx[y_bin] * ny[x_bin], prob = ifelse(possible > 0, ws / possible, 0))
  max(obs$prob, na.rm = TRUE)
}

.jpwMaxForCandidates <- function(res, bins, weight_mode) {
  if (is.null(res)) return(NA_real_)
  node_attr <- res$node_attributes
  interactions_df <- res$interactions
  bests <- c()
  for (model_name in names(res$candidates_full_list)) {
    cand <- res$candidates_full_list[[model_name]]
    pass <- tryCatch({ cand %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id) }, error = function(e) character(0))
    if (length(pass) < 2) next
    sorted_genes <- node_attr %>% dplyr::filter(Gene %in% pass) %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    n_genes <- length(sorted_genes); if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
    if (nrow(sym) < 1) next
    coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
    row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
    xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
    w <- getJpNormalizedWeights(sym, node_attr, weight_mode)
    if (length(w) == 0) next
    obs <- data.frame(x_bin = xb, y_bin = yb, w = w) %>% group_by(y_bin, x_bin) %>% summarise(ws = sum(w), .groups = "drop")
    genes_idx <- 1:n_genes
    cbin <- .binIndex(genes_idx / n_genes, bins)
    rbin <- .binIndex(genes_idx / n_genes, bins)
    nx <- as.integer(tabulate(cbin, nbins = bins))
    ny <- as.integer(tabulate(rbin, nbins = bins))
    if (nrow(obs) == 0) next
    obs <- obs %>% mutate(possible = nx[y_bin] * ny[x_bin], prob = ifelse(possible > 0, ws / possible, 0))
    bests <- c(bests, max(obs$prob, na.rm = TRUE))
  }
  if (length(bests) == 0) NA_real_ else max(bests, na.rm = TRUE)
}

.jpwProbValsForResult <- function(res, bins, weight_mode) {
  if (is.null(res)) return(numeric(0))
  node_attr <- res$node_attributes
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes); if (n_genes < 2) return(numeric(0))
  gene_map <- setNames(1:n_genes, sorted_genes)
  interactions_df <- res$interactions
  filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) return(numeric(0))
  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
  w <- getJpNormalizedWeights(sym, node_attr, weight_mode)
  if (length(w) == 0) return(numeric(0))
  obs <- data.frame(x_bin = xb, y_bin = yb, w = w) %>% group_by(y_bin, x_bin) %>% summarise(ws = sum(w), .groups = "drop")
  genes_idx <- 1:n_genes
  cbin <- .binIndex(genes_idx / n_genes, bins)
  rbin <- .binIndex(genes_idx / n_genes, bins)
  nx <- as.integer(tabulate(cbin, nbins = bins))
  ny <- as.integer(tabulate(rbin, nbins = bins))
  if (nrow(obs) == 0) return(numeric(0))
  prob <- with(obs, ifelse(nx[y_bin] * ny[x_bin] > 0, ws / (nx[y_bin] * ny[x_bin]), 0))
  prob[is.finite(prob) & prob > 0]
}

.jpwProbValsForCandidates <- function(res, bins, weight_mode) {
  if (is.null(res)) return(numeric(0))
  node_attr <- res$node_attributes
  interactions_df <- res$interactions
  vals <- c()
  for (model_name in names(res$candidates_full_list)) {
    cand <- res$candidates_full_list[[model_name]]
    pass <- tryCatch({ cand %>% dplyr::filter(PassThreshold == "Pass") %>% dplyr::pull(protein_id) }, error = function(e) character(0))
    if (length(pass) < 2) next
    sorted_genes <- node_attr %>% dplyr::filter(Gene %in% pass) %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
    n_genes <- length(sorted_genes); if (n_genes < 2) next
    gene_map <- setNames(1:n_genes, sorted_genes)
    filtered <- interactions_df %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
    sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
    if (nrow(sym) < 1) next
    coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
    row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
    xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)
    w <- getJpNormalizedWeights(sym, node_attr, weight_mode)
    if (length(w) == 0) next
    obs <- data.frame(x_bin = xb, y_bin = yb, w = w) %>% group_by(y_bin, x_bin) %>% summarise(ws = sum(w), .groups = "drop")
    genes_idx <- 1:n_genes
    cbin <- .binIndex(genes_idx / n_genes, bins)
    rbin <- .binIndex(genes_idx / n_genes, bins)
    nx <- as.integer(tabulate(cbin, nbins = bins))
    ny <- as.integer(tabulate(rbin, nbins = bins))
    if (nrow(obs) == 0) next
    prob <- with(obs, ifelse(nx[y_bin] * ny[x_bin] > 0, ws / (nx[y_bin] * ny[x_bin]), 0))
    vals <- c(vals, prob[is.finite(prob) & prob > 0])
  }
  vals
}

# Calculate weighted JP all/candidates upper/lower limits
jpw_max <- setNames(rep(0, length(jp_weight_modes)), jp_weight_modes)
jpw_cand_max <- setNames(rep(0, length(jp_weight_modes)), jp_weight_modes)
jpw_min <- setNames(rep(0, length(jp_weight_modes)), jp_weight_modes)
jpw_cand_min <- setNames(rep(0, length(jp_weight_modes)), jp_weight_modes)

for (mode in jp_weight_modes) {
  # Upper limit
  vals <- c(); vals_c <- c()
  for (res in all_results) {
    if (is.null(res)) next
    z <- suppressWarnings(.jpwMaxForResult(res, jp_bins, mode))
    if (is.finite(z)) vals <- c(vals, z)
    zc <- suppressWarnings(.jpwMaxForCandidates(res, jp_bins, mode))
    if (is.finite(zc)) vals_c <- c(vals_c, zc)
  }
  if (length(vals) > 0) {
    jpw_max[[mode]] <- as.numeric(stats::quantile(vals, probs = max(0, min(1, jpw_max_quantile)), na.rm = TRUE))
  }
  if (length(vals_c) > 0) {
    jpw_cand_max[[mode]] <- as.numeric(stats::quantile(vals_c, probs = max(0, min(1, jpw_candidates_max_quantile)), na.rm = TRUE))
  }

  # Lower limit (non-zero cells)
  all_vals <- c(); cand_vals <- c()
  for (res in all_results) {
    if (is.null(res)) next
    all_vals <- c(all_vals, suppressWarnings(.jpwProbValsForResult(res, jp_bins, mode)))
    cand_vals <- c(cand_vals, suppressWarnings(.jpwProbValsForCandidates(res, jp_bins, mode)))
  }
  if (length(all_vals) > 0 && is.finite(jpw_min_quantile) && jpw_min_quantile > 0) {
    jpw_min[[mode]] <- as.numeric(stats::quantile(all_vals, probs = max(0, min(1, jpw_min_quantile)), na.rm = TRUE))
  }
  if (length(cand_vals) > 0 && is.finite(jpw_candidates_min_quantile) && jpw_candidates_min_quantile > 0) {
    jpw_cand_min[[mode]] <- as.numeric(stats::quantile(cand_vals, probs = max(0, min(1, jpw_candidates_min_quantile)), na.rm = TRUE))
  }
}

plot_scales_config$jpw_max_prob_by_mode <- jpw_max
plot_scales_config$jpw_candidates_max_prob_by_mode <- jpw_cand_max
plot_scales_config$jpw_min_prob_by_mode <- jpw_min
plot_scales_config$jpw_candidates_min_prob_by_mode <- jpw_cand_min



# -------------------------------------------------------------------
# 4. Plotting Function Definitions (Revised)
# -------------------------------------------------------------------

# 4a. ROC Curves and Youden Index Plots (Fully Implemented)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_model_performance <- function(result, plot_dir_roc, plot_dir_youden, cfg) {
  if (is.null(result) || length(result$auc_scores) == 0) return()
  
  df_name <- result$df_name; group_name <- result$group_name
  plot_data_full <- result$node_attributes
  roc_plots <- list(); youden_plots <- list()
  
  # Store plotting data for all models
  all_roc_data <- list()
  all_youden_data <- list()
  
  models_to_run <- result$candidates_full_list
  for (model_name in names(models_to_run)) {
    scored_df <- models_to_run[[model_name]]
    roc_obj <- roc(response = plot_data_full$plot_group, predictor = scored_df$model_score,
                   levels = c("Other", "TP"), quiet = TRUE, direction = "<")
    auc_val <- auc(roc_obj)
    
    # Extract ROC curve data
    roc_coords <- coords(roc_obj, "all", ret = c("threshold", "specificity", "sensitivity"))
    roc_coords$model <- model_name
    roc_coords$AUC <- as.numeric(auc_val)
    roc_coords$df_name <- df_name
    roc_coords$group_name <- group_name
    all_roc_data[[model_name]] <- roc_coords
    
    # ROC Plot
    p_roc <- ggroc(roc_obj, color = cfg$roc_line_color, size = 1) +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
      ggtitle(model_name, subtitle = paste0("AUC = ", round(auc_val, 3))) +
      theme(plot.title = element_text(size=10), plot.subtitle = element_text(size=9))
    roc_plots[[model_name]] <- p_roc
    
    # --- YOUDEN INDEX PLOT CORRECTION ---
    # 1. Get sensitivity and specificity instead of 'youden' directly
    coords_df <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
    
    # 2. Manually calculate the Youden index based on its definition
    coords_df$youden <- coords_df$sensitivity + coords_df$specificity - 1
    
    # Find the threshold that maximizes our correctly calculated Youden index
    best_idx <- which.max(coords_df$youden)
    best_thresh <- coords_df$threshold[best_idx]
    
    # Save Youden data
    coords_df$model <- model_name
    coords_df$best_threshold <- best_thresh
    coords_df$df_name <- df_name
    coords_df$group_name <- group_name
    all_youden_data[[model_name]] <- coords_df
    
    # 3. Plot the correctly calculated Youden index
    p_youden <- ggplot(coords_df, aes(x = threshold, y = youden)) +
      geom_line(color = cfg$youden_line_color, size = 1) +
      geom_vline(xintercept = best_thresh, linetype = "dashed", color = "red") +
      # Use the configured y-axis range, which now correctly fits the data
      coord_cartesian(ylim = cfg$youden_y_axis_range) +
      ggtitle(model_name, subtitle = paste("Best Threshold =", round(best_thresh, 3))) +
      labs(x = "Threshold", y = "Youden Index (Sens + Spec - 1)") +
      theme(plot.title = element_text(size=10), plot.subtitle = element_text(size=9))
    youden_plots[[model_name]] <- p_youden
  }
  
  if (length(roc_plots) > 0) {
    roc_combined <- wrap_plots(roc_plots, ncol = 4)
    ggsave(file.path(plot_dir_roc, paste0(df_name, "_", group_name, "_ROC.pdf")),
           roc_combined, width = 16, height = ceiling(length(roc_plots)/4) * 4)
    youden_combined <- wrap_plots(youden_plots, ncol = 4)
    ggsave(file.path(plot_dir_youden, paste0(df_name, "_", group_name, "_Youden.pdf")),
           youden_combined, width = 16, height = ceiling(length(youden_plots)/4) * 4)
    
    # Export ROC and Youden plotting data
    roc_data_combined <- bind_rows(all_roc_data)
    write.csv(roc_data_combined, 
              file.path(plot_dir_roc, paste0(df_name, "_", group_name, "_ROC_PlotData.csv")),
              row.names = FALSE)
    
    youden_data_combined <- bind_rows(all_youden_data)
    write.csv(youden_data_combined,
              file.path(plot_dir_youden, paste0(df_name, "_", group_name, "_Youden_PlotData.csv")),
              row.names = FALSE)
    
    cat(paste("  -> ROC and Youden plotting data exported:", df_name, "-", group_name, "\n"))
  }
}

# 4b. Scatterplots with Decision Boundaries (Fully Implemented)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_decision_boundaries <- function(result, scales_config, plot_dir, cfg) {
  if (is.null(result) || length(result$auc_scores) == 0) return()
  
  df_name <- result$df_name; group_name <- result$group_name
  logFC_col <- group_info[[group_name]]$logFC
  plot_data_full <- result$node_attributes
  thresholds <-sapply(result$candidates_full_list, function(x) x$Model_best_threshold[1])
  all_plots <- list()
  
  get_lim <- function(override, fallback) { if (!is.null(override)) override else fallback }
  
  # Define all base plots with controlled axes
  p_base_ab_lfc <- ggplot(plot_data_full, aes(x = mean_group_medianNorm, y = .data[[logFC_col]])) +
    geom_point(aes(color = plot_group), alpha = 0.5) + scale_color_manual(values = TP_vs_other_colors, name="Group") +
    coord_cartesian(xlim = get_lim(cfg$axis_range_abundance, scales_config$abundance_range), ylim = get_lim(cfg$axis_range_logfc, scales_config$logfc_range)) + 
    labs(x = "Abundance", y = "logFC")
  
  p_base_ab_ppi_g <- ggplot(plot_data_full, aes(x = mean_group_medianNorm, y = log_ppi_score_global)) +
    geom_point(aes(color = plot_group), alpha = 0.5) + scale_color_manual(values = TP_vs_other_colors, name="Group") +
    coord_cartesian(xlim = get_lim(cfg$axis_range_abundance, scales_config$abundance_range), ylim = get_lim(cfg$axis_range_ppi_global, scales_config$ppi_global_range)) + 
    labs(x = "Abundance", y = "log(Global PPI)")
  
  p_base_ab_ppi_l <- ggplot(plot_data_full, aes(x = mean_group_medianNorm, y = log_ppi_score_local)) +
    geom_point(aes(color = plot_group), alpha = 0.5) + scale_color_manual(values = TP_vs_other_colors, name="Group") +
    coord_cartesian(xlim = get_lim(cfg$axis_range_abundance, scales_config$abundance_range), ylim = get_lim(cfg$axis_range_ppi_local, scales_config$ppi_local_range)) + 
    labs(x = "Abundance", y = "log(Local PPI)")
  
  p_base_lfc_ppi_g <- ggplot(plot_data_full, aes(x = .data[[logFC_col]], y = log_ppi_score_global)) +
    geom_point(aes(color = plot_group), alpha = 0.5) + scale_color_manual(values = TP_vs_other_colors, name="Group") +
    coord_cartesian(xlim = get_lim(cfg$axis_range_logfc, scales_config$logfc_range), ylim = get_lim(cfg$axis_range_ppi_global, scales_config$ppi_global_range)) + 
    labs(x = "logFC", y = "log(Global PPI)")
  
  p_base_lfc_ppi_l <- ggplot(plot_data_full, aes(x = .data[[logFC_col]], y = log_ppi_score_local)) +
    geom_point(aes(color = plot_group), alpha = 0.5) + scale_color_manual(values = TP_vs_other_colors, name="Group") +
    coord_cartesian(xlim = get_lim(cfg$axis_range_logfc, scales_config$logfc_range), ylim = get_lim(cfg$axis_range_ppi_local, scales_config$ppi_local_range)) + 
    labs(x = "logFC", y = "log(Local PPI)")
  
  # 1D Models
  all_plots[['1D_Abundance']] <- p_base_ab_lfc + geom_vline(xintercept = thresholds["1D_Abundance"], color = "red", linetype="dashed") + ggtitle("1D Abundance Boundary")
  all_plots[['1D_logFC']] <- p_base_ab_lfc + geom_hline(yintercept = thresholds["1D_logFC"], color = "blue", linetype="dashed") + ggtitle("1D logFC Boundary")
  all_plots[['1D_PPI_Global']] <- p_base_ab_ppi_g + geom_hline(yintercept = thresholds["1D_PPI_Global"], color = "purple", linetype="dashed") + ggtitle("1D Global PPI Boundary")
  all_plots[['1D_PPI_Local']] <- p_base_ab_ppi_l + geom_hline(yintercept = thresholds["1D_PPI_Local"], color = "orange", linetype="dashed") + ggtitle("1D Local PPI Boundary")
  
  # Helper function to get GLM boundary parameters
  get_boundary <- function(model_formula, threshold) {
    model <- glm(as.formula(model_formula), data = plot_data_full, family = "binomial")
    if (is.na(threshold) || threshold >= 1) { logit_thresh <- 1e6 } 
    else if (threshold <= 0) { logit_thresh <- -1e6 }
    else { logit_thresh <- log(threshold / (1 - threshold)) }
    list(intercept = (logit_thresh - coef(model)[1]) / coef(model)[3], slope = -coef(model)[2] / coef(model)[3])
  }
  
  # 2D Models with tryCatch for robustness
  tryCatch({
    b <- get_boundary(paste("plot_group ~ mean_group_medianNorm +", logFC_col), thresholds["2D_Abundance_logFC"])
    all_plots[['2D_Abundance_logFC']] <- p_base_ab_lfc + geom_abline(intercept = b$intercept, slope = b$slope, color = "darkgreen", size=1) + ggtitle("2D Abundance + logFC Boundary")
  }, error = function(e){cat("Could not plot 2D_Abundance_logFC boundary for", df_name, group_name, "\n")})
  
  tryCatch({
    b <- get_boundary("plot_group ~ mean_group_medianNorm + log_ppi_score_global", thresholds["2D_Abundance_PPI_Global"])
    all_plots[['2D_Abundance_PPI_Global']] <- p_base_ab_ppi_g + geom_abline(intercept = b$intercept, slope = b$slope, color = "darkgreen", size=1) + ggtitle("2D Abundance + Global PPI Boundary")
  }, error = function(e){cat("Could not plot 2D_Abundance_PPI_Global boundary for", df_name, group_name, "\n")})
  
  tryCatch({
    b <- get_boundary("plot_group ~ mean_group_medianNorm + log_ppi_score_local", thresholds["2D_Abundance_PPI_Local"])
    all_plots[['2D_Abundance_PPI_Local']] <- p_base_ab_ppi_l + geom_abline(intercept = b$intercept, slope = b$slope, color = "darkgreen", size=1) + ggtitle("2D Abundance + Local PPI Boundary")
  }, error = function(e){cat("Could not plot 2D_Abundance_PPI_Local boundary for", df_name, group_name, "\n")})
  
  tryCatch({
    b <- get_boundary(paste("plot_group ~", logFC_col, "+ log_ppi_score_global"), thresholds["2D_logFC_PPI_Global"])
    all_plots[['2D_logFC_PPI_Global']] <- p_base_lfc_ppi_g + geom_abline(intercept = b$intercept, slope = b$slope, color = "darkgreen", size=1) + ggtitle("2D logFC + Global PPI Boundary")
  }, error = function(e){cat("Could not plot 2D_logFC_PPI_Global boundary for", df_name, group_name, "\n")})
  
  tryCatch({
    b <- get_boundary(paste("plot_group ~", logFC_col, "+ log_ppi_score_local"), thresholds["2D_logFC_PPI_Local"])
    all_plots[['2D_logFC_PPI_Local']] <- p_base_lfc_ppi_l + geom_abline(intercept = b$intercept, slope = b$slope, color = "darkgreen", size=1) + ggtitle("2D logFC + Local PPI Boundary")
  }, error = function(e){cat("Could not plot 2D_logFC_PPI_Local boundary for", df_name, group_name, "\n")})
  
  
  if(length(all_plots) %% 2 != 0) { all_plots[['spacer']] <- plot_spacer() }
  
  combined_plot <- wrap_plots(all_plots, ncol=2, guides = "collect") + 
    plot_annotation(title=paste("Decision Boundaries for", df_name, "-", group_name)) & 
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, paste0(df_name, "_", group_name, "_Boundaries.pdf")),
         combined_plot, width = 12, height = ceiling(length(all_plots)/2) * 5)
}



# 4c. 3D-like Scatter plots (Fully Implemented)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_3d_scatter <- function(result, scales_config, cfg, plot_dir) {
  if (is.null(result) || length(result$auc_scores) == 0) return()
  
  df_name <- result$df_name; group_name <- result$group_name; logFC_col <- group_info[[group_name]]$logFC
  plot_data_full <- result$node_attributes
  
  # Get decision boundary thresholds
  thresholds <- sapply(result$candidates_full_list, function(x) x$Model_best_threshold[1])
  
  # Use configured axis ranges, consistent with decision boundaries
  get_lim <- function(override, fallback) { if (!is.null(override)) override else fallback }
  xlim_range <- get_lim(cfg$axis_range_abundance, scales_config$abundance_range)
  ylim_range <- get_lim(cfg$axis_range_logfc, scales_config$logfc_range)
  
  # Calculate decision boundary parameters for 2D models
  boundary_params <- NULL
  tryCatch({
    model <- glm(as.formula(paste("plot_group ~ mean_group_medianNorm +", logFC_col)), 
                 data = plot_data_full, family = "binomial")
    threshold <- thresholds["2D_Abundance_logFC"]
    if (!is.na(threshold) && threshold > 0 && threshold < 1) {
      logit_thresh <- log(threshold / (1 - threshold))
      boundary_params <- list(
        intercept = (logit_thresh - coef(model)[1]) / coef(model)[3],
        slope = -coef(model)[2] / coef(model)[3]
      )
    }
  }, error = function(e) { cat("Warning: Could not calculate 2D boundary for", df_name, group_name, "\n") })
  
  # Plot first figure: colored by Global PPI
  p1 <- ggplot(plot_data_full, aes(x = mean_group_medianNorm, y = .data[[logFC_col]])) +
    geom_point(aes(color = log_ppi_score_global)) +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    scale_color_gradientn(colors = cfg$scatter_3d_color_palette, # Configurable palette
                          name = "log(Global PPI + 1)", 
                          limits = scales_config$ppi_global_range) +
    labs(x = "Abundance", y = "logFC", title = "Color by Global PPI")
  
  # Add decision boundary lines
  if (!is.na(thresholds["1D_Abundance"])) {
    p1 <- p1 + geom_vline(xintercept = thresholds["1D_Abundance"], linetype = "dashed", color = "red", size = 0.8)
  }
  if (!is.na(thresholds["1D_logFC"])) {
    p1 <- p1 + geom_hline(yintercept = thresholds["1D_logFC"], linetype = "dashed", color = "blue", size = 0.8)
  }
  if (!is.null(boundary_params)) {
    p1 <- p1 + geom_abline(intercept = boundary_params$intercept, slope = boundary_params$slope, 
                           linetype = "dashed", color = "darkgreen", size = 0.8)
  }
  
  # Plot second figure: colored by Local PPI
  p2 <- ggplot(plot_data_full, aes(x = mean_group_medianNorm, y = .data[[logFC_col]])) +
    geom_point(aes(color = log_ppi_score_local)) +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    scale_color_gradientn(colors = cfg$scatter_3d_color_palette, # Configurable palette
                          name = "log(Local PPI + 1)", 
                          limits = scales_config$ppi_local_range) +
    labs(x = "Abundance", y = "logFC", title = "Color by Local PPI")
  
  # Add decision boundary lines
  if (!is.na(thresholds["1D_Abundance"])) {
    p2 <- p2 + geom_vline(xintercept = thresholds["1D_Abundance"], linetype = "dashed", color = "red", size = 0.8)
  }
  if (!is.na(thresholds["1D_logFC"])) {
    p2 <- p2 + geom_hline(yintercept = thresholds["1D_logFC"], linetype = "dashed", color = "blue", size = 0.8)
  }
  if (!is.null(boundary_params)) {
    p2 <- p2 + geom_abline(intercept = boundary_params$intercept, slope = boundary_params$slope, 
                           linetype = "dashed", color = "darkgreen", size = 0.8)
  }
  
  all_plots <- p1 | p2
  ggsave(file.path(plot_dir, paste0(df_name, "_", group_name, "_3D_Scatters.pdf")),
         all_plots, width = 14, height = 6)
}


# 4d. PPI Heatmap (Added Safety Check and Symmetry)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ppi_heatmap <- function(result, cfg, plot_dir){
  if (is.null(result)) return()
  df_name <- result$df_name; group_name <- result$group_name
  node_attr <- result$node_attributes; interactions_df <- result$interactions
  sorted_df <- node_attr %>% arrange(desc(mean_group_medianNorm))
  sorted_genes <- sorted_df$Gene
  gene_map <- setNames(1:length(sorted_genes), sorted_genes)
  
  # Filter interactions for current group
  filtered_interactions <- interactions_df %>% 
    dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  
  if(nrow(filtered_interactions) == 0) {
    cat("INFO: Skipping PPI heatmap for", df_name, group_name, "due to no interactions found within the group.\n")
    return()
  }
  
  # Symmetrize interactions: add (B,A) for each (A,B) to ensure matrix symmetry along diagonal
  symmetric_interactions <- bind_rows(
    filtered_interactions,
    filtered_interactions %>% 
      rename(protein1 = protein2, protein2 = protein1, combined_score = combined_score)
  ) %>%
  distinct(protein1, protein2, .keep_all = TRUE)  # Remove duplicates, keep first occurrence score
  
  plot_df <- symmetric_interactions %>%
    mutate(x = gene_map[protein1], y = gene_map[protein2], score = combined_score) %>%
    dplyr::select(x, y, score)
  
  # Export plotting data matrix for inspection
  output_csv_path <- file.path(plot_dir, paste0(df_name, "_", group_name, "_PPI_Heatmap_PlotData.csv"))
  write.csv(plot_df, output_csv_path, row.names = FALSE)
  cat(paste("  -> PPI Heatmap plotting data exported:", output_csv_path, "\n"))
  
  p_main <- ggplot(plot_df, aes(x = x, y = y, fill = score)) +
    geom_tile() +
    scale_fill_gradientn(colors = cfg$heatmap_color_palette, name = "STRING Score", limits = cfg$heatmap_score_range) +
    coord_fixed() + theme_void() + theme(legend.position = "right")
  
  # Compatible with both plot_group and group_Localization column names
  group_col_name <- if ("plot_group" %in% colnames(sorted_df)) {
    "plot_group"
  } else if ("group_Localization" %in% colnames(sorted_df)) {
    "group_Localization"
  } else {
    NULL
  }
  
  annot_df <- sorted_df %>% mutate(index = 1:n())
  
  if (!is.null(group_col_name)) {
    annot_df <- annot_df %>% 
      mutate(group_label = factor(.data[[group_col_name]], levels = c("Other", "TP")))
    p_top <- ggplot(annot_df, aes(x = index, y = 1, fill = group_label)) + geom_tile() + scale_fill_manual(values = TP_vs_other_colors) + theme_void() + theme(legend.position = "none")
    p_right <- ggplot(annot_df, aes(x = 1, y = index, fill = group_label)) + geom_tile() + scale_fill_manual(values = TP_vs_other_colors) + theme_void() + theme(legend.position = "none")
  } else {
    p_top <- plot_spacer()
    p_right <- plot_spacer()
  }
  
  p_final <- (p_top | plot_spacer()) / (p_main | p_right) + 
    plot_layout(widths = c(8, 1), heights = c(1, 8), guides = 'collect') +
    plot_annotation(title = paste(df_name, "-", group_name), subtitle = "STRING Interaction Heatmap")
  
  ggsave(file.path(plot_dir, paste0(df_name, "_", group_name, "_PPI_Heatmap.pdf")), p_final, width = 9, height = 9)
}

# 4f. Candidate Proteins PPI Heatmaps (based on protein list after PassThreshold) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_candidate_ppi_heatmaps <- function(result, cfg, plot_dir){
  if (is.null(result) || length(result$candidates_full_list) == 0) return()
  
  df_name <- result$df_name
  group_name <- result$group_name
  interactions_df <- result$interactions
  
  cat(paste("  Generating candidate protein PPI heatmaps for", df_name, "-", group_name, "...\n"))
  
  # Generate heatmap for each model
  for (model_name in names(result$candidates_full_list)) {
    tryCatch({
      # Get candidate proteins that passed threshold for this model
      candidates_df <- result$candidates_full_list[[model_name]]
      pass_proteins <- candidates_df %>%
        filter(PassThreshold == "Pass") %>%
        pull(protein_id)
      
      if (length(pass_proteins) < 5) {
        cat(paste("    ⚠️  Model", model_name, "has too few candidate proteins (<5), skipping heatmap generation\n"))
        next
      }
      
      cat(paste("    - Model", model_name, ": ", length(pass_proteins), "candidate proteins\n"))
      
      # Sort candidate proteins by abundance
      sorted_candidates <- result$node_attributes %>%
        filter(Gene %in% pass_proteins) %>%
        arrange(desc(mean_group_medianNorm))
      
      sorted_genes <- sorted_candidates$Gene
      gene_map <- setNames(1:length(sorted_genes), sorted_genes)
      
      # Filter interactions between candidate proteins
      filtered_interactions <- interactions_df %>%
        filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
      
      if (nrow(filtered_interactions) == 0) {
        cat(paste("    ⚠️  No interactions found between candidate proteins, skipping\n"))
        next
      }
      
      # Symmetrize interactions
      symmetric_interactions <- bind_rows(
        filtered_interactions,
        filtered_interactions %>% 
          rename(protein1 = protein2, protein2 = protein1)
      ) %>%
      distinct(protein1, protein2, .keep_all = TRUE)
      
      plot_df <- symmetric_interactions %>%
        mutate(x = gene_map[protein1], y = gene_map[protein2], score = combined_score) %>%
        dplyr::select(x, y, score)
      
      # Export plotting data
      output_csv_path <- file.path(plot_dir, 
                                    paste0(df_name, "_", group_name, "_", model_name, 
                                           "_CandidatePPI_PlotData.csv"))
      write.csv(plot_df, output_csv_path, row.names = FALSE)
      
      # Plot heatmap
      p_main <- ggplot(plot_df, aes(x = x, y = y, fill = score)) +
        geom_tile() +
        scale_fill_gradientn(colors = cfg$heatmap_color_palette, 
                            name = "STRING Score", 
                            limits = cfg$heatmap_score_range) +
        coord_fixed() + 
        theme_void() + 
        theme(legend.position = "right")
      
      # Add annotation bars
      # Compatible with both plot_group and group_Localization column names
      group_col <- if ("plot_group" %in% colnames(sorted_candidates)) {
        "plot_group"
      } else if ("group_Localization" %in% colnames(sorted_candidates)) {
        "group_Localization"
      } else {
        cat("    ⚠️  Grouping column not found, skipping annotation bars\n")
        NULL
      }
      
      if (!is.null(group_col)) {
        annot_df <- sorted_candidates %>% 
          mutate(index = 1:n(),
                 group_label = factor(.data[[group_col]], levels = c("Other", "TP")))
        
        p_top <- ggplot(annot_df, aes(x = index, y = 1, fill = group_label)) + 
          geom_tile() + 
          scale_fill_manual(values = TP_vs_other_colors) + 
          theme_void() + 
          theme(legend.position = "none")
        
        p_right <- ggplot(annot_df, aes(x = 1, y = index, fill = group_label)) + 
          geom_tile() + 
          scale_fill_manual(values = TP_vs_other_colors) + 
          theme_void() + 
          theme(legend.position = "none")
      } else {
        # If no grouping column, use blank spacer
        p_top <- plot_spacer()
        p_right <- plot_spacer()
      }
      
      p_final <- (p_top | plot_spacer()) / (p_main | p_right) + 
        plot_layout(widths = c(8, 1), heights = c(1, 8), guides = 'collect') +
        plot_annotation(
          title = paste(df_name, "-", group_name, "-", model_name),
          subtitle = paste("Candidate Proteins PPI Heatmap (n=", length(pass_proteins), ")")
        )
      
      output_pdf_path <- file.path(plot_dir, 
                                    paste0(df_name, "_", group_name, "_", model_name, 
                                           "_CandidatePPI_Heatmap.pdf"))
      ggsave(output_pdf_path, p_final, width = 9, height = 9)
      
      cat(paste("    ✓ Saved:", basename(output_pdf_path), "\n"))
      
    }, error = function(e) {
      cat(paste("    ❌ Model", model_name, "heatmap generation failed:", e$message, "\n"))
    })
  }
}

# 4h. Grid Tools and Core Matrix Calculation ####

#' @title Calculate PM Matrix with Lower-Left Corner Accumulation
#' @description Two-dimensional probability mass (PM) matrix: bin interaction points by normalized rank coordinates, where grid weight sum divided by total weight sum equals probability mass; and calculate lower-left corner multi-threshold accumulation.
#' @param node_attr Node attribute data frame (containing Gene, mean_group_medianNorm, etc.)
#' @param interactions Interaction data frame (protein1, protein2, combined_score)
#' @param bins Integer, number of grid cells (recommended 100)
#' @param cfg Global configuration (reads pm_lower_left_percents)
#' @param weight_mode Weight mode: none/combined_score_tp/combined_score/tp
#' @param candidate_genes Optional character vector, only count candidate subset
#' @return List: matrix_df(long format, x_bin, y_bin, prob), summary_df(lower-left accumulation), bins
computePmMatrixWithSummary <- function(node_attr, interactions, bins, cfg, weight_mode = "none", candidate_genes = NULL) {
  if (!is.null(candidate_genes)) {
    node_attr <- node_attr %>% dplyr::filter(Gene %in% candidate_genes)
  }
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes)
  if (n_genes < 2) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), prob = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }
  gene_map <- setNames(1:n_genes, sorted_genes)

  filtered <- interactions %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), prob = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }

  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)

  w <- getPmWeightsFromMode(sym, node_attr, weight_mode)
  if (length(w) == 0 || sum(w) <= 0) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), prob = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }

  df <- data.frame(x_bin = xb, y_bin = yb, w = w)
  agg <- df %>% group_by(y_bin, x_bin) %>% summarise(v = sum(w), .groups = "drop")
  # Fill full grid, missing cells set to 0, ensure background is non-empty
  full_grid <- expand.grid(x_bin = 1:bins, y_bin = 1:bins)
  agg <- full_grid %>% left_join(agg, by = c("x_bin", "y_bin")) %>% mutate(v = ifelse(is.na(v), 0, v))
  total <- sum(agg$v); if (total <= 0) total <- 1
  agg <- agg %>% mutate(prob = v / total)

  ll <- lapply(cfg$pm_lower_left_percents, function(p) {
    max_bin <- max(1L, floor(bins * p))
    s <- agg %>% dplyr::filter(x_bin <= max_bin & y_bin <= max_bin) %>% summarise(s = sum(prob)) %>% pull(s)
    data.frame(percent = p, lower_left_prob = ifelse(is.na(s), 0, s))
  }) %>% bind_rows()

  list(matrix_df = agg %>% dplyr::select(x_bin, y_bin, prob), summary_df = ll, bins = bins)
}

#' @title Calculate JP Matrix with Lower-Left Corner Accumulation (unweighted, directed)
#' @description Calculate probability matrix of observed PPI pairs to possible pairs ratio by binning, and provide lower-left corner accumulation.
#' @param node_attr Node attributes
#' @param interactions Interactions
#' @param bins Number of grid cells
#' @param candidate_genes Optional, candidate subset
#' @return List: matrix_df(long format, x_bin, y_bin, probability), summary_df(lower-left accumulation)
computeJpMatrixWithSummary <- function(node_attr, interactions, bins, candidate_genes = NULL) {
  if (!is.null(candidate_genes)) {
    node_attr <- node_attr %>% dplyr::filter(Gene %in% candidate_genes)
  }
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes)
  if (n_genes < 2) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), probability = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }
  gene_map <- setNames(1:n_genes, sorted_genes)

  filtered <- interactions %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), probability = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }

  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)

  # Observation count
  obs <- data.frame(x_bin = xb, y_bin = yb) %>% group_by(y_bin, x_bin) %>% summarise(cnt = n(), .groups = "drop")
  full_grid <- expand.grid(x_bin = 1:bins, y_bin = 1:bins)
  obs <- full_grid %>% left_join(obs, by = c("x_bin", "y_bin")) %>% mutate(cnt = ifelse(is.na(cnt), 0L, cnt))
  # Bin protein counts on axes (directed)
  genes_idx <- 1:n_genes
  cbin <- .binIndex(genes_idx / n_genes, bins)
  rbin <- .binIndex(genes_idx / n_genes, bins)
  nx <- as.integer(tabulate(cbin, nbins = bins))
  ny <- as.integer(tabulate(rbin, nbins = bins))

  obs <- obs %>% mutate(possible = nx[y_bin] * ny[x_bin],
                        probability = ifelse(possible > 0, cnt / possible, 0))

  ll <- lapply(pm_lower_left_percents, function(p) {
    max_bin <- max(1L, floor(bins * p))
    s <- obs %>% dplyr::filter(x_bin <= max_bin & y_bin <= max_bin) %>% summarise(s = sum(probability)) %>% pull(s)
    data.frame(percent = p, lower_left_prob = ifelse(is.na(s), 0, s))
  }) %>% bind_rows()

  list(matrix_df = obs %>% dplyr::select(x_bin, y_bin, probability), summary_df = ll, bins = bins)
}

#' @title Calculate Weighted JP Matrix with Lower-Left Corner Accumulation (directed)
#' @description Numerator = weighted sum within grid (weights normalized and powered); denominator = possible pairs (nx*ny); output probability matrix and lower-left accumulation.
#' @param node_attr Node attributes
#' @param interactions Interactions
#' @param bins Number of grid cells
#' @param weight_mode Weight mode: combined_score_tp/combined_score/tp
#' @param candidate_genes Optional candidate subset
computeJpWeightedMatrixWithSummary <- function(node_attr, interactions, bins, weight_mode = "combined_score_tp", candidate_genes = NULL) {
  if (!is.null(candidate_genes)) {
    node_attr <- node_attr %>% dplyr::filter(Gene %in% candidate_genes)
  }
  sorted_genes <- node_attr %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  n_genes <- length(sorted_genes)
  if (n_genes < 2) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), probability = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }
  gene_map <- setNames(1:n_genes, sorted_genes)

  filtered <- interactions %>% dplyr::filter(protein1 %in% sorted_genes & protein2 %in% sorted_genes)
  sym <- bind_rows(filtered, filtered %>% rename(protein1 = protein2, protein2 = protein1)) %>% dplyr::distinct()
  if (nrow(sym) < 1) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), probability = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }

  coords <- sym %>% mutate(row = gene_map[protein1], col = gene_map[protein2]) %>% dplyr::select(row, col) %>% as.matrix()
  row01 <- coords[,1] / n_genes; col01 <- coords[,2] / n_genes
  xb <- .binIndex(col01, bins); yb <- .binIndex(row01, bins)

  w <- getJpNormalizedWeights(sym, node_attr, weight_mode)
  if (length(w) == 0) {
    return(list(
      matrix_df = data.frame(x_bin = integer(), y_bin = integer(), probability = numeric()),
      summary_df = data.frame(percent = numeric(), lower_left_prob = numeric()),
      bins = bins
    ))
  }

  obs <- data.frame(x_bin = xb, y_bin = yb, w = w) %>% group_by(y_bin, x_bin) %>% summarise(ws = sum(w), .groups = "drop")
  full_grid <- expand.grid(x_bin = 1:bins, y_bin = 1:bins)
  obs <- full_grid %>% left_join(obs, by = c("x_bin", "y_bin")) %>% mutate(ws = ifelse(is.na(ws), 0, ws))

  genes_idx <- 1:n_genes
  cbin <- .binIndex(genes_idx / n_genes, bins)
  rbin <- .binIndex(genes_idx / n_genes, bins)
  nx <- as.integer(tabulate(cbin, nbins = bins))
  ny <- as.integer(tabulate(rbin, nbins = bins))

  obs <- obs %>% mutate(possible = nx[y_bin] * ny[x_bin], probability = ifelse(possible > 0, ws / possible, 0))

  ll <- lapply(pm_lower_left_percents, function(p) {
    max_bin <- max(1L, floor(bins * p))
    s <- obs %>% dplyr::filter(x_bin <= max_bin & y_bin <= max_bin) %>% summarise(s = sum(probability)) %>% pull(s)
    data.frame(percent = p, lower_left_prob = ifelse(is.na(s), 0, s))
  }) %>% bind_rows()

  list(matrix_df = obs %>% dplyr::select(x_bin, y_bin, probability), summary_df = ll, bins = bins)
}

# 4i. PM/JP Plotting Functions (All) ####

#' @title JP Heatmap (All, Unweighted, Directed)
plotJointProbHeatmap <- function(result, scales_config, cfg, plot_dir) {
  if (is.null(result)) return()
  df_name <- result$df_name; group_name <- result$group_name
  bins <- cfg$jp_grid_bins_fixed
  node_attr <- result$node_attributes; interactions <- result$interactions

  obj <- computeJpMatrixWithSummary(node_attr, interactions, bins, candidate_genes = NULL)
  mat_df <- obj$matrix_df; ll_df <- obj$summary_df
  csv_mat <- file.path(plot_dir, paste0(df_name, "_", group_name, "_JP_Probability_Matrix.csv"))
  write.csv(mat_df, csv_mat, row.names = FALSE)
  csv_ll <- file.path(plot_dir, paste0(df_name, "_", group_name, "_JP_LowerLeftSummary.csv"))
  write.csv(ll_df, csv_ll, row.names = FALSE)

  vmax <- scales_config$jp_max_prob
  vmin <- 0
  if (!is.null(scales_config$jp_min_prob)) {
    vmin <- max(0, as.numeric(scales_config$jp_min_prob))
  }
  p <- ggplot(mat_df, aes(x = x_bin, y = y_bin)) +
    geom_raster(aes(fill = probability)) +
    scale_fill_gradientn(colors = cfg$density_color_palette, limits = c(vmin, vmax), oob = scales::squish, name = "Probability") +
    coord_fixed() +
    labs(title = paste(df_name, "-", group_name), subtitle = "Joint Probability Heatmap (directed)",
         x = "Column Bins", y = "Row Bins")
  pdf_path <- file.path(plot_dir, paste0(df_name, "_", group_name, "_JP.pdf"))
  ggsave(pdf_path, p, width = 9, height = 8)
}

# 4l. JP Weighted Heatmap (All) ####

#' @title JP Weighted Heatmap (All, directed)
plotJointProbWeightedHeatmap <- function(result, scales_config, cfg, plot_dir, weight_mode = "combined_score_tp") {
  if (is.null(result)) return()
  df_name <- result$df_name; group_name <- result$group_name
  bins <- cfg$jp_grid_bins_fixed
  node_attr <- result$node_attributes; interactions <- result$interactions

  obj <- computeJpWeightedMatrixWithSummary(node_attr, interactions, bins, weight_mode = weight_mode, candidate_genes = NULL)
  mat_df <- obj$matrix_df; ll_df <- obj$summary_df
  mode_tag <- weight_mode
  csv_mat <- file.path(plot_dir, paste0(df_name, "_", group_name, "_JPw_", mode_tag, "_Probability_Matrix.csv"))
  write.csv(mat_df, csv_mat, row.names = FALSE)
  csv_ll <- file.path(plot_dir, paste0(df_name, "_", group_name, "_JPw_", mode_tag, "_LowerLeftSummary.csv"))
  write.csv(ll_df, csv_ll, row.names = FALSE)

  vmax <- scales_config$jpw_max_prob_by_mode[[weight_mode]]
  vmin <- 0
  if (!is.null(scales_config$jpw_min_prob_by_mode) && !is.null(scales_config$jpw_min_prob_by_mode[[weight_mode]])) {
    vmin <- max(0, as.numeric(scales_config$jpw_min_prob_by_mode[[weight_mode]]))
  }
  p <- ggplot(mat_df, aes(x = x_bin, y = y_bin)) +
    geom_raster(aes(fill = probability)) +
    scale_fill_gradientn(colors = cfg$density_color_palette, limits = c(vmin, vmax), oob = scales::squish, name = "Weighted Probability") +
    coord_fixed() +
    labs(title = paste(df_name, "-", group_name), subtitle = paste0("Weighted JP (", weight_mode, ")"),
         x = "Column Bins", y = "Row Bins")
  pdf_path <- file.path(plot_dir, paste0(df_name, "_", group_name, "_JPw_", mode_tag, ".pdf"))
  ggsave(pdf_path, p, width = 9, height = 8)
}


# -------------------------------------------------------------------
# 5. Main Execution Loop
# -------------------------------------------------------------------
# (The main loop is updated to pass the new 'cfg' object to functions)
cat("\nStarting to generate all visualization plots...\n")

# Create the config object to pass to functions
config_params <- list(
  roc_line_color = roc_line_color,
  youden_line_color = youden_line_color,
  youden_y_axis_range = youden_y_axis_range,
  scatter_3d_color_palette = scatter_3d_color_palette,
  heatmap_color_palette = heatmap_color_palette,
  heatmap_score_range = heatmap_score_range,
  density_kmeans_k = density_kmeans_k,
  density_color_palette = density_color_palette,
  enable_weighted_density = enable_weighted_density,
  weighted_density_sample_size = weighted_density_sample_size,
  weighted_density_suffix = weighted_density_suffix,
  density_coord_mode = density_coord_mode,
  density_bandwidth = density_bandwidth,
  density_grid_size = density_grid_size,
  axis_range_abundance = axis_range_abundance,
  axis_range_logfc = axis_range_logfc,
  axis_range_ppi_global = axis_range_ppi_global,
  axis_range_ppi_local = axis_range_ppi_local,
  weighted_density_mode = weighted_density_mode,
  tp_weight_both = tp_weight_both,
  tp_weight_single = tp_weight_single,
  tp_weight_none = tp_weight_none,
  # PM/JP Configuration
  pm_enable = pm_enable,
  pm_grid_mode = pm_grid_mode,
  pm_grid_bins_fixed = pm_grid_bins_fixed,
  pm_lower_left_percents = pm_lower_left_percents,
  pm_weight_modes = pm_weight_modes,
  pm_weight_power = pm_weight_power,
  pm_min_quantile = pm_min_quantile,
  pm_candidates_min_quantile = pm_candidates_min_quantile,
  pm_max_quantile = pm_max_quantile,
  pm_candidates_max_quantile = pm_candidates_max_quantile,
  jp_enable = jp_enable,
  jp_pair_mode = jp_pair_mode,
  jp_grid_mode = jp_grid_mode,
  jp_grid_bins_fixed = jp_grid_bins_fixed,
  jp_min_quantile = jp_min_quantile,
  jp_candidates_min_quantile = jp_candidates_min_quantile,
  jp_max_quantile = jp_max_quantile,
  jp_candidates_max_quantile = jp_candidates_max_quantile,
  # JP weighted
  jp_weight_enable = jp_weight_enable,
  jp_weight_modes = jp_weight_modes,
  jp_weight_power = jp_weight_power,
  jpw_min_quantile = jpw_min_quantile,
  jpw_candidates_min_quantile = jpw_candidates_min_quantile,
  jpw_max_quantile = jpw_max_quantile,
  jpw_candidates_max_quantile = jpw_candidates_max_quantile
)

for (i in seq_along(all_results)) {
  result <- all_results[[i]]
  if (is.null(result)) {
    cat(paste("Skipping plotting for task", i, "due to missing results.\n"))
    next
  }
  cat(paste("Generating plots for", result$df_name, "-", result$group_name, "...\n"))
  
  tryCatch({
    plot_model_performance(result, file.path(output_dir_plots, "1_ROC_Curves"), file.path(output_dir_plots, "2_Youden_Index"), cfg = config_params)
    plot_decision_boundaries(result, plot_scales_config, file.path(output_dir_plots, "3_Decision_Boundaries"), cfg = config_params)
    plot_3d_scatter(result, plot_scales_config, cfg = config_params, file.path(output_dir_plots, "4_3D_Like_Scatters"))
    plot_ppi_heatmap(result, cfg = config_params, file.path(output_dir_plots, "5_PPI_Heatmaps"))
    plot_candidate_ppi_heatmaps(result, cfg = config_params, file.path(output_dir_plots, "7_Candidate_PPI_Heatmaps"))

    # 13: JP All (Unweighted)
    if (isTRUE(config_params$jp_enable)) {
      plotJointProbHeatmap(result, plot_scales_config, cfg = config_params, plot_dir = file.path(output_dir_plots, "13_JointProb_Heatmaps"))
    }

    # 15: JP All (Weighted, Multi-mode)
    if (isTRUE(config_params$jp_weight_enable)) {
      for (mode in config_params$jp_weight_modes) {
        plotJointProbWeightedHeatmap(result, plot_scales_config, cfg = config_params, plot_dir = file.path(output_dir_plots, "15_JointProb_Heatmaps_Weighted"), weight_mode = mode)
      }
    }
  }, error = function(e){
    cat("ERROR: Plotting failed for", result$df_name, "-", result$group_name, "\nMessage:", conditionMessage(e), "\n")
  })
}

cat("\nAll plots have been generated successfully.\n")



# Step29 Phase3 Different Models FinalList and Comparision ----------------

# ===================================================================
#
### Phase 3: Candidate Protein Set Filtering, Summarization and Visualization Comparison (v3.0 - Advanced Configuration)
#
# ===================================================================

# -------------------------------------------------------------------
# 1. Environment Setup and Data Loading
# -------------------------------------------------------------------

# Ensure required R packages are loaded
# install.packages(c("dplyr", "tidyr", "ggplot2", "openxlsx", "purrr", "ggrepel"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(purrr)
library(ggrepel)
library(rlang)

# Define functional error handler
my_error_handler <- function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
  return(NULL) # Return NULL for subsequent processing
}

# --- Load Phase 1 analysis results ---
tryCatch({
  load("Step29_Analysis_Phase1_DataExport/analysis_phase1_output.RData")
  cat("Phase 1 analysis results 'analysis_phase1_output.RData' loaded successfully.\n")
}, error = function(e) {
  stop("Failed to load Phase 1 RData file. Please ensure the file path 'Step29_Analysis_Phase1_DataExport/analysis_phase1_output.RData' is correct.")
})

# --- Create output directory ---
output_dir_phase3 <- "Step29_Analysis_Phase3_Comparison"
dir.create(output_dir_phase3, showWarnings = FALSE)

# ===================================================================
#
### Section A: User Configuration Area
#
# ===================================================================

# 1. Set thresholds for traditional filtering methods
PARAM_FDR <- 0.05
PARAM_FC <- 0.5

# 2. Select models to include in plots
#    All model names can be obtained from all_results[[1]]$auc_scores
#    "logFC>0.5 & FDR<0.05" and "Control" are reserved names for traditional methods and control groups
#    Note: These names must match the model names in all_results[[i]]$candidates_full_list
#    Available model names include:
#    "1D_Abundance", "1D_logFC", "1D_PPI_Global","1D_PPI_Local","2D_Abundance_logFC",
#    "2D_Abundance_PPI_Global", "2D_Abundance_PPI_Local","2D_logFC_PPI_Global",
#    "2D_logFC_PPI_Local", "3D_with_Global_PPI", "3D_with_Local_PPI",
#    "logFC>0.5 & FDR<0.05","Control"

# names(all_results[[1]]$candidates_full_list)

# Comment out models that are not needed #
PARAM_MODELS_TO_PLOT <- c(
  "1D_Abundance", 
  "1D_logFC", 
  #"1D_PPI_Global", 
  #"1D_PPI_Local",
  "2D_Abundance_logFC", 
  #"2D_Abundance_PPI_Global", 
  #"2D_Abundance_PPI_Local",
  #"2D_logFC_PPI_Global", 
  #"2D_logFC_PPI_Local", 
  #"3D_with_Global_PPI",
  #"3D_with_Local_PPI", 
  "logFC>0.5 & FDR<0.05", 
  "Control"
)

# 3. Select column for grouping and color filling
#    Available values: "MultiBait_Localization", "HaloMap_Localization", "GO_Localization", "group_Localization"

PARAM_GROUPING_COLUMN <- "group_Localization"

# 4. Select category for annotation on plots (must be a value from PARAM_GROUPING_COLUMN)
PARAM_ANNOTATION_CATEGORY <- "TP"

# 5. [To be confirmed and modified by user] Define color scheme
#    Note: This will be checked in the next step. After seeing the category information output by the script, ensure that the naming here exactly matches the categories
PARAM_LOCALIZATION_COLORS <- c(
  "TP" = "#e41a1c",
  "Other" = "#bdbdbd"
  # Example: If using MultiBait_Localization, you need to provide the complete color scheme as follows
  # "Cytosol&SGs"="#DB6968", "Nuclear&Cytosol&SGs"="#EFCA72",
  # "Nuclear&SGs"="#4D97CD", "SGs&Other"="#B379B4", "Other"="grey89",
  # "Mitochondrion"="grey89", "Nuclear"="grey89", "Nuclear_Cytosol"="grey89",
  # "Cytosol"="grey89"
)

# ===================================================================
#
### Section B: Data Preprocessing and Configuration Check
#
# ===================================================================


# --- 1. Data Preprocessing (Logic Optimization) ---
# a. Add group_Localization column to all_results object
for (i in seq_along(all_results)) {
  if (!is.null(all_results[[i]])) {
    if ("plot_group" %in% colnames(all_results[[i]]$node_attributes)) {
      all_results[[i]]$node_attributes <- all_results[[i]]$node_attributes %>%
        rename(group_Localization = plot_group)
    }
    group_localization_mapping <- all_results[[i]]$node_attributes %>% dplyr::select(Gene, group_Localization)
    if (nrow(group_localization_mapping) > 0) {
      for (model_name in names(all_results[[i]]$candidates_full_list)) {
        all_results[[i]]$candidates_full_list[[model_name]] <- all_results[[i]]$candidates_full_list[[model_name]] %>%
          dplyr::select(-any_of("group_Localization")) %>%
          left_join(group_localization_mapping, by = c("protein_id" = "Gene"))
      }
    }
  }
}

names(all_results[[1]]$node_attributes)
unique(all_results[[1]]$node_attributes$group_Localization)

cat("Preprocessing 1/2 completed: Added 'group_Localization' column to 'all_results'.\n")

# b. Add group_Localization column to each data frame in original data ForStep19
for (i in seq_along(ForStep19)) {
  if ("MultiBait_Localization" %in% names(ForStep19[[i]])) {
    ForStep19[[i]] <- ForStep19[[i]] %>%
      mutate(group_Localization = if_else(
        MultiBait_Localization %in% myTP_vector, 
        "TP", "Other"
      ))
  }
}
cat("Preprocessing 2/2 completed: Added 'group_Localization' column to all data frames in 'ForStep19'.\n")
names(ForStep19$noMBR_Local_QNorm_New_B)
unique(ForStep19$noMBR_Local_QNorm_New_B$group_Localization)


# --- 2. Dynamic Color Configuration Check ---
cat("\n--- Checking color configuration ---\n")
detected_levels_from_results <- unlist(lapply(all_results, function(res) {
  if(!is.null(res) && PARAM_GROUPING_COLUMN %in% colnames(res$node_attributes)) {
    unique(as.character(res$node_attributes[[PARAM_GROUPING_COLUMN]]))
  }
}))
levels_from_forstep19 <- unlist(lapply(ForStep19, function(df) {
  if (PARAM_GROUPING_COLUMN %in% colnames(df)) {
    unique(as.character(df[[PARAM_GROUPING_COLUMN]]))
  }
}))
detected_levels <- unique(c(detected_levels_from_results, levels_from_forstep19))
detected_levels <- detected_levels[!is.na(detected_levels)]

cat(paste("The following categories were detected in your selected grouping column '", PARAM_GROUPING_COLUMN, "':\n", sep=""))
cat(paste(detected_levels, collapse = ", "), "\n")

defined_colors <- names(PARAM_LOCALIZATION_COLORS)
if (!all(detected_levels %in% defined_colors)) {
  missing_levels <- setdiff(detected_levels, defined_colors)
  stop(paste("Color configuration error! The following categories exist in the data but do not have corresponding colors defined in 'PARAM_LOCALIZATION_COLORS':\n",
             paste(missing_levels, collapse = ", ")))
}
if (!all(defined_colors %in% detected_levels)) {
  extra_colors <- setdiff(defined_colors, detected_levels)
  warning(paste("Color configuration warning: The following colors are defined but no corresponding categories were found in the data:\n",
                paste(extra_colors, collapse = ", ")))
}
cat("Color configuration check passed.\n\n")


# ===================================================================
#
### Section C: Core Functions and Analysis Loop
#
# ===================================================================
#### Add a new function process_replicates for convenience in calculations ####

#' @title Quality control and summarization of replicate sample data by group
#' @description Automatically performs quality control on columns within a specified range to ensure they belong to the same group of replicate samples, then performs summarization calculations.
#' @param data A data frame.
#' @param group_col Column name for grouping (bare column name or string).
#' @param RepRange A numeric range specifying the positions of quantitative columns to be summarized, e.g., `2:4`.
#'
#' @return If QC passes, returns a new data frame with summarization calculations. If QC fails, throws an error.

process_replicates <- function(data, group_col, RepRange) {
  # --- 1. Quality Control (QC) Step ---
  # Check if RepRange exceeds data frame boundaries
  if (max(RepRange) > ncol(data) || min(RepRange) < 1) {
    stop("QC Error: RepRange exceeds the column range of the data frame.")
  }
  # Get column names within the specified range
  quant_cols <- colnames(data)[RepRange]
  # QC Check 1: Do all column names contain "LFQ"
  if (!all(grepl("LFQ", quant_cols))) {
    # Find column names that don't contain LFQ to provide clearer error message
    failed_cols <- quant_cols[!grepl("LFQ", quant_cols)]
    stop(paste("QC Error: Some column names in the specified range do not contain 'LFQ'. Failed columns:", 
               paste(failed_cols, collapse = ", ")))
  }
  # QC Check 2: Are the strings before "LFQ" all the same
  # Use sub() to extract the part before "LFQ"
  prefixes <- sub("_LFQ.*", "", quant_cols)
  
  if (length(unique(prefixes)) > 1) {
    stop(paste("QC Error: Sample prefixes (part before LFQ) of column names in the specified range are inconsistent. Detected prefixes:",
               paste(unique(prefixes), collapse = ", ")))
  }
  # QC passed
  message(paste("QC passed. Sample prefix:", unique(prefixes)[1]))
  # --- 2. Core Calculation Step ---
  summarised_data <- data %>%
    group_by({{ group_col }}) %>%
    summarise(
      count = n(),
      across(all_of(quant_cols), ~mean(.x, na.rm = TRUE), .names = "Mean_{.col}"),
      across(all_of(quant_cols), ~sum(2^.x, na.rm = TRUE), .names = "Sum_{.col}")
    ) %>%
    ungroup() %>%
    mutate(
      Percent = count / sum(count),
      across(starts_with("Sum_"), ~ .x / sum(.x), .names = "{.col}_Percent")
    ) %>%
    rowwise() %>%
    mutate(
      MeanSum = mean(c_across(ends_with("_Percent")), na.rm = TRUE)
    ) %>%
    ungroup()
  return(summarised_data)
}

final_plotting_data_list <- list()
detailed_candidate_list_all_sources <- list()

for (df_idx in seq_along(ForStep19)) {
  original_data <- ForStep19[[df_idx]]
  data_source_name <- names(ForStep19)[df_idx]
  cat(paste("\n========================================================\n"))
  cat(paste("Processing data source:", data_source_name, "\n"))
  cat(paste("========================================================\n\n"))
  
  all_summaries_list <- list()
  detailed_candidate_list_current_source <- list()
  
  method_name_trad <- paste0("logFC>", PARAM_FC, " & FDR<", PARAM_FDR)
  if (method_name_trad %in% PARAM_MODELS_TO_PLOT) {
    cat("--- (A) Processing traditional logFC/FDR filtering method...\n")
    for (group in names(group_info)) {
      tryCatch({
        # Use actual column names defined in group_info, not fixed NES format
        logfc_col <- group_info[[group]]$logFC; fdr_col <- group_info[[group]]$logFC_FDR; lfg_cols <- group_info[[group]]$samples
        if (!all(c(logfc_col, fdr_col) %in% colnames(original_data))) {
          cat(paste("    Warning: Experimental group", group, "lacks required logFC or FDR column, skipping traditional filtering.\n"))
          next
        }
        
        # Add debugging information
        cat(paste("    - Experimental group:", group, "\n"))
        cat(paste("      logFC column:", logfc_col, "\n"))
        cat(paste("      FDR column:", fdr_col, "\n"))
        
        candidates_step1 <- original_data %>%
          filter(.data[[logfc_col]] > PARAM_FC & .data[[fdr_col]] < PARAM_FDR)
        cat(paste("      First round filtering (logFC >", PARAM_FC, "& FDR <", PARAM_FDR, ") yielded", nrow(candidates_step1), "candidate proteins\n"))
        
        candidates_trad_df <- original_data %>%
          filter(.data[[logfc_col]] > PARAM_FC & .data[[fdr_col]] < PARAM_FDR) %>%
          filter(rowSums(!is.na(dplyr::select(., one_of(intersect(lfg_cols, names(.)))))) >= PARAM_MIN_LFQ_NONNA)
        
        if (nrow(candidates_trad_df) > 0) {
          list_name <- paste(group, "trad", sep="_")
          detailed_candidate_list_current_source[[list_name]] <- candidates_trad_df
          
          rep_range_dynamic <- which(colnames(candidates_trad_df) %in% lfg_cols)
          summary_trad <- process_replicates(candidates_trad_df, !!sym(PARAM_GROUPING_COLUMN), rep_range_dynamic)
          summary_trad$DataSource <- data_source_name; summary_trad$Group <- group; summary_trad$Method <- method_name_trad
          all_summaries_list[[list_name]] <- summary_trad
        }
      }, error = my_error_handler)
    }
  }
  
  if ("Control" %in% PARAM_MODELS_TO_PLOT) {
    cat("--- (B) Processing spatial control groups...\n")
    loc_cols <- colnames(original_data)[endsWith(colnames(original_data), "Localization")]
    for (control_name in names(spatial_control_groups)) {
      tryCatch({
        control_cfg <- spatial_control_groups[[control_name]]
        pattern <- control_cfg$pattern
        display_name <- control_cfg$display_name
        if (is.null(pattern) || pattern == "") {
          cat(paste0("    - Control group ", control_name, " has no column name matching pattern provided, skipping.\n"))
          next
        }
        if (is.null(display_name) || is.na(display_name)) {
          display_name <- control_name
        }
        control_cols <- colnames(original_data)[grepl(pattern, colnames(original_data))]
        if (length(control_cols) == 0) {
          cat(paste0("    - Control group ", control_name, ": Matching pattern '", pattern, "' not found in column names.\n"))
          next
        }
        candidates_control <- original_data %>%
          dplyr::select(Gene, one_of(control_cols), one_of(loc_cols), one_of("group_Localization")) %>%
          na.omit()
        if (nrow(candidates_control) == 0) {
          cat(paste0("    - Control group ", control_name, ": Matching columns exist but all are NA, skipping.\n"))
          next
        }
        detailed_name <- paste0(control_name, "_Control")
        detailed_candidate_list_current_source[[detailed_name]] <- candidates_control
        rep_range <- which(colnames(candidates_control) %in% control_cols)
        summary_control <- process_replicates(candidates_control, !!sym(PARAM_GROUPING_COLUMN), rep_range)
        summary_control$DataSource <- data_source_name
        summary_control$Group <- display_name
        summary_control$Method <- "Control"
        all_summaries_list[[paste0("control_", control_name)]] <- summary_control
      }, error = my_error_handler)
    }
  }
  
  cat("--- (C) Processing filtering results for all AUC models...\n")
  results_for_df <- all_results[sapply(all_results, function(res) !is.null(res) && res$df_name == data_source_name)]
  if (length(results_for_df) > 0) {
    for (res in results_for_df) {
      group_name <- res$group_name
      lfg_cols <- group_info[[group_name]]$samples
      cat(paste("  - Processing experimental group:", group_name, "\n"))
      
      for (model_name in names(res$candidates_full_list)) {
        if (model_name %in% PARAM_MODELS_TO_PLOT) {
          tryCatch({
            model_candidates_df <- res$candidates_full_list[[model_name]]
            pass_genes <- model_candidates_df %>% filter(PassThreshold == "Pass") %>% pull(protein_id)
            
            if (length(pass_genes) > 0) {
              full_candidates_info <- res$node_attributes %>% filter(Gene %in% pass_genes)
              
              if(nrow(full_candidates_info) > 0) {
                list_name <- paste(group_name, model_name, sep="_")
                detailed_candidate_list_current_source[[list_name]] <- full_candidates_info
                
                rep_range_dynamic <- which(colnames(full_candidates_info) %in% lfg_cols)
                summary_model <- process_replicates(full_candidates_info, !!sym(PARAM_GROUPING_COLUMN), rep_range_dynamic)
                summary_model$DataSource<-data_source_name; summary_model$Group<-group_name; summary_model$Method<-model_name
                all_summaries_list[[list_name]] <- summary_model
              }
            }
          }, error = my_error_handler)
        }
      }
    }
  }
  
  if (length(all_summaries_list) > 0) {
    final_plotting_data <- bind_rows(all_summaries_list) %>%
      dplyr::select(DataSource, Group, Method, !!sym(PARAM_GROUPING_COLUMN), count, Percent, MeanSum)
    final_plotting_data_list[[data_source_name]] <- final_plotting_data
    cat(paste("\n--- All results for data source '", data_source_name, "' have been summarized. ---\n", sep=""))
  }
  if (length(detailed_candidate_list_current_source) > 0) {
    detailed_candidate_list_all_sources[[data_source_name]] <- detailed_candidate_list_current_source
  }
}

# ===================================================================
#
### Section D: Visualization and File Export
#
# ===================================================================
# ===================================================================
#
### Section D: Visualization and File Export (Corrected Version)
#
# ===================================================================
cat("\n--- Starting plot generation and data export ---\n")
for (data_source_name in names(detailed_candidate_list_all_sources)) {
  wb <- createWorkbook()
  sheet_data_list <- detailed_candidate_list_all_sources[[data_source_name]]
  for(sheet_name in names(sheet_data_list)){
    safe_sheet_name <- substr(sheet_name, 1, 31)
    addWorksheet(wb, safe_sheet_name)
    writeData(wb, safe_sheet_name, sheet_data_list[[sheet_name]])
  }
  output_xlsx_path <- file.path(output_dir_phase3, paste0("Detailed_Candidates_", data_source_name, ".xlsx"))
  tryCatch({
    saveWorkbook(wb, output_xlsx_path, overwrite = TRUE)
    cat(paste("Detailed protein list saved to:", output_xlsx_path, "\n"))
  }, error = function(e){cat("Excel file save failed:", conditionMessage(e), "\n")})
}
 
for (data_source_name in names(final_plotting_data_list)) {
  plot_df <- final_plotting_data_list[[data_source_name]]
  if (is.null(plot_df) || nrow(plot_df) == 0) {next}
  
  output_csv_path <- file.path(output_dir_phase3, paste0("PlottingData_Comparison_", data_source_name, ".csv"))
  write.csv(plot_df, output_csv_path, row.names = FALSE)
  
  # Prepare plot annotation data
  label_df <- plot_df %>%
    filter(!!sym(PARAM_GROUPING_COLUMN) == PARAM_ANNOTATION_CATEGORY) %>%
    group_by(DataSource, Group, Method) %>%
    filter(n() > 0) %>%
    summarise(
      label_text = paste0(
        PARAM_ANNOTATION_CATEGORY, ": ", count, "\n",
        "(", round(Percent * 100), "%)"
      ),
      .groups = 'drop'
    ) %>%
    distinct()
  
  # Plot count percentage chart
  p_count <- ggplot(plot_df, aes(x = Method, y = Percent, fill = !!sym(PARAM_GROUPING_COLUMN))) +
    geom_col(position = "stack") +
    # Correction: add x = Method in aes()
    geom_text_repel(
      data = label_df, 
      aes(x = Method, label = label_text, y = 1.0), # <--- Key correction point
      inherit.aes = FALSE,
      nudge_y = 0.05, direction = "x", size = 3.5, 
      min.segment.length = 0, box.padding = 0.5
    ) +
    facet_wrap(~ Group, scales = "free_x", ncol = 4) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = PARAM_LOCALIZATION_COLORS, name = "Localization", drop = FALSE) +
    labs(
      title = paste("Localization by Protein Count Percentage:", data_source_name),
      subtitle = "Comparison of different filtering methods across experimental groups",
      x = "Filtering Method", y = "Percentage of Proteins"
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 10), legend.position = "bottom")
  
  output_pdf_count_path <- file.path(output_dir_phase3, paste0("Plot_CountPercent_Comparison_", data_source_name, ".pdf"))
  ggsave(output_pdf_count_path, p_count, width = 20, height = 12)
  cat(paste("Count percentage chart saved to:", output_pdf_count_path, "\n"))
  
  # Prepare annotation for abundance chart
  label_df_abundance <- plot_df %>%
    filter(!!sym(PARAM_GROUPING_COLUMN) == PARAM_ANNOTATION_CATEGORY) %>%
    group_by(DataSource, Group, Method) %>%
    filter(n() > 0) %>%
    summarise(
      label_text = paste0(
        PARAM_ANNOTATION_CATEGORY, ": ", count, "\n",
        "(", round(MeanSum * 100), "%)"
      ),
      .groups = 'drop'
    ) %>%
    distinct()
  
  # Plot abundance percentage chart
  p_abundance <- ggplot(plot_df, aes(x = Method, y = MeanSum, fill = !!sym(PARAM_GROUPING_COLUMN))) +
    geom_col(position = "stack") +
    # Correction: add x = Method in aes()
    geom_text_repel(
      data = label_df_abundance, 
      aes(x = Method, label = label_text, y = 1.0), # <--- Key correction point
      inherit.aes = FALSE,
      nudge_y = 0.05, direction = "x", size = 3.5,
      min.segment.length = 0, box.padding = 0.5
    ) +
    facet_wrap(~ Group, scales = "free_x", ncol = 4) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = PARAM_LOCALIZATION_COLORS, name = "Localization", drop = FALSE) +
    labs(
      title = paste("Localization by Protein Abundance Percentage:", data_source_name),
      subtitle = "Comparison of different filtering methods across experimental groups",
      x = "Filtering Method", y = "Percentage of Total Abundance"
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 10), legend.position = "bottom")
  
  output_pdf_abundance_path <- file.path(output_dir_phase3, paste0("Plot_AbundancePercent_Comparison_", data_source_name, ".pdf"))
  ggsave(output_pdf_abundance_path, p_abundance, width = 20, height = 12)
  cat(paste("Abundance percentage chart saved to:", output_pdf_abundance_path, "\n"))
}

cat("\n--- Phase 3 analysis completed. ---\n")


# ===================================================================
#
### Phase 3B: MultiBait_Localization Stacked Bar Charts (Figure 7 Supplement) ####
#
# ===================================================================

cat("\n--- Starting MultiBait_Localization stacked bar charts (Figure 7) generation ---\n")

# MultiBait colors and stacking order are configured in MULTIBAIT_* interfaces at the beginning of the script

# Reprocess data using MultiBait_Localization
final_plotting_data_multibait_list <- list()

for (df_idx in seq_along(ForStep19)) {
  original_data <- ForStep19[[df_idx]]
  data_source_name <- names(ForStep19)[df_idx]
  cat(paste("\nProcessing data source (MultiBait):", data_source_name, "\n"))
  
  all_summaries_list <- list()
  
  # (A) Traditional method
  method_name_trad <- paste0("logFC>", PARAM_FC, " & FDR<", PARAM_FDR)
  if (method_name_trad %in% PARAM_MODELS_TO_PLOT) {
    for (group in names(group_info)) {
      tryCatch({
        logfc_col <- group_info[[group]]$logFC
        fdr_col <- group_info[[group]]$logFC_FDR
        lfg_cols <- group_info[[group]]$samples
        
        if (!all(c(logfc_col, fdr_col) %in% colnames(original_data))) {
          next
        }
        
        candidates_trad_df <- original_data %>%
          filter(.data[[logfc_col]] > PARAM_FC & .data[[fdr_col]] < PARAM_FDR) %>%
          filter(rowSums(!is.na(dplyr::select(., one_of(intersect(lfg_cols, names(.)))))) >= PARAM_MIN_LFQ_NONNA)
        
        if (nrow(candidates_trad_df) > 0 && "MultiBait_Localization" %in% colnames(candidates_trad_df)) {
          rep_range_dynamic <- which(colnames(candidates_trad_df) %in% lfg_cols)
          summary_trad <- process_replicates(candidates_trad_df, MultiBait_Localization, rep_range_dynamic)
          summary_trad$DataSource <- data_source_name
          summary_trad$Group <- group
          summary_trad$Method <- method_name_trad
          all_summaries_list[[paste(group, "trad", sep="_")]] <- summary_trad
        }
      }, error = my_error_handler)
    }
  }
  
  # (B) Control group
  if ("Control" %in% PARAM_MODELS_TO_PLOT) {
    if (!"MultiBait_Localization" %in% colnames(original_data)) {
      cat("  - Warning: Current data source lacks MultiBait_Localization column, skipping control group processing.\n")
    } else {
      for (control_name in names(spatial_control_groups)) {
        tryCatch({
          control_cfg <- spatial_control_groups[[control_name]]
          pattern <- control_cfg$pattern
          display_name <- control_cfg$display_name
          if (is.null(pattern) || pattern == "") {
            next
          }
          if (is.null(display_name) || is.na(display_name)) {
            display_name <- control_name
          }
          control_cols <- colnames(original_data)[grepl(pattern, colnames(original_data))]
          if (length(control_cols) == 0) {
            next
          }
          candidates_control <- original_data %>%
            dplyr::select(Gene, one_of(control_cols), MultiBait_Localization) %>%
            na.omit()
          if (nrow(candidates_control) == 0) {
            next
          }
          rep_range <- which(colnames(candidates_control) %in% control_cols)
          summary_control <- process_replicates(candidates_control, MultiBait_Localization, rep_range)
          summary_control$DataSource <- data_source_name
          summary_control$Group <- display_name
          summary_control$Method <- "Control"
          all_summaries_list[[paste0("control_", control_name)]] <- summary_control
        }, error = my_error_handler)
      }
    }
  }
  
  # (C) AUC models
  results_for_df <- all_results[sapply(all_results, function(res) !is.null(res) && res$df_name == data_source_name)]
  if (length(results_for_df) > 0) {
    for (res in results_for_df) {
      group_name <- res$group_name
      lfg_cols <- group_info[[group_name]]$samples
      
      for (model_name in names(res$candidates_full_list)) {
        if (model_name %in% PARAM_MODELS_TO_PLOT) {
          tryCatch({
            model_candidates_df <- res$candidates_full_list[[model_name]]
            pass_genes <- model_candidates_df %>% filter(PassThreshold == "Pass") %>% pull(protein_id)
            
            if (length(pass_genes) > 0) {
              full_candidates_info <- res$node_attributes %>% filter(Gene %in% pass_genes)
              
              if (nrow(full_candidates_info) > 0 && "MultiBait_Localization" %in% colnames(full_candidates_info)) {
                rep_range_dynamic <- which(colnames(full_candidates_info) %in% lfg_cols)
                summary_model <- process_replicates(full_candidates_info, MultiBait_Localization, rep_range_dynamic)
                summary_model$DataSource <- data_source_name
                summary_model$Group <- group_name
                summary_model$Method <- model_name
                all_summaries_list[[paste(group_name, model_name, sep="_")]] <- summary_model
              }
            }
          }, error = my_error_handler)
        }
      }
    }
  }
  
  if (length(all_summaries_list) > 0) {
    final_plotting_data <- bind_rows(all_summaries_list) %>%
      dplyr::select(DataSource, Group, Method, MultiBait_Localization, count, Percent, MeanSum)
    final_plotting_data_multibait_list[[data_source_name]] <- final_plotting_data
  }
}

# Generate MultiBait_Localization stacked bar charts
output_dir_phase3b <- file.path(output_dir_phase3, "MultiBait_Plots")
dir.create(output_dir_phase3b, showWarnings = FALSE)

for (data_source_name in names(final_plotting_data_multibait_list)) {
  plot_df <- final_plotting_data_multibait_list[[data_source_name]]
  if (is.null(plot_df) || nrow(plot_df) == 0) { next }
  
  # Set factor order to control stacking order
  plot_df$MultiBait_Localization <- factor(
    plot_df$MultiBait_Localization,
    levels = MULTIBAIT_LEVEL_ORDER
  )
  
  # Export CSV data
  output_csv_path <- file.path(output_dir_phase3b, 
                                paste0("PlottingData_MultiBait_", data_source_name, ".csv"))
  write.csv(plot_df, output_csv_path, row.names = FALSE)
  cat(paste("MultiBait plotting data saved:", output_csv_path, "\n"))
  
  # Prepare annotation data (showing counts of main SGs-related categories)
  label_df <- plot_df %>%
    filter(MultiBait_Localization %in% MULTIBAIT_SG_CATEGORIES) %>%
    group_by(DataSource, Group, Method, MultiBait_Localization) %>%
    summarise(
      count = first(count),
      Percent = first(Percent),
      .groups = 'drop'
    ) %>%
    group_by(DataSource, Group, Method) %>%
    summarise(
      label_text = paste0(
        paste(MultiBait_Localization, ": ", count, collapse = "\n")
      ),
      .groups = 'drop'
    )
  
  # Plot count percentage chart
  p_count <- ggplot(plot_df, aes(x = Method, y = Percent, fill = MultiBait_Localization)) +
    geom_col(position = "stack") +
    facet_wrap(~ Group, scales = "free_x", ncol = 4) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = MULTIBAIT_COLORS, name = "Localization", drop = FALSE) +
    labs(
      title = paste("MultiBait Localization by Count Percentage:", data_source_name),
      subtitle = "Comparison of different filtering methods across experimental groups",
      x = "Filtering Method", y = "Percentage of Proteins"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 80, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    )
  
  output_pdf_count_path <- file.path(output_dir_phase3b, 
                                      paste0("Plot_MultiBait_CountPercent_", data_source_name, ".pdf"))
  ggsave(output_pdf_count_path, p_count, width = 20, height = 12)
  cat(paste("MultiBait count percentage chart saved:", output_pdf_count_path, "\n"))
  
  # Plot abundance percentage chart
  p_abundance <- ggplot(plot_df, aes(x = Method, y = MeanSum, fill = MultiBait_Localization)) +
    geom_col(position = "stack") +
    facet_wrap(~ Group, scales = "free_x", ncol = 4) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = MULTIBAIT_COLORS, name = "Localization", drop = FALSE) +
    labs(
      title = paste("MultiBait Localization by Abundance Percentage:", data_source_name),
      subtitle = "Comparison of different filtering methods across experimental groups",
      x = "Filtering Method", y = "Percentage of Total Abundance"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 80, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    )
  
  output_pdf_abundance_path <- file.path(output_dir_phase3b, 
                                         paste0("Plot_MultiBait_AbundancePercent_", data_source_name, ".pdf"))
  ggsave(output_pdf_abundance_path, p_abundance, width = 20, height = 12)
  cat(paste("MultiBait abundance percentage chart saved:", output_pdf_abundance_path, "\n"))
}

cat("\n--- MultiBait_Localization stacked bar charts (Figure 7) generation completed ---\n")


# ===================================================================
#
### Phase 3C: TP Percentage Curves (Figure 9) ####
#
# ===================================================================

cat("\n--- Starting TP Percentage Curves (Figure 9) generation ---\n")

# Create output directory
output_dir_phase3c <- file.path(output_dir_phase3, "TP_Percent_Curves")
dir.create(output_dir_phase3c, showWarnings = FALSE)

# Generate TP percentage curves for all data sources (grouped by model, comparing different proximity labeling methods)
for (df_idx in seq_along(ForStep19)) {
  original_data <- ForStep19[[df_idx]]
  data_source_name <- names(ForStep19)[df_idx]
  
  cat(paste("\nProcessing data source (TP percentage curves):", data_source_name, "\n"))
  
  # Store data for all models × experimental groups
  all_model_data <- list()
  
  # Iterate through each experimental group to collect data
  for (group in names(group_info)) {
    cat(paste("  - Collecting experimental group data:", group, "\n"))
    
    # Get logFC column and sample columns for this experimental group
    logfc_col <- group_info[[group]]$logFC
    fdr_col <- group_info[[group]]$logFC_FDR
    lfg_cols <- group_info[[group]]$samples
    
    # (A) TP percentage curve for traditional method
    method_name_trad <- paste0("logFC>", PARAM_FC, " & FDR<", PARAM_FDR)
    if (method_name_trad %in% PARAM_MODELS_TO_PLOT) {
      tryCatch({
        if (all(c(logfc_col, fdr_col) %in% colnames(original_data))) {
          # Sort and label TP/FP, while filtering out SOD1 protein
          candidates_sorted <- original_data %>%
            filter(.data[[logfc_col]] > PARAM_FC & .data[[fdr_col]] < PARAM_FDR) %>%
            filter(rowSums(!is.na(dplyr::select(., one_of(intersect(lfg_cols, names(.)))))) >= PARAM_MIN_LFQ_NONNA_TPR) %>%
            filter(Gene != "SOD1") %>%  # Filter out SOD1 protein
            arrange(desc(.data[[logfc_col]])) %>%  # Sort by logFC in descending order
            mutate(
              rank_pct = (row_number() / n()) * 100,
              is_TP = MultiBait_Localization %in% myTP_vector2
            )
          
          # Calculate cumulative TP and FP lists
          n_rows <- nrow(candidates_sorted)
          tp_lists <- character(n_rows)
          fp_lists <- character(n_rows)
          
          for (i in 1:n_rows) {
            current_genes <- candidates_sorted$Gene[1:i]
            current_is_tp <- candidates_sorted$is_TP[1:i]
            tp_lists[i] <- paste(current_genes[current_is_tp], collapse = ";")
            fp_lists[i] <- paste(current_genes[!current_is_tp], collapse = ";")
          }
          
          candidates_trad_df <- candidates_sorted %>%
            mutate(
              cumsum_TP = cumsum(is_TP),
              cumsum_total = row_number(),
              TP_Percent = (cumsum_TP / cumsum_total) * 100,
              TP_List = tp_lists,
              FP_List = fp_lists
            ) %>%
            dplyr::select(rank_pct, TP_Percent, TP_List, FP_List)
          
          if (nrow(candidates_trad_df) > 0) {
            candidates_trad_df$Method <- method_name_trad
            candidates_trad_df$Group <- group
            
            # Store in all_model_data (by model name)
            if (!method_name_trad %in% names(all_model_data)) {
              all_model_data[[method_name_trad]] <- list()
            }
            all_model_data[[method_name_trad]][[group]] <- candidates_trad_df
          }
        }
      }, error = my_error_handler)
    }
    
    # (B) TP percentage curve for AUC models
    results_for_group <- all_results[sapply(all_results, function(res) 
      !is.null(res) && res$df_name == data_source_name && res$group_name == group)]
    
    if (length(results_for_group) > 0) {
      res <- results_for_group[[1]]
      
      for (model_name in names(res$candidates_full_list)) {
        if (model_name %in% PARAM_MODELS_TO_PLOT) {
          tryCatch({
            model_candidates_df <- res$candidates_full_list[[model_name]]
            pass_genes <- model_candidates_df %>% 
              filter(PassThreshold == "Pass") %>% 
              pull(protein_id)
            
            if (length(pass_genes) > 0) {
              # Get complete information, first filter out SOD1 protein, then match pass_genes
              full_candidates_info <- res$node_attributes %>% 
                filter(Gene != "SOD1") %>%  # Filter out SOD1 protein first
                filter(Gene %in% pass_genes)
              
              # Select sorting column based on model type
              if (grepl("^1D_Abundance", model_name)) {
                sort_col <- "mean_group_medianNorm"
              } else if (grepl("^1D_logFC", model_name)) {
                sort_col <- logfc_col
              } else if (grepl("^1D_PPI", model_name)) {
                if (grepl("Global", model_name)) {
                  sort_col <- "log_ppi_score_global"
                } else {
                  sort_col <- "log_ppi_score_local"
                }
              } else {
                # 2D and 3D models use model_score
                model_score <- model_candidates_df %>% 
                  filter(PassThreshold == "Pass") %>%
                  dplyr::select(protein_id, model_score)
                
                full_candidates_info <- full_candidates_info %>%
                  left_join(model_score, by = c("Gene" = "protein_id"))
                sort_col <- "model_score"
              }
              
              # Calculate TP percentage curve
              # Sort and label TP/FP first
              candidates_sorted <- full_candidates_info %>%
                arrange(desc(.data[[sort_col]])) %>%
                mutate(
                  rank_pct = (row_number() / n()) * 100,
                  is_TP = MultiBait_Localization %in% myTP_vector2
                )
              
              # Calculate cumulative TP and FP lists
              n_rows <- nrow(candidates_sorted)
              tp_lists <- character(n_rows)
              fp_lists <- character(n_rows)
              
              for (i in 1:n_rows) {
                current_genes <- candidates_sorted$Gene[1:i]
                current_is_tp <- candidates_sorted$is_TP[1:i]
                tp_lists[i] <- paste(current_genes[current_is_tp], collapse = ";")
                fp_lists[i] <- paste(current_genes[!current_is_tp], collapse = ";")
              }
              
              tpr_df <- candidates_sorted %>%
                mutate(
                  cumsum_TP = cumsum(is_TP),
                  cumsum_total = row_number(),
                  TP_Percent = (cumsum_TP / cumsum_total) * 100,
                  TP_List = tp_lists,
                  FP_List = fp_lists
                ) %>%
                dplyr::select(rank_pct, TP_Percent, TP_List, FP_List)
              
              if (nrow(tpr_df) > 0) {
                tpr_df$Method <- model_name
                tpr_df$Group <- group
                
                # Store in all_model_data (by model name)
                if (!model_name %in% names(all_model_data)) {
                  all_model_data[[model_name]] <- list()
                }
                all_model_data[[model_name]][[group]] <- tpr_df
              }
            }
          }, error = my_error_handler)
        }
      }
    }
  }
  
  # Generate plots for each model (comparing different proximity labeling methods)
  cat(paste("\n  Generating TP percentage curves for each model...\n"))
  
  for (model_name in names(all_model_data)) {
    model_groups_data <- all_model_data[[model_name]]
    
    if (length(model_groups_data) > 0) {
      # Combine data from all experimental groups for this model
      model_combined <- bind_rows(model_groups_data)
      
      # Clean file name: remove Windows-disallowed characters < > : " / \ | ? *
      safe_model_name <- gsub("[<>:\"/\\|?*]", "_", model_name)
      safe_model_name <- gsub("&", "and", safe_model_name)
      
      # Export CSV data
      output_csv_path <- file.path(output_dir_phase3c, 
                                    paste0("TP_Percent_Curve_", data_source_name, "_", safe_model_name, ".csv"))
      write.csv(model_combined, output_csv_path, row.names = FALSE)
      
      # Plot TP percentage curve (different colors represent different experimental groups)
      p_tp <- ggplot(model_combined, aes(x = rank_pct, y = TP_Percent, color = Group)) +
        geom_line(size = 1) +
        geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", size = 0.5) +
        scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
        labs(
          title = paste("TP Percentage Comparison:", data_source_name),
          subtitle = paste("Model:", model_name, "| TP defined by MultiBait_Localization containing 'SG'"),
          x = "Rank Percentile (Top X%)",
          y = "TP Percentage in Top-X% (%)",
          color = "Proximity Labeling\nMethod"
        ) +
        theme_bw(base_size = 12) +
        theme(
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9)
        )
      
      output_pdf_path <- file.path(output_dir_phase3c, 
                                    paste0("TP_Percent_Curve_", data_source_name, "_", safe_model_name, ".pdf"))
      ggsave(output_pdf_path, p_tp, width = 10, height = 6)
      
      cat(paste("    ✓", model_name, "curve saved\n"))
    }
  }
}

cat("\n--- TP Percentage Curves (Figure 9) generation completed ---\n")


save.image(file = "Phase1toPhase3.RData")
print("Workspace successfully saved to Phase1toPhase3.RData")
