# Military Medical Academy SGs Sample Analysis
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
library(purrr)
library(openxlsx)
library(preprocessCore)
library(limma)
library(colourpicker)
library(RColorBrewer)
library(writexl)


#2. Working Directory and Data Loading ####
setwd("D:/Bioinfomatics/MS/CodeNorm_Nuc_Batch1")
#source("../202501017 3DAUC models base.R")
## Replace with any environment that can complete the main pipeline up to FC/Scatter plot
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
myTP_vector %in% unique(ForStep16[[1]]$MultiBait_Localization) #myTP_Vector needs to be redefined in the downstream process_group function as global variables cannot be inherited
ForStep19=ForStep16

str(ForStep19)


#Step29 Multi-Model Complex Requirements (Version) #####

# ===================================================================
#
###    Configuration Interface - User Optional Settings ####
#
# ===================================================================

#1. Group Selection Configuration ####
## Available groups: c("K69A1B3", "K69C3","C3")
## Set to NULL to use all groups, or specify a vector of specific groups
names(group_info)
selected_groups <- names(group_info)


#2. Model Selection Configuration ####
## Available model types:
## 1D models: "1D_Abundance", "1D_logFC", "1D_PPI_Global", "1D_PPI_Local"
## 2D models: "2D_Abundance_logFC", "2D_Abundance_PPI_Global", "2D_Abundance_PPI_Local", "2D_logFC_PPI_Global", "2D_logFC_PPI_Local"
## 3D models: "3D_with_Global_PPI", "3D_with_Local_PPI"
## Traditional method: "logFC>0.5 & FDR<0.05"
## Set to NULL to run all models, or specify a vector of specific models
selected_models <- c(
  "logFC>0.5 & FDR<0.05",
  "1D_logFC",
  "1D_Abundance",
  "2D_Abundance_logFC")  # Optional: c("1D_Abundance", "2D_Abundance_logFC") etc.


#3. Output and Visualization Control ####
## Whether to generate various files and plots
enable_cytoscape_export <- TRUE    # Whether to export Cytoscape files
enable_visualization <- TRUE       # Whether to generate visualization plots
enable_robustness_analysis <- TRUE # Whether to perform robustness analysis
enable_group_comparison <- TRUE    # Whether to perform inter-group comparison

## Visualization subplot control
enable_plot_auc_summary <- TRUE                # Figure 1: AUC_summary grouped forest plot
enable_plot_auc_summary_vs_logFC <- TRUE       # Figure 1b: Grouped forest plot vs 1D_logFC
enable_plot_group_comparison <- TRUE           # Figure 2: Cross-proximity method comparison (bar plot + heatmap)
enable_plot_intersection_dumbbell <- TRUE      # Figure 3: Sample intersection fairness dumbbell plot
enable_plot_robustness_scatter <- TRUE         # Figure 4: Performance-stability scatter plot
enable_plot_controlled_comparison <- TRUE      # Figure 5: Controlled variable comparison gain plot
enable_plot_localization_distribution <- FALSE  # Figure 6: Candidate protein localization distribution stacked bar plot
enable_plot_controlled_comparison_vs_logFC <- TRUE  # Figure 7: Controlled variable gain plot vs 1D_logFC

#4. Parallel Computing Configuration ####
## Parallel core number setting (recommended to keep 12, using half of 24 total cores)
parallel_cores <- 20

#4.1 Controlled Variable Comparison Plot Configuration ####
## Figures 5a-5c: Specify models to display
models_for_controlled_comparison <- c(
  "1D_Abundance",
  "1D_logFC", 
  "2D_Abundance_logFC"
  # Traditional model will automatically serve as reference, no need to list here
)

## Figures 5a-5c: Axis range settings
controlled_comparison_xlim <- c(-0.15, 0.5)  # x-axis range for Specificity/Sensitivity gain
controlled_comparison_ylim <- c(-0.15, 0.5)  # y-axis range for Sensitivity gain (for figure 5c)

## Figures 5a-5c: Bootstrap resampling number setting
n_bootstrap <- 500  # Optional: 200(quick test), 500(routine), 1000(publication), 2000+(strict)
                     # Note: 1000 × 6 groups × 3 models ≈ 36,000 resamplings, requires longer time

#5. Validation and Initialization Configuration ####
cat("=== Configuration Validation Started ===\n")

# Validate group selection
available_groups <- names(group_info)
if (is.null(selected_groups)) {
  selected_groups <- available_groups
  cat("✓ Using all available groups:", paste(selected_groups, collapse = ", "), "\n")
} else {
  invalid_groups <- setdiff(selected_groups, available_groups)
  if (length(invalid_groups) > 0) {
    stop("❌ Error: Invalid group selection - ", paste(invalid_groups, collapse = ", "), 
         "\n   Available groups: ", paste(available_groups, collapse = ", "))
  }
  cat("✓ Selected groups:", paste(selected_groups, collapse = ", "), "\n")
}

# Validate model selection
available_models <- c(
  "1D_Abundance", "1D_logFC", "1D_PPI_Global", "1D_PPI_Local",
  "2D_Abundance_logFC", "2D_Abundance_PPI_Global", "2D_Abundance_PPI_Local", 
  "2D_logFC_PPI_Global", "2D_logFC_PPI_Local",
  "3D_with_Global_PPI", "3D_with_Local_PPI",
  "logFC>0.5 & FDR<0.05"
)
if (is.null(selected_models)) {
  selected_models <- setdiff(available_models, "logFC>0.5 & FDR<0.05")  # Default: exclude traditional method
  cat("✓ Using all available models (excluding traditional method):", length(selected_models), "models\n")
} else {
  invalid_models <- setdiff(selected_models, available_models)
  if (length(invalid_models) > 0) {
    stop("❌ Error: Invalid model selection - ", paste(invalid_models, collapse = ", "), 
         "\n   Available models: ", paste(available_models, collapse = ", "))
  }
  cat("✓ Selected models:", length(selected_models), "models -", paste(selected_models, collapse = ", "), "\n")
}

# Validate controlled variable comparison configuration
if (!is.null(models_for_controlled_comparison)) {
  invalid_cc_models <- setdiff(models_for_controlled_comparison, available_models)
  if (length(invalid_cc_models) > 0) {
    warning("Invalid models in controlled variable comparison plot configuration: ", paste(invalid_cc_models, collapse = ", "))
    models_for_controlled_comparison <- intersect(models_for_controlled_comparison, available_models)
  }
  cat("✓ Controlled variable comparison plot display models:", paste(models_for_controlled_comparison, collapse = ", "), "\n")
  cat("✓ Axis range - X-axis:", paste(controlled_comparison_xlim, collapse = " to "), "\n")
  cat("✓ Axis range - Y-axis:", paste(controlled_comparison_ylim, collapse = " to "), "\n")
  if (exists("n_bootstrap")) {
    cat("✓ Bootstrap resampling number:", n_bootstrap, "times\n")
  }
}

# Display other configurations
cat("✓ Cytoscape export:", ifelse(enable_cytoscape_export, "Enabled", "Disabled"), "\n")
cat("✓ Visualization generation:", ifelse(enable_visualization, "Enabled", "Disabled"), "\n")
plot_toggle_status <- c(
  "Figure 1 AUC Grouped Forest Plot" = enable_plot_auc_summary,
  "Figure 1b AUC vs 1D_logFC" = enable_plot_auc_summary_vs_logFC,
  "Figure 2 Inter-group Comparison" = enable_plot_group_comparison,
  "Figure 3 Sample Intersection Dumbbell Plot" = enable_plot_intersection_dumbbell,
  "Figure 4 Performance-Stability Scatter Plot" = enable_plot_robustness_scatter,
  "Figure 5 Controlled Variable Gain" = enable_plot_controlled_comparison,
  "Figure 6 Candidate Protein Localization Distribution" = enable_plot_localization_distribution,
  "Figure 7 Controlled Variable Gain vs 1D_logFC" = enable_plot_controlled_comparison_vs_logFC
)
cat("  Visualization subplot control:\n")
for (toggle_name in names(plot_toggle_status)) {
  cat("   -", toggle_name, ":", ifelse(plot_toggle_status[[toggle_name]], "Enabled", "Disabled"), "\n")
}
cat("✓ Robustness analysis:", ifelse(enable_robustness_analysis, "Enabled", "Disabled"), "\n")
cat("✓ Inter-group comparison:", ifelse(enable_group_comparison, "Enabled", "Disabled"), "\n")
cat("✓ Parallel cores:", parallel_cores, "\n")

cat("=== Configuration Validation Complete ===\n")
cat("=== Usage Instructions ===\n")
cat("To modify configuration, edit the above variables:\n")
cat("1. selected_groups: Select specific groups, e.g., c(\"E7A2B4\", \"C2\")\n")
cat("2. selected_models: Select specific models, e.g., c(\"1D_Abundance\", \"2D_Abundance_logFC\")\n")
cat("3. enable_* variables: Control whether to generate specific types of output files and plots\n")
cat("4. parallel_cores: Adjust the number of cores for parallel computing\n")
cat("=== Starting Analysis ===\n\n")

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
# Using 20 cores as requested
plan(multisession, workers = 20)
cat(paste("Parallel computing started, will use", nbrOfWorkers(), "cores.\n"))

# --- Key Parameters ---
# This parameter can be adjusted for the local PPI calculation
local_ppi_window_size <- 10

# Create output directories
output_dir <- "Step29_ModelComparision"
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
LOCAL_STRING_DIR <- file.path(getwd(), "Reference", "StringDb")
LOCAL_STRING_ALIAS_FILE <- file.path(LOCAL_STRING_DIR, "9606.protein.aliases.v12.0.txt")
LOCAL_STRING_LINKS_FILE <- file.path(LOCAL_STRING_DIR, "9606.protein.links.detailed.v12.0.txt")
SCORE_THRESHOLD <- 400  # Consistent with the original online version

# Validate file paths
cat(paste("STRING file directory:", LOCAL_STRING_DIR, "\n"))
cat(paste("Alias file:", LOCAL_STRING_ALIAS_FILE, "\n"))
cat(paste("Interaction file:", LOCAL_STRING_LINKS_FILE, "\n"))

# === Step 1: Load STRING Alias Mapping Table (load only once) ===
cat("Loading STRING alias mapping table...\n")
if (!file.exists(LOCAL_STRING_ALIAS_FILE)) {
  stop(paste("Error: Mapping file not found", LOCAL_STRING_ALIAS_FILE))
}

string_aliases <- read_tsv(LOCAL_STRING_ALIAS_FILE, 
                           col_names = c("STRING_id", "alias", "source"),
                           col_types = "ccc",
                           skip = 1,
                           show_col_types = FALSE)

# Filter aliases from gene name sources (retain multiple sources to improve mapping rate)
string_gene_map <- string_aliases %>%
  filter(source %in% c("BioMart_HUGO", "Ensembl_gene_name", "BLAST_UniProt_GN", 
                       "Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)")) %>%
  select(STRING_id, gene = alias) %>%
  distinct()

cat(paste("  ✅ Loaded", nrow(string_gene_map), "gene mapping records\n"))

# === Step 2: Load STRING Interaction Network (load only once) ===
cat("Loading STRING interaction network...\n")
if (!file.exists(LOCAL_STRING_LINKS_FILE)) {
  stop(paste("Error: Interaction file not found", LOCAL_STRING_LINKS_FILE))
}

# STRING files use space separator, need to specify explicitly
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

# === Step 3: Process PPI for Each Data Source ===
ppi_data <- list(dfs = ForStep19, interactions = list(), ppi_global = list())

for(df_name in names(ForStep19)) {
  cat(paste("Processing PPI information for data source:", df_name, "\n"))
  
  tryCatch({
    # Get protein list
    protein_list <- unique(ForStep19[[df_name]]$Gene)
    cat(paste("  - Data source contains", length(protein_list), "proteins\n"))
    
    # Map gene names to STRING IDs
    mapped_proteins <- string_gene_map %>%
      filter(gene %in% protein_list) %>%
      group_by(gene) %>%
      slice(1) %>%  # Keep only the first STRING ID for each gene
      ungroup()
    
    mapping_rate <- round(nrow(mapped_proteins) / length(protein_list) * 100, 1)
    cat(paste("  - Successfully mapped", nrow(mapped_proteins), "/", length(protein_list), 
              "proteins (", mapping_rate, "%)\n"))
    
    if (nrow(mapped_proteins) == 0) {
      cat("  ⚠️  Warning: No proteins were successfully mapped to STRING IDs, skipping this data source\n")
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
    
    # Convert STRING IDs back to gene names and scale scores to 0-1
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
    
    # Calculate global PPI scores (total interaction strength for each protein)
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
  
  # Define myTP_vector within function to avoid global variable dependency
  myTP_vector=c("Nucleolus","Other&Nucleolus","Nuclear&Nucleolus","Nuclear&Cytosol&Nucleolus","Cytosol&Nucleolus")
  
  # Define local_ppi_window_size within function
  local_ppi_window_size <- 10
  
  # --- Unpack data for the current task ---
  # Add safety checks for data availability
  if (is.null(ppi_data$dfs[[df_name]])) {
    warning(paste("No data found for df_name:", df_name))
    return(NULL)
  }
  
  df <- ppi_data$dfs[[df_name]]
  interactions_df <- ppi_data$interactions[[df_name]]
  ppi_feature_global <- ppi_data$ppi_global[[df_name]]
  
  # Handle missing interactions data gracefully
  if (is.null(interactions_df)) {
    warning(paste("No PPI interactions data for:", df_name))
    interactions_df <- data.frame(protein1 = character(0), protein2 = character(0), combined_score = numeric(0))
  }
  
  if (is.null(ppi_feature_global)) {
    warning(paste("No PPI global features for:", df_name))
    ppi_feature_global <- data.frame(Gene = character(0), ppi_score_global = numeric(0))
  }
  
  # --- Data Preparation ---
  # Check if group exists in group_info
  if (is.null(group_info[[group]])) {
    warning(paste("No group info found for group:", group))
    return(NULL)
  }
  
  samples <- group_info[[group]]$samples
  logFC_col <- group_info[[group]]$logFC
  fdr_col <- group_info[[group]]$logFC_FDR
  
  # Verify required columns exist
  if (is.null(samples) || is.null(logFC_col) || is.null(fdr_col)) {
    warning(paste("Missing required columns for group:", group))
    return(NULL)
  }
  
  # Check if required columns exist in dataframe
  missing_cols <- setdiff(c(samples, logFC_col, fdr_col, "MultiBait_Localization", "Gene"), names(df))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in dataframe for", df_name, "-", group, ":", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
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
  
  # --- Local PPI Feature Calculation ---
  sorted_genes <- plot_data %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  
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
  auc_summary_df <- NULL
  roc_objects <- list()
  
  # Ensure there are at least two groups to compare for modeling
  if (n_distinct(plot_data_full$plot_group) < 2 || any(table(plot_data_full$plot_group) < 5)) {
    warning(paste("Skipping modeling for", df_name, "-", group, "due to insufficient data in one or both classes."))
    return(NULL)
  }
  
  tryCatch({
    # Define all available models
    all_models_available <- list(
      "1D_Abundance" = list(type="1D", predictor="mean_group_medianNorm"),
      "1D_logFC" = list(type="1D", predictor=logFC_col),
      "1D_PPI_Global" = list(type="1D", predictor="log_ppi_score_global"),
      "1D_PPI_Local" = list(type="1D", predictor="log_ppi_score_local"),
      "2D_Abundance_logFC" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col)),
      "2D_Abundance_PPI_Global" = list(type="GLM", formula="plot_group ~ mean_group_medianNorm + log_ppi_score_global"),
      "2D_Abundance_PPI_Local" = list(type="GLM", formula="plot_group ~ mean_group_medianNorm + log_ppi_score_local"),
      "2D_logFC_PPI_Global" = list(type="GLM", formula=paste("plot_group ~", logFC_col, "+ log_ppi_score_global")),
      "2D_logFC_PPI_Local" = list(type="GLM", formula=paste("plot_group ~", logFC_col, "+ log_ppi_score_local")),
      "3D_with_Global_PPI" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col, "+ log_ppi_score_global")),
      "3D_with_Local_PPI" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col, "+ log_ppi_score_local"))
    )
    
    # Filter models based on user selection
    models_to_run <- all_models_available[intersect(names(all_models_available), selected_models)]
    
    # Use the pre-initialized roc_objects list

    for(model_name in names(models_to_run)) {
      model_info <- models_to_run[[model_name]]
      predictor_values <- NULL
      
      # Train model and get predictor values
      if(model_info$type == "1D") {
        model <- glm(as.formula(paste("plot_group ~", model_info$predictor)), data = plot_data_full, family = "binomial")
        predictor_values <- plot_data_full[[model_info$predictor]]
        model_predictions[[paste0("prob_", model_name)]] <- predict(model, type = "response")
      } else { # GLM for 2D/3D
        model <- glm(as.formula(model_info$formula), data = plot_data_full, family = "binomial")
        predictor_values <- predict(model, type = "response")
        model_predictions[[paste0("prob_", model_name)]] <- predictor_values
      }
      
      # Perform ROC analysis and store the roc object
      roc_objects[[model_name]] <- roc(response = plot_data_full$plot_group, predictor = predictor_values, 
                                       levels = c("Other", "TP"), quiet = TRUE, direction = "<")
      
      # Get best threshold using Youden index
      best_coords <- coords(roc_objects[[model_name]], "best", best.method = "youden", ret = "threshold")
      threshold_value <- best_coords$threshold
      
      # Score all proteins based on the model and create the detailed output table
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

    # --- Add Traditional and Control Methods to ROC objects list ---
    # Only calculate if user selected traditional method
    if ("logFC>0.5 & FDR<0.05" %in% selected_models) {
      # Handle NA values in logFC and FDR columns gracefully
      traditional_predictor <- ifelse(
        !is.na(plot_data_full[[logFC_col]]) & !is.na(plot_data_full[[fdr_col]]) & 
        plot_data_full[[logFC_col]] > 0.5 & plot_data_full[[fdr_col]] < 0.05, 
        1, 0
      )
      
      # Only create ROC if there's variation in the predictor
      if (length(unique(traditional_predictor)) > 1) {
        roc_objects[["logFC>0.5 & FDR<0.05"]] <- roc(response = plot_data_full$plot_group, predictor = traditional_predictor, levels = c("Other", "TP"), quiet = TRUE, direction = "<")
      } else {
        warning(paste("Traditional method has no variation for", df_name, "-", group, ", skipping ROC creation"))
      }
    }
    
    set.seed(123)
    random_predictor <- runif(nrow(plot_data_full))
    roc_objects[["Control"]] <- roc(response = plot_data_full$plot_group, predictor = random_predictor, levels = c("Other", "TP"), quiet = TRUE, direction = "<")

    # --- Calculate AUC, CIs, and P-values vs. Traditional Method ---
    # Get AUC scores and CIs into a single data frame
    auc_summary_df <- purrr::map_df(roc_objects, ~ {
      ci_vals <- pROC::ci.auc(.x, method = "delong")
      tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
    }, .id = "ModelName")

    # Get p-values vs. the traditional model
    traditional_roc <- roc_objects[["logFC>0.5 & FDR<0.05"]]
    p_values_vs_trad <- purrr::map_dbl(roc_objects, ~ {
      # Don't compare the traditional model to itself or if traditional_roc is NULL
      if (identical(.x, traditional_roc) || is.null(traditional_roc)) return(NA)
      tryCatch({
        test_result <- roc.test(.x, traditional_roc, method = "delong")
        test_result$p.value
      }, error = function(e) { 
        warning(paste("ROC test failed:", e$message))
        NA 
      })
    })
    
    # Add p-values to the summary dataframe
    auc_summary_df <- auc_summary_df %>%
      mutate(p_value_vs_traditional = p_values_vs_trad)
    
    # --- Calculate P-values vs. 1D_logFC Model (for Figure 1b) ---
    logFC_roc <- roc_objects[["1D_logFC"]]
    p_values_vs_logFC <- purrr::map_dbl(roc_objects, ~ {
      # Don't compare 1D_logFC to itself or if logFC_roc is NULL
      if (identical(.x, logFC_roc) || is.null(logFC_roc)) return(NA)
      tryCatch({
        test_result <- roc.test(.x, logFC_roc, method = "delong")
        test_result$p.value
      }, error = function(e) { 
        warning(paste("ROC test vs 1D_logFC failed:", e$message))
        NA 
      })
    })
    
    # Add p-values vs 1D_logFC to the summary dataframe
    auc_summary_df <- auc_summary_df %>%
      mutate(p_value_vs_1D_logFC = p_values_vs_logFC)

    # --- New Analysis: Fixed Sensitivity/Specificity Comparison ---
    # Controlled variable comparison: fix Sensitivity to compare Specificity, or fix Specificity to compare Sensitivity
    controlled_comparison_results <- NULL
    
    if (!is.null(traditional_roc)) {
      # Get performance metrics of traditional method
      traditional_coords <- coords(traditional_roc, x = 1, input = "threshold", 
                                   ret = c("sensitivity", "specificity"))
      traditional_sens <- traditional_coords$sensitivity
      traditional_spec <- traditional_coords$specificity
      
      # Perform controlled variable comparison for each continuous model
      controlled_comparison_list <- list()
      
      for (model_name in setdiff(names(roc_objects), c("logFC>0.5 & FDR<0.05", "Control"))) {
        model_roc <- roc_objects[[model_name]]
        
        tryCatch({
          # 1. Fix Sensitivity, compare Specificity
          # Find the point on ROC curve closest to traditional method's sensitivity
          model_coords_all <- coords(model_roc, x = "all", ret = c("threshold", "sensitivity", "specificity"))
          
          # Find the point closest to traditional method's sensitivity
          sens_diff <- abs(model_coords_all$sensitivity - traditional_sens)
          closest_sens_idx <- which.min(sens_diff)
          model_spec_at_fixed_sens <- model_coords_all$specificity[closest_sens_idx]
          actual_sens_used <- model_coords_all$sensitivity[closest_sens_idx]
          
          # 2. Fix Specificity, compare Sensitivity
          # Find the point closest to traditional method's specificity
          spec_diff <- abs(model_coords_all$specificity - traditional_spec)
          closest_spec_idx <- which.min(spec_diff)
          model_sens_at_fixed_spec <- model_coords_all$sensitivity[closest_spec_idx]
          actual_spec_used <- model_coords_all$specificity[closest_spec_idx]
          
          # Calculate gains
          specificity_gain <- model_spec_at_fixed_sens - traditional_spec
          sensitivity_gain <- model_sens_at_fixed_spec - traditional_sens
          
          # Use bootstrap method to calculate confidence intervals and p-values
          # Bootstrap for specificity comparison (fixed sensitivity)
          set.seed(123)
          n_boot <- if(exists("n_bootstrap")) n_bootstrap else 1000  # Use global bootstrap configuration, default 1000
          boot_spec_diff <- numeric(n_boot)
          
          # Get original data
          response_vec <- as.numeric(model_roc$response) - 1  # Convert to 0/1
          predictor_traditional <- as.numeric(traditional_roc$predictor)
          predictor_model <- as.numeric(model_roc$predictor)
          
          for (b in 1:n_boot) {
            # Bootstrap sample
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_trad <- predictor_traditional[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            # Calculate traditional method specificity
            boot_trad_pred_binary <- ifelse(boot_pred_trad >= 1, 1, 0)
            boot_trad_spec <- sum(boot_response == 0 & boot_trad_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            # For model, find threshold that gives similar sensitivity to traditional
            boot_trad_sens <- sum(boot_response == 1 & boot_trad_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            # Find model threshold
            thresholds <- sort(unique(boot_pred_model))
            best_thresh <- thresholds[1]
            min_sens_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
              sens_diff_boot <- abs(boot_model_sens - boot_trad_sens)
              if (sens_diff_boot < min_sens_diff) {
                min_sens_diff <- sens_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            boot_spec_diff[b] <- boot_model_spec - boot_trad_spec
          }
          
          # Calculate CI and p-value for specificity gain
          spec_gain_ci_lower <- quantile(boot_spec_diff, 0.025, na.rm = TRUE)
          spec_gain_ci_upper <- quantile(boot_spec_diff, 0.975, na.rm = TRUE)
          # Use conservative P-value estimation: (x+1)/(n+1) to avoid ambiguity of P=0
          # This is a simplified version of Agresti-Coull adjustment
          n_valid_boot <- sum(!is.na(boot_spec_diff))
          n_unfavorable <- sum(boot_spec_diff <= 0, na.rm = TRUE)
          spec_gain_p_value <- (n_unfavorable + 1) / (n_valid_boot + 1)
          
          # Bootstrap for sensitivity comparison (fixed specificity)
          boot_sens_diff <- numeric(n_boot)
          
          for (b in 1:n_boot) {
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_trad <- predictor_traditional[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            # Calculate traditional method sensitivity
            boot_trad_pred_binary <- ifelse(boot_pred_trad >= 1, 1, 0)
            boot_trad_sens <- sum(boot_response == 1 & boot_trad_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            # For model, find threshold that gives similar specificity to traditional
            boot_trad_spec <- sum(boot_response == 0 & boot_trad_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            # Find model threshold
            thresholds <- sort(unique(boot_pred_model), decreasing = TRUE)
            best_thresh <- thresholds[1]
            min_spec_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
              spec_diff_boot <- abs(boot_model_spec - boot_trad_spec)
              if (spec_diff_boot < min_spec_diff) {
                min_spec_diff <- spec_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            boot_sens_diff[b] <- boot_model_sens - boot_trad_sens
          }
          
          # Calculate CI and p-value for sensitivity gain
          sens_gain_ci_lower <- quantile(boot_sens_diff, 0.025, na.rm = TRUE)
          sens_gain_ci_upper <- quantile(boot_sens_diff, 0.975, na.rm = TRUE)
          # Use conservative P-value estimation: (x+1)/(n+1) to avoid ambiguity of P=0
          n_valid_boot_sens <- sum(!is.na(boot_sens_diff))
          n_unfavorable_sens <- sum(boot_sens_diff <= 0, na.rm = TRUE)
          sens_gain_p_value <- (n_unfavorable_sens + 1) / (n_valid_boot_sens + 1)
          
          # Store results
          controlled_comparison_list[[model_name]] <- data.frame(
            ModelName = model_name,
            # Traditional method metrics
            Traditional_Sensitivity = traditional_sens,
            Traditional_Specificity = traditional_spec,
            # Fixed Sensitivity comparison
            Fixed_Sensitivity = actual_sens_used,
            Model_Specificity_at_Fixed_Sens = model_spec_at_fixed_sens,
            Specificity_Gain = specificity_gain,
            Specificity_Gain_CI_Lower = spec_gain_ci_lower,
            Specificity_Gain_CI_Upper = spec_gain_ci_upper,
            Specificity_Gain_P_Value = spec_gain_p_value,
            # Fixed Specificity comparison
            Fixed_Specificity = actual_spec_used,
            Model_Sensitivity_at_Fixed_Spec = model_sens_at_fixed_spec,
            Sensitivity_Gain = sensitivity_gain,
            Sensitivity_Gain_CI_Lower = sens_gain_ci_lower,
            Sensitivity_Gain_CI_Upper = sens_gain_ci_upper,
            Sensitivity_Gain_P_Value = sens_gain_p_value,
            stringsAsFactors = FALSE
          )
          
        }, error = function(e) {
          warning(paste("Controlled comparison failed for", model_name, ":", e$message))
        })
      }
      
      # Merge controlled variable comparison results for all models
      if (length(controlled_comparison_list) > 0) {
        controlled_comparison_results <- bind_rows(controlled_comparison_list)
      }
    }
    
    # --- New Analysis: Comparison vs 1D_logFC Model (for Figure 7a-7c) ---
    # Controlled variable comparison vs 1D_logFC model
    controlled_comparison_vs_logFC <- NULL
    
    logFC_roc <- roc_objects[["1D_logFC"]]
    
    if (!is.null(logFC_roc)) {
      # Get performance metrics of 1D_logFC model at Youden optimal threshold
      logFC_best_coords <- coords(logFC_roc, "best", best.method = "youden", 
                                   ret = c("threshold", "sensitivity", "specificity"))
      logFC_sens <- logFC_best_coords$sensitivity
      logFC_spec <- logFC_best_coords$specificity
      
      # Perform controlled variable comparison for each other continuous model
      controlled_comparison_vs_logFC_list <- list()
      
      for (model_name in setdiff(names(roc_objects), c("1D_logFC", "logFC>0.5 & FDR<0.05", "Control"))) {
        model_roc <- roc_objects[[model_name]]
        
        tryCatch({
          # 1. Fix Sensitivity, compare Specificity
          model_coords_all <- coords(model_roc, x = "all", ret = c("threshold", "sensitivity", "specificity"))
          
          sens_diff <- abs(model_coords_all$sensitivity - logFC_sens)
          closest_sens_idx <- which.min(sens_diff)
          model_spec_at_fixed_sens <- model_coords_all$specificity[closest_sens_idx]
          actual_sens_used <- model_coords_all$sensitivity[closest_sens_idx]
          
          # 2. Fix Specificity, compare Sensitivity
          spec_diff <- abs(model_coords_all$specificity - logFC_spec)
          closest_spec_idx <- which.min(spec_diff)
          model_sens_at_fixed_spec <- model_coords_all$sensitivity[closest_spec_idx]
          actual_spec_used <- model_coords_all$specificity[closest_spec_idx]
          
          # Calculate gains
          specificity_gain <- model_spec_at_fixed_sens - logFC_spec
          sensitivity_gain <- model_sens_at_fixed_spec - logFC_sens
          
          # Use bootstrap method to calculate confidence intervals and p-values
          set.seed(123)
          n_boot <- if(exists("n_bootstrap")) n_bootstrap else 1000
          boot_spec_diff <- numeric(n_boot)
          
          # Get original data
          response_vec <- as.numeric(model_roc$response) - 1
          predictor_logFC <- as.numeric(logFC_roc$predictor)
          predictor_model <- as.numeric(model_roc$predictor)
          
          # Bootstrap for specificity comparison (fixed sensitivity)
          for (b in 1:n_boot) {
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_logFC <- predictor_logFC[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            # Calculate logFC model specificity at its best threshold
            logFC_best_thresh <- logFC_best_coords$threshold
            boot_logFC_pred_binary <- ifelse(boot_pred_logFC >= logFC_best_thresh, 1, 0)
            boot_logFC_spec <- sum(boot_response == 0 & boot_logFC_pred_binary == 0) / max(1, sum(boot_response == 0))
            boot_logFC_sens <- sum(boot_response == 1 & boot_logFC_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            # Find model threshold that gives similar sensitivity
            thresholds <- sort(unique(boot_pred_model))
            best_thresh <- thresholds[1]
            min_sens_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
              sens_diff_boot <- abs(boot_model_sens - boot_logFC_sens)
              if (sens_diff_boot < min_sens_diff) {
                min_sens_diff <- sens_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            boot_spec_diff[b] <- boot_model_spec - boot_logFC_spec
          }
          
          spec_gain_ci_lower <- quantile(boot_spec_diff, 0.025, na.rm = TRUE)
          spec_gain_ci_upper <- quantile(boot_spec_diff, 0.975, na.rm = TRUE)
          # Use conservative P-value estimation: (x+1)/(n+1)
          n_valid_boot <- sum(!is.na(boot_spec_diff))
          n_unfavorable <- sum(boot_spec_diff <= 0, na.rm = TRUE)
          spec_gain_p_value <- (n_unfavorable + 1) / (n_valid_boot + 1)
          
          # Bootstrap for sensitivity comparison (fixed specificity)
          boot_sens_diff <- numeric(n_boot)
          
          for (b in 1:n_boot) {
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_logFC <- predictor_logFC[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            boot_logFC_pred_binary <- ifelse(boot_pred_logFC >= logFC_best_thresh, 1, 0)
            boot_logFC_sens <- sum(boot_response == 1 & boot_logFC_pred_binary == 1) / max(1, sum(boot_response == 1))
            boot_logFC_spec <- sum(boot_response == 0 & boot_logFC_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            # Find model threshold that gives similar specificity
            thresholds <- sort(unique(boot_pred_model), decreasing = TRUE)
            best_thresh <- thresholds[1]
            min_spec_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
              spec_diff_boot <- abs(boot_model_spec - boot_logFC_spec)
              if (spec_diff_boot < min_spec_diff) {
                min_spec_diff <- spec_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            boot_sens_diff[b] <- boot_model_sens - boot_logFC_sens
          }
          
          sens_gain_ci_lower <- quantile(boot_sens_diff, 0.025, na.rm = TRUE)
          sens_gain_ci_upper <- quantile(boot_sens_diff, 0.975, na.rm = TRUE)
          # Use conservative P-value estimation: (x+1)/(n+1)
          n_valid_boot_sens <- sum(!is.na(boot_sens_diff))
          n_unfavorable_sens <- sum(boot_sens_diff <= 0, na.rm = TRUE)
          sens_gain_p_value <- (n_unfavorable_sens + 1) / (n_valid_boot_sens + 1)
          
          # Store results
          controlled_comparison_vs_logFC_list[[model_name]] <- data.frame(
            ModelName = model_name,
            # 1D_logFC model metrics
            Baseline_Sensitivity = logFC_sens,
            Baseline_Specificity = logFC_spec,
            # Fixed Sensitivity comparison
            Fixed_Sensitivity = actual_sens_used,
            Model_Specificity_at_Fixed_Sens = model_spec_at_fixed_sens,
            Specificity_Gain = specificity_gain,
            Specificity_Gain_CI_Lower = spec_gain_ci_lower,
            Specificity_Gain_CI_Upper = spec_gain_ci_upper,
            Specificity_Gain_P_Value = spec_gain_p_value,
            # Fixed Specificity comparison
            Fixed_Specificity = actual_spec_used,
            Model_Sensitivity_at_Fixed_Spec = model_sens_at_fixed_spec,
            Sensitivity_Gain = sensitivity_gain,
            Sensitivity_Gain_CI_Lower = sens_gain_ci_lower,
            Sensitivity_Gain_CI_Upper = sens_gain_ci_upper,
            Sensitivity_Gain_P_Value = sens_gain_p_value,
            stringsAsFactors = FALSE
          )
          
        }, error = function(e) {
          warning(paste("Controlled comparison vs 1D_logFC failed for", model_name, ":", e$message))
        })
      }
      
      # Merge comparison results for all models vs 1D_logFC
      if (length(controlled_comparison_vs_logFC_list) > 0) {
        controlled_comparison_vs_logFC <- bind_rows(controlled_comparison_vs_logFC_list)
      }
    }

    # Add traditional method to the scored list (for consistency in Phase 3) only if ROC was created
    # Note: Threshold is binary (1), so candidates are simply those that pass
    if ("logFC>0.5 & FDR<0.05" %in% names(roc_objects) && "logFC>0.5 & FDR<0.05" %in% selected_models) {
      all_proteins_scored_list[["logFC>0.5 & FDR<0.05"]] <- plot_data_full %>%
          mutate(
            model_score = traditional_predictor,
            Model_best_threshold = 1,
            PassThreshold = if_else(model_score >= 1, "Pass", "Fail"),
            isCandidate = if_else(plot_group == "Other" & PassThreshold == "Pass", "TP_candidate", NA_character_)
          ) %>%
          select(
            protein_id = Gene, group = plot_group, contains("Localization"),
            abundance = mean_group_medianNorm, logFC = all_of(logFC_col),
            ppi_score_global, ppi_score_local, log_ppi_score_global, log_ppi_score_local,
            model_score, Model_best_threshold, PassThreshold, isCandidate
          )
    }

  }, error = function(e) { 
    warning(paste("Model evaluation failed for", df_name, "-", group, ":", e$message))
    # Return empty objects to prevent downstream errors
    auc_summary_df <<- data.frame()
    roc_objects <<- list()
    all_proteins_scored_list <<- list()
  })
  
  # --- Package all results for return ---
  return(list(
    df_name = df_name, 
    group_name = group,
    auc_summary_df = if(is.null(auc_summary_df)) data.frame() else auc_summary_df, # AUC, CIs, and p-value vs traditional
    roc_objects = roc_objects, # Return full ROC objects for cross-group comparison
    node_attributes = plot_data_full,
    model_predictions = model_predictions, 
    candidates_full_list = all_proteins_scored_list,
    interactions = interactions_df,
    controlled_comparison = if(is.null(controlled_comparison_results)) data.frame() else controlled_comparison_results, # Controlled variable comparison vs traditional
    controlled_comparison_vs_logFC = if(is.null(controlled_comparison_vs_logFC)) data.frame() else controlled_comparison_vs_logFC # Controlled variable comparison vs 1D_logFC
  ))
}


# -------------------------------------------------------------------
# 5. Execute Parallel Analysis
# -------------------------------------------------------------------
cat("Starting parallel analysis for all combinations...\n")

# Enable progress bar tracking
handlers(global = TRUE)
handlers("progress")

# Create a grid of tasks based on selected groups
tasks <- expand_grid(df_name = names(ForStep19), group = selected_groups)

with_progress({
  # Setup the progressor
  p <- progressr::progressor(steps = nrow(tasks))
  
  # Use future_pmap to iterate over all combinations in parallel
  all_results <- future_pmap(list(
    df_name = tasks$df_name,
    group = tasks$group
  ), function(df_name, group) {
    # Run the main worker function
    result <- tryCatch({
      process_group(df_name, group, ppi_data)
    }, error = function(e) {
      warning(paste("Task failed for", df_name, "-", group, ":", e$message))
      return(NULL)
    })
    
    # Update progress only if progressor is still active
    tryCatch({
      if (!is.null(p)) {
        p(message = paste("Completed", df_name, "-", group))
      }
    }, error = function(e) {
      # Silently ignore progress update errors - this is expected when progressor completes
    })
    
    return(result)
  }, .options = furrr_options(seed = TRUE)) # Use seed for reproducibility
  
  # Assign meaningful names to the entries in all_results
  names(all_results) <- paste(tasks$df_name, tasks$group, sep = "_")
})


cat("\nAll parallel computing tasks have completed.\n")


# -------------------------------------------------------------------
# 6. Consolidate Results and Export Files
# -------------------------------------------------------------------
cat("Consolidating results and exporting files...\n")

# --- Consolidate AUC Summaries and Apply Multiple-Test Correction ---
# Extract the auc_summary_df from each task result
all_auc_summaries <- purrr::map_df(all_results, ~ {
  if (is.null(.x) || !is.list(.x) || is.null(.x$auc_summary_df)) {
    return(NULL)
  }
  .x$auc_summary_df
}, .id = "Analysis_Group")

# Apply BH correction to p-values within each analysis group
final_auc_summary <- all_auc_summaries %>% 
  group_by(Analysis_Group) %>% 
  mutate(
    q_value_vs_traditional = p.adjust(p_value_vs_traditional, method = "BH")
  ) %>% 
  ungroup() %>%
  select(Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, p_value_vs_traditional, q_value_vs_traditional, everything())

# --- Export Main AUC Summary ---
writexl::write_xlsx(list(AUC_Summary = final_auc_summary), path = file.path(output_dir, "AUC_summary.xlsx"))

# --- Consolidate and Export Candidate Lists and Cytoscape Files ---
final_candidates_list <- list()

# Get the task list again for matching failed tasks
tasks_check <- expand_grid(df_name = names(ForStep19), group = names(group_info))

for(i in seq_along(all_results)) {
  result <- all_results[[i]]
  
  # Check for failed tasks
  if (is.null(result) || !is.list(result) || is.null(result$candidates_full_list)) {
    failed_task <- tasks_check[i, ]
    warning(paste("Task", failed_task$df_name, "-", failed_task$group, "results incomplete, some export files may be missing."))
    next
  }
  
  df_name <- result$df_name
  group <- result$group_name

  # Consolidate scored protein list
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
  
  # Export files for Cytoscape (if enabled)
  if (enable_cytoscape_export) {
    cytoscape_path_prefix <- file.path(output_dir, "Cytoscape_Files", paste0(df_name, "_", group))
    dir.create(dirname(cytoscape_path_prefix), recursive = TRUE, showWarnings = FALSE)
    logFC_col <- group_info[[group]]$logFC
    
    # Node attributes file
    node_attributes_for_export <- result$node_attributes %>% 
      left_join(result$model_predictions, by = c("Gene" = "protein_id")) %>% 
      select(
        protein_id = Gene, 
        group = plot_group,
        contains("Localization"),
        logFC = all_of(logFC_col),
        everything()
      )
    write.csv(node_attributes_for_export, 
              paste0(cytoscape_path_prefix, "_node_attributes.csv"), 
              row.names = FALSE,quote=FALSE)
    
    # Edge attributes file
    write.csv(result$interactions, paste0(cytoscape_path_prefix, "_edge_attributes.csv"), 
              row.names = FALSE,quote=FALSE)
  }
}
#install.packages("writexl")
library("writexl")

# Combine and export candidate list
final_candidates_df <- bind_rows(final_candidates_list)
writexl::write_xlsx(list(Candidates = final_candidates_df), path = file.path(output_dir, "detailed_candidates_by_model.xlsx"))

cat("\nAll result files successfully exported to directory:", output_dir, "\n")

# --- Consolidate and Export Controlled Comparison Results ---
cat("\nConsolidating controlled variable comparison results (fixed Sensitivity/Specificity)...\n")

all_controlled_comparisons <- purrr::map_df(all_results, ~ {
  if (is.null(.x) || !is.list(.x) || is.null(.x$controlled_comparison)) {
    return(NULL)
  }
  .x$controlled_comparison
}, .id = "Analysis_Group")

if (nrow(all_controlled_comparisons) > 0) {
  # Add Df_Name column
  group_to_df_map_controlled <- purrr::map_df(all_results, ~tibble::tibble(Df_Name = .x$df_name), .id = "Analysis_Group") %>% dplyr::distinct()
  
  final_controlled_comparison <- all_controlled_comparisons %>%
    dplyr::left_join(group_to_df_map_controlled, by = "Analysis_Group") %>%
    select(Df_Name, Analysis_Group, everything()) %>%
    # Apply BH multiple testing correction
    group_by(Analysis_Group) %>%
    mutate(
      Specificity_Gain_Q_Value = p.adjust(Specificity_Gain_P_Value, method = "BH"),
      Sensitivity_Gain_Q_Value = p.adjust(Sensitivity_Gain_P_Value, method = "BH"),
      # Add P-value theoretical lower limit note (based on bootstrap number)
      P_Value_Lower_Limit = 1 / (n_bootstrap + 1),
      # Add significance markers
      Specificity_Significant = if_else(Specificity_Gain_Q_Value < 0.05, "Yes", "No"),
      Sensitivity_Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "Yes", "No")
    ) %>%
    ungroup()
  
  # Export to Excel
  writexl::write_xlsx(
    list(Controlled_Comparison = final_controlled_comparison), 
    path = file.path(output_dir, "Controlled_Variable_Comparison.xlsx")
  )
  
  cat("Controlled variable comparison results exported to Controlled_Variable_Comparison.xlsx\n")
  cat("  - Fixed Sensitivity compare Specificity:", sum(!is.na(final_controlled_comparison$Specificity_Gain)), "comparisons\n")
  cat("  - Fixed Specificity compare Sensitivity:", sum(!is.na(final_controlled_comparison$Sensitivity_Gain)), "comparisons\n")
} else {
  cat("Warning: Could not extract any controlled variable comparison information from analysis results.\n")
}

# --- Consolidate and Export Controlled Comparison vs 1D_logFC Results ---
cat("\nConsolidating controlled variable comparison results vs 1D_logFC model...\n")

all_controlled_comparisons_vs_logFC <- purrr::map_df(all_results, ~ {
  if (is.null(.x) || !is.list(.x) || is.null(.x$controlled_comparison_vs_logFC)) {
    return(NULL)
  }
  .x$controlled_comparison_vs_logFC
}, .id = "Analysis_Group")

if (nrow(all_controlled_comparisons_vs_logFC) > 0) {
  # Add Df_Name column
  group_to_df_map_logFC <- purrr::map_df(all_results, ~tibble::tibble(Df_Name = .x$df_name), .id = "Analysis_Group") %>% dplyr::distinct()
  
  final_controlled_comparison_vs_logFC <- all_controlled_comparisons_vs_logFC %>%
    dplyr::left_join(group_to_df_map_logFC, by = "Analysis_Group") %>%
    select(Df_Name, Analysis_Group, everything()) %>%
    # Apply BH multiple testing correction
    group_by(Analysis_Group) %>%
    mutate(
      Specificity_Gain_Q_Value = p.adjust(Specificity_Gain_P_Value, method = "BH"),
      Sensitivity_Gain_Q_Value = p.adjust(Sensitivity_Gain_P_Value, method = "BH"),
      # Add P-value theoretical lower limit note (based on bootstrap number)
      P_Value_Lower_Limit = 1 / (n_bootstrap + 1),
      # Add significance markers
      Specificity_Significant = if_else(Specificity_Gain_Q_Value < 0.05, "Yes", "No"),
      Sensitivity_Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "Yes", "No")
    ) %>%
    ungroup()
  
  # Export to Excel
  writexl::write_xlsx(
    list(Controlled_Comparison_vs_1D_logFC = final_controlled_comparison_vs_logFC), 
    path = file.path(output_dir, "Controlled_Variable_Comparison_vs_1D_logFC.xlsx")
  )
  
  cat("Controlled variable comparison results vs 1D_logFC model exported to Controlled_Variable_Comparison_vs_1D_logFC.xlsx\n")
  cat("  - Fixed Sensitivity compare Specificity:", sum(!is.na(final_controlled_comparison_vs_logFC$Specificity_Gain)), "comparisons\n")
  cat("  - Fixed Specificity compare Sensitivity:", sum(!is.na(final_controlled_comparison_vs_logFC$Sensitivity_Gain)), "comparisons\n")
} else {
  cat("Warning: Could not extract any controlled variable comparison information vs 1D_logFC from analysis results.\n")
}


# ===================================================================
#
###    Phase 1.1: Data Export and Robustness Checks #####
#
# ===================================================================

cat("\n\n--- Starting Phase 1.1: Data Export and Robustness Checks ---\n")

# --- 0. Consolidate, Correct, and Export All Model AUCs ---
cat("\n--- Consolidating, correcting, and exporting all model AUC results... ---\n")

# Extract the auc_summary dataframe from each result, adding the group name
all_auc_summaries_raw <- purrr::map_dfr(all_results, ~ .x$auc_summary_df, .id = "Analysis_Group")

# Create a mapping from Analysis_Group (list names) to Df_Name (field in list)
# This is crucial for linking results back to their original dataset (e.g., APEX, BioID)
group_to_df_map <- purrr::map_df(all_results, ~tibble::tibble(Df_Name = .x$df_name), .id = "Analysis_Group") %>% dplyr::distinct()

# Check if the summary data frame is not empty
if (nrow(all_auc_summaries_raw) > 0) {
  
  # Detect whether it's a single comparison case
  total_comparisons <- length(selected_groups) * length(setdiff(selected_models, "Control"))
  single_comparison_flag <- (total_comparisons == 1)
  
  if (single_comparison_flag) {
    cat("✓ Detected single comparison case - only 1 group using 1 model, p-values don't need multiple testing correction\n")
  } else {
    cat("✓ Detected multiple comparison case - total", total_comparisons, "comparisons, will apply BH multiple testing correction\n")
  }
  
  # Apply BH correction to p-values within each analysis group to get q-values
  # This is the definitive, final summary table for AUCs.
  final_auc_summary <- all_auc_summaries_raw %>% 
    group_by(Analysis_Group) %>% 
    mutate(
      q_value_vs_traditional = if (single_comparison_flag) {
        # For single comparison, q-value equals p-value
        p_value_vs_traditional
      } else {
        # For multiple comparisons, apply BH correction
        p.adjust(p_value_vs_traditional, method = "BH")
      }
    ) %>%
    ungroup() %>% 
    # Add the Df_Name column by joining the map
    dplyr::left_join(group_to_df_map, by = "Analysis_Group") %>%
    # Reorder columns for clarity, bringing Df_Name to the front
    select(Df_Name, Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, p_value_vs_traditional, q_value_vs_traditional, everything()) %>%
    # Add marker column indicating whether it's a single comparison
    mutate(single_comparison = single_comparison_flag)

  # Export the final, corrected summary to its own Excel file
  auc_summary_path <- file.path(output_dir, "AUC_summary.xlsx")
  writexl::write_xlsx(list("AUC_Summary" = final_auc_summary), path = auc_summary_path)
  cat("All model AUC summary results (with multiple testing correction) exported to AUC_summary.xlsx\n")
  
  # --- Export AUC Summary with comparison to 1D_logFC (for Figure 1b) ---
  cat("\n--- Creating AUC summary results compared to 1D_logFC... ---\n")
  
  # Check if p_value_vs_1D_logFC column exists
  if ("p_value_vs_1D_logFC" %in% names(final_auc_summary)) {
    
    # Create summary table for 1D_logFC comparison
    final_auc_summary_b <- final_auc_summary %>%
      # Keep only columns related to 1D_logFC comparison
      select(Df_Name, Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, p_value_vs_1D_logFC) %>%
      # Filter out 1D_logFC itself and Control
      filter(!ModelName %in% c("1D_logFC", "Control")) %>%
      # Apply BH correction to p-values vs 1D_logFC
      group_by(Analysis_Group) %>%
      mutate(
        q_value_vs_1D_logFC = p.adjust(p_value_vs_1D_logFC, method = "BH"),
        # Add significance marker
        is_significant_vs_1D_logFC = if_else(
          !is.na(q_value_vs_1D_logFC) & q_value_vs_1D_logFC < 0.05,
          "Yes",
          "No"
        )
      ) %>%
      ungroup() %>%
      # Reorder columns
      select(Df_Name, Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, 
             p_value_vs_1D_logFC, q_value_vs_1D_logFC, is_significant_vs_1D_logFC, everything())
    
    # Export AUC_Summary_b
    auc_summary_b_path <- file.path(output_dir, "AUC_Summary_b.xlsx")
    writexl::write_xlsx(list("AUC_Summary_vs_1D_logFC" = final_auc_summary_b), path = auc_summary_b_path)
    cat("AUC summary results compared to 1D_logFC model exported to AUC_Summary_b.xlsx\n")
    
  } else {
    cat("Warning: p_value_vs_1D_logFC column not found, skipping AUC_Summary_b.xlsx generation.\n")
  }
  
} else {
  cat("Warning: Could not extract any AUC summary information from analysis results, skipping AUC_summary.xlsx generation.\n")
}


# --- 1. Type 2 Robustness (Model Stability) Analysis ---
if (enable_robustness_analysis) {
  cat("\nExecuting Type 2 robustness analysis (model stability)...\n")

  # Calculate Mean and StdDev of AUCs for each group across all models
  # This now uses the `final_auc_summary` table which contains q-values
  group_robustness_summary <- final_auc_summary %>%
    # Exclude the control and traditional methods from robustness calculation
    filter(!ModelName %in% c("Control", "logFC>0.5 & FDR<0.05")) %>%
    group_by(Df_Name, Analysis_Group) %>%
    summarise(
      Mean_AUC = mean(AUC, na.rm = TRUE),
      AUC_SD = sd(AUC, na.rm = TRUE),
      Num_Models_Tested = n(),
      .groups = 'drop'  # Avoid nested grouping issues
    ) %>%
    arrange(desc(Mean_AUC), AUC_SD) # Sort by best and most stable

  # Export the robustness summary
  writexl::write_xlsx(
    list(Group_Robustness_Summary = group_robustness_summary),
    path = file.path(output_dir, "Group_Robustness_Summary.xlsx")
  )

  cat("\nModel stability results exported to Group_Robustness_Summary.xlsx\n")
} else {
  cat("\nRobustness analysis disabled, skipping related calculations.\n")
}


# --- New Section: Cross-Group (Inter-Method) Comparison ---
if (enable_group_comparison) {
  cat("\nStarting cross-proximity method comparison analysis...\n")

# Define the baseline model for this comparison
baseline_model_for_comparison <- "2D_Abundance_logFC"

# Group results by the original dataframe name (e.g., "ForStep19_APEX", "ForStep19_BioID2")
# Ensure Phase 1.2 data copies exist when running this block alone
if (!exists("all_results_p1_2")) {
  all_results_p1_2 <- rlang::duplicate(all_results, shallow = FALSE)
}
if (!exists("group_info_p1_2") && exists("group_info")) {
  group_info_p1_2 <- rlang::duplicate(group_info, shallow = FALSE)
}

results_by_df <- all_results_p1_2 %>% 
  purrr::keep(~!is.null(.x)) %>% # Filter out failed runs
  split(., purrr::map_chr(., "df_name"))

# List to store the comparison results for each dataframe
all_group_comparisons <- list()

for (df_name in names(results_by_df)) {
  
  # Extract the ROC objects for the baseline model from each group
  roc_objects_for_comparison <- results_by_df[[df_name]] %>% 
    purrr::map(~.x$roc_objects[[baseline_model_for_comparison]]) %>% 
    # Name the list elements by their group name for clarity
    setNames(purrr::map_chr(results_by_df[[df_name]], "group_name")) %>% 
    # Remove any NULL entries if a model failed for a specific group
    purrr::compact()

  if (length(roc_objects_for_comparison) < 2) {
    warning(paste("Skipping group comparison for", df_name, "as it has fewer than 2 valid groups for the model", baseline_model_for_comparison))
    next
  }
  
  # --- Perform pairwise DeLong's test ---
  group_names <- names(roc_objects_for_comparison)
  p_value_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                           dimnames = list(group_names, group_names))
  
  group_pairs <- combn(group_names, 2, simplify = FALSE)
  
  p_values <- purrr::map_dbl(group_pairs, ~ {
    tryCatch({
      test_result <- roc.test(roc_objects_for_comparison[[.x[1]]], roc_objects_for_comparison[[.x[2]]], method = "delong")
      test_result$p.value
    }, error = function(e) { NA })
  })
  
  for (i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    p_value_matrix[pair[1], pair[2]] <- p_values[i]
    p_value_matrix[pair[2], pair[1]] <- p_values[i]
  }
  diag(p_value_matrix) <- 1.0

  # --- Consolidate results into a summary table ---
  # Get AUC and CIs for each group
  auc_ci_summary <- purrr::map_df(roc_objects_for_comparison, ~{
    ci_vals <- pROC::ci.auc(.x, method = "delong")
    tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
  }, .id = "GroupName")
  
  # Convert p-value matrix to a long format and apply BH correction
  p_values_long_temp <- as.data.frame(p_value_matrix) %>% 
    tibble::rownames_to_column("Group1") %>% 
    tidyr::pivot_longer(-Group1, names_to = "Group2", values_to = "p_value") %>% 
    filter(Group1 < Group2) # Avoid duplicates and self-comparisons
  
  p_values_long <- p_values_long_temp %>%
    mutate(
      q_value = if (nrow(p_values_long_temp) == 1) {
        # For single inter-group comparison, q-value equals p-value
        p_value
      } else {
        # For multiple inter-group comparisons, apply BH correction
        p.adjust(p_value, method = "BH")
      },
      single_group_comparison = (nrow(p_values_long_temp) == 1)
    )

  # Store results
  all_group_comparisons[[df_name]] <- list(
    auc_ci_summary = auc_ci_summary,
    pairwise_tests = p_values_long
  )
}

# --- Export Cross-Group Comparison Results ---
if (length(all_group_comparisons) > 0) {
  # Create a single Excel file with sheets for each type of result
  export_list_group_comp <- list()
  
  # Add AUC/CI summaries for each dataset
  auc_sheets <- purrr::map(all_group_comparisons, "auc_ci_summary")
  names(auc_sheets) <- paste("C_AUC", names(auc_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, auc_sheets)
  
  # Add pairwise test results for each dataset
  pairwise_sheets <- purrr::map(all_group_comparisons, "pairwise_tests")
  names(pairwise_sheets) <- paste("C_PVal", names(pairwise_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, pairwise_sheets)
  
  writexl::write_xlsx(export_list_group_comp, path = file.path(output_dir, "Group_Comparison_Summary.xlsx"))
  cat("\nCross-proximity method comparison results exported to Group_Comparison_Summary.xlsx\n")
} else {
  cat("\nCould not generate any cross-proximity method comparison results.\n")
}

} else {
  cat("\nInter-group comparison analysis disabled, skipping related calculations.\n")
}

# --- New Section: Type 1 Robustness (Sample Intersection Fairness) Analysis ---
if (enable_robustness_analysis) {
  cat("\nStarting Type 1 robustness analysis (sample intersection fairness)...\n")

# Define the baseline model again for clarity
baseline_model_for_intersection <- "2D_Abundance_logFC"

# List to store the intersection comparison results
all_intersection_comparisons <- list()

for (df_name in names(results_by_df)) {
  
  current_groups_results <- results_by_df[[df_name]]
  
  # --- Find intersecting proteins for the current df_name ---
  list_of_protein_ids <- current_groups_results %>% 
    purrr::map(~.x$node_attributes$Gene)
  
  intersecting_proteins <- Reduce(intersect, list_of_protein_ids)
  
  cat(paste("\nFor", df_name, ", found", length(intersecting_proteins), "proteins common to all groups for intersection analysis.\n"))
  
  if (length(intersecting_proteins) < 20) { # Heuristic check
      warning(paste("Skipping intersection analysis for", df_name, "due to very small protein intersection (", length(intersecting_proteins), ")"))
      next
  }

  # --- Re-run model on intersected data for each group ---
  roc_objects_intersect <- list()
  
  for(result in current_groups_results) {
    group_name <- result$group_name
    
    # Filter the original full data to the intersection
    intersected_data <- result$node_attributes %>% 
      filter(Gene %in% intersecting_proteins)
    
    # Check if there are enough samples in both classes
    if (n_distinct(intersected_data$plot_group) < 2 || any(table(intersected_data$plot_group) < 5)) {
      warning(paste("Skipping group", group_name, "in intersection analysis for", df_name, "due to insufficient data post-intersection."))
      next
    }
    
    # Re-fit the baseline model on the intersected data
    logFC_col <- group_info[[group_name]]$logFC
    formula_str <- paste("plot_group ~ mean_group_medianNorm +", logFC_col)
    
    tryCatch({
      model_intersect <- glm(as.formula(formula_str), data = intersected_data, family = "binomial")
      predictor_intersect <- predict(model_intersect, type = "response")
      
      # Calculate and store the new ROC object
      roc_objects_intersect[[group_name]] <- roc(
        response = intersected_data$plot_group, 
        predictor = predictor_intersect, 
        levels = c("Other", "TP"), 
        quiet = TRUE, 
        direction = "<"
      )
    }, error = function(e) {
      warning(paste("Model re-fitting failed for group", group_name, "in intersection analysis for", df_name, ":", e$message))
    })
  }
  
  if (length(roc_objects_intersect) < 2) {
    warning(paste("Skipping intersection comparison for", df_name, "as it has fewer than 2 valid groups after re-fitting."))
    next
  }
  
  # --- Perform pairwise DeLong's test on intersection ROCs ---
  group_names_intersect <- names(roc_objects_intersect)
  p_matrix_intersect <- matrix(NA, nrow = length(group_names_intersect), ncol = length(group_names_intersect),
                               dimnames = list(group_names_intersect, group_names_intersect))
  
  group_pairs_intersect <- combn(group_names_intersect, 2, simplify = FALSE)
  
  p_values_intersect <- purrr::map_dbl(group_pairs_intersect, ~ {
    tryCatch({
      roc.test(roc_objects_intersect[[.x[1]]], roc_objects_intersect[[.x[2]]], method = "delong")$p.value
    }, error = function(e) { NA })
  })
  
  for (i in seq_along(group_pairs_intersect)) {
    pair <- group_pairs_intersect[[i]]
    p_matrix_intersect[pair[1], pair[2]] <- p_values_intersect[i]
    p_matrix_intersect[pair[2], pair[1]] <- p_values_intersect[i]
  }
  diag(p_matrix_intersect) <- 1.0
  
  # --- Consolidate and store results ---
  auc_ci_intersect <- purrr::map_df(roc_objects_intersect, ~{
    ci_vals <- pROC::ci.auc(.x, method = "delong")
    tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
  }, .id = "GroupName")
  
  p_values_long_intersect <- as.data.frame(p_matrix_intersect) %>% 
    tibble::rownames_to_column("Group1") %>% 
    tidyr::pivot_longer(-Group1, names_to = "Group2", values_to = "p_value") %>% 
    filter(Group1 < Group2) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
    
  all_intersection_comparisons[[df_name]] <- list(
    auc_ci_summary = auc_ci_intersect,
    pairwise_tests = p_values_long_intersect
  )
}

  # --- Export Intersection Analysis Results ---
if (length(all_intersection_comparisons) > 0) {
  export_list_intersect <- list()
  
  auc_sheets_intersect <- purrr::map(all_intersection_comparisons, "auc_ci_summary")
  names(auc_sheets_intersect) <- paste("I_AUC", names(auc_sheets_intersect), sep = "_")
  export_list_intersect <- c(export_list_intersect, auc_sheets_intersect)
  
  pairwise_sheets_intersect <- purrr::map(all_intersection_comparisons, "pairwise_tests")
  names(pairwise_sheets_intersect) <- paste("I_PVal", names(pairwise_sheets_intersect), sep = "_")
  export_list_intersect <- c(export_list_intersect, pairwise_sheets_intersect)
  
  writexl::write_xlsx(export_list_intersect, path = file.path(output_dir, "Group_Comparison_Intersection_Summary.xlsx"))
  cat("\nSample intersection fairness analysis results exported to Group_Comparison_Intersection_Summary.xlsx\n")
} else {
  cat("\nCould not generate any sample intersection fairness analysis results.\n")
}

} else {
  cat("\nRobustness analysis disabled, skipping sample intersection fairness analysis.\n")
}

# ===================================================================
#
###    Phase 1.2: Visualization (Isolated Data) #####
#
# ===================================================================

if (enable_visualization) {
  cat("\n\n--- Starting Phase 1.2: Result Visualization ---\n")

  # To prevent the analysis in this section from accidentally modifying the original
  # `all_results` and `group_info` objects, which are needed for Phase 2
  # and Phase 3, we create deep copies of them specifically for this section.
  all_results_p1_2 <- rlang::duplicate(all_results, shallow = FALSE)
  group_info_p1_2 <- rlang::duplicate(group_info, shallow = FALSE)


# Create a directory for plots
plots_dir <- file.path(output_dir, "Plots")
dir.create(plots_dir, showWarnings = FALSE)

# --- 1. Visualization for AUC_summary.xlsx: Faceted Forest Plot ---
if (enable_plot_auc_summary) {
  cat("\nCreating grouped forest plot for AUC_summary.xlsx...\n")

  # Check if the final summary data frame exists and is not empty
  if (exists("final_auc_summary") && nrow(final_auc_summary) > 0) {

    # Prepare data for plotting directly from the final_auc_summary dataframe
    plot_data_auc <- final_auc_summary %>%
      # Exclude the random control as it clutters the plot
      filter(ModelName != "Control") %>%
      # Create a more descriptive label for faceting
      mutate(Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n")) %>%
      # For each group, find the AUC of the traditional method to use as a reference line
      group_by(Analysis_Group) %>%
      mutate(
        traditional_auc = tryCatch(AUC[ModelName == "logFC>0.5 & FDR<0.05"][1], error = function(e) NA),
        # The q-value is already calculated, so we just create the significance flag
        is_significant = if_else(
          !is.na(q_value_vs_traditional) & q_value_vs_traditional < 0.05, 
          "Significant (q < 0.05)", 
          "Not Significant"
        )
      ) %>%
      ungroup() %>%
      # Make ModelName a factor to control plotting order
      mutate(ModelName = factor(ModelName, levels = rev(c("logFC>0.5 & FDR<0.05", setdiff(unique(ModelName), "logFC>0.5 & FDR<0.05")))))

    # Generate the plot
    forest_plot <- ggplot(plot_data_auc, aes(y = ModelName, x = AUC, xmin = AUC_Lower_CI, xmax = AUC_Upper_CI)) +
      # Add reference line for the traditional method's AUC
      geom_vline(aes(xintercept = traditional_auc), linetype = "dashed", color = "grey50", na.rm = TRUE) +
      # Add points and error bars, colored by significance
      geom_point(aes(color = is_significant), size = 2.5) +
      geom_errorbarh(aes(color = is_significant), height = 0.2) +
      # Facet by the new combined label
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      # Customize colors and labels
      scale_color_manual(values = c("Not Significant" = "black", "Significant (q < 0.05)" = "#d62728")) +
      labs(
        title = "Model Performance Comparison within Each Analysis Group",
        subtitle = "AUCs compared to the traditional method (dashed line). Red models are significantly better (q < 0.05).",
        x = "Area Under Curve (AUC) with 95% Confidence Interval",
        y = "Model",
        color = "Significance vs. Traditional"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      # Set x-axis limits for better visualization
      coord_cartesian(xlim = c(0.4, 1.0))

    # Save the plot
    ggsave(
      filename = file.path(plots_dir, "1_AUC_Faceted_Forest_Plot.pdf"),
      plot = forest_plot,
      width = 16,
      height = 12,
      bg = "white"
    )

    cat("Grouped forest plot saved to:", file.path(plots_dir, "1_AUC_Faceted_Forest_Plot.pdf"), "\n")

  } else {
    cat("Warning: AUC_summary.xlsx not found, skipping grouped forest plot generation.\n")
  }
} else {
  cat("\nAUC_summary grouped forest plot generation disabled, skipping this plot.\n")
}


# --- 1b. Visualization for AUC_Summary_b.xlsx: Faceted Forest Plot vs 1D_logFC ---
if (enable_plot_auc_summary_vs_logFC) {
  cat("\nCreating grouped forest plot (vs 1D_logFC) for AUC_Summary_b.xlsx...\n")

  # Check if the AUC_Summary_b data exists
  auc_summary_b_path <- file.path(output_dir, "AUC_Summary_b.xlsx")
  if (file.exists(auc_summary_b_path)) {
    auc_summary_b_data <- readxl::read_excel(auc_summary_b_path)
    
    if (nrow(auc_summary_b_data) > 0) {
      
      # Prepare data for plotting
      # First create 1D_logFC AUC reference table
      logFC_auc_reference <- final_auc_summary %>%
        filter(ModelName == "1D_logFC") %>%
        select(Analysis_Group, logFC_auc = AUC)
      
      plot_data_auc_b <- auc_summary_b_data %>%
        # Create a more descriptive label for faceting
        mutate(Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n")) %>%
        # Join with 1D_logFC AUC as reference line
        left_join(logFC_auc_reference, by = "Analysis_Group") %>%
        # Create significance flag based on q_value_vs_1D_logFC
        mutate(
          is_significant = if_else(
            !is.na(q_value_vs_1D_logFC) & q_value_vs_1D_logFC < 0.05,
            "Significant (q < 0.05)",
            "Not Significant"
          )
        ) %>%
        # Make ModelName a factor to control plotting order
        mutate(ModelName = factor(ModelName, levels = rev(c("logFC>0.5 & FDR<0.05", "1D_Abundance", "2D_Abundance_logFC", setdiff(unique(ModelName), c("logFC>0.5 & FDR<0.05", "1D_Abundance", "2D_Abundance_logFC"))))))
      
      # Generate the forest plot
      forest_plot_b <- ggplot(plot_data_auc_b, aes(y = ModelName, x = AUC, xmin = AUC_Lower_CI, xmax = AUC_Upper_CI)) +
        # Add reference line for 1D_logFC's AUC
        geom_vline(aes(xintercept = logFC_auc), linetype = "dashed", color = "grey50", na.rm = TRUE) +
        # Add points and error bars, colored by significance
        geom_point(aes(color = is_significant), size = 2.5) +
        geom_errorbarh(aes(color = is_significant), height = 0.2) +
        # Facet by the combined label
        facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
        # Customize colors and labels
        scale_color_manual(values = c("Not Significant" = "black", "Significant (q < 0.05)" = "#d62728")) +
        labs(
          title = "Model Performance Comparison vs 1D_logFC within Each Analysis Group",
          subtitle = "AUCs compared to the 1D_logFC model (dashed line). Red models are significantly better (q < 0.05).",
          x = "Area Under Curve (AUC) with 95% Confidence Interval",
          y = "Model",
          color = "Significance vs. 1D_logFC"
        ) +
        theme_bw(base_size = 12) +
        theme(
          legend.position = "bottom",
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 7, face = "bold"),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        ) +
        # Set x-axis limits for better visualization
        coord_cartesian(xlim = c(0.4, 1.0))
      
      # Save the plot
      ggsave(
        filename = file.path(plots_dir, "1b_AUC_Faceted_Forest_Plot_vs_1D_logFC.pdf"),
        plot = forest_plot_b,
        width = 16,
        height = 12,
        bg = "white"
      )
      
      cat("Grouped forest plot vs 1D_logFC saved to:", file.path(plots_dir, "1b_AUC_Faceted_Forest_Plot_vs_1D_logFC.pdf"), "\n")
      
    } else {
      cat("Warning: AUC_Summary_b.xlsx is empty, skipping figure 1b generation.\n")
    }
  } else {
    cat("Warning: AUC_Summary_b.xlsx not found, skipping figure 1b generation.\n")
  }
} else {
  cat("\nGrouped forest plot vs 1D_logFC generation disabled, skipping this plot.\n")
}


# --- 2. Visualization for Group_Comparison_Summary.xlsx: Bar Plots and Heatmaps ---
group_comp_path <- file.path(output_dir, "Group_Comparison_Summary.xlsx")
if (enable_plot_group_comparison) {
  cat("\nCreating visualization plots for Group_Comparison_Summary.xlsx...\n")

  if (file.exists(group_comp_path)) {
    
    # Get all sheet names from the Excel file
    sheet_names <- readxl::excel_sheets(group_comp_path)
    
    # Identify the base names for each analysis (e.g., "ForStep19_APEX")
    base_names <- unique(sub("^(C_AUC_|C_PVal_)", "", sheet_names))

    for (base_name in base_names) {
      cat(paste("  - Processing group comparison:", base_name, "\n"))
      
      auc_sheet_name <- paste0("C_AUC_", base_name)
      pval_sheet_name <- paste0("C_PVal_", base_name)
      
      # --- a) Ordered Bar Plot with CIs ---
      if (auc_sheet_name %in% sheet_names) {
        auc_data <- readxl::read_excel(group_comp_path, sheet = auc_sheet_name)
        
        if (nrow(auc_data) > 0) {
          # Order according to selected_groups, ensuring only actual existing groups are included
          groups_from_data <- unique(auc_data$GroupName)
          ordered_groups <- intersect(selected_groups, groups_from_data)
          missing_groups <- setdiff(groups_from_data, selected_groups)
          final_group_order <- c(ordered_groups, missing_groups)
          
          # Convert GroupName to ordered factor
          auc_data <- auc_data %>%
            mutate(GroupName = factor(GroupName, levels = final_group_order))
          
          bar_plot <- ggplot(auc_data, aes(x = GroupName, y = AUC)) +
            geom_col(fill = "#0072B2", width = 0.7) +
            geom_errorbar(aes(ymin = AUC_Lower_CI, ymax = AUC_Upper_CI), width = 0.2, color = "gray30") +
            labs(
              title = paste("Cross-Group AUC Comparison for:", base_name),
              subtitle = paste("Using baseline model '2D_Abundance_logFC' - Groups ordered by selection"),
              x = "Proximity Labeling Method (Group)",
              y = "AUC (with 95% CI)"
            ) +
            theme_minimal(base_size = 12) +
            theme(
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5, size = 3.5, color = "black") + 
            ylim(NA, max(auc_data$AUC_Upper_CI, na.rm = TRUE) * 1.05) # Adjust y-axis limit

          ggsave(
            filename = file.path(plots_dir, paste0("2a_Group_Comparison_Bar_Plot_", base_name, ".pdf")),
            plot = bar_plot,
            width = 10, height = 8, bg = "white"
          )
        }
      }
      
      # --- b) Pairwise Comparison Heatmap ---
      if (pval_sheet_name %in% sheet_names) {
        pval_data <- readxl::read_excel(group_comp_path, sheet = pval_sheet_name)
        
        if (nrow(pval_data) > 0) {
          # Detect whether it's a single comparison
          is_single_comparison <- nrow(pval_data) == 1 || 
                                 (exists("single_group_comparison") && any(pval_data$single_group_comparison == TRUE, na.rm = TRUE))
          
          # Determine whether to use p-value or q-value, and corresponding labels
          value_column <- if (is_single_comparison) "p_value" else "q_value"
          value_label <- if (is_single_comparison) "p-value" else "q-value (BH-adj.)"
          subtitle_text <- if (is_single_comparison) {
            "Single comparison - using p-values (no multiple testing correction needed)"
          } else {
            "Multiple comparisons - using BH-adjusted q-values"
          }
          
          # Use selected_groups order, ensuring only actual existing groups are included
          all_groups_from_data <- union(pval_data$Group1, pval_data$Group2)
          ordered_groups <- intersect(selected_groups, all_groups_from_data)
          # If there are groups in data but not in selected_groups, append to the end
          missing_groups <- setdiff(all_groups_from_data, selected_groups)
          all_groups <- c(ordered_groups, missing_groups)
          
          full_pval_data <- pval_data %>% 
            # Add the reverse pairs for a symmetric matrix
            bind_rows(rename(., Group1 = Group2, Group2 = Group1)) %>% 
            # Add diagonal comparisons
            bind_rows(tibble(Group1 = all_groups, Group2 = all_groups, 
                            p_value = NA, q_value = NA, single_group_comparison = NA)) %>%
            # Convert groups to ordered factors according to selected_groups order
            mutate(
              Group1 = factor(Group1, levels = all_groups),
              Group2 = factor(Group2, levels = all_groups)
            )

          # Dynamically select fill value
          heatmap_plot <- ggplot(full_pval_data, aes(x = Group1, y = Group2, fill = .data[[value_column]])) +
            geom_tile(color = "white", size = 0.5) +
            # Use a color scale that highlights significant values
            scale_fill_gradientn(
              colors = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF"),
              na.value = "grey80",
              limits = c(0, 1),
              name = value_label
            ) +
            geom_text(aes(label = ifelse(is.na(.data[[value_column]]), "", sprintf("%.3f", .data[[value_column]]))), 
                     size = 3, color = "white") +
            coord_fixed() + # Ensure tiles are square
            # Set axis order, cancel default alphabetical sorting
            scale_x_discrete(limits = all_groups) +
            scale_y_discrete(limits = rev(all_groups)) + # Reverse y-axis to match matrix display convention
            labs(
              title = paste("Pairwise DeLong Test", value_label, "for:", base_name),
              subtitle = subtitle_text,
              x = "", y = ""
            ) +
            theme_minimal(base_size = 12) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right"
            )

          ggsave(
            filename = file.path(plots_dir, paste0("2b_Group_Comparison_Heatmap_", base_name, ".pdf")),
            plot = heatmap_plot,
            width = 10, height = 9, bg = "white"
          )
        }
      }
    }
    cat("Cross-proximity method comparison visualization plots saved.\n")
  } else {
    cat("Warning: Group_Comparison_Summary.xlsx not found, skipping related plot generation.\n")
  }
} else {
  cat("\nCross-proximity method comparison plot generation disabled, skipping this plot.\n")
}


# --- 3. Visualization for Intersection Analysis: Dumbbell Plot ---
intersection_comp_path <- file.path(output_dir, "Group_Comparison_Intersection_Summary.xlsx")
if (enable_plot_intersection_dumbbell) {
  cat("\nCreating dumbbell plot for sample intersection fairness analysis...\n")

  # We need both files to create the comparison plot
  if (file.exists(group_comp_path) && file.exists(intersection_comp_path)) {

    # Get sheet names from both files
    union_sheets <- readxl::excel_sheets(group_comp_path)
    intersection_sheets <- readxl::excel_sheets(intersection_comp_path)
    
    # Identify the base names for which both analyses were successful
    union_base_names <- unique(sub("^C_AUC_", "", union_sheets[grep("^C_AUC_", union_sheets)]))
    intersection_base_names <- unique(sub("^I_AUC_", "", intersection_sheets[grep("^I_AUC_", intersection_sheets)]))
    common_base_names <- intersect(union_base_names, intersection_base_names)

    for (base_name in common_base_names) {
      cat(paste("  - Creating dumbbell plot for:", base_name, "\n"))
      
      # Load data from both analyses
      auc_union <- readxl::read_excel(group_comp_path, sheet = paste0("C_AUC_", base_name)) %>%
        mutate(AnalysisType = "Union (All Proteins)")
      
      auc_intersection <- readxl::read_excel(intersection_comp_path, sheet = paste0("I_AUC_", base_name)) %>%
        mutate(AnalysisType = "Intersection (Common Proteins)")
        
      # Combine the data
      combined_auc_data <- bind_rows(auc_union, auc_intersection)

      if (nrow(combined_auc_data) > 0) {
        # Align groups according to selected_groups order and keep only those that exist
        groups_from_data <- unique(combined_auc_data$GroupName)
        ordered_groups <- intersect(selected_groups, groups_from_data)
        missing_groups <- setdiff(groups_from_data, selected_groups)
        final_group_order <- c(ordered_groups, missing_groups)
        
        # Convert GroupName into an ordered factor
        combined_auc_data <- combined_auc_data %>%
          mutate(GroupName = factor(GroupName, levels = final_group_order))
        
        dumbbell_plot <- ggplot(combined_auc_data, aes(x = AUC, y = GroupName)) +
          # The line connecting the two points
          geom_line(aes(group = GroupName), color = "gray", linewidth = 1.5) +
          # The points for each analysis type
          geom_point(aes(color = AnalysisType), size = 4) +
          # Customize colors
          scale_color_manual(values = c("Union (All Proteins)" = "#1f77b4", "Intersection (Common Proteins)" = "#ff7f0e")) +
          labs(
            title = paste("Robustness to Protein Set for:", base_name),
            subtitle = "Comparing AUC on all available proteins (Union) vs. only common proteins (Intersection) - Groups ordered by selection",
            x = "AUC Score",
            y = "Proximity Labeling Method (Group)",
            color = "Analysis Set"
          ) +
          theme_bw(base_size = 12) +
          theme(
            legend.position = "top",
            panel.grid.major.y = element_line(linetype = "dashed")
          )

        ggsave(
          filename = file.path(plots_dir, paste0("3_Intersection_Dumbbell_Plot_", base_name, ".pdf")),
          plot = dumbbell_plot,
          width = 11, height = 8, bg = "white"
        )
      }
    }
    cat("Sample intersection fairness analysis dumbbell plots saved.\n")
  } else {
    cat("Warning: Group_Comparison_Summary.xlsx or Group_Comparison_Intersection_Summary.xlsx not found, skipping dumbbell plot generation.\n")
  }
} else {
  cat("\nSample intersection fairness dumbbell plot generation disabled, skipping this plot.\n")
}


# --- 4. Visualization for Model Choice Robustness: Performance vs. Stability Plot ---
robustness_path <- file.path(output_dir, "Group_Robustness_Summary.xlsx")
if (enable_plot_robustness_scatter) {
  cat("\nCreating performance-stability scatter plot for model choice robustness...\n")

  if (file.exists(robustness_path)) {
    robustness_data <- readxl::read_excel(robustness_path)

    if (nrow(robustness_data) > 0) {
      # Create a new column for combined labels
      robustness_data <- robustness_data %>%
        mutate(Combined_Label = paste(Df_Name, Analysis_Group, sep = " - "))

      perf_stability_plot <- ggplot(robustness_data, aes(x = AUC_SD, y = Mean_AUC)) +
        geom_point(aes(color = Combined_Label), size = 5, alpha = 0.8) +
        # Add text labels with a bit of offset
        ggrepel::geom_text_repel(aes(label = Combined_Label), box.padding = 0.5, max.overlaps = Inf) +
        labs(
          title = "Performance vs. Stability of Proximity Labeling Methods",
          subtitle = "Across all defined 1D, 2D, and 3D models",
          x = "AUC Standard Deviation (Stability, lower is better)",
          y = "Mean AUC (Performance, higher is better)",
          caption = "Each point represents a method's average performance and variability across different models."
        ) +
        scale_color_viridis_d() + # Use a nice color scale
        theme_bw(base_size = 12) +
        theme(legend.position = "none") # Hide legend as labels are on the plot

      ggsave(
        filename = file.path(plots_dir, "4_Performance_vs_Stability_Scatter_Plot.pdf"),
        plot = perf_stability_plot,
        width = 12, height = 10, bg = "white"
      )
      cat("Performance-stability scatter plot saved.\n")
    } else {
      cat("Warning: Group_Robustness_Summary.xlsx is empty, skipping scatter plot generation.\n")
    }
  } else {
    cat("Warning: Group_Robustness_Summary.xlsx not found, skipping scatter plot generation.\n")
  }
} else {
  cat("\nModel choice robustness scatter plot generation disabled, skipping this plot.\n")
}



# --- 5. Visualization for Controlled Variable Comparison ---
controlled_comp_path <- file.path(output_dir, "Controlled_Variable_Comparison.xlsx")
if (enable_plot_controlled_comparison) {
  cat("\nCreating gain plots for controlled variable comparison...\n")

if (file.exists(controlled_comp_path)) {
  controlled_comp_data <- readxl::read_excel(controlled_comp_path)
  
  if (nrow(controlled_comp_data) > 0) {
    
    # --- 5a. Fixed Sensitivity - Specificity Gain Plot ---
    cat("  - Creating Specificity gain plot at fixed Sensitivity...\n")
    
    # Filter specified models
    spec_gain_data <- controlled_comp_data %>%
      filter(ModelName %in% models_for_controlled_comparison) %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Specificity_Gain_Q_Value < 0.05, "Significant (q<0.05)", "Not Significant"),
        ModelName = factor(ModelName, levels = rev(models_for_controlled_comparison))
      )
    
    spec_gain_plot <- ggplot(spec_gain_data, aes(y = ModelName, x = Specificity_Gain, 
                                                   xmin = Specificity_Gain_CI_Lower, 
                                                   xmax = Specificity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("Not Significant" = "black", "Significant (q<0.05)" = "#d62728")) +
      labs(
        title = "Controlled Variable Comparison: Specificity Gain at Fixed Sensitivity",
        subtitle = "Specificity gain of continuous models at traditional method's Sensitivity level (positive values indicate improvement)",
        x = "Specificity Gain (relative to traditional method)",
        y = "Model",
        color = "Statistical Significance",
        caption = "Error bars represent 95% Bootstrap confidence intervals"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "5a_Controlled_Comparison_Specificity_Gain.pdf"),
      plot = spec_gain_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 5b. Fixed Specificity - Sensitivity Gain Plot ---
    cat("  - Creating Sensitivity gain plot at fixed Specificity...\n")
    
    # Filter specified models
    sens_gain_data <- controlled_comp_data %>%
      filter(ModelName %in% models_for_controlled_comparison) %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "Significant (q<0.05)", "Not Significant"),
        ModelName = factor(ModelName, levels = rev(models_for_controlled_comparison))
      )
    
    sens_gain_plot <- ggplot(sens_gain_data, aes(y = ModelName, x = Sensitivity_Gain, 
                                                   xmin = Sensitivity_Gain_CI_Lower, 
                                                   xmax = Sensitivity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("Not Significant" = "black", "Significant (q<0.05)" = "#d62728")) +
      labs(
        title = "Controlled Variable Comparison: Sensitivity Gain at Fixed Specificity",
        subtitle = "Sensitivity gain of continuous models at traditional method's Specificity level (positive values indicate improvement)",
        x = "Sensitivity Gain (relative to traditional method)",
        y = "Model",
        color = "Statistical Significance",
        caption = "Error bars represent 95% Bootstrap confidence intervals"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "5b_Controlled_Comparison_Sensitivity_Gain.pdf"),
      plot = sens_gain_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 5c. Combined Gain Comparison (Scatter Plot) ---
    cat("  - Creating combined gain comparison scatter plot...\n")
    
    # Filter specified models
    combined_gain_data <- controlled_comp_data %>%
      filter(ModelName %in% models_for_controlled_comparison) %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = " - "),
        Both_Significant = if_else(
          Specificity_Gain_Q_Value < 0.05 & Sensitivity_Gain_Q_Value < 0.05,
          "Both Significant", 
          if_else(
            Specificity_Gain_Q_Value < 0.05 | Sensitivity_Gain_Q_Value < 0.05,
            "Partially Significant",
            "Neither Significant"
          )
        ),
        ModelName = factor(ModelName, levels = models_for_controlled_comparison)
      )
    
    combined_gain_plot <- ggplot(combined_gain_data, 
                                  aes(x = Specificity_Gain, y = Sensitivity_Gain,
                                      xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper,
                                      ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper,
                                      color = Both_Significant, shape = ModelName)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      # Add horizontal error bar (CI for Specificity gain)
      geom_errorbarh(aes(xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper),
                     height = 0, alpha = 0.5, size = 0.5) +
      # Add vertical error bar (CI for Sensitivity gain)
      geom_errorbar(aes(ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper),
                    width = 0, alpha = 0.5, size = 0.5) +
      # Points on top of error bars
      geom_point(size = 4, alpha = 0.8) +
      facet_wrap(~ Combined_Label, scales = "free") +
      scale_color_manual(
        values = c("Both Significant" = "#2ca02c", "Partially Significant" = "#ff7f0e", "Neither Significant" = "grey60")
      ) +
      labs(
        title = "Controlled Variable Comparison: Combined Comparison of Sensitivity and Specificity Gains",
        subtitle = "Performance gains relative to traditional method (logFC>0.5 & FDR<0.05) (with 95% CI error bars)",
        x = "Specificity Gain (at fixed Sensitivity)",
        y = "Sensitivity Gain (at fixed Specificity)",
        color = "Significance Status",
        shape = "Model Type",
        caption = "Upper-right quadrant indicates both metrics are better than traditional method; error bars represent 95% Bootstrap confidence intervals"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim, ylim = controlled_comparison_ylim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(breaks = seq(controlled_comparison_ylim[1], controlled_comparison_ylim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      filename = file.path(plots_dir, "5c_Controlled_Comparison_Combined_Gain.pdf"),
      plot = combined_gain_plot,
      width = 16, height = 12, bg = "white"
    )
    
    cat("Controlled variable comparison visualization plots saved.\n")
  } else {
    cat("Warning: Controlled_Variable_Comparison.xlsx is empty, skipping related plot generation.\n")
  }
} else {
  cat("Warning: Controlled_Variable_Comparison.xlsx not found, skipping controlled variable comparison plot generation.\n")
}
} else {
  cat("\nControlled variable comparison gain plot generation disabled, skipping this plot.\n")
}


# --- 7. Visualization for Controlled Variable Comparison vs 1D_logFC ---
controlled_comp_vs_logFC_path <- file.path(output_dir, "Controlled_Variable_Comparison_vs_1D_logFC.xlsx")
if (enable_plot_controlled_comparison_vs_logFC) {
  cat("\nCreating gain plots for controlled variable comparison vs 1D_logFC model...\n")

if (file.exists(controlled_comp_vs_logFC_path)) {
  controlled_comp_vs_logFC_data <- readxl::read_excel(controlled_comp_vs_logFC_path)
  
  if (nrow(controlled_comp_vs_logFC_data) > 0) {
    
    # --- 7a. Fixed Sensitivity - Specificity Gain Plot (vs 1D_logFC) ---
    cat("  - Creating Specificity gain plot at fixed Sensitivity (vs 1D_logFC)...\n")
    
    spec_gain_logFC_data <- controlled_comp_vs_logFC_data %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Specificity_Gain_Q_Value < 0.05, "Significant (q<0.05)", "Not Significant"),
        ModelName = factor(ModelName, levels = rev(unique(ModelName)))
      )
    
    spec_gain_logFC_plot <- ggplot(spec_gain_logFC_data, aes(y = ModelName, x = Specificity_Gain, 
                                                   xmin = Specificity_Gain_CI_Lower, 
                                                   xmax = Specificity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("Not Significant" = "black", "Significant (q<0.05)" = "#d62728")) +
      labs(
        title = "Controlled Variable Comparison (vs 1D_logFC): Specificity Gain at Fixed Sensitivity",
        subtitle = "Specificity gain of other models at 1D_logFC model's Sensitivity level (positive values indicate improvement)",
        x = "Specificity Gain (relative to 1D_logFC model)",
        y = "Model",
        color = "Statistical Significance",
        caption = "Error bars represent 95% Bootstrap confidence intervals"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "7a_Controlled_Comparison_vs_1D_logFC_Specificity_Gain.pdf"),
      plot = spec_gain_logFC_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 7b. Fixed Specificity - Sensitivity Gain Plot (vs 1D_logFC) ---
    cat("  - Creating Sensitivity gain plot at fixed Specificity (vs 1D_logFC)...\n")
    
    sens_gain_logFC_data <- controlled_comp_vs_logFC_data %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "Significant (q<0.05)", "Not Significant"),
        ModelName = factor(ModelName, levels = rev(unique(ModelName)))
      )
    
    sens_gain_logFC_plot <- ggplot(sens_gain_logFC_data, aes(y = ModelName, x = Sensitivity_Gain, 
                                                   xmin = Sensitivity_Gain_CI_Lower, 
                                                   xmax = Sensitivity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("Not Significant" = "black", "Significant (q<0.05)" = "#d62728")) +
      labs(
        title = "Controlled Variable Comparison (vs 1D_logFC): Sensitivity Gain at Fixed Specificity",
        subtitle = "Sensitivity gain of other models at 1D_logFC model's Specificity level (positive values indicate improvement)",
        x = "Sensitivity Gain (relative to 1D_logFC model)",
        y = "Model",
        color = "Statistical Significance",
        caption = "Error bars represent 95% Bootstrap confidence intervals"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "7b_Controlled_Comparison_vs_1D_logFC_Sensitivity_Gain.pdf"),
      plot = sens_gain_logFC_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 7c. Combined Gain Comparison (Scatter Plot vs 1D_logFC) ---
    cat("  - Creating combined gain comparison scatter plot (vs 1D_logFC)...\n")
    
    combined_gain_logFC_data <- controlled_comp_vs_logFC_data %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = " - "),
        Both_Significant = if_else(
          Specificity_Gain_Q_Value < 0.05 & Sensitivity_Gain_Q_Value < 0.05,
          "Both Significant", 
          if_else(
            Specificity_Gain_Q_Value < 0.05 | Sensitivity_Gain_Q_Value < 0.05,
            "Partially Significant",
            "Neither Significant"
          )
        ),
        ModelName = factor(ModelName, levels = unique(ModelName))
      )
    
    combined_gain_logFC_plot <- ggplot(combined_gain_logFC_data, 
                                  aes(x = Specificity_Gain, y = Sensitivity_Gain,
                                      xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper,
                                      ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper,
                                      color = Both_Significant, shape = ModelName)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      # Add horizontal error bar (CI for Specificity gain)
      geom_errorbarh(aes(xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper),
                     height = 0, alpha = 0.5, size = 0.5) +
      # Add vertical error bar (CI for Sensitivity gain)
      geom_errorbar(aes(ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper),
                    width = 0, alpha = 0.5, size = 0.5) +
      # Points on top of error bars
      geom_point(size = 4, alpha = 0.8) +
      facet_wrap(~ Combined_Label, scales = "free") +
      scale_color_manual(
        values = c("Both Significant" = "#2ca02c", "Partially Significant" = "#ff7f0e", "Neither Significant" = "grey60")
      ) +
      labs(
        title = "Controlled Variable Comparison (vs 1D_logFC): Combined Comparison of Sensitivity and Specificity Gains",
        subtitle = "Performance gains relative to 1D_logFC model (with 95% CI error bars)",
        x = "Specificity Gain (at fixed Sensitivity)",
        y = "Sensitivity Gain (at fixed Specificity)",
        color = "Significance Status",
        shape = "Model Type",
        caption = "Upper-right quadrant indicates both metrics are better than 1D_logFC model; error bars represent 95% Bootstrap confidence intervals"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim, ylim = controlled_comparison_ylim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(breaks = seq(controlled_comparison_ylim[1], controlled_comparison_ylim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      filename = file.path(plots_dir, "7c_Controlled_Comparison_vs_1D_logFC_Combined_Gain.pdf"),
      plot = combined_gain_logFC_plot,
      width = 16, height = 12, bg = "white"
    )
    
    cat("Controlled variable comparison visualization plots vs 1D_logFC model saved.\n")
  } else {
    cat("Warning: Controlled_Variable_Comparison_vs_1D_logFC.xlsx is empty, skipping related plot generation.\n")
  }
} else {
  cat("Warning: Controlled_Variable_Comparison_vs_1D_logFC.xlsx not found, skipping controlled variable comparison plot generation vs 1D_logFC.\n")
}
} else {
  cat("\nControlled variable gain plot generation vs 1D_logFC disabled, skipping this plot.\n")
}


# --- 6. Visualization for Candidate Protein Localization Distribution ---
candidates_path <- file.path(output_dir, "detailed_candidates_by_model.xlsx")
if (enable_plot_localization_distribution) {
  cat("\nCreating stacked bar plots for candidate protein localization distribution...\n")

# Read candidate protein data

if (file.exists(candidates_path)) {
  candidates_data <- readxl::read_excel(candidates_path)
  
  if (nrow(candidates_data) > 0) {
    
    # --- 6a. Prepare data: Statistics for each model+group localization distribution ---
    cat("  - Calculating candidate protein localization distribution statistics...\n")
    
    # Keep only candidates that passed threshold (PassThreshold == "Pass")
    candidates_passed <- candidates_data %>%
      filter(PassThreshold == "Pass")
    
    # Count proteins for each model, group, and localization category
    localization_summary <- candidates_passed %>%
      group_by(DataAnalysisMethod, ExperimentGroup, Model, MultiBait_Localization) %>%
      summarise(
        Count = n(),
        .groups = "drop"
      ) %>%
      # Calculate percentage within each model+group combination
      group_by(DataAnalysisMethod, ExperimentGroup, Model) %>%
      mutate(
        Total_Count = sum(Count),
        Percentage = Count / Total_Count * 100,
        # Create combined label for plotting
        Model_Group = paste(Model, ExperimentGroup, sep = " - ")
      ) %>%
      ungroup()
    
    # Export detailed statistics
    localization_summary_export <- localization_summary %>%
      select(
        DataAnalysisMethod,
        ExperimentGroup, 
        Model,
        MultiBait_Localization,
        Count,
        Total_Count,
        Percentage
      ) %>%
      arrange(DataAnalysisMethod, ExperimentGroup, Model, desc(Count))
    
    writexl::write_xlsx(
      list(Localization_Distribution = localization_summary_export),
      path = file.path(output_dir, "Candidate_Localization_Distribution.xlsx")
    )
    cat("  - Localization distribution statistics exported to Candidate_Localization_Distribution.xlsx\n")
    
    # --- 6b. Create stacked bar plots for each DataAnalysisMethod ---
    
    # Define consistent color scheme (based on example plots)
    localization_colors <- c(
      "Nuclear&SGs" = "#4472C4",           # Blue
      "Nuclear&Cytosol&SGs" = "#70AD47",  # Cyan
      "SGs&Other" = "#FFC000",             # Yellow
      "Cytosol&SGs" = "#C00000"            # Red
    )
    
    # Get all unique DataAnalysisMethod
    data_methods <- unique(localization_summary$DataAnalysisMethod)
    
    for (data_method in data_methods) {
      cat(paste("  - Creating localization distribution plot for", data_method, "...\n"))
      
      # Filter data for current data method
      plot_data_loc <- localization_summary %>%
        filter(DataAnalysisMethod == data_method) %>%
        # Arrange according to selected groups order
        mutate(
          ExperimentGroup = factor(ExperimentGroup, levels = selected_groups),
          # Ensure localization categories are ordered (from bottom to top)
          # Order: Red Cytosol&SGs (bottom) → Yellow SGs&Other → Green Nuclear&Cytosol&SGs → Blue Nuclear&SGs (top)
          MultiBait_Localization = factor(
            MultiBait_Localization,
            levels = c("Cytosol&SGs", "SGs&Other", "Nuclear&Cytosol&SGs", "Nuclear&SGs")
          )
        )
      
      # Create ordered labels for x-axis
      # Format: Model (Group)
      plot_data_loc <- plot_data_loc %>%
        mutate(
          X_Label = paste0(Model, "\n(", ExperimentGroup, ")"),
          # Create sort key: first by ExperimentGroup then by Model
          Sort_Key = paste(
            sprintf("%02d", as.numeric(ExperimentGroup)),
            Model
          )
        ) %>%
        arrange(Sort_Key) %>%
        mutate(
          X_Label = factor(X_Label, levels = unique(X_Label))
        )
      
      # Create stacked bar plot
      loc_barplot <- ggplot(plot_data_loc, aes(x = X_Label, y = Percentage, fill = MultiBait_Localization)) +
        geom_col(width = 0.8, color = "white", size = 0.3) +
        # Add protein count annotation on top of bars
        geom_text(
          data = plot_data_loc %>% 
            group_by(X_Label) %>% 
            summarise(
              Total_Count = first(Total_Count),
              Y_Position = 101,
              .groups = "drop"
            ),
          aes(x = X_Label, y = Y_Position, label = Total_Count),
          inherit.aes = FALSE,
          vjust = 0,
          size = 3,
          fontface = "bold"
        ) +
        # Display percentage in each stacked segment (only show >5%)
        geom_text(
          aes(label = ifelse(Percentage > 5, paste0(round(Percentage), "%"), "")),
          position = position_stack(vjust = 0.5),
          size = 2.5,
          color = "white",
          fontface = "bold"
        ) +
        scale_fill_manual(
          values = localization_colors,
          name = "Localization",
          breaks = c("Nuclear&SGs", "Nuclear&Cytosol&SGs", "SGs&Other", "Cytosol&SGs")
        ) +
        labs(
          title = paste("Candidate Protein Localization Distribution -", data_method),
          subtitle = paste("Localization distribution of threshold-passing", myTP_vector[1], "related proteins by model and experiment group"),
          x = "Model (Experiment Group)",
          y = "Relative Abundance (%)",
          caption = "Numbers on top of bars indicate total candidate proteins for that combination"
        ) +
        scale_y_continuous(
          limits = c(0, 105),
          breaks = seq(0, 100, 25),
          expand = c(0, 0)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        )
      
      # Save plot
      ggsave(
        filename = file.path(plots_dir, paste0("6_Localization_Distribution_", data_method, ".pdf")),
        plot = loc_barplot,
        width = max(14, length(unique(plot_data_loc$X_Label)) * 0.8),
        height = 10,
        bg = "white"
      )
    }
    
    # --- 6c. Create comparison plot: Compare localization preferences across different models ---
    cat("  - Creating model-to-model localization distribution comparison plot...\n")
    
    # Calculate average localization distribution for each model across all groups
    model_avg_distribution <- localization_summary %>%
      group_by(DataAnalysisMethod, Model, MultiBait_Localization) %>%
      summarise(
        Avg_Percentage = mean(Percentage, na.rm = TRUE),
        Total_Proteins = sum(Count),
        .groups = "drop"
      ) %>%
      mutate(
        # Ensure localization categories are ordered (from bottom to top)
        # Order: Red Cytosol&SGs (bottom) → Yellow SGs&Other → Green Nuclear&Cytosol&SGs → Blue Nuclear&SGs (top)
        MultiBait_Localization = factor(
          MultiBait_Localization,
          levels = c("Cytosol&SGs", "SGs&Other", "Nuclear&Cytosol&SGs", "Nuclear&SGs")
        )
      )
    
    for (data_method in data_methods) {
      plot_data_model_comp <- model_avg_distribution %>%
        filter(DataAnalysisMethod == data_method) %>%
        # Sort by model dimension
        mutate(
          Model = factor(Model, levels = c(
            "logFC>0.5 & FDR<0.05",
            setdiff(unique(Model), "logFC>0.5 & FDR<0.05")
          ))
        )
      
      model_comp_plot <- ggplot(plot_data_model_comp, 
                                 aes(x = Model, y = Avg_Percentage, fill = MultiBait_Localization)) +
        geom_col(width = 0.7, color = "white", size = 0.3) +
        geom_text(
          data = plot_data_model_comp %>% 
            group_by(Model) %>% 
            summarise(Total_Proteins = sum(Total_Proteins), .groups = "drop"),
          aes(x = Model, y = 102, label = Total_Proteins),
          inherit.aes = FALSE,
          vjust = 0,
          size = 3.5,
          fontface = "bold"
        ) +
        geom_text(
          aes(label = ifelse(Avg_Percentage > 5, paste0(round(Avg_Percentage), "%"), "")),
          position = position_stack(vjust = 0.5),
          size = 3,
          color = "white",
          fontface = "bold"
        ) +
        scale_fill_manual(
          values = localization_colors,
          name = "Localization",
          breaks = c("Nuclear&SGs", "Nuclear&Cytosol&SGs", "SGs&Other", "Cytosol&SGs")
        ) +
        labs(
          title = paste("Localization Preference Comparison Between Models -", data_method),
          subtitle = "Average localization distribution across all experimental groups for each model",
          x = "Model",
          y = "Average Relative Abundance (%)",
          caption = "Numbers indicate the total candidate proteins for that model across all groups"
        ) +
        scale_y_continuous(
          limits = c(0, 110),
          breaks = seq(0, 100, 25),
          expand = c(0, 0)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        )
      
      ggsave(
        filename = file.path(plots_dir, paste0("6_Model_Localization_Comparison_", data_method, ".pdf")),
        plot = model_comp_plot,
        width = 12,
        height = 8,
        bg = "white"
      )
    }
    
    cat("Candidate protein localization distribution visualization completed.\n")
    
  } else {
    cat("Warning: detailed_candidates_by_model.xlsx is empty, skipping localization distribution analysis.\n")
  }
} else {
  cat("Warning: detailed_candidates_by_model.xlsx not found, skipping localization distribution analysis.\n")
}
} else {
  cat("\nCandidate protein localization distribution plot generation disabled, skipping this plot.\n")
}


cat("\n--- Visualization Phase Complete ---\n")

} else {
  cat("\nVisualization generation disabled, skipping plot creation.\n")
}


# --- New Section: Cross-Group (Inter-Method) Comparison ---
if (enable_group_comparison) {
  cat("\nStarting cross-proximity method comparison analysis (duplicate section)...\n")

# Define the baseline model for this comparison
baseline_model_for_comparison <- "2D_Abundance_logFC"

# --- FIX: Create a dedicated deep copy for this section to avoid side effects ---
all_results_for_comparison <- rlang::duplicate(all_results, shallow = FALSE)

# Group results by the original dataframe name (e.g., "ForStep19_APEX", "ForStep19_BioID2")
results_by_df <- all_results_for_comparison %>% 
  purrr::keep(~!is.null(.x)) %>% # Filter out failed runs
  split(., purrr::map_chr(., "df_name"))

# List to store the comparison results for each dataframe
all_group_comparisons <- list()

for (df_name in names(results_by_df)) {
  
  # Extract the ROC objects for the baseline model from each group
  roc_objects_for_comparison <- results_by_df[[df_name]] %>% 
    purrr::map(~.x$roc_objects[[baseline_model_for_comparison]]) %>% 
    # Name the list elements by their group name for clarity
    setNames(purrr::map_chr(results_by_df[[df_name]], "group_name")) %>% 
    # Remove any NULL entries if a model failed for a specific group
    purrr::compact()

  if (length(roc_objects_for_comparison) < 2) {
    warning(paste("Skipping group comparison for", df_name, "as it has fewer than 2 valid groups for the model", baseline_model_for_comparison))
    next
  }
  
  # --- Perform pairwise DeLong's test ---
  group_names <- names(roc_objects_for_comparison)
  p_value_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                           dimnames = list(group_names, group_names))
  
  group_pairs <- combn(group_names, 2, simplify = FALSE)
  
  p_values <- purrr::map_dbl(group_pairs, ~ {
    tryCatch({
      test_result <- roc.test(roc_objects_for_comparison[[.x[1]]], roc_objects_for_comparison[[.x[2]]], method = "delong")
      test_result$p.value
    }, error = function(e) { NA })
  })
  
  for (i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    p_value_matrix[pair[1], pair[2]] <- p_values[i]
    p_value_matrix[pair[2], pair[1]] <- p_values[i]
  }
  diag(p_value_matrix) <- 1.0

  # --- Consolidate results into a summary table ---
  # Get AUC and CIs for each group
  auc_ci_summary <- purrr::map_df(roc_objects_for_comparison, ~{
    ci_vals <- pROC::ci.auc(.x, method = "delong")
    tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
  }, .id = "GroupName")
  
  # Convert p-value matrix to a long format and apply BH correction
  p_values_long_temp <- as.data.frame(p_value_matrix) %>% 
    tibble::rownames_to_column("Group1") %>% 
    tidyr::pivot_longer(-Group1, names_to = "Group2", values_to = "p_value") %>% 
    filter(Group1 < Group2) # Avoid duplicates and self-comparisons
  
  p_values_long <- p_values_long_temp %>%
    mutate(
      q_value = if (nrow(p_values_long_temp) == 1) {
        # For single inter-group comparison, q-value equals p-value
        p_value
      } else {
        # For multiple inter-group comparisons, apply BH correction
        p.adjust(p_value, method = "BH")
      },
      single_group_comparison = (nrow(p_values_long_temp) == 1)
    )

  # Store results
  all_group_comparisons[[df_name]] <- list(
    auc_ci_summary = auc_ci_summary,
    pairwise_tests = p_values_long
  )
}

# --- Export Cross-Group Comparison Results ---
if (length(all_group_comparisons) > 0) {
  # Create a single Excel file with sheets for each type of result
  export_list_group_comp <- list()
  
  # Add AUC/CI summaries for each dataset
  auc_sheets <- purrr::map(all_group_comparisons, "auc_ci_summary")
  names(auc_sheets) <- paste("AUCs", names(auc_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, auc_sheets)
  
  # Add pairwise test results for each dataset
  pairwise_sheets <- purrr::map(all_group_comparisons, "pairwise_tests")
  names(pairwise_sheets) <- paste("P-Values", names(pairwise_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, pairwise_sheets)
  
  writexl::write_xlsx(export_list_group_comp, path = file.path(output_dir, "Group_Comparison_Summary.xlsx"))
  cat("\nCross-proximity method comparison results exported to Group_Comparison_Summary.xlsx\n")
} else {
  cat("\nCould not generate any cross-proximity method comparison results.\n")
}

} else {
  cat("\nInter-group comparison analysis disabled, skipping duplicate inter-group comparison analysis.\n")
}

# -------------------------------------------------------------------
# 7. Shutdown
# -------------------------------------------------------------------

# --- (REMOVED) Save Phase 1 Results to RData ---
# cat("\nSaving the key Phase 1 results to an RData file...\n")
# save(
#   all_results, 
#   all_group_comparisons, 
#   final_auc_summary,
#   file = file.path(output_dir, "analysis_phase1_output.RData")
# )
# cat("Phase 1 results saved to", file.path(output_dir, "analysis_phase1_output.RData"), "\n")

# Gracefully terminate the parallel backend
plan(sequential)
cat("\nParallel session terminated. Script execution completed.\n")


save.image(file = "Add_Phase1_2_alldata.RData")
print("Work environment successfully saved to Add_Phase1_2_alldata.RData")
