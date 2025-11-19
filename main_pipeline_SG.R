# ============================================================================
# Mass Spectrometry Analysis Pipeline - Main Control Script
# ============================================================================
# Instructions:
# 1. Run modules sequentially, check output after each module completes
# 2. Data is passed between modules via RData files
# 3. Ensure Reference/ and Rawdata/ directories exist before first run
# ============================================================================

# ============================================================================
# Global Configuration
# ============================================================================

# Working directory (modify according to actual situation)
setwd("D:/Bioinfomatics/MS/CodeNorm_SG_Batch2")

# Load required packages
required_packages <- c(
  "readr", "tidyverse", "readxl", "purrr", "ggplot2", 
  "ggsci", "ggthemes", "gridExtra", "pheatmap", "pROC", 
  "openxlsx", "preprocessCore", "limma", "scales", "ggrepel"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("✓ Packages loaded successfully\n")

# Directory configuration
dir_config <- list(
  root = getwd(),
  module = file.path(getwd(), "Module"),
  dev = file.path(getwd(), "Dev"),
  output = file.path(getwd(), "Output"),
  reference = file.path(getwd(), "Reference"),
  rawdata = file.path(getwd(), "Rawdata")
)

# Data reading configuration
data_config <- list(
  # TSV file name matching pattern (regular expression)
  file_pattern = "SG_batch2*\\.tsv$"
  # Can be modified according to actual file names, for example:
  # file_pattern = ".*\\.tsv$"  # All TSV files
  # file_pattern = "sg_.*\\.tsv$"  # TSV files starting with sg_
)

# ============================================================================
# Module 1: Environment Initialization
# ============================================================================
cat("\n========================================\n")
cat("Module 1: Environment Initialization\n")
cat("========================================\n")

source(file.path(dir_config$module, "module01_setup.R"))
config <- module01_setup(dir_config)

# Save Module01 environment (contains all data)
save.image(file = "Module01_workspace.RData")
cat("✓ Saved: Module01_workspace.RData (working directory)\n")
cat("  Contains: dir_config, config\n")

cat("✓ Module 1 completed\n")

# ============================================================================
# Module 2: Data Reading and Grouping Table
# ============================================================================
cat("\n========================================\n")
cat("Module 2: Data Reading and Grouping Table\n")
cat("========================================\n")

source(file.path(dir_config$module, "module02_data_import.R"))

# Step 1: Generate template
result <- module02_read_and_generate_template(dir_config, file_pattern = data_config$file_pattern)

# Step 2: Wait for user input
cat("\n⚠ Please complete the following steps before continuing:\n")
cat("  1. Open Module02_sampleGroup_template.csv\n")
cat("  2. Fill in all required columns (bioGroup is required)\n")
cat("  3. Optional: Modify Order column to adjust data column order\n")
cat("  4. Save Module02_sampleGroup_template.csv\n")
cat("  5. Press Enter to continue (program will automatically generate Module02_sampleGroup.csv)...\n")
readline()

# Step 3: Process user-filled data
result <- module02_process_samplegroup(dir_config, result$data)

# Assign returned data to global variables
data_raw <- result$data
sampleGroup <- result$sampleGroup

# Save Module02 environment (contains all Module01 data + Module02 data)
save.image(file = "Module02_workspace.RData")
cat("✓ Saved: Module02_workspace.RData (working directory)\n")
cat("  Contains: All Module01 data + data_raw + sampleGroup\n")

cat("✓ Module 2 completed\n")

# ============================================================================
# Module 3: Annotation System
# ============================================================================
cat("\n========================================\n")
cat("Module 3: Annotation System\n")
cat("========================================\n")

source(file.path(dir_config$module, "module03_annotation.R"))

# Select annotation mode:
#   - "SG": Use HaloMap/GO/MultiBait SGs presets
#   - "Nucleolus": Use HPA Nucleolus + CLL nucleolus reference
annotation_mode <- "SG"

# Advanced: If custom TP annotation columns are needed, provide a list here (will override annotation_mode when set)
custom_annotations_override <- NULL
# Example:
# custom_annotations_override <- list(
#   list(
#     column_name = "Custom_Localization",
#     TP_source = "GO_SGs",
#     TP_column = "Gene",
#     TP_label = "SGs"
#   )
# )

# Optional: Append additional annotation columns after selecting SG/Nucleolus (will not override existing configuration)
additional_annotations <- list()
# Example:
# additional_annotations <- list(
#   list(
#     column_name = "MyTurbo_Localization",
#     TP_source = "MyTurboSet",
#     TP_column = "Gene",
#     TP_label = "TurboTP"
#   )
# )

# Optional: Declare custom TP reference files (can be placed in Reference directory or provide absolute path)
custom_tp_sources <- list()
# Example:
# custom_tp_sources <- list(
#   list(
#     source_name = "MyNucleolusSet",
#     file = "MyNucleolus.csv",   # Default relative to Reference directory
#     file_type = "csv",          # Can be omitted, automatically determined by extension
#     gene_column = "Gene"        # Gene column name in the file
#   )
# )

result <- module03_annotation(
  dir_config,
  data_raw,
  sampleGroup,
  custom_annotations = custom_annotations_override,
  annotation_mode = annotation_mode,
  custom_tp_sources = custom_tp_sources,
  additional_annotations = additional_annotations
)

# If only basic annotations (Cytosol/Nuclear/Mitochondrion) are needed, keep custom_annotations_override as NULL,
# and select "SG"/"Nucleolus" mode within module03_annotation

# Assign returned data to global variables
data_annotated <- result$data_annotated
annotation_references <- result$annotation_references

# Save Module03 environment (contains all Module01-02 data + Module03 data)
save.image(file = "Module03_workspace.RData")
cat("✓ Saved: Module03_workspace.RData (working directory)\n")
cat("  Contains: All Module01-02 data + data_annotated + annotation_references\n")

cat("✓ Module 3 completed\n")

# ============================================================================
# Module 4: Standardization
# ============================================================================
cat("\n========================================\n")
cat("Module 4: Standardization\n")
cat("========================================\n")

source(file.path(dir_config$module, "module04_standardization.R"))

# Configure standardization types
# Optional values: "noNorm", "Global_QNorm", "Global_MNorm", "Local_QNorm", "Local_MNorm"
# Note: Must include at least one global standardization (Global_QNorm or Global_MNorm)
norm_types <- c("noNorm","Global_QNorm", "Local_QNorm")

result <- module04_standardization(dir_config, data_annotated, sampleGroup, norm_types)

# Assign returned data to global variables
standardized_data_list <- result$standardized_data_list
norm_types_used <- result$norm_types_used

# Save Module04 environment (contains all Module01-03 data + Module04 data)
save.image(file = "Module04_workspace.RData")
cat("✓ Saved: Module04_workspace.RData (working directory)\n")
cat("  Contains: All Module01-03 data + standardized_data_list + norm_types_used\n")

cat("✓ Module 4 completed\n")

# ============================================================================
# Module 5: Missing Value Imputation
# ============================================================================
cat("\n========================================\n")
cat("Module 5: Missing Value Imputation\n")
cat("========================================\n")

source(file.path(dir_config$module, "module05_imputation.R"))

# Configure imputation parameters
# impute_cat_mean: Whether to use mean imputation for Cat group with n_valid=2 (default FALSE)
# random_seed: Random seed for ensuring reproducibility of Perseus imputation
impute_cat_mean <- FALSE  # Can be changed to TRUE to enable Cat group imputation
random_seed <- 123        # Can be changed to other values to alter random results

result <- module05_imputation(dir_config, standardized_data_list, sampleGroup, 
                              impute_cat_mean, random_seed)

# Assign returned data to global variables
imputed_data_list <- result$imputed_data_list
imputation_params <- result$imputation_params

# Save Module05 environment (contains all Module01-04 data + Module05 data)
save.image(file = "Module05_workspace.RData")
cat("✓ Saved: Module05_workspace.RData (working directory)\n")
cat("  Contains: All Module01-04 data + imputed_data_list + imputation_params\n")

cat("✓ Module 5 completed\n")

# ============================================================================
# Module 6: Heatmap System
# ============================================================================
cat("\n========================================\n")
cat("Module 6: Heatmap System\n")
cat("========================================\n")

source(file.path(dir_config$module, "module06_heatmap.R"))

# Configure heatmap parameters
# selected_versions: Select standardized versions for plotting heatmaps (supports vector)
# heatmap_types: Types of heatmaps to plot
# correlation_config: Correlation heatmap configuration
# localization_columns: Annotation columns for classification (NULL means auto-detect)
# color_params: Heatmap color parameters (for AllLocalization and by_localization)

selected_versions <- c("noNorm_Imputed","Global_QNorm_Imputed", "Local_QNorm_Imputed")  # Supports multiple versions
#selected_versions <- c("Global_QNorm_Imputed", "Local_QNorm_Imputed")  # Supports multiple versions
## Optional values: "noNorm_Imputed", "Global_QNorm_Imputed", "Global_MNorm_Imputed", "Local_QNorm_Imputed", "Local_MNorm _Imputed"
heatmap_types <- c("all", "correlation")
#heatmap_types <- c("all", "correlation", "by_localization")
# Correlation heatmap configuration
correlation_config <- list(
  # Correlation range and colors
  corr_min = 0.8,
  corr_max = 1,
  corr_center = 0.9,
  col_low = "#4D97CD",
  col_center = "white",
  col_high = "#DB6968",
  
  # Method 1: Use exclusion rules (exclude unwanted samples)
  exclude_context = NULL,       # Example: c("Control") excludes all Control groups
  exclude_pltype = NULL,         # Example: c("PL") excludes all PL groups
  exclude_catalytic = c("NoCat"),      # Example: c("NoCat") excludes all NoCat groups
  
  # Method 2: Directly specify included bioGroups (takes priority over exclusion rules)
  include_biogroups = NULL       # Example: c("Light_Experiment", "H2O2_Experiment")
)

# Example of using exclusion rules:
# correlation_config$exclude_context <- c("Control")  # Keep only Experiment groups
# correlation_config$exclude_pltype <- c("PL")        # Exclude PL groups

# Example of directly specifying bioGroups:
# correlation_config$include_biogroups <- c("Light_Experiment", "H2O2_Experiment")

localization_columns <- NULL  # NULL means auto-detect all Localization columns

# Color parameters for AllLocalization and by_localization heatmaps
color_params <- list(
  custom_min = -4,
  custom_max = 4,
  custom_center = 0,
  col_low = "#4D97CD",
  col_center = "white",
  col_high = "#DB6968"
)

result <- module06_heatmap(dir_config, imputed_data_list, sampleGroup,
                          selected_versions, heatmap_types, 
                          correlation_config, localization_columns, color_params)

# Assign returned information to global variables
heatmap_info <- result

# Save Module06 environment (contains all Module01-05 data + Module06 data)
save.image(file = "Module06_workspace.RData")
cat("✓ Saved: Module06_workspace.RData (working directory)\n")
cat("  Contains: All Module01-05 data + heatmap_info\n")

cat("✓ Module 6 completed\n")

# ============================================================================
# Module 7: First Differential Analysis
# ============================================================================
cat("\n========================================\n")
cat("Module 7: First Differential Analysis\n")
cat("========================================\n")

source(file.path(dir_config$module, "module07_diff_analysis1.R"))

# Configure differential analysis parameters
# Module 7 automatically builds comparison groups based on sampleGroup$FirstROCgroup
# FirstROCgroup logic:
#   - Groups containing the same elements belong to the same group (e.g., A, A/B, A&B all belong to group A)
#   - Within the same group, bioGroups with Context=Experiment vs bioGroups with Context=Control
#
# selected_versions: Data versions to analyze (NULL means all versions)
#   - Optional: c("noNorm_Imputed", "Local_QNorm_Imputed", "Global_QNorm_Imputed", etc.)

# Select analysis versions (NULL means all versions)
selected_versions_diff <- c("Local_QNorm_Imputed")

result <- module07_diff_analysis1(dir_config, imputed_data_list, sampleGroup,
                                   selected_versions_diff)

# Assign returned information to global variables
diff_results1 <- result$diff_results1
comparisons_used <- result$comparisons_used
comparison_info <- result$comparison_info

# Save Module07 environment (contains all Module01-06 data + Module07 data)
save.image(file = "Module07_workspace.RData")
cat("✓ Saved: Module07_workspace.RData (working directory)\n")
cat("  Contains: All Module01-06 data + diff_results1, comparisons_used, comparison_info\n")

cat("✓ Module 7 completed\n")

# ============================================================================
# Module 8: First ROC Analysis
# ============================================================================
cat("\n========================================\n")
cat("Module 8: First ROC Analysis\n")
cat("========================================\n")

source(file.path(dir_config$module, "module08_roc_analysis1.R"))

# Configure ROC analysis parameters
# SubMito conversion: Replace Mitochondrion in annotations with MitoCarta3 sub-localizations (MIM, Matrix, MOM, IMS)
# ROC analysis: Use pROC package to calculate ROC curves and optimal thresholds based on logFC and annotation columns
#
# selected_versions: Data versions to analyze (NULL means all versions)
#   - Optional: c("noNorm_Imputed", "Local_QNorm_Imputed", "Global_QNorm_Imputed", etc.)
# roc_annotation_column: Annotation column for ROC analysis (default "GO_Localization")
# tp_label: True Positive label (default "SGs")
# fp_label: False Positive label (default "Matrix")
# enable_submito: Whether to enable SubMito conversion (default TRUE)
# submito_annotation_columns: Annotation columns for SubMito conversion (NULL means auto-detect all Localization columns)
# min_tp: Minimum TP threshold (default 0.3, used for filtering optimal thresholds)

# Select analysis versions (NULL means all versions)
selected_versions_roc <- c("Local_QNorm_Imputed")

# ROC analysis parameters
roc_annotation_column <- "GO_Localization"  # Annotation column for ROC analysis
tp_label <- "SGs"                            # True Positive label
fp_label <- "Matrix"                         # False Positive label (after SubMito conversion)
enable_submito <- TRUE                       # Enable SubMito conversion
submito_annotation_columns <- NULL           # Auto-detect all Localization columns
min_tp <- 0.3                                # Minimum TP threshold

result <- module08_roc_analysis1(dir_config, diff_results1, annotation_references,
                                  selected_versions_roc,
                                  roc_annotation_column,
                                  tp_label,
                                  fp_label,
                                  enable_submito,
                                  submito_annotation_columns,
                                  min_tp)

# Assign returned information to global variables
roc_results1 <- result$roc_results
roc_thresholds1 <- result$all_thresholds
data_with_submito <- result$data_with_submito
expr_fdr_df_list <- result$expr_fdr_df_list

# Save Module08 environment (contains all Module01-07 data + Module08 data)
save.image(file = "Module08_workspace.RData")
cat("✓ Saved: Module08_workspace.RData (working directory)\n")
cat("  Contains: All Module01-07 data + roc_results1, roc_thresholds1, data_with_submito, expr_fdr_df_list\n")

cat("✓ Module 8 completed\n")

# ============================================================================
# Module 9: Background Subtraction
# ============================================================================
cat("\n========================================\n")
cat("Module 9: Background Subtraction\n")
cat("========================================\n")

source(file.path(dir_config$module, "module09_background_subtraction.R"))

# Configure background subtraction parameters
# Automatically generates two types of results: Method A and Method B:
# - Method A: Do not use FDR filtering for specified comparisons
# - Method B: Use FDR filtering for all comparisons
# - bioGroups with Context="Experiment": Use threshold filtering
# - bioGroups with Context="Spatial": Do not use threshold filtering
# - bioGroups with Context="Control": No operation
#
# selected_versions: Data versions to analyze (NULL means all versions)
# fdr_threshold: FDR threshold (default 0.05)
# no_fdr_comparisons: Vector of comparisons that do not use FDR filtering in Method A (default NULL, i.e., process as Method B)
# use_roc_threshold: Whether to use ROC threshold (TRUE) or fixed FC threshold (FALSE)
# fixed_fc_threshold: Fixed FC threshold to use if not using ROC threshold
# min_valid_lfq: Minimum number of valid LFQ values (default 2, used for catalytic group filtering)
# annotation_column: Annotation column for grouping and plotting (default "GO_Localization")
# plot_all_annotations: Whether to show all annotations (TRUE) or only TP (FALSE)
# tp_label: TP label (default "SGs")
# tp_color: TP color (default "#DB6968")

# Select analysis versions (NULL means all versions)
selected_versions_bg <- c("Local_QNorm_Imputed")

# Background subtraction parameters
fdr_threshold <- 0.05
# Method A: Specify vector of comparisons that do not use FDR filtering (if NULL, Method A is equivalent to Method B)
no_fdr_comparisons <- c("K69A1B3_Light_vs_A1B3_Light","K69C3_Light_vs_C3_Light")  # Can be set to NULL to use default Method B
use_roc_threshold <- TRUE
fixed_fc_threshold <- NULL  # If use_roc_threshold=FALSE, set fixed threshold such as 2.0
min_valid_lfq <- 2
annotation_column <- "GO_Localization"
plot_all_annotations <- TRUE  # FALSE shows only TP, TRUE shows all annotations
tp_label <- "SGs"
tp_color <- "#DB6968"

result <- module09_background_subtraction(
  dir_config, sampleGroup, diff_results1, roc_thresholds1, comparison_info,
  data_with_submito,
  expr_fdr_df_list,
  selected_versions_bg, fdr_threshold, no_fdr_comparisons,
  use_roc_threshold, fixed_fc_threshold, min_valid_lfq,
  annotation_column, plot_all_annotations, tp_label, tp_color
)

# Assign returned information to global variables
filtered_data_A <- result$filtered_data_A
filtered_data_B <- result$filtered_data_B
merged_data_A <- result$merged_data_A
merged_data_B <- result$merged_data_B

# Save Module09 environment (contains all Module01-08 data + Module09 data)
save.image(file = "Module09_workspace.RData")
cat("✓ Saved: Module09_workspace.RData (working directory)\n")
cat("  Contains: All Module01-08 data + filtered_data_A, filtered_data_B, merged_data_A, merged_data_B\n")

cat("✓ Module 9 completed\n")

# ============================================================================
# Module 10: Data Replacement
# ============================================================================
cat("\n========================================\n")
cat("Module 10: Data Replacement\n")
cat("========================================\n")

source(file.path(dir_config$module, "module10_data_replacement.R"))

# Use Module 5's imputed_data_list, prefer Global_QNorm_Imputed (or Global_MNorm_Imputed if not available)
# Perform point replacement on Module 9's merged data (AfterROC_merged)
# Replacement logic: Keep all NAs, replace all non-NA values with Global standardization results
result <- module10_data_replacement(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  merged_data_A = merged_data_A,        # Use merged data (not filtered_data)
  merged_data_B = merged_data_B,        # Use merged data (not filtered_data)
  imputed_data_list = imputed_data_list,
  selected_versions = selected_versions_bg,
  prefer_version = "Global_QNorm_Imputed",
  fallback_version = "Global_MNorm_Imputed",
  test_n = 20
)

# Assign returned data to global variables
replaced_data_A <- result$replaced_data_A
replaced_data_B <- result$replaced_data_B

# Save Module10 environment (contains all Module01-09 data + Module10 data)
save.image(file = "Module10_workspace.RData")
cat("✓ Saved: Module10_workspace.RData (working directory)\n")
cat("  Contains: All Module01-09 data + replaced_data_A, replaced_data_B\n")

cat("✓ Module 10 completed\n")


load("Module10_workspace.RData")
# ============================================================================
# Module 11: Second Differential Analysis
# ============================================================================
cat("\n========================================\n")
cat("Module 11: Second Differential Analysis\n")
cat("========================================\n")

source(file.path(dir_config$module, "module11_diff_analysis2.R"))

# Perform second differential analysis on Module 10's data replacement results
# Only select Experiment and Spatial samples, build Exp vs Exp and Exp vs Spatial comparisons
result <- module11_diff_analysis2(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  replaced_data_A = replaced_data_A,
  replaced_data_B = replaced_data_B,
  selected_versions = selected_versions_bg
)

# Assign returned data to global variables
FDR_combined_df_list_2nd <- result$FDR_combined_df_list_2nd
sample_info_2nd <- result$sample_info
comparisons_list_2nd <- result$comparisons_list
bioGroups_selected_2nd <- result$bioGroups_selected
bioGroup_info_2nd <- result$bioGroup_info

# Save Module11 environment (contains all Module01-10 data + Module11 data)
save.image(file = "Module11_workspace.RData")
cat("✓ Saved: Module11_workspace.RData (working directory)\n")
cat("  Contains: All Module01-10 data + FDR_combined_df_list_2nd, etc.\n")

cat("✓ Module 11 completed\n")

# ============================================================================
# Module 12: Second ROC Analysis
# ============================================================================
cat("\n========================================\n")
cat("Module 12: Second ROC Analysis\n")
cat("========================================\n")

source(file.path(dir_config$module, "module12_second_roc.R"))

# Prepare expression matrix list (merge Method A/B and clean names)
expr_data_list_2nd <- list()
if (length(replaced_data_A) > 0) {
  expr_data_list_2nd <- c(
    expr_data_list_2nd,
    setNames(replaced_data_A, paste0(names(replaced_data_A), "_A"))
  )
}
if (length(replaced_data_B) > 0) {
  expr_data_list_2nd <- c(
    expr_data_list_2nd,
    setNames(replaced_data_B, paste0(names(replaced_data_B), "_B"))
  )
}

name_cleanup <- function(x) stringr::str_replace(x, "_Imputed", "_New")
names(expr_data_list_2nd) <- name_cleanup(names(expr_data_list_2nd))
names(FDR_combined_df_list_2nd) <- name_cleanup(names(FDR_combined_df_list_2nd))

names(FDR_combined_df_list_2nd) # View available versions

# Configure second ROC analysis interface (no SubMito)
# second_roc_versions: Specify versions to enter second ROC (names must match FDR list)
#   Example: c("noNorm_New_A","Local_QNorm_New_A","noNorm_New_B","Local_QNorm_New_B")
#   Set to NULL to automatically include all versions
second_roc_versions <- c(
  "Local_QNorm_New_A"
)
second_roc_annotation_column <- "GO_Localization"
second_roc_tp_label <- "SGs"
second_roc_fp_label <- "Cytosol"
second_roc_min_tp <- 0.3

second_roc_order <- second_roc_versions
if (is.null(second_roc_order)) {
  second_roc_order <- names(FDR_combined_df_list_2nd)
}
second_roc_order <- second_roc_order[second_roc_order %in% names(FDR_combined_df_list_2nd)]

result <- module12_second_roc(
  dir_config = dir_config,
  FDR_combined_df_list_2nd = FDR_combined_df_list_2nd,
  expr_data_list = expr_data_list_2nd,
  comparisons_list = if (exists("comparisons_list_2nd")) comparisons_list_2nd else NULL,
  annotation_column = second_roc_annotation_column,
  tp_label = second_roc_tp_label,
  fp_label = second_roc_fp_label,
  desired_order = if (length(second_roc_order) > 0) second_roc_order else names(FDR_combined_df_list_2nd),
  min_tp = second_roc_min_tp,
  sample_columns = if (exists("sample_info_2nd")) sample_info_2nd$SampleName else NULL,
  youden_ylim = c(0, 0.8)
)

FDR_combined_df_list_2nd <- result$FDR_combined_df_list_2nd
Expr_FDR_df_list_2nd <- result$Expr_FDR_df_list_2nd
all_desired_thresholds_2nd <- result$all_desired_thresholds_2nd

save.image(file = "Module12_workspace.RData")
cat("✓ Saved: Module12_workspace.RData (working directory)\n")
cat("  Contains: All Module01-11 data + second ROC results\n")

cat("✓ Module 12 completed\n")

# ============================================================================
# Module 13: Sub Annotation (SG / Nucleolus)
# ============================================================================
cat("\n========================================\n")
cat("Module 13: Sub Annotation\n")
cat("========================================\n")

source(file.path(dir_config$module, "module13_subsg_annotation.R"))

# Configure Sub annotation mode and column (subSG defaults to MultiBait, subNucleolus defaults to Sub_HPA)
subsg_mode <- "subSG"  # Optional: "subSG" or "subNucleolus"
subsg_target_column <- if (subsg_mode == "subNucleolus") "Sub_HPA_Localization" else "MultiBait_Localization"
subsg_levels_order <- NULL  # Use default order, can assign character vector for custom order
names(Expr_FDR_df_list_2nd)
subsg_forstep16_versions <- c(
  "Local_QNorm_New_A"
)

result <- module13_subsg_annotation(
  dir_config = dir_config,
  expr_fdr_df_list_2nd = Expr_FDR_df_list_2nd,
  annotation_references = annotation_references,
  mode = subsg_mode,
  target_column = subsg_target_column,
  levels_order = subsg_levels_order,
  forstep16_versions = subsg_forstep16_versions
)

Expr_FDR_df_list_2nd_SubSGs <- result$Expr_FDR_df_list_2nd_SubSGs
ForStep16 <- result$ForStep16
subsg_output_files <- result$output_files
names(Expr_FDR_df_list_2nd_SubSGs)
save.image(file = "Module13_workspace.RData")
cat("✓ Saved: Module13_workspace.RData (working directory)\n")
cat("  Contains: All Module01-12 data + Expr_FDR_df_list_2nd_SubSGs\n")

cat("✓ Module 13 completed\n")
load("Module13_workspace.RData")
# ============================================================================
# Module 14: Volcano Plots
# ============================================================================
cat("\n========================================\n")
cat("Module 14: Volcano Plots\n")
cat("========================================\n")

source(file.path(dir_config$module, "module14_volcano_plots.R"))

if (!exists("Expr_FDR_df_list_2nd_SubSGs")) {
  stop("✗ Error: Expr_FDR_df_list_2nd_SubSGs not detected, please complete Module 13 first")
}
if (length(Expr_FDR_df_list_2nd_SubSGs) == 0) {
  stop("✗ Error: Expr_FDR_df_list_2nd_SubSGs is empty, cannot plot volcano plots")
}
cat("✓ Detected Expr_FDR_df_list_2nd_SubSGs, can proceed to configure Module 14 volcano plot parameters\n")

module14_comparison_metadata <- NULL
if (exists("comparisons_list_2nd") && exists("bioGroup_info_2nd")) {
  module14_comparison_metadata <- module14_build_comparison_metadata(
    comparisons_list = comparisons_list_2nd,
    bioGroup_info = bioGroup_info_2nd
  )
}
if (is.null(module14_comparison_metadata)) {
  cat("⚠ Note: Failed to obtain comparison category information (comparisons_list_2nd or bioGroup_info_2nd missing)\n")
} else {
  cat("✓ Built comparison categories based on Module11 information (Exp_vs_Exp / Exp_vs_Spatial)\n")
}

# Configure volcano plot interface
# - volcano_versions: Versions to enter volcano plots (NULL = use all)
# - volcano_annotation_column: Select annotation column for coloring
# - volcano_annotation_color_map: Annotation -> color mapping, can be added/deleted as needed
# - volcano_use_threshold_colors: TRUE means color by threshold, ignoring annotation colors
# - volcano_logfc_threshold / volcano_fdr_threshold: Thresholds
# - volcano_label_mode: c("with","without") means output both with/without label PDFs
# - volcano_label_annotations: Which annotation categories are allowed to show labels
# - volcano_comparison_sets: Define logFC & FDR column combinations
volcano_versions <- c("Local_QNorm_New_A") # Optional: names(Expr_FDR_df_list_2nd_SubSGs), NULL=all
volcano_annotation_column <- "HaloMap_Localization" # Optional: *_Localization columns printed above
volcano_annotation_color_map <- c(
  "SGs" = "red",
  "Mitochondrion" = "black",
  "Nuclear" = "blue"
) # Unlisted annotations will automatically fall back to gray
volcano_use_threshold_colors <- FALSE # TRUE=color by threshold and override annotation colors
volcano_logfc_threshold <- 0.5       # Threshold logic: |logFC| ≥ this value
volcano_fdr_threshold <- 0.05        # Threshold logic: FDR ≤ this value
volcano_threshold_color_above <- "red"   # Threshold coloring: meets condition
volcano_threshold_color_below <- "grey70" # Threshold coloring: does not meet condition
volcano_label_mode <- c("without") # Allowed "with" or "without", can fill both
volcano_label_annotations <- c("SGs", "Mitochondrion") # Only these annotations are allowed to show labels
volcano_label_size <- 3
volcano_label_max_overlaps <- 15
# Optional: c("Exp_vs_Exp","Exp_vs_Spatial"), or select a single category, NULL=all
volcano_comparison_categories <- c("Exp_vs_Spatial")
volcano_xlim_override <- c(-3, 6)
volcano_ylim_override <- c(0, 15)
# Custom example (if need to replace default values):
# volcano_comparison_sets <- list(
#   list(
#     name = "Custom1",
#     fc_cols = c("K69A1B3_vs_K69C3_logFC", "K69A1B3_vs_C3_logFC"),
#     fdr_cols = c("K69A1B3_vs_K69C3_adj.P.Val", "K69A1B3_vs_C3_adj.P.Val"),
#     xlim = c(-2, 6),
#     ylim = c(0, 13)
#   ),
#   list(
#     name = "Custom2",
#     fc_cols = c("K69C3_vs_K20_logFC"),
#     fdr_cols = c("K69C3_vs_K20_adj.P.Val")
#   )
# )
volcano_comparison_sets <- module14_default_comparison_sets(module14_comparison_metadata)

module14_versions <- if (is.null(volcano_versions)) NULL else volcano_versions

volcano_result <- module14_volcano_plots(
  dir_config = dir_config,
  expr_fdr_df_list_2nd_subsgs = Expr_FDR_df_list_2nd_SubSGs,
  versions = module14_versions,
  annotation_column = volcano_annotation_column,
  annotation_color_map = volcano_annotation_color_map,
  use_threshold_colors = volcano_use_threshold_colors,
  logfc_threshold = volcano_logfc_threshold,
  fdr_threshold = volcano_fdr_threshold,
  threshold_color_above = volcano_threshold_color_above,
  threshold_color_below = volcano_threshold_color_below,
  label_mode = volcano_label_mode,
  label_annotations = volcano_label_annotations,
  label_size = volcano_label_size,
  label_max_overlaps = volcano_label_max_overlaps,
  comparison_sets = volcano_comparison_sets,
  xlim_override = volcano_xlim_override,
  ylim_override = volcano_ylim_override,
  comparison_metadata = module14_comparison_metadata,
  comparison_categories = volcano_comparison_categories
)

save.image(file = "Module14_workspace.RData")
cat("✓ Saved: Module14_workspace.RData (working directory)\n")
cat("  Contains: All Module01-13 data + volcano plot PDFs\n")

cat("✓ Module 14 completed\n")

# ============================================================================
cat("\n========================================\n")
cat("All modules completed!\n")
cat("========================================\n")
load("Module14_workspace.RData")
forstep16_r_object <- file.path(dir_config$output, "Module14_ForStep16.rds")
saveRDS(ForStep16, forstep16_r_object)
cat(sprintf("✓ Exported ForStep16 R object: %s\n", forstep16_r_object))

forstep16_excel <- file.path(dir_config$output, "Module14_ForStep16.xlsx")
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  warning("⚠ Note: Missing openxlsx package, cannot export ForStep16 Excel, please install openxlsx first")
} else {
  wb_forstep16 <- openxlsx::createWorkbook()
  existing_names <- character()
  for (sheet_name in names(ForStep16)) {
    df <- ForStep16[[sheet_name]]
    if (!is.data.frame(df)) {
      df <- as.data.frame(df)
    }
    sanitized <- gsub("[^A-Za-z0-9]+", "_", sheet_name)
    if (nchar(sanitized) == 0) sanitized <- "Sheet"
    base_name <- substr(sanitized, 1, 31)
    unique_name <- base_name
    suffix <- 1
    while (unique_name %in% existing_names) {
      suffix_label <- paste0("_", suffix)
      max_base_len <- max(31 - nchar(suffix_label), 1)
      unique_name <- paste0(substr(base_name, 1, max_base_len), suffix_label)
      suffix <- suffix + 1
    }
    existing_names <- c(existing_names, unique_name)
    openxlsx::addWorksheet(wb_forstep16, unique_name)
    openxlsx::writeData(wb_forstep16, unique_name, df)
  }
  openxlsx::saveWorkbook(wb_forstep16, forstep16_excel, overwrite = TRUE)
  cat(sprintf("✓ Exported ForStep16 multi-sheet Excel: %s\n", forstep16_excel))
}
# Next, you can run Base_model.R script to build 1D/2D filtering models and output final protein list

