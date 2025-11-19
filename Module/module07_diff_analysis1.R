# ============================================================================
# Module 7: First Differential Analysis
# ============================================================================
# Functions:
#   1. Automatically construct comparison groups based on sampleGroup$FirstROCgroup
#   2. Perform differential expression analysis using limma package
#   3. Apply eBayes empirical Bayes moderation
#   4. Extract logFC and adj.P.Val (FDR correction)
#   5. Output merged differential analysis results and expression matrix
#
# FirstROCgroup Logic:
#   - Groups containing the same elements belong to the same group (e.g., A, A/B, A&B all belong to group A)
#   - Within the same group, bioGroup with Context=Experiment vs bioGroup with Context=Control
#   - Example:
#     * FirstROCgroup=A, Context=Experiment: K69A1B3_Light
#     * FirstROCgroup=A, Context=Control: K69A1B3_noLight, A1B3_Light
#     * FirstROCgroup=A&B, Context=Control: K69_Light (belongs to both group A and B)
#     * Comparisons: K69A1B3_Light vs K69A1B3_noLight
#                    K69A1B3_Light vs K69_Light
#                    K69A1B3_Light vs A1B3_Light
#
# Input:
#   - dir_config: Directory configuration (must contain output)
#   - imputed_data_list: Imputed data list output from Module05
#   - sampleGroup: Sample grouping information (must contain FinalName, bioGroup, FirstROCgroup, Context)
#   - selected_versions: Data versions to analyze (default NULL means all versions)
#
# Output:
#   - Module07_workspace.RData: Contains all data from Module01-06 + diff_results1
#   - CSV files (Output directory):
#     - Module07_<version>_DiffAnalysis.csv - Merged differential analysis results
#     - Module07_<version>_ExprMatrix.csv - Expression matrix (Gene × Samples)
#     - Module07_<version>_LogFC_Wide.csv - LogFC wide table (Gene × Comparisons)
#   - XLSX files (Output directory):
#     - Module07_<version>_DiffAnalysis_Full.xlsx - Full topTable results (multiple sheets)
#
# Dependencies: limma, dplyr, tidyr, openxlsx

module07_diff_analysis1 <- function(dir_config,
                                     imputed_data_list,
                                     sampleGroup,
                                     selected_versions = NULL) {
  
  # Load required packages
  if (!require("limma", quietly = TRUE)) {
    stop("limma package required: BiocManager::install('limma')")
  }
  if (!require("dplyr", quietly = TRUE)) {
    stop("dplyr package required: install.packages('dplyr')")
  }
  if (!require("tidyr", quietly = TRUE)) {
    stop("tidyr package required: install.packages('tidyr')")
  }
  if (!require("openxlsx", quietly = TRUE)) {
    stop("openxlsx package required: install.packages('openxlsx')")
  }
  
  cat("\n========================================\n")
  cat("Module 7: First Differential Analysis\n")
  cat("========================================\n")
  
  # Validate input
  if (!("output" %in% names(dir_config))) {
    stop("dir_config must contain 'output' path")
  }
  
  if (is.null(imputed_data_list) || length(imputed_data_list) == 0) {
    stop("imputed_data_list cannot be empty")
  }
  
  # Validate required columns in sampleGroup
  required_cols <- c("FinalName", "bioGroup", "FirstROCgroup", "Context")
  missing_cols <- setdiff(required_cols, colnames(sampleGroup))
  if (length(missing_cols) > 0) {
    stop(paste("sampleGroup missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # If no version specified, use all versions
  if (is.null(selected_versions)) {
    selected_versions <- names(imputed_data_list)
  }
  
  # Validate that selected versions exist
  missing_versions <- setdiff(selected_versions, names(imputed_data_list))
  if (length(missing_versions) > 0) {
    stop(paste("The following versions do not exist:", paste(missing_versions, collapse = ", ")))
  }
  
  cat(sprintf("  - Analysis versions: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("  - Number of samples: %d\n", nrow(sampleGroup)))
  
  # ========================================
  # Step 1: Parse FirstROCgroup and construct comparison groups
  # ========================================
  cat("\nStep 1: Parse FirstROCgroup and construct comparison groups...\n")
  
  # Filter out samples with NA or empty FirstROCgroup
  sampleGroup_filtered <- sampleGroup %>%
    filter(!is.na(FirstROCgroup) & FirstROCgroup != "")
  
  cat(sprintf("  - Number of samples participating in first differential analysis: %d\n", nrow(sampleGroup_filtered)))
  
  # Parse FirstROCgroup to extract all independent elements
  # Example: "A&B" or "A/B" -> c("A", "B")
  parse_roc_group <- function(roc_str) {
    if (is.na(roc_str) || roc_str == "") return(character(0))
    # Separators may be &, /, or others
    elements <- unlist(strsplit(roc_str, "[&/,;|]+"))
    # Remove whitespace
    elements <- trimws(elements)
    # Remove empty elements
    elements <- elements[elements != ""]
    return(unique(elements))
  }
  
  # Extract FirstROCgroup elements for each sample
  sampleGroup_filtered$ROCgroups <- lapply(sampleGroup_filtered$FirstROCgroup, parse_roc_group)
  
  # Get all unique ROC group elements
  all_roc_elements <- unique(unlist(sampleGroup_filtered$ROCgroups))
  cat(sprintf("  - Detected ROC groups: %s\n", paste(all_roc_elements, collapse = ", ")))
  
  # Construct comparisons for each ROC group
  comparisons <- list()
  comparison_info <- list()
  
  for (roc_element in all_roc_elements) {
    cat(sprintf("\n  Processing ROC group: %s\n", roc_element))
    
    # Find all bioGroups belonging to this ROC group
    belongs_to_group <- sapply(sampleGroup_filtered$ROCgroups, function(x) roc_element %in% x)
    group_samples <- sampleGroup_filtered[belongs_to_group, ]
    
    cat(sprintf("    - Number of samples in this group: %d\n", nrow(group_samples)))
    cat(sprintf("    - bioGroups in this group: %s\n", paste(unique(group_samples$bioGroup), collapse = ", ")))
    
    # Find Experiment groups (Context = Experiment)
    exp_groups <- group_samples %>%
      filter(Context == "Experiment") %>%
      pull(bioGroup) %>%
      unique()
    
    # Find Control groups (Context = Control)
    ctrl_groups <- group_samples %>%
      filter(Context == "Control") %>%
      pull(bioGroup) %>%
      unique()
    
    cat(sprintf("    - Experiment groups: %s\n", 
                ifelse(length(exp_groups) > 0, paste(exp_groups, collapse = ", "), "None")))
    cat(sprintf("    - Control groups: %s\n", 
                ifelse(length(ctrl_groups) > 0, paste(ctrl_groups, collapse = ", "), "None")))
    
    # Construct comparisons: each Experiment group vs each Control group
    if (length(exp_groups) > 0 && length(ctrl_groups) > 0) {
      for (exp_group in exp_groups) {
        for (ctrl_group in ctrl_groups) {
          comp <- c(exp_group, ctrl_group)
          comparisons[[length(comparisons) + 1]] <- comp
          
          # Save comparison information
          comparison_info[[length(comparison_info) + 1]] <- list(
            roc_group = roc_element,
            exp_group = exp_group,
            ctrl_group = ctrl_group,
            comparison_name = paste0(make.names(exp_group), "_vs_", make.names(ctrl_group))
          )
          
          cat(sprintf("    - Added comparison: %s vs %s\n", exp_group, ctrl_group))
        }
      }
    } else {
      cat(sprintf("    - Warning: ROC group %s lacks Experiment or Control groups, skipping\n", roc_element))
    }
  }
  
  if (length(comparisons) == 0) {
    stop("Failed to construct any comparison groups, please check FirstROCgroup and Context columns in sampleGroup")
  }
  
  cat(sprintf("\n✓ Total of %d comparison groups constructed\n", length(comparisons)))
  
  # Store differential analysis results for all versions
  diff_results1 <- list()
  
  # ========================================
  # Step 2: Perform differential analysis for each version
  # ========================================
  for (version in selected_versions) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Analyzing version: %s\n", version))
    cat(sprintf("========================================\n"))
    
    # Get data
    data_imputed <- imputed_data_list[[version]]
    
    # Identify data columns and annotation columns
    data_cols <- intersect(sampleGroup$FinalName, colnames(data_imputed))
    anno_cols <- setdiff(colnames(data_imputed), data_cols)
    
    if (length(data_cols) == 0) {
      warning(sprintf("Version %s has no matching data columns, skipping", version))
      next
    }
    
    cat(sprintf("  - Number of data columns: %d\n", length(data_cols)))
    cat(sprintf("  - Number of annotation columns: %d\n", length(anno_cols)))
    
    # Extract expression matrix (keep only data columns)
    expr_matrix <- as.matrix(data_imputed[, data_cols, drop = FALSE])
    
    # Set row names as Gene (assuming Gene column exists)
    if ("Gene" %in% anno_cols) {
      rownames(expr_matrix) <- data_imputed$Gene
      gene_col <- data_imputed$Gene
    } else {
      warning("No Gene column in data, using row numbers as row names")
      rownames(expr_matrix) <- paste0("Protein_", 1:nrow(expr_matrix))
      gene_col <- rownames(expr_matrix)
    }
    
    # Output expression matrix
    expr_matrix_df <- data.frame(Gene = gene_col, expr_matrix, check.names = FALSE)
    expr_matrix_file <- file.path(dir_config$output, 
                                   paste0("Module07_", version, "_ExprMatrix.csv"))
    write.csv(expr_matrix_df, expr_matrix_file, row.names = FALSE)
    cat(sprintf("✓ Saved expression matrix: %s\n", basename(expr_matrix_file)))
    
    # Ensure sampleGroup order matches expr_matrix column order
    sampleGroup_ordered <- sampleGroup[match(data_cols, sampleGroup$FinalName), ]
    
    # Construct design matrix
    cat("  - Constructing design matrix...\n")
    Group <- factor(sampleGroup_ordered$bioGroup)
    design <- model.matrix(~ 0 + Group)
    colnames(design) <- levels(Group)
    
    # Fit linear model
    cat("  - Fitting linear model...\n")
    fit <- limma::lmFit(expr_matrix, design)
    
    # Construct contrast matrix
    cat("  - Constructing contrast matrix...\n")
    
    # Get bioGroups actually present in design matrix
    available_groups <- colnames(design)
    cat(sprintf("  - bioGroups in design matrix: %s\n", paste(available_groups, collapse = ", ")))
    
    # Filter out comparison groups not in design matrix
    valid_comparisons <- list()
    contrast_strings <- character()
    contrast_names <- character()
    
    for (i in seq_along(comparisons)) {
      comp <- comparisons[[i]]
      group1 <- comp[1]
      group2 <- comp[2]
      
      # Check if both groups are in design matrix
      group1_safe <- make.names(group1)
      group2_safe <- make.names(group2)
      
      if (!(group1_safe %in% available_groups)) {
        cat(sprintf("  - Warning: Skipping comparison %s vs %s (%s not in data)\n", group1, group2, group1))
        next
      }
      if (!(group2_safe %in% available_groups)) {
        cat(sprintf("  - Warning: Skipping comparison %s vs %s (%s not in data)\n", group1, group2, group2))
        next
      }
      
      # Create comparison name (GroupA_vs_GroupB)
      contrast_name <- paste0(group1_safe, "_vs_", group2_safe)
      contrast_names <- c(contrast_names, contrast_name)
      
      # Create contrast string (GroupA - GroupB)
      contrast_strings <- c(contrast_strings, paste0(group1_safe, " - ", group2_safe))
      
      # Save valid comparison
      valid_comparisons[[length(valid_comparisons) + 1]] <- comp
    }
    
    if (length(contrast_names) == 0) {
      warning(sprintf("Version %s has no valid comparison groups, skipping", version))
      next
    }
    
    names(contrast_strings) <- contrast_names
    
    # Construct contrast matrix
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_strings,
      levels = design
    )
    
    cat(sprintf("  - Number of valid comparison groups: %d\n", length(contrast_names)))
    for (name in contrast_names) {
      cat(sprintf("    - %s\n", name))
    }
    
    # Calculate contrasts and apply empirical Bayes moderation
    cat("  - Applying empirical Bayes moderation...\n")
    fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
    fit_ebayes <- limma::eBayes(fit_contrasts)
    
    # Extract results
    cat("  - Extracting differential analysis results...\n")
    
    # Get actual column names of contrast matrix (names generated by limma)
    actual_coef_names <- colnames(contrast_matrix)
    cat(sprintf("  - Contrast matrix column names: %s\n", paste(actual_coef_names, collapse = ", ")))
    
    # Store results for each comparison
    FDR_test_list <- list()
    Raw_FDR_test_list <- list()
    LogFC_list <- list()  # For constructing LogFC wide table
    successful_comparisons_list <- list()  # Save successful comparisons
    
    for (i in seq_along(contrast_names)) {
      contrast_name <- contrast_names[i]
      actual_coef_name <- actual_coef_names[i]  # Use actual column name from contrast matrix
      
      # Extract topTable results (using actual coefficient name)
      tryCatch({
        top_table <- limma::topTable(fit_ebayes, 
                                     coef = actual_coef_name,  # Use actual coefficient name
                                     n = Inf, 
                                     adjust.method = "BH")
        
        # Add Gene column
        top_table$Gene <- rownames(top_table)
        
        # Save full results (using friendly name as key)
        Raw_FDR_test_list[[contrast_name]] <- top_table
        
        # Extract logFC and adj.P.Val
        result_df <- top_table %>%
          select(Gene, logFC, adj.P.Val) %>%
          rename(!!paste0(contrast_name, "_logFC") := logFC,
                 !!paste0(contrast_name, "_adj.P.Val") := adj.P.Val)
        
        FDR_test_list[[contrast_name]] <- result_df
        
        # Extract logFC for wide table
        logfc_df <- top_table %>%
          select(Gene, logFC) %>%
          rename(!!contrast_name := logFC)
        
        LogFC_list[[contrast_name]] <- logfc_df
        
        # Save successful comparison
        successful_comparisons_list[[contrast_name]] <- valid_comparisons[[i]]
        
        cat(sprintf("  ✓ Completed comparison: %s\n", contrast_name))
        
      }, error = function(e) {
        cat(sprintf("  - Error: Comparison %s failed - %s\n", contrast_name, e$message))
      })
    }
    
    # Check if there are any successful comparisons
    if (length(FDR_test_list) == 0) {
      warning(sprintf("Version %s has no successful comparisons, skipping", version))
      next
    }
    
    # Merge results from all comparisons
    cat(sprintf("  - Merging comparison results (%d/%d successful comparisons)...\n", 
                length(FDR_test_list), length(contrast_names)))
    FDR_combined_df <- Reduce(function(x, y) full_join(x, y, by = "Gene"), FDR_test_list)
    
    # Add annotation columns (if they exist)
    if (length(anno_cols) > 0) {
      anno_df <- data_imputed %>% select(all_of(anno_cols))
      FDR_combined_df <- FDR_combined_df %>%
        left_join(anno_df, by = "Gene")
    }
    
    # Save merged results as CSV
    csv_file <- file.path(dir_config$output, 
                          paste0("Module07_", version, "_DiffAnalysis.csv"))
    write.csv(FDR_combined_df, csv_file, row.names = FALSE)
    cat(sprintf("✓ Saved merged results: %s\n", basename(csv_file)))
    
    # Merge LogFC wide table
    if (length(LogFC_list) > 0) {
      LogFC_wide <- Reduce(function(x, y) full_join(x, y, by = "Gene"), LogFC_list)
    } else {
      LogFC_wide <- data.frame(Gene = gene_col)
    }
    logfc_wide_file <- file.path(dir_config$output, 
                                  paste0("Module07_", version, "_LogFC_Wide.csv"))
    write.csv(LogFC_wide, logfc_wide_file, row.names = FALSE)
    cat(sprintf("✓ Saved LogFC wide table: %s\n", basename(logfc_wide_file)))
    
    # Save full results as XLSX (multiple sheets)
    xlsx_file <- file.path(dir_config$output, 
                           paste0("Module07_", version, "_DiffAnalysis_Full.xlsx"))
    wb <- openxlsx::createWorkbook()
    
    # Add merged results sheet
    openxlsx::addWorksheet(wb, "Combined")
    openxlsx::writeData(wb, "Combined", FDR_combined_df)
    
    # Add LogFC wide table sheet
    openxlsx::addWorksheet(wb, "LogFC_Wide")
    openxlsx::writeData(wb, "LogFC_Wide", LogFC_wide)
    
    # Add expression matrix sheet
    openxlsx::addWorksheet(wb, "ExprMatrix")
    openxlsx::writeData(wb, "ExprMatrix", expr_matrix_df)
    
    # Add full results sheet for each comparison
    for (contrast_name in names(Raw_FDR_test_list)) {
      sheet_name <- substr(contrast_name, 1, 31)  # Excel sheet name limit is 31 characters
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, Raw_FDR_test_list[[contrast_name]])
    }
    
    openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
    cat(sprintf("✓ Saved full results: %s\n", basename(xlsx_file)))
    
    # Store results
    diff_results1[[version]] <- list(
      combined = FDR_combined_df,
      logfc_wide = LogFC_wide,
      expr_matrix = expr_matrix_df,
      raw_results = Raw_FDR_test_list,
      comparisons = successful_comparisons_list,  # Only save successful comparisons
      comparison_info = comparison_info,
      contrast_names = names(Raw_FDR_test_list)  # Only save successful comparison names
    )
  }
  
  cat("\n✓ Module 7 completed\n")
  cat(sprintf("  - Number of analyzed versions: %d\n", length(diff_results1)))
  cat(sprintf("  - Total comparison groups: %d\n", length(comparisons)))
  
  # Return results
  return(list(
    diff_results1 = diff_results1,
    comparisons_used = comparisons,  # Return all constructed comparisons (may not all be valid in some versions)
    comparison_info = comparison_info
  ))
}
