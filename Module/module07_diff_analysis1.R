# ============================================================================
# Module 7: First Differential Analysis
# ============================================================================
# Features:
#   1. Automatically build comparison groups based on sampleGroup$FirstROCgroup
#   2. Perform differential expression with limma
#   3. Apply eBayes empirical Bayes adjustment
#   4. Extract logFC and adj.P.Val (FDR)
#   5. Output combined differential results and expression matrices
#
# FirstROCgroup logic:
#   - Groups with overlapping elements are considered the same (e.g., A, A/B, A&B all belong to A)
#   - Within a group: Context=Experiment bioGroup vs Context=Control bioGroup
#   - Examples:
#     * FirstROCgroup=A, Context=Experiment: K69A1B3_Light
#     * FirstROCgroup=A, Context=Control: K69A1B3_noLight, A1B3_Light
#     * FirstROCgroup=A&B, Context=Control: K69_Light (belongs to both A and B)
#     * Comparisons: K69A1B3_Light vs K69A1B3_noLight
#                   K69A1B3_Light vs K69_Light
#                   K69A1B3_Light vs A1B3_Light
#
# Input:
#   - dir_config: directory config (must include output)
#   - imputed_data_list: imputed data list from Module 05
#   - sampleGroup: sample grouping info (must include FinalName, bioGroup, FirstROCgroup, Context)
#   - selected_versions: data versions to analyze (NULL means all)
#
# Output:
#   - Module07_workspace.RData: Module01-06 data + diff_results1
#   - CSV files (Output directory):
#     - Module07_<version>_DiffAnalysis.csv - combined differential results
#     - Module07_<version>_ExprMatrix.csv - expression matrix (Gene × Samples)
#     - Module07_<version>_LogFC_Wide.csv - LogFC wide table (Gene × Comparisons)
#   - XLSX files (Output directory):
#     - Module07_<version>_DiffAnalysis_Full.xlsx - full topTable results (multi-sheet)
#
# Dependencies: limma, dplyr, tidyr, openxlsx

module07_diff_analysis1 <- function(dir_config,
                                     imputed_data_list,
                                     sampleGroup,
                                     selected_versions = NULL) {
  
  # Load required packages
  if (!require("limma", quietly = TRUE)) {
    stop("Need limma package: BiocManager::install('limma')")
  }
  if (!require("dplyr", quietly = TRUE)) {
    stop("Need dplyr package: install.packages('dplyr')")
  }
  if (!require("tidyr", quietly = TRUE)) {
    stop("Need tidyr package: install.packages('tidyr')")
  }
  if (!require("openxlsx", quietly = TRUE)) {
    stop("Need openxlsx package: install.packages('openxlsx')")
  }
  
  cat("\n========================================\n")
  cat("Module 7: First Differential Analysis\n")
  cat("========================================\n")
  
  # Validate input
  if (!("output" %in% names(dir_config))) {
    stop("dir_config must include 'output' path")
  }
  
  if (is.null(imputed_data_list) || length(imputed_data_list) == 0) {
    stop("imputed_data_list cannot be empty")
  }
  
  # Validate required sampleGroup columns
  required_cols <- c("FinalName", "bioGroup", "FirstROCgroup", "Context")
  missing_cols <- setdiff(required_cols, colnames(sampleGroup))
  if (length(missing_cols) > 0) {
    stop(paste("sampleGroup missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Use all versions if none specified
  if (is.null(selected_versions)) {
    selected_versions <- names(imputed_data_list)
  }
  
  # Validate selected versions exist
  missing_versions <- setdiff(selected_versions, names(imputed_data_list))
  if (length(missing_versions) > 0) {
    stop(paste("Versions not found:", paste(missing_versions, collapse = ", ")))
  }
  
  cat(sprintf("  - Versions to analyze: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("  - Sample count: %d\n", nrow(sampleGroup)))
  
  # ========================================
  # Step 1: Parse FirstROCgroup and build comparisons
  # ========================================
  cat("\nStep 1: Parse FirstROCgroup and build comparisons...\n")
  
  # Filter out samples with NA/empty FirstROCgroup
  sampleGroup_filtered <- sampleGroup %>%
    filter(!is.na(FirstROCgroup) & FirstROCgroup != "")
  
  cat(sprintf("  - Samples in first differential analysis: %d\n", nrow(sampleGroup_filtered)))
  
  # Parse FirstROCgroup to extract distinct elements
  # e.g., "A&B" or "A/B" -> c("A", "B")
  parse_roc_group <- function(roc_str) {
    if (is.na(roc_str) || roc_str == "") return(character(0))
    # Separators may be &, /, or others
    elements <- unlist(strsplit(roc_str, "[&/,;|]+"))
    # Trim spaces
    elements <- trimws(elements)
    # Drop empties
    elements <- elements[elements != ""]
    return(unique(elements))
  }
  
  # Extract FirstROCgroup elements per sample
  sampleGroup_filtered$ROCgroups <- lapply(sampleGroup_filtered$FirstROCgroup, parse_roc_group)
  
  # Collect all unique ROC group elements
  all_roc_elements <- unique(unlist(sampleGroup_filtered$ROCgroups))
  cat(sprintf("  - Detected ROC groups: %s\n", paste(all_roc_elements, collapse = ", ")))
  
  # Build comparisons for each ROC group
  comparisons <- list()
  comparison_info <- list()
  
  for (roc_element in all_roc_elements) {
    cat(sprintf("\n  Processing ROC group: %s\n", roc_element))
    
    # Find all bioGroup in this ROC group
    belongs_to_group <- sapply(sampleGroup_filtered$ROCgroups, function(x) roc_element %in% x)
    group_samples <- sampleGroup_filtered[belongs_to_group, ]
    
    cat(sprintf("    - Samples in group: %d\n", nrow(group_samples)))
    cat(sprintf("    - bioGroup in group: %s\n", paste(unique(group_samples$bioGroup), collapse = ", ")))
    
    # Experiment groups (Context == Experiment)
    exp_groups <- group_samples %>%
      filter(Context == "Experiment") %>%
      pull(bioGroup) %>%
      unique()
    
    # Control groups (Context == Control)
    ctrl_groups <- group_samples %>%
      filter(Context == "Control") %>%
      pull(bioGroup) %>%
      unique()
    
    cat(sprintf("    - Experiment groups: %s\n", 
                ifelse(length(exp_groups) > 0, paste(exp_groups, collapse = ", "), "None")))
    cat(sprintf("    - Control groups: %s\n", 
                ifelse(length(ctrl_groups) > 0, paste(ctrl_groups, collapse = ", "), "None")))
    
    # Build comparisons: each Experiment vs each Control
    if (length(exp_groups) > 0 && length(ctrl_groups) > 0) {
      for (exp_group in exp_groups) {
        for (ctrl_group in ctrl_groups) {
          comp <- c(exp_group, ctrl_group)
          comparisons[[length(comparisons) + 1]] <- comp
          
          # Save comparison info
          comparison_info[[length(comparison_info) + 1]] <- list(
            roc_group = roc_element,
            exp_group = exp_group,
            ctrl_group = ctrl_group,
            comparison_name = paste0(make.names(exp_group), "_vs_", make.names(ctrl_group))
          )
          
          cat(sprintf("    - Add comparison: %s vs %s\n", exp_group, ctrl_group))
        }
      }
    } else {
      cat(sprintf("    - Warning: ROC group %s lacks Experiment or Control; skipped\n", roc_element))
    }
  }
  
  if (length(comparisons) == 0) {
    stop("No comparisons constructed; check sampleGroup FirstROCgroup and Context columns")
  }
  
  cat(sprintf("\n✓ Built %d comparison groups\n", length(comparisons)))
  
  # Store differential analysis results for all versions
  diff_results1 <- list()
  
  # ========================================
  # Step 2: Run differential analysis per version
  # ========================================
  for (version in selected_versions) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Analyze version: %s\n", version))
    cat(sprintf("========================================\n"))
    
    # Get data
    data_imputed <- imputed_data_list[[version]]
    
    # Identify data/annotation columns
    data_cols <- intersect(sampleGroup$FinalName, colnames(data_imputed))
    anno_cols <- setdiff(colnames(data_imputed), data_cols)
    
    if (length(data_cols) == 0) {
      warning(sprintf("Version %s has no matching data columns; skipped", version))
      next
    }
    
    cat(sprintf("  - Data columns: %d\n", length(data_cols)))
    cat(sprintf("  - Annotation columns: %d\n", length(anno_cols)))
    
    # Extract expression matrix (data columns only)
    expr_matrix <- as.matrix(data_imputed[, data_cols, drop = FALSE])
    
    # Set rownames to Gene (if present)
    if ("Gene" %in% anno_cols) {
      rownames(expr_matrix) <- data_imputed$Gene
      gene_col <- data_imputed$Gene
    } else {
      warning("Gene column not found; using row numbers as names")
      rownames(expr_matrix) <- paste0("Protein_", 1:nrow(expr_matrix))
      gene_col <- rownames(expr_matrix)
    }
    
    # Output expression matrix
    expr_matrix_df <- data.frame(Gene = gene_col, expr_matrix, check.names = FALSE)
    expr_matrix_file <- file.path(dir_config$output, 
                                   paste0("Module07_", version, "_ExprMatrix.csv"))
    write.csv(expr_matrix_df, expr_matrix_file, row.names = FALSE)
    cat(sprintf("✓ Saved expression matrix: %s\n", basename(expr_matrix_file)))
    
    # Ensure sampleGroup order matches expr_matrix columns
    sampleGroup_ordered <- sampleGroup[match(data_cols, sampleGroup$FinalName), ]
    
    # Build design matrix
    cat("  - Building design matrix...\n")
    Group <- factor(sampleGroup_ordered$bioGroup)
    design <- model.matrix(~ 0 + Group)
    colnames(design) <- levels(Group)
    
    # Fit linear model
    cat("  - Fitting linear model...\n")
    fit <- limma::lmFit(expr_matrix, design)
    
    # Build contrast matrix
    cat("  - Building contrast matrix...\n")
    
    # bioGroup present in design matrix
    available_groups <- colnames(design)
    cat(sprintf("  - bioGroup in design matrix: %s\n", paste(available_groups, collapse = ", ")))
    
    # Filter out comparisons not in design
    valid_comparisons <- list()
    contrast_strings <- character()
    contrast_names <- character()
    
    for (i in seq_along(comparisons)) {
      comp <- comparisons[[i]]
      group1 <- comp[1]
      group2 <- comp[2]
      
      # Require both groups to be present in design
      group1_safe <- make.names(group1)
      group2_safe <- make.names(group2)
      
      if (!(group1_safe %in% available_groups)) {
        cat(sprintf("  - Warning: skip %s vs %s (%s not in data)\n", group1, group2, group1))
        next
      }
      if (!(group2_safe %in% available_groups)) {
        cat(sprintf("  - Warning: skip %s vs %s (%s not in data)\n", group1, group2, group2))
        next
      }
      
      # Create contrast name (GroupA_vs_GroupB)
      contrast_name <- paste0(group1_safe, "_vs_", group2_safe)
      contrast_names <- c(contrast_names, contrast_name)
      
      # Create contrast string (GroupA - GroupB)
      contrast_strings <- c(contrast_strings, paste0(group1_safe, " - ", group2_safe))
      
      # Save valid comparison
      valid_comparisons[[length(valid_comparisons) + 1]] <- comp
    }
    
    if (length(contrast_names) == 0) {
      warning(sprintf("Version %s has no valid comparisons; skipped", version))
      next
    }
    
    names(contrast_strings) <- contrast_names
    
    # Build contrast matrix
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_strings,
      levels = design
    )
    
    cat(sprintf("  - Valid contrasts: %d\n", length(contrast_names)))
    for (name in contrast_names) {
      cat(sprintf("    - %s\n", name))
    }
    
    # Compute contrasts and apply empirical Bayes
    cat("  - Empirical Bayes adjustment...\n")
    fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
    fit_ebayes <- limma::eBayes(fit_contrasts)
    
    # Extract results
    cat("  - Extract differential results...\n")
    
    # Actual coefficient names (from limma contrasts)
    actual_coef_names <- colnames(contrast_matrix)
    cat(sprintf("  - Contrast matrix columns: %s\n", paste(actual_coef_names, collapse = ", ")))
    
    # Store per-contrast results
    FDR_test_list <- list()
    Raw_FDR_test_list <- list()
    LogFC_list <- list()  # for building LogFC wide table
    successful_comparisons_list <- list()  # successful contrasts
    
    for (i in seq_along(contrast_names)) {
      contrast_name <- contrast_names[i]
      actual_coef_name <- actual_coef_names[i]  # use actual contrast column
      
      # Extract topTable using actual coefficient name
      tryCatch({
        top_table <- limma::topTable(fit_ebayes, 
                                     coef = actual_coef_name,
                                     n = Inf, 
                                     adjust.method = "BH")
        
        # Add Gene column
        top_table$Gene <- rownames(top_table)
        
        # Save full result (friendly name as key)
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
        
        # Store successful contrasts
        successful_comparisons_list[[contrast_name]] <- valid_comparisons[[i]]
        
        cat(sprintf("  ✓ Completed contrast: %s\n", contrast_name))
        
      }, error = function(e) {
        cat(sprintf("  - Error: contrast %s failed - %s\n", contrast_name, e$message))
      })
    }
    
    # Ensure at least one successful contrast
    if (length(FDR_test_list) == 0) {
      warning(sprintf("Version %s has no successful contrasts; skipped", version))
      next
    }
    
    # Merge all contrast results
    cat(sprintf("  - Merging contrast results (success %d/%d)...\n", 
                length(FDR_test_list), length(contrast_names)))
    FDR_combined_df <- Reduce(function(x, y) full_join(x, y, by = "Gene"), FDR_test_list)
    
    # Add annotation columns (if present)
    if (length(anno_cols) > 0) {
      anno_df <- data_imputed %>% select(all_of(anno_cols))
      FDR_combined_df <- FDR_combined_df %>%
        left_join(anno_df, by = "Gene")
    }
    
    # Save merged results to CSV
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
    
    # Save full results as XLSX (multi-sheet)
    xlsx_file <- file.path(dir_config$output, 
                           paste0("Module07_", version, "_DiffAnalysis_Full.xlsx"))
    wb <- openxlsx::createWorkbook()
    
    # Combined sheet
    openxlsx::addWorksheet(wb, "Combined")
    openxlsx::writeData(wb, "Combined", FDR_combined_df)
    
    # LogFC wide sheet
    openxlsx::addWorksheet(wb, "LogFC_Wide")
    openxlsx::writeData(wb, "LogFC_Wide", LogFC_wide)
    
    # Expression matrix sheet
    openxlsx::addWorksheet(wb, "ExprMatrix")
    openxlsx::writeData(wb, "ExprMatrix", expr_matrix_df)
    
    # Per-contrast full result sheets
    for (contrast_name in names(Raw_FDR_test_list)) {
      sheet_name <- substr(contrast_name, 1, 31)  # Excel sheet name limit 31 chars
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
      comparisons = successful_comparisons_list,  # store successful contrasts only
      comparison_info = comparison_info,
      contrast_names = names(Raw_FDR_test_list)  # store successful contrast names only
    )
  }
  
  cat("\n✓ Module 7 complete\n")
  cat(sprintf("  - Versions analyzed: %d\n", length(diff_results1)))
  cat(sprintf("  - Total comparison groups: %d\n", length(comparisons)))
  
  # Return result
  return(list(
    diff_results1 = diff_results1,
    comparisons_used = comparisons,  # all constructed contrasts (some may be invalid per version)
    comparison_info = comparison_info
  ))
}
