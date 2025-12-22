#' Module 09: Background Subtraction
#'
#' Filter data using Module 8 ROC thresholds and generate Method A and Method B results:
#' - Method A: skip FDR filter for specified comparisons
#' - Method B: apply FDR filter to all comparisons
#' - bioGroup with Context "Experiment": apply thresholds
#' - bioGroup with Context "Spatial": no threshold; select columns and drop NA
#' - bioGroup with Context "Control": no operation
#' - Produce stats per bioGroup (count, mean, sum, etc.)
#' - Plot stacked bar charts (count ratio and abundance ratio)
#'
#' @param dir_config Directory config list (includes output, etc.)
#' @param sampleGroup Sample grouping info (must include Context)
#' @param diff_results Differential results list (from Module 7)
#' @param roc_thresholds ROC threshold list (from Module 8)
#' @param comparison_info Comparison info (from Module 7; mapping comparison to bioGroup)
#' @param data_with_submito SubMito-transformed data list from Module 8 (optional)
#' @param selected_versions Versions to analyze (NULL means all)
#' @param fdr_threshold FDR threshold (default 0.05)
#' @param no_fdr_comparisons Comparisons without FDR in Method A (default NULL = same as Method B)
#' @param use_roc_threshold Whether to use ROC threshold (default TRUE)
#' @param fixed_fc_threshold Fixed FC threshold if not using ROC (default NULL)
#' @param min_valid_lfq Minimum valid LFQ count (default 2, for catalytic filtering)
#' @param annotation_column Annotation column for grouping/plotting (default "GO_Localization")
#' @param plot_all_annotations Show all annotations (TRUE) or TP only (FALSE)
#' @param tp_label TP label (default "SGs")
#' @param tp_color TP color (default "#DB6968")
#'
#' @return List with filtered_data_A, filtered_data_B, merged_data_A, merged_data_B

module09_background_subtraction <- function(dir_config,
                                             sampleGroup,
                                             diff_results,
                                             roc_thresholds,
                                             comparison_info,
                                             data_with_submito = NULL,
                                             expr_fdr_df_list = NULL,
                                             selected_versions = NULL,
                                             fdr_threshold = 0.05,
                                             no_fdr_comparisons = NULL,
                                             use_roc_threshold = TRUE,
                                             fixed_fc_threshold = NULL,
                                             min_valid_lfq = 2,
                                             annotation_column = "GO_Localization",
                                             plot_all_annotations = FALSE,
                                             tp_label = "SGs",
                                             tp_color = "#DB6968") {
  
  cat("\n=== Module 09: Background Subtraction ===\n")
  
  # Load required packages
  if (!require("dplyr")) stop("dplyr package required")
  if (!require("tidyr")) stop("tidyr package required")
  if (!require("purrr")) stop("purrr package required")
  if (!require("ggplot2")) stop("ggplot2 package required")
  if (!require("openxlsx")) stop("openxlsx package required")
  
  # Parameter validation
  if (!use_roc_threshold && is.null(fixed_fc_threshold)) {
    stop("If not using ROC threshold, fixed_fc_threshold must be provided")
  }
  
  if (!"Context" %in% colnames(sampleGroup)) {
    stop("sampleGroup must contain Context column")
  }
  
  # Determine versions to process
  if (is.null(selected_versions)) {
    selected_versions <- names(diff_results)
  }
  
  # Initialize result lists
  filtered_data_A <- list()  # Method A result
  filtered_data_B <- list()  # Method B result
  merged_data_A <- list()    # Method A merged data
  merged_data_B <- list()    # Method B merged data
  
  # Process each version
  for (version in selected_versions) {
    cat(sprintf("\nProcessing version: %s\n", version))
    cat(sprintf("  FDR threshold: %.2f\n", fdr_threshold))
    cat(sprintf("  Min valid LFQ: %d\n", min_valid_lfq))
    
    if (!version %in% names(diff_results)) {
      cat(sprintf("  - Warning: version %s not found in diff_results; skipped\n", version))
      next
    }
    
    if (!version %in% names(roc_thresholds)) {
      cat(sprintf("  - Warning: version %s not found in roc_thresholds; skipped\n", version))
      next
    }
    
    # Get data and thresholds
    thresholds <- roc_thresholds[[version]]
    
    # Data priority: Expr_FDR_df (Module 8) > SubMito-transformed (Module 8) > Module 7 combined
    if (!is.null(expr_fdr_df_list) && version %in% names(expr_fdr_df_list)) {
      cat("  - Using Module 8 Expr_FDR_df (expression + logFC/FDR + annotations)\n")
      combined_data <- expr_fdr_df_list[[version]]
    } else if (!is.null(data_with_submito) && version %in% names(data_with_submito)) {
      cat("  - Using Module 8 SubMito-transformed data\n")
      combined_data <- data_with_submito[[version]]
    } else {
      cat("  - Note: Module 8 SubMito data not found; using Module 7 combined results\n")
      combined_data <- diff_results[[version]]$combined
    }
    
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      cat("  - Warning: data empty; skipped\n")
      next
    }
    
    # If no Expr_FDR_df, try attaching LFQ columns from ExprMatrix
    if (is.null(expr_fdr_df_list) || !(version %in% names(expr_fdr_df_list))) {
      expr_matrix <- diff_results[[version]]$expr_matrix
      if (!is.null(expr_matrix) && "Gene" %in% colnames(expr_matrix)) {
        expr_lfq_cols <- intersect(sampleGroup$FinalName, colnames(expr_matrix))
        expr_df <- expr_matrix %>%
          select(Gene, all_of(expr_lfq_cols)) %>%
          distinct(Gene, .keep_all = TRUE)
        
        if (ncol(expr_df) > 1) {
          combined_data <- combined_data %>%
            left_join(expr_df, by = "Gene")
          cat(sprintf("  - Attached %d LFQ columns\n", ncol(expr_df) - 1))
        } else {
          cat("  - Warning: ExprMatrix has no matching LFQ columns\n")
        }
      } else {
        cat("  - Warning: ExprMatrix not found; cannot attach LFQ columns\n")
      }
    }
    
    # Identify all bioGroups and their Context
    biogroups <- unique(sampleGroup$bioGroup)
    cat(sprintf("  Found %d bioGroups\n", length(biogroups)))
    
    # Group bioGroups by Context (Experiment first, Spatial next)
    exp_groups <- c()
    spatial_groups <- c()
    control_groups <- c()
    
    for (bg in biogroups) {
      bg_context <- unique(sampleGroup$Context[sampleGroup$bioGroup == bg])
      if (length(bg_context) > 1) {
        cat(sprintf("  - Warning: bioGroup %s has multiple Context; using first: %s\n", bg, bg_context[1]))
        bg_context <- bg_context[1]
      }
      
      if (bg_context == "Experiment") {
        exp_groups <- c(exp_groups, bg)
      } else if (bg_context == "Spatial") {
        spatial_groups <- c(spatial_groups, bg)
      } else if (bg_context == "Control") {
        control_groups <- c(control_groups, bg)
      }
    }
    
    cat(sprintf("  Experiment groups: %d\n", length(exp_groups)))
    cat(sprintf("  Spatial groups: %d\n", length(spatial_groups)))
    cat(sprintf("  Control groups: %d\n", length(control_groups)))
    
    # ========== Method A: some comparisons without FDR ==========
    cat("\n  === Method A ===\n")
    if (is.null(no_fdr_comparisons) || length(no_fdr_comparisons) == 0) {
      cat("  - Note: no_fdr_comparisons not set; Method A equals Method B\n")
    } else {
      cat(sprintf("  - Comparisons without FDR filter: %s\n", paste(no_fdr_comparisons, collapse=", ")))
    }
    
    result_A <- process_one_method(
      combined_data, thresholds, sampleGroup, comparison_info,
      exp_groups, spatial_groups, control_groups,
      fdr_threshold, no_fdr_comparisons, use_roc_threshold, 
      fixed_fc_threshold, min_valid_lfq, annotation_column,
      method_name = "A"
    )
    
    # ========== Method B: all comparisons use FDR ==========
    cat("\n  === Method B ===\n")
    cat("  - All comparisons use FDR filter\n")
    
    result_B <- process_one_method(
      combined_data, thresholds, sampleGroup, comparison_info,
      exp_groups, spatial_groups, control_groups,
      fdr_threshold, NULL,  # Method B does not specify no_fdr_comparisons
      use_roc_threshold, fixed_fc_threshold, min_valid_lfq, annotation_column,
      method_name = "B"
    )
    
    # Save results
    filtered_data_A[[version]] <- result_A$filtered_list
    filtered_data_B[[version]] <- result_B$filtered_list
    merged_data_A[[version]] <- result_A$merged_data
    merged_data_B[[version]] <- result_B$merged_data
    
    # Output Excel files
    output_results(result_A$filtered_list, version, "A", dir_config$output)
    output_results(result_B$filtered_list, version, "B", dir_config$output)
    # Extra output: summarise and merged workbooks in source order
    output_summarise_workbook(result_A$filtered_list, version, "A", dir_config$output, exp_groups, spatial_groups)
    output_summarise_workbook(result_B$filtered_list, version, "B", dir_config$output, exp_groups, spatial_groups)
    output_merged_workbook(result_A$merged_data, version, "A", dir_config$output)
    output_merged_workbook(result_B$merged_data, version, "B", dir_config$output)
    
    # Plot stacked bar charts
    plot_stacked_bar(result_A$filtered_list, result_A$threshold_info,
                     version, "A", fdr_threshold, annotation_column,
                     plot_all_annotations, tp_label, tp_color, dir_config$output)
    
    plot_stacked_bar(result_B$filtered_list, result_B$threshold_info,
                     version, "B", fdr_threshold, annotation_column,
                     plot_all_annotations, tp_label, tp_color, dir_config$output)
  }
  
  cat("\n✓ Module 09 complete\n")
  
  return(list(
    filtered_data_A = filtered_data_A,
    filtered_data_B = filtered_data_B,
    merged_data_A = merged_data_A,
    merged_data_B = merged_data_B
  ))
}


#' Process a single method (A or B)
#'
#' @return List containing filtered_list, merged_data, threshold_info

process_one_method <- function(combined_data, thresholds, sampleGroup, comparison_info,
                                exp_groups, spatial_groups, control_groups,
                                fdr_threshold, no_fdr_comparisons,
                                use_roc_threshold, fixed_fc_threshold, min_valid_lfq,
                                annotation_column, method_name) {
  
  filtered_list <- list()  # store all filtered data and summaries
  data_list <- list()      # for merging filtered data only
  threshold_info <- list() # thresholds used per bioGroup
  
  # Process in order: Exp groups then Spatial groups
  all_groups <- c(exp_groups, spatial_groups)
  
  for (bg in all_groups) {
    cat(sprintf("    Processing bioGroup: %s\n", bg))
    
    # Get Context for this bioGroup
    bg_context <- unique(sampleGroup$Context[sampleGroup$bioGroup == bg])[1]
    
    # Get LFQ columns for this bioGroup
    lfq_cols <- sampleGroup$FinalName[sampleGroup$bioGroup == bg]
    lfq_cols <- lfq_cols[lfq_cols %in% colnames(combined_data)]
    
    if (length(lfq_cols) == 0) {
      cat(sprintf("      - Warning: LFQ column not found; skipped\n"))
      next
    }
    
    # Spatial groups: no threshold filtering
    if (bg_context == "Spatial") {
      cat(sprintf("      - Spatial groups: no threshold filter\n"))
      
      # Select columns: Gene, LFQ, all annotations
      annotation_cols <- grep("_Localization$", colnames(combined_data), value = TRUE)
      select_cols <- c("Gene", lfq_cols, annotation_cols)
      select_cols <- select_cols[select_cols %in% colnames(combined_data)]
      
      filtered_data <- combined_data %>%
        select(all_of(select_cols)) %>%
        na.omit()
      
      threshold_info[[bg]] <- "No threshold (Spatial)"
      
    } else if (bg_context == "Experiment") {
      # Experiment groups: apply thresholds
      cat(sprintf("      - Experiment groups: apply thresholds\n"))
      
      # Find comparisons involving this bioGroup
      bg_comparisons_idx <- vapply(comparison_info, function(comp) {
        comp$exp_group == bg
      }, logical(1))
      candidate_info <- comparison_info[bg_comparisons_idx]
      
      if (length(candidate_info) == 0) {
        cat(sprintf("      - Warning: comparison not found; skipped\n"))
        next
      }
      
      bg_comparisons <- vapply(candidate_info, function(comp) {
        if (!is.null(comp$comparison_name) && comp$comparison_name != "") {
          comp$comparison_name
        } else if (!is.null(names(comp))) {
          names(comp)[1]
        } else {
          NA_character_
        }
      }, character(1))
      bg_comparisons <- bg_comparisons[!is.na(bg_comparisons)]
      
      if (length(bg_comparisons) == 0) {
        cat(sprintf("      - Warning: comparison not found; skipped\n"))
        next
      }
      
      cat(sprintf("      - Found %d comparisons\n", length(bg_comparisons)))
      
      # Build filter conditions
      filtered_data <- combined_data
      threshold_text <- c()
      
      for (comp in bg_comparisons) {
        logfc_col <- paste0(comp, "_logFC")
        fdr_col <- paste0(comp, "_adj.P.Val")
        
        if (!logfc_col %in% colnames(filtered_data) || !fdr_col %in% colnames(filtered_data)) {
          cat(sprintf("        - Warning: column for comparison %s missing; skipped\n", comp))
          next
        }
        
        # Determine threshold to use
        if (use_roc_threshold) {
          # Use ROC threshold
          threshold <- thresholds[comp]
          if (is.na(threshold)) {
            cat(sprintf("        - Warning: ROC threshold for %s not found; skipped\n", comp))
            next
          }
        } else {
          # Use fixed FC threshold
          threshold <- fixed_fc_threshold
        }
        
        # Determine FDR threshold
        if (!is.null(no_fdr_comparisons) && comp %in% no_fdr_comparisons) {
          fdr_cutoff <- 1  # No FDR filter
          fdr_text <- "no FDR filter"
        } else {
          fdr_cutoff <- fdr_threshold
          fdr_text <- sprintf("FDR<%.2f", fdr_cutoff)
        }
        
        # Apply filters
        filtered_data <- filtered_data %>%
          filter(!!sym(logfc_col) > threshold,
                 !!sym(fdr_col) < fdr_cutoff)
        
        threshold_text <- c(threshold_text, 
                            sprintf("%s: logFC>%.2f, %s", comp, threshold, fdr_text))
      }
      
      # Select columns: Gene, LFQ, logFC/FDR, all annotations
      annotation_cols <- grep("_Localization$", colnames(combined_data), value = TRUE)
      logfc_cols <- paste0(bg_comparisons, "_logFC")
      fdr_cols <- paste0(bg_comparisons, "_adj.P.Val")
      select_cols <- c("Gene", lfq_cols, logfc_cols, fdr_cols, annotation_cols)
      select_cols <- select_cols[select_cols %in% colnames(filtered_data)]
      
      filtered_data <- filtered_data %>%
        select(all_of(select_cols))
      
      # Catalytic valid-value filter (LFQ columns with at least min_valid_lfq non-NA values)
      if (length(lfq_cols) >= min_valid_lfq) {
        filtered_data <- filtered_data %>%
          filter(rowSums(!is.na(select(., all_of(lfq_cols)))) >= min_valid_lfq)
      }
      
      threshold_info[[bg]] <- paste(threshold_text, collapse = "; ")
      
    } else {
      # Control groups: no action
      cat(sprintf("      - Control groups: skip\n"))
      next
    }
    
    cat(sprintf("      - After filtering retained %d proteins\n", nrow(filtered_data)))
    
    # Check minimum values (validate filtering)
    if (bg_context == "Experiment") {
      for (comp in bg_comparisons) {
        logfc_col <- paste0(comp, "_logFC")
        if (logfc_col %in% colnames(filtered_data)) {
          min_logfc <- min(filtered_data[[logfc_col]], na.rm = TRUE)
          cat(sprintf("      - %s min logFC: %.2f\n", comp, min_logfc))
        }
      }
    }
    
    # Save filtered data
    filtered_list[[bg]] <- filtered_data
    data_list[[bg]] <- filtered_data
    
    # Generate summary (grouped by annotation_column)
    if (annotation_column %in% colnames(filtered_data)) {
      summarise_data <- filtered_data %>%
        group_by(!!sym(annotation_column)) %>%
        summarise(
          count = n(),
          across(
            all_of(lfq_cols),
            list(
              mean = ~mean(.x, na.rm = TRUE),
              sum = ~sum(2^.x, na.rm = TRUE)  # back-transform to linear then sum
            ),
            .names = "{.col}_{.fn}"
          ),
          .groups = "drop"
        )
      
      # Compute percentage
      summarise_data <- summarise_data %>%
        mutate(Percent = count / sum(count))
      
      # Compute abundance percent (sum percent per LFQ column, then average)
      sum_cols <- grep("_sum$", colnames(summarise_data), value = TRUE)
      if (length(sum_cols) > 0) {
        for (col in sum_cols) {
          percent_col <- paste0(col, "_percent")
          summarise_data[[percent_col]] <- summarise_data[[col]] / sum(summarise_data[[col]], na.rm = TRUE)
        }
        
        # MeanSum (average of all sum_percent)
        percent_cols <- grep("_sum_percent$", colnames(summarise_data), value = TRUE)
        if (length(percent_cols) > 0) {
          summarise_data <- summarise_data %>%
            mutate(MeanSum = rowMeans(select(., all_of(percent_cols)), na.rm = TRUE))
        }
      }
      
      # Save summary
      summarise_name <- paste0(bg, "_Summarise")
      filtered_list[[summarise_name]] <- summarise_data
      cat(sprintf("      - Generated summary: %d groups\n", nrow(summarise_data)))
    }
  }
  
  # Merge all filtered data (AfterROC_merged logic)
  merged_data <- NULL
  if (length(data_list) > 0) {
    # Extract all annotation columns (from original combined_data)
    annotation_cols <- grep("_Localization$", colnames(combined_data), value = TRUE)
    anno_data <- combined_data %>% 
      select(Gene, all_of(annotation_cols)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    # full_join all filtered data
    merged_data <- reduce(data_list, full_join, by = "Gene")
    
    # Keep Gene and LFQ columns only
    lfq_cols_all <- c()
    for (bg in names(data_list)) {
      bg_lfq <- sampleGroup$FinalName[sampleGroup$bioGroup == bg]
      lfq_cols_all <- c(lfq_cols_all, bg_lfq)
    }
    lfq_cols_all <- lfq_cols_all[lfq_cols_all %in% colnames(merged_data)]
    
    merged_data <- merged_data %>%
      select(Gene, all_of(lfq_cols_all))
    
    # Add annotation columns
    merged_data <- left_join(merged_data, anno_data, by = "Gene")
    
    cat(sprintf("    - Merged data: %d proteins, %d columns\n", nrow(merged_data), ncol(merged_data)))
  }
  
  return(list(
    filtered_list = filtered_list,
    merged_data = merged_data,
    threshold_info = threshold_info
  ))
}


#' Output Excel files
#'
#' @param filtered_list Filtered data list (includes summaries)
#' @param version Data version
#' @param method Method name (A or B)
#' @param output_dir Output directory

output_results <- function(filtered_list, version, method, output_dir) {
  if (length(filtered_list) == 0) {
    cat(sprintf("    - Method %s: no data to output\n", method))
    return(invisible(NULL))
  }
  
  output_xlsx <- file.path(output_dir, 
                           sprintf("Module09_%s_%s_FilteredData.xlsx", method, version))
  wb <- createWorkbook()
  
  # Follow source order: Exp data + summaries then Spatial data + summaries
  for (name in names(filtered_list)) {
    sheet_name <- substr(name, 1, 31)  # Excel limit
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, filtered_list[[name]])
  }
  
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  cat(sprintf("    ✓ Method %s: Saved %s\n", method, basename(output_xlsx)))
}


#' Extra output: Summarise workbook ordered Exp→Spatial, data→summary
#'
#' @param filtered_list Filtered data and summaries
#' @param version Version
#' @param method Method (A/B)
#' @param output_dir Output directory
#' @param exp_groups Experiment bioGroup vector
#' @param spatial_groups Spatial bioGroup vector
output_summarise_workbook <- function(filtered_list, version, method, output_dir, exp_groups, spatial_groups) {
  if (length(filtered_list) == 0) {
    return(invisible(NULL))
  }
  ordered_names <- c()
  # 1) All "data" sheets: Exp → Spatial
  for (bg in exp_groups) {
    if (bg %in% names(filtered_list)) ordered_names <- c(ordered_names, bg)
  }
  for (bg in spatial_groups) {
    if (bg %in% names(filtered_list)) ordered_names <- c(ordered_names, bg)
  }
  # 2) Then all "Summarise" sheets: Exp → Spatial
  for (bg in exp_groups) {
    sum_name <- paste0(bg, "_Summarise")
    if (sum_name %in% names(filtered_list)) ordered_names <- c(ordered_names, sum_name)
  }
  for (bg in spatial_groups) {
    sum_name <- paste0(bg, "_Summarise")
    if (sum_name %in% names(filtered_list)) ordered_names <- c(ordered_names, sum_name)
  }
  # Fallback: append remaining if missed
  remaining <- setdiff(names(filtered_list), ordered_names)
  ordered_names <- c(ordered_names, remaining)
  
  output_xlsx <- file.path(output_dir,
                           sprintf("Module09_%s_%s_AfterROC_data_summarise.xlsx", method, version))
  wb <- createWorkbook()
  for (name in ordered_names) {
    sheet_name <- substr(name, 1, 31)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, filtered_list[[name]])
  }
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  cat(sprintf("    ✓ Method %s: Saved %s\n", method, basename(output_xlsx)))
}


#' Extra output: merged data workbook
#'
#' @param merged_data Merged data frame
#' @param version Version
#' @param method Method (A/B)
#' @param output_dir Output directory
output_merged_workbook <- function(merged_data, version, method, output_dir) {
  if (is.null(merged_data) || nrow(merged_data) == 0) {
    return(invisible(NULL))
  }
  output_xlsx <- file.path(output_dir,
                           sprintf("Module09_%s_%s_AfterROC_data_merged.xlsx", method, version))
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = substr(version, 1, 31))
  writeData(wb, sheet = substr(version, 1, 31), merged_data)
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  cat(sprintf("    ✓ Method %s: Saved %s\n", method, basename(output_xlsx)))
}


#' Plot stacked bar charts
#'
#' @param filtered_list Filtered data list (includes summaries)
#' @param threshold_info Threshold info list
#' @param version Data version
#' @param method Method name
#' @param fdr_threshold FDR threshold
#' @param annotation_column Annotation column name
#' @param plot_all_annotations Show all annotations
#' @param tp_label TP label
#' @param tp_color TP color
#' @param output_dir Output directory

plot_stacked_bar <- function(filtered_list, threshold_info, version, method,
                             fdr_threshold, annotation_column,
                             plot_all_annotations, tp_label, tp_color, output_dir) {
  
  # Extract summary data (names containing "_Summarise")
  summarise_names <- grep("_Summarise$", names(filtered_list), value = TRUE)
  
  if (length(summarise_names) == 0) {
    cat(sprintf("    - Method %s: no summary data for plotting\n", method))
    return(invisible(NULL))
  }
  
  # Prepare data for plotting
  plot_data_list <- list()
  
  for (name in summarise_names) {
    bg <- sub("_Summarise$", "", name)
    data <- filtered_list[[name]]
    
    if (!annotation_column %in% colnames(data)) {
      next
    }
    
    # Add bioGroup column
    data <- data %>%
      mutate(bioGroup = bg,
             bioGroup_factor = factor(bg, levels = sub("_Summarise$", "", summarise_names)))
    
    plot_data_list[[bg]] <- data
  }
  
  if (length(plot_data_list) == 0) {
    cat(sprintf("    - Method %s: no valid data to plot\n", method))
    return(invisible(NULL))
  }
  
  # Combine data
  merged_data <- bind_rows(plot_data_list)
  
  # Prepare threshold text (for plot subtitle)
  threshold_text <- sapply(names(threshold_info), function(bg) {
    paste0(bg, ": ", threshold_info[[bg]])
  })
  threshold_subtitle <- paste(threshold_text, collapse = "\n")
  
  # Extract TP data (for labeling)
  tp_data <- merged_data %>%
    filter(!!sym(annotation_column) == tp_label)
  
  # Set colors
  if (plot_all_annotations) {
    # Show all annotations: TP in red, others with palette
    all_annotations <- unique(merged_data[[annotation_column]])
    fill_levels <- c(tp_label, setdiff(all_annotations, tp_label))
    color_values <- setNames(
      c(tp_color, scales::hue_pal()(length(all_annotations) - 1)),
      fill_levels
    )
    
    # Ensure TP at bottom
    merged_data[[annotation_column]] <- factor(
      merged_data[[annotation_column]],
      levels = fill_levels
    )
  } else {
    # Show TP only; others grey
    all_annotations <- unique(merged_data[[annotation_column]])
    fill_levels <- c(tp_label, setdiff(all_annotations, tp_label))
    color_values <- setNames(
      c(tp_color, rep("grey", length(all_annotations) - 1)),
      fill_levels
    )
    
    # Ensure TP at bottom
    merged_data[[annotation_column]] <- factor(
      merged_data[[annotation_column]],
      levels = fill_levels
    )
  }
  
  # Explicit stack order: TP at bottom (smaller order lower)
  merged_data <- merged_data %>%
    mutate(.stack_order = ifelse(!!sym(annotation_column) == tp_label, 0L, 1L))
  
  # Plot count proportion
  output_file <- file.path(output_dir, 
                           sprintf("Module09_%s_%s_CountPercent_StackBar.pdf", 
                                   method, version))
  pdf(output_file, width = 10, height = 7)
  
  p <- ggplot(merged_data, aes(x = bioGroup_factor, y = Percent, fill = !!sym(annotation_column), order = .stack_order)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_text(data = tp_data,
              aes(x = bioGroup_factor, y = 1,
                  label = sprintf("%s:\n%d\n(%.1f%%)", 
                                  tp_label, count, Percent * 100)),
              inherit.aes = FALSE,
              vjust = -0.25,
              size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(title = sprintf("%s - Count Percent (Method %s, FDR<%.2f)", 
                         version, method, fdr_threshold),
         subtitle = threshold_subtitle,
         y = "Percent",
         x = "bioGroup") +
    scale_fill_manual(values = color_values, breaks = fill_levels, limits = fill_levels) +
    guides(fill = guide_legend(title = sprintf("Localization\n(%s)", annotation_column))) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0, size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  print(p)
  dev.off()
  cat(sprintf("    ✓ Method %s: Saved %s\n", method, basename(output_file)))
  
  # Plot abundance proportion (if MeanSum column present)
  if ("MeanSum" %in% colnames(merged_data)) {
    output_file <- file.path(output_dir, 
                             sprintf("Module09_%s_%s_AbundancePercent_StackBar.pdf", 
                                     method, version))
    pdf(output_file, width = 10, height = 7)
    
    tp_data_abundance <- merged_data %>%
      filter(!!sym(annotation_column) == tp_label)
    
    p <- ggplot(merged_data, aes(x = bioGroup_factor, y = MeanSum, fill = !!sym(annotation_column), order = .stack_order)) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(data = tp_data_abundance,
                aes(x = bioGroup_factor, y = 1,
                    label = sprintf("%s:\n%d\n(%.1f%%)", 
                                    tp_label, count, MeanSum * 100)),
                inherit.aes = FALSE,
                vjust = -0.25,
                size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(title = sprintf("%s - Abundance Percent (Method %s, FDR<%.2f)", 
                           version, method, fdr_threshold),
           subtitle = threshold_subtitle,
           y = "Percent Abundance",
           x = "bioGroup") +
      scale_fill_manual(values = color_values, breaks = fill_levels, limits = fill_levels) +
      guides(fill = guide_legend(title = sprintf("Localization\n(%s)", annotation_column))) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    print(p)
    dev.off()
    cat(sprintf("    ✓ Method %s: Saved %s\n", method, basename(output_file)))
  }
  
  return(invisible(NULL))
}
