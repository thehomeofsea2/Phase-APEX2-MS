# ============================================================================
# Module 8: First ROC Analysis
# ============================================================================
# Features:
# 1. SubMito localization conversion: replace Mitochondrion with MitoCarta3 sublocalizations (MIM, Matrix, MOM, IMS)
# 2. ROC analysis using the pROC package
# 3. Output ROC plots, Youden Index plots, ROC data tables, and optimal threshold tables
#
# Input:
# - Module07_workspace.RData (contains diff_results1, comparisons_used, etc.)
# - MitoCarta3.0 annotation data (for SubMito conversion)
#
# Output:
# - Module08_workspace.RData (working directory)
# - Output/Module08_{version}_ROC_curves.pdf
# - Output/Module08_{version}_Youden_Index.pdf
# - Output/Module08_{version}_ROC_data.xlsx
# - Output/Module08_{version}_thresholds.xlsx
# - Output/Module08_{version}_data_with_SubMito.csv
#
# ============================================================================

library(dplyr)
library(pROC)
library(openxlsx)
library(stringr)

#' Apply SubMito localization conversion
#'
#' @param data Data frame with annotation columns
#' @param mitocarta_anno MitoCarta annotation data
#' @param annotation_columns Annotation columns to convert
#' @param enable_submito Enable SubMito conversion (default TRUE)
#' @return Transformed data frame
apply_submito_transformation <- function(data, 
                                          mitocarta_anno, 
                                          annotation_columns = NULL,
                                          enable_submito = TRUE) {
  
  if (!enable_submito) {
    cat("  - SubMito conversion disabled; skipping\n")
    return(data)
  }
  
  # Auto-detect annotation columns if not provided
  if (is.null(annotation_columns)) {
    # Find columns containing 'Localization'
    annotation_columns <- grep("Localization", colnames(data), value = TRUE)
    if (length(annotation_columns) == 0) {
      cat("  - Warning: no annotation columns found; skipping SubMito conversion\n")
      return(data)
    }
  }
  
  cat(sprintf("  - Detected %d annotation columns to convert\n", length(annotation_columns)))
  
  # Prepare mitochondrial sublocalization lookup
  mito_lookup <- mitocarta_anno %>%
    select(Gene = Symbol, mito_subClass = MitoCarta3.0_SubMitoLocalization) %>%
    filter(mito_subClass %in% c("MIM", "Matrix", "MOM", "IMS"))
  
  cat(sprintf("  - MitoCarta3 sublocalization categories: %s\n", 
              paste(unique(mito_lookup$mito_subClass), collapse = ", ")))
  
  # Add sublocalization via left_join
  data_transformed <- data %>%
    left_join(mito_lookup, by = "Gene")
  
  # Convert each annotation column
  for (col in annotation_columns) {
    if (!(col %in% colnames(data_transformed))) {
      cat(sprintf("  - Warning: column %s not found; skipping\n", col))
      next
    }
    
    # Only replace rows with value 'Mitochondrion'; keep other annotations (Cytosol, Nuclear, SGs, Other, etc.)
    # Keeps annotation priority from previous modules
    data_transformed[[col]] <- ifelse(
      data_transformed[[col]] == "Mitochondrion" & !is.na(data_transformed$mito_subClass),
      data_transformed$mito_subClass,
      data_transformed[[col]]
    )
    
    cat(sprintf("  - Converted column: %s (only Mitochondrion -> SubMito)\n", col))
  }
  
  # Remove helper column
  data_transformed <- data_transformed %>%
    select(-mito_subClass)
  
  cat("  ✓ SubMito conversion complete\n")
  return(data_transformed)
}


#' Perform ROC analysis
#'
#' @param data Data frame with logFC and annotation columns
#' @param logfc_columns Vector of logFC column names
#' @param annotation_column Annotation column used for ROC
#' @param tp_label True Positive label (e.g., 'SGs')
#' @param fp_label False Positive label (e.g., 'Matrix')
#' @param direction ROC direction (default '<')
#' @return List with ROC objects and data frames
perform_roc_analysis <- function(data,
                                  logfc_columns,
                                  annotation_column,
                                  tp_label,
                                  fp_label,
                                  direction = "<") {
  
  cat(sprintf("  - Run ROC analysis: %s vs %s\n", tp_label, fp_label))
  cat(sprintf("  - Annotation column: %s\n", annotation_column))
  cat(sprintf("  - LogFC columns: %d\n", length(logfc_columns)))
  
  # Filter data: keep TP and FP only
  data_for_roc <- data %>%
    filter(!!sym(annotation_column) %in% c(tp_label, fp_label))
  
  n_tp <- sum(data_for_roc[[annotation_column]] == tp_label)
  n_fp <- sum(data_for_roc[[annotation_column]] == fp_label)
  
  cat(sprintf("  - %s count: %d\n", tp_label, n_tp))
  cat(sprintf("  - %s count: %d\n", fp_label, n_fp))
  
  if (n_tp < 5 || n_fp < 5) {
    cat("  - Warning: insufficient samples (TP or FP < 5); skip ROC analysis\n")
    return(list(roc_objects = list(), roc_dfs = list()))
  }
  
  # Run ROC analysis for each logFC column
  roc_objects <- list()
  roc_dfs <- list()
  
  for (i in seq_along(logfc_columns)) {
    logfc_col <- logfc_columns[i]
    
    if (!(logfc_col %in% colnames(data_for_roc))) {
      cat(sprintf("  - Warning: column %s not found; skipping\n", logfc_col))
      next
    }
    
    # Remove NA values
    data_clean <- data_for_roc %>%
      filter(!is.na(!!sym(logfc_col)))
    
    if (nrow(data_clean) < 10) {
      cat(sprintf("  - Warning: %s lacks sufficient data; skipping\n", logfc_col))
      next
    }
    
    tryCatch({
      # Run ROC analysis
      roc_obj <- roc(
        response = data_clean[[annotation_column]],
        predictor = data_clean[[logfc_col]],
        levels = c(fp_label, tp_label),
        direction = direction,
        quiet = TRUE
      )
      
      roc_name <- paste0("roc", i)
      roc_objects[[roc_name]] <- roc_obj
      
      # Extract ROC data
      roc_data <- data.frame(
        Threshold = roc_obj$thresholds,
        TP = roc_obj$sensitivities,
        FP = 1 - roc_obj$specificities,
        TP_FP = roc_obj$sensitivities - (1 - roc_obj$specificities)
      )
      
      # Rename columns using logFC column name
      colnames(roc_data) <- c(
        paste0(logfc_col, "_Threshold"),
        paste0(logfc_col, "_TP"),
        paste0(logfc_col, "_FP"),
        paste0(logfc_col, "_TP_FP")
      )
      
      roc_dfs[[roc_name]] <- roc_data
      
      cat(sprintf("  ✓ Completed: %s (AUC = %.3f)\n", logfc_col, auc(roc_obj)))
      
    }, error = function(e) {
      cat(sprintf("  - Error: %s - %s\n", logfc_col, e$message))
    })
  }
  
  return(list(roc_objects = roc_objects, roc_dfs = roc_dfs))
}


#' Compute optimal thresholds
#'
#' @param roc_dfs List of ROC data frames
#' @param min_tp Minimum TP threshold (default 0.3)
#' @return Threshold vector
calculate_optimal_thresholds <- function(roc_dfs, min_tp = 0.3) {
  
  threshold_vec <- numeric(length(roc_dfs))
  comparison_names <- character(length(roc_dfs))
  
  for (i in seq_along(roc_dfs)) {
    df <- roc_dfs[[i]]
    
    # Fetch column names dynamically
    tp_fp_col <- grep("TP_FP", colnames(df), value = TRUE)[1]
    thres_col <- grep("Threshold", colnames(df), value = TRUE)[1]
    tp_col <- grep("_TP$", colnames(df), value = TRUE)[1]
    
    # Extract comparison name (remove suffix)
    comparison_names[i] <- sub("_logFC_TP_FP$", "", tp_fp_col)
    
    # Find row with max TP_FP
    max_tp_fp <- max(df[[tp_fp_col]], na.rm = TRUE)
    candidate_rows <- df %>%
      filter(!!sym(tp_fp_col) == max_tp_fp)
    
    # Check if TP in those rows are all < min_tp
    if (all(candidate_rows[[tp_col]] < min_tp)) {
      threshold_vec[i] <- 0  # Does not meet requirement; set to 0
    } else {
      # For TP >= min_tp, choose smallest threshold
      valid_candidates <- candidate_rows %>%
        filter(!!sym(tp_col) >= min_tp)
      threshold_vec[i] <- min(valid_candidates[[thres_col]], na.rm = TRUE)
    }
  }
  
  names(threshold_vec) <- comparison_names
  return(threshold_vec)
}


#' Plot ROC curves
#'
#' @param roc_objects List of ROC objects
#' @param roc_dfs List of ROC data frames (for column names)
#' @param output_file Output PDF file path
#' @param main_title Main title
plot_roc_curves <- function(roc_objects, roc_dfs, output_file, main_title = "") {
  
  n_plots <- length(roc_objects)
  if (n_plots == 0) {
    cat("  - Warning: no ROC curves to plot\n")
    return(invisible(NULL))
  }
  
  # Compute layout
  n_cols <- min(4, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  pdf(output_file, width = 12, height = 6 * n_rows / 2)
  par(mfrow = c(n_rows, n_cols))
  
  for (i in seq_along(roc_objects)) {
    roc_obj <- roc_objects[[i]]
    
    # Extract comparison name from column (same as Youden)
    roc_name <- names(roc_objects)[i]  # e.g., "roc1"
    df <- roc_dfs[[roc_name]]
    
    # Get TP_FP column and strip suffix for comparison name
    tp_fp_col <- grep("_TP_FP$", colnames(df), value = TRUE)[1]
    # Remove '_logFC_TP_FP' suffix for comparison name
    predictor_name <- sub("_logFC_TP_FP$", "", tp_fp_col)
    
    plot(roc_obj,
         xlab = "FPR (1 - Specificity)",
         ylab = "TPR (Sensitivity)",
         print.auc = TRUE,
         auc.polygon = TRUE,
         max.auc.polygon = TRUE,
         max.auc.polygon.col = "white",
         auc.polygon.col = "#DB6B691A",
         grid = c(0.2, 0.2),
         grid.col = c("black", "black"),
         print.thres = TRUE,
         print.thres.col = "black",
         print.thres.cex = 1,
         print.auc.cex = 1,
         print.auc.col = "#DB6968",
         cex.lab = 1.5,
         cex.axis = 1,
         col = "#DB6968",
         legacy.axes = TRUE
    )
    
    # Add title separately (pROC plot lacks 'main')
    title(main = predictor_name, cex.main = 1.5)
  }
  
  dev.off()
  cat(sprintf("  ✓ Saved: %s\n", output_file))
}


#' Plot Youden Index
#'
#' @param roc_dfs List of ROC data frames
#' @param output_file Output PDF file path
plot_youden_index <- function(roc_dfs, output_file) {
  
  n_plots <- length(roc_dfs)
  if (n_plots == 0) {
    cat("  - Warning: no Youden Index to plot\n")
    return(invisible(NULL))
  }
  
  # Compute layout
  n_cols <- min(4, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  pdf(output_file, width = 12, height = 6 * n_rows / 2)
  par(mfrow = c(n_rows, n_cols))
  
  for (i in seq_along(roc_dfs)) {
    df <- roc_dfs[[i]]
    
    # Fetch column names dynamically
    tp_fp_col <- grep("_TP_FP$", colnames(df), value = TRUE)
    thres_col <- grep("_Threshold$", colnames(df), value = TRUE)
    
    # Extract clean title (remove '_logFC_TP_FP' suffix for comparison name)
    main_title <- sub("_logFC_TP_FP$", "", tp_fp_col)
    
    # Find optimal point
    optimal_point <- df[which.max(df[[tp_fp_col]]), ]
    
    plot(
      x = df[[thres_col]],
      y = df[[tp_fp_col]],
      type = "l",
      col = "#DB6968",
      lwd = 2,
      main = main_title,
      xlab = "Log2 FC",
      ylab = "TPR - FPR (Youden's Index)",
      cex.main = 1.5,
      cex.lab = 1.5,
      cex.axis = 1,
      ylim = c(0, max(0.7, max(df[[tp_fp_col]], na.rm = TRUE) * 1.1))
    )
    
    # Add vertical line and optimal point
    abline(v = optimal_point[[thres_col]], col = "#8B96AD", lty = 2)
    points(
      x = optimal_point[[thres_col]],
      y = optimal_point[[tp_fp_col]],
      col = "black",
      pch = 18,
      cex = 2
    )
    
    # Add annotation
    text(
      x = optimal_point[[thres_col]],
      y = optimal_point[[tp_fp_col]],
      labels = sprintf("(%.2f, %.2f)", 
                       optimal_point[[thres_col]], 
                       optimal_point[[tp_fp_col]]),
      pos = 4,
      col = "black",
      cex = 1.2
    )
  }
  
  dev.off()
  cat(sprintf("  ✓ Saved: %s\n", output_file))
}


#' Module 8: First ROC Analysis
#'
#' @param dir_config Directory config
#' @param diff_results1 Differential results list (from Module 7)
#' @param annotation_references Annotation references (from Module 3)
#' @param selected_versions Versions to analyze (NULL for all)
#' @param roc_annotation_column Annotation column for ROC (default "GO_Localization")
#' @param tp_label True Positive label (default 'SGs')
#' @param fp_label False Positive label (default 'Matrix')
#' @param enable_submito Enable SubMito conversion (default TRUE)
#' @param submito_annotation_columns Annotation columns for SubMito conversion (NULL auto-detect)
#' @param min_tp Minimum TP threshold (default 0.3)
#' @return List containing ROC analysis results
module08_roc_analysis1 <- function(dir_config,
                                    diff_results1,
                                    annotation_references,
                                    selected_versions = NULL,
                                    roc_annotation_column = "GO_Localization",
                                    tp_label = "SGs",
                                    fp_label = "Matrix",
                                    enable_submito = TRUE,
                                    submito_annotation_columns = NULL,
                                    min_tp = 0.3) {
  
  cat("\n========================================\n")
  cat("Module 8: First ROC Analysis\n")
  cat("========================================\n")
  
  # Validate input
  if (!all(c("reference", "output") %in% names(dir_config))) {
    stop("Error: dir_config must include 'reference' and 'output' paths")
  }
  
  if (length(diff_results1) == 0) {
    stop("Error: diff_results1 is empty; run Module 7 first")
  }
  
  # Select versions to analyze
  if (is.null(selected_versions)) {
    selected_versions <- names(diff_results1)
  } else {
    # Validate version names
    missing_versions <- setdiff(selected_versions, names(diff_results1))
    if (length(missing_versions) > 0) {
      stop(sprintf("Error: version not found: %s", paste(missing_versions, collapse = ", ")))
    }
  }
  
  cat(sprintf("Versions to analyze: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("ROCAnnotation column: %s\n", roc_annotation_column))
  cat(sprintf("TP label: %s\n", tp_label))
  cat(sprintf("FP label: %s\n", fp_label))
  cat(sprintf("SubMito conversion: %s\n", ifelse(enable_submito, "enabled", "disabled")))
  
  # Fetch MitoCarta annotations
  mitocarta_anno <- annotation_references$MitoCarta
  if (is.null(mitocarta_anno) && enable_submito) {
    cat("  - Warning: MitoCarta annotation not found; SubMito conversion disabled\n")
    enable_submito <- FALSE
  }
  
  # Store results
  roc_results <- list()
  all_thresholds <- list()
  data_with_submito_list <- list()
  expr_fdr_df_list <- list()
  
  # Run ROC analysis for each version
  for (version in selected_versions) {
    cat(sprintf("\nProcessing version: %s\n", version))
    cat("----------------------------------------\n")
    
    # Get differential data (Gene + logFC + adj.P.Val + annotations)
    diff_data <- diff_results1[[version]]$combined
    
    if (is.null(diff_data) || nrow(diff_data) == 0) {
      cat("  - Warning: data empty; skipping\n")
      next
    }
    
    # Apply SubMito conversion
    cat("Step 1: SubMito localization conversion\n")
    data_transformed <- apply_submito_transformation(
      data = diff_data,
      mitocarta_anno = mitocarta_anno,
      annotation_columns = submito_annotation_columns,
      enable_submito = enable_submito
    )
    
    # Save transformed data
    data_with_submito_list[[version]] <- data_transformed
    output_csv <- file.path(dir_config$output, 
                            paste0("Module08_", version, "_data_with_SubMito.csv"))
    write.csv(data_transformed, output_csv, row.names = FALSE)
    cat(sprintf("  ✓ Saved: %s\n", basename(output_csv)))
    
    # Generate Expr_FDR_df (expression + logFC/FDR + annotations) for later background subtraction
    expr_matrix <- diff_results1[[version]]$expr_matrix
    if (!is.null(expr_matrix) && "Gene" %in% colnames(expr_matrix)) {
      expr_df <- expr_matrix %>%
        select(Gene, everything()) %>%
        distinct(Gene, .keep_all = TRUE)
      # Build right table (original annotations + FC/FDR, no SubMito)
      anno_cols <- grep("_Localization$", colnames(diff_data), value = TRUE)
      fdr_cols <- grep("(_logFC$|_adj\\.P\\.Val$)", colnames(diff_data), value = TRUE)
      right_core <- diff_data %>% select(Gene, all_of(anno_cols), all_of(fdr_cols))
      # Final table: Gene + expression + original annotations + FC/FDR (no SubMito)
      expr_fdr_df <- expr_df %>% left_join(right_core, by = "Gene")
      expr_fdr_df_list[[version]] <- expr_fdr_df
      cat("  ✓ Built: Expr_FDR_df (expression + original annotations + FC/FDR, no SubMito)\n")
    } else {
      cat("  - Warning: expr_matrix missing or Gene column absent; cannot build Expr_FDR_df\n")
    }
    
    # Check whether annotation column exists
    if (!(roc_annotation_column %in% colnames(data_transformed))) {
      cat(sprintf("  - Warning: annotation column %s not found; skipping ROC analysis\n", roc_annotation_column))
      next
    }
    
    # Get logFC columns
    logfc_columns <- grep("_logFC$", colnames(data_transformed), value = TRUE)
    if (length(logfc_columns) == 0) {
      cat("  - Warning: no logFC columns found; skipping ROC analysis\n")
      next
    }
    
    # Run ROC analysis
    cat("\nStep 2: ROC analysis\n")
    roc_analysis <- perform_roc_analysis(
      data = data_transformed,
      logfc_columns = logfc_columns,
      annotation_column = roc_annotation_column,
      tp_label = tp_label,
      fp_label = fp_label,
      direction = "<"
    )
    
    if (length(roc_analysis$roc_objects) == 0) {
      cat("  - Warning: ROC analysis failed; skipping\n")
      next
    }
    
    # Compute optimal thresholds
    cat("\nStep 3: Compute optimal thresholds\n")
    thresholds <- calculate_optimal_thresholds(roc_analysis$roc_dfs, min_tp = min_tp)
    all_thresholds[[version]] <- thresholds
    cat(sprintf("  - Thresholds: %s\n", paste(round(thresholds, 3), collapse = ", ")))
    
    # Save ROC data
    cat("\nStep 4: Save ROC data\n")
    output_xlsx <- file.path(dir_config$output, 
                             paste0("Module08_", version, "_ROC_data.xlsx"))
    wb <- createWorkbook()
    for (i in seq_along(roc_analysis$roc_dfs)) {
      sheet_name <- names(roc_analysis$roc_dfs)[i]
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, roc_analysis$roc_dfs[[i]])
    }
    saveWorkbook(wb, output_xlsx, overwrite = TRUE)
    cat(sprintf("  ✓ Saved: %s\n", basename(output_xlsx)))
    
    # Plot ROC curves
    cat("\nStep 5: Plot ROC curves\n")
    output_roc_pdf <- file.path(dir_config$output, 
                                 paste0("Module08_", version, "_ROC_curves.pdf"))
    plot_roc_curves(roc_analysis$roc_objects, roc_analysis$roc_dfs, output_roc_pdf, version)
    
    # Plot Youden Index
    cat("\nStep 6: Plot Youden Index\n")
    output_youden_pdf <- file.path(dir_config$output, 
                                    paste0("Module08_", version, "_Youden_Index.pdf"))
    plot_youden_index(roc_analysis$roc_dfs, output_youden_pdf)
    
    # Save results
    roc_results[[version]] <- list(
      roc_objects = roc_analysis$roc_objects,
      roc_dfs = roc_analysis$roc_dfs,
      thresholds = thresholds,
      data_transformed = data_transformed
    )
    
    cat(sprintf("\n✓ Version %s complete\n", version))
  }
  
  # Save all thresholds
  if (length(all_thresholds) > 0) {
    cat("\nSave all optimal thresholds...\n")
    output_thresholds_xlsx <- file.path(dir_config$output, 
                                         "Module08_all_thresholds.xlsx")
    wb <- createWorkbook()
    for (version in names(all_thresholds)) {
      sheet_name <- substr(version, 1, 31)  # Excel limit
      addWorksheet(wb, sheetName = sheet_name)
      
      # Get ROC data frame for this version to extract comparison names
      version_roc_dfs <- roc_results[[version]]$roc_dfs
      comparison_names <- sapply(version_roc_dfs, function(df) {
        # Extract comparison names from column names (same as Youden plot)
        tp_fp_col <- grep("_TP_FP$", colnames(df), value = TRUE)[1]
        # Remove '_logFC_TP_FP' suffix for comparison name
        sub("_logFC_TP_FP$", "", tp_fp_col)
      })
      
      threshold_df <- data.frame(
        Comparison = comparison_names,
        Threshold = all_thresholds[[version]],
        stringsAsFactors = FALSE
      )
      writeData(wb, sheet = sheet_name, threshold_df)
    }
    saveWorkbook(wb, output_thresholds_xlsx, overwrite = TRUE)
    cat(sprintf("  ✓ Saved: %s\n", basename(output_thresholds_xlsx)))
  }
  
  # Save Expr_FDR_df_list (one sheet per version)
  if (length(expr_fdr_df_list) > 0) {
    cat("\nSave Expr_FDR_df_list...\n")
    expr_fdr_xlsx <- file.path(dir_config$output, "Module08_Expr_FDR_df_list.xlsx")
    wb_expr <- createWorkbook()
    for (version in names(expr_fdr_df_list)) {
      sheet_name <- substr(version, 1, 31)
      addWorksheet(wb_expr, sheetName = sheet_name)
      writeData(wb_expr, sheet = sheet_name, expr_fdr_df_list[[version]])
    }
    saveWorkbook(wb_expr, expr_fdr_xlsx, overwrite = TRUE)
    cat(sprintf("  ✓ Saved: %s\n", basename(expr_fdr_xlsx)))
  }
  
  cat("\n========================================\n")
  cat("Module 8 complete\n")
  cat("========================================\n")
  
  return(list(
    roc_results = roc_results,
    all_thresholds = all_thresholds,
    data_with_submito = data_with_submito_list,
    expr_fdr_df_list = expr_fdr_df_list
  ))
}
