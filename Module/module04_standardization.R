# Module 04: Data Standardization
# Functions: Log2 transform, global normalization (QNorm/QNorm_Fix/QNorm_Limma/VSN/RefNorm/MNorm), local normalization (Local QNorm/QNorm_Fix/QNorm_Limma/VSN/RefNorm/MNorm)
# Author: CodeNorm Pipeline
# Date: 2024

# NA-safe quantile normalization: fill with column medians, normalize, then restore NA
quantile_normalize_preserve_na <- function(mat) {
  na_mask <- is.na(mat)
  col_medians <- apply(mat, 2, function(x) {
    med <- median(x, na.rm = TRUE)
    if (is.na(med)) 0 else med
  })
  
  filled_mat <- mat
  for (i in seq_len(ncol(filled_mat))) {
    filled_mat[is.na(filled_mat[, i]), i] <- col_medians[i]
  }
  
  normalized <- normalize.quantiles(filled_mat)
  colnames(normalized) <- colnames(mat)
  rownames(normalized) <- rownames(mat)
  
  normalized[na_mask] <- NA
  normalized
}

# Align target column to reference distribution (RefNorm helper)
align_to_reference <- function(vec, ref_vec) {
  if (all(is.na(vec))) return(vec)
  ref_sorted <- sort(ref_vec, na.last = NA)
  if (length(ref_sorted) == 0) return(vec)
  ref_probs <- seq(0, 1, length.out = length(ref_sorted))
  
  idx <- which(!is.na(vec))
  if (length(idx) == 0) return(vec)
  
  vec_non_na <- vec[idx]
  vec_order <- order(vec_non_na)
  
  target_vals <- approx(ref_probs, ref_sorted, xout = seq(0, 1, length.out = length(idx)), rule = 2)$y
  
  new_vals <- vec_non_na
  new_vals[vec_order] <- target_vals
  
  vec[idx] <- new_vals
  vec
}

# Select reference sample (RefNorm)
choose_reference_sample <- function(cols, data_matrix, sampleGroup, context_only_control = FALSE) {
  ctx_map <- setNames(sampleGroup$Context, sampleGroup$FinalName)
  contexts <- ctx_map[cols]
  valid_counts <- apply(data_matrix[, cols, drop = FALSE], 2, function(x) sum(!is.na(x)))
  
  pick_by_max <- function(pool_cols) {
    if (length(pool_cols) == 0) return(NULL)
    pool_counts <- valid_counts[pool_cols]
    pool_cols[which.max(pool_counts)]
  }
  
  if (context_only_control) {
    med_count <- median(valid_counts, na.rm = TRUE)
    if (is.na(med_count)) med_count <- 0
    diffs <- abs(valid_counts - med_count)
    ref_col <- cols[which.min(diffs)]
    return(ref_col)
  }
  
  exp_cols <- cols[contexts == "Experiment"]
  ref_col <- pick_by_max(exp_cols)
  if (is.null(ref_col)) {
    ref_col <- pick_by_max(cols)
  }
  ref_col
}

# VSN transform while keeping all-NA rows in place (avoid "only NA elements" warning)
vsn_transform_preserve_na <- function(mat) {
  if (!require(vsn, quietly = TRUE)) {
    stop("❌ vsn package required: BiocManager::install('vsn')")
  }
  row_all_na <- apply(mat, 1, function(x) all(is.na(x)))
  mat_use <- mat[!row_all_na, , drop = FALSE]
  if (nrow(mat_use) == 0) return(mat)
  min_val <- suppressWarnings(min(mat_use, na.rm = TRUE))
  if (is.finite(min_val) && min_val <= 0) {
    shift <- abs(min_val) + 1e-3
    mat_use <- mat_use + shift
  }
  fit_vsn <- vsn::vsn2(mat_use)
  transformed <- predict(fit_vsn, newdata = mat_use)
  out <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  out[!row_all_na, ] <- transformed
  colnames(out) <- colnames(mat)
  rownames(out) <- rownames(mat)
  out
}

#' Module 04: Data Standardization
#' 
#' @description
#' Apply log2 transformation and multiple normalization approaches, including global and local methods.
#' 
#' @param dir_config Directory config list (from Module 1)
#' @param data_annotated Annotated data frame (from Module 3)
#' @param sampleGroup Sample grouping information (from Module 2)
#' @param norm_types Character vector specifying normalization types to run.
#'   Options: "noNorm", "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
#'            "Global_VSN", "Global_RefNorm", "Global_MNorm",
#'            "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
#'            "Local_VSN", "Local_RefNorm", "Local_MNorm"
#'   Default: c("noNorm", "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
#'            "Global_VSN", "Global_RefNorm", "Global_MNorm",
#'            "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
#'            "Local_VSN", "Local_RefNorm", "Local_MNorm")
#'   Note: must include at least one global normalization (Global_QNorm / Global_QNorm_Fix / Global_QNorm_Limma / Global_VSN / Global_RefNorm / Global_MNorm)
#' 
#' @return List containing:
#'   - standardized_data_list: list of all normalized versions
#'   - norm_types_used: normalization types actually used
#' 
#' @details
#' - Log2 transform: apply log2 to numeric columns
#' - Global normalization: normalize across all samples
#' - Local normalization: normalize within each bioGroup then merge
#' - Outputs boxplots comparing normalization methods
#' 
#' @export
module04_standardization <- function(dir_config, 
                                     data_annotated, 
                                     sampleGroup,
                                     norm_types = c("noNorm",
                                                   "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
                                                   "Global_VSN", "Global_RefNorm", "Global_MNorm",
                                                   "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
                                                   "Local_VSN", "Local_RefNorm", "Local_MNorm")) {
  
  cat("\n=== Module 04: Data Standardization ===\n")
  
  # 1. Input validation ####
  cat("\n[1] Validate input data...\n")
  
  if (!all(c("reference", "output") %in% names(dir_config))) {
    stop("❌ dir_config must include reference and output paths")
  }
  
  if (!is.data.frame(data_annotated)) {
    stop("❌ data_annotated must be a data frame")
  }
  
  if (!"Gene" %in% colnames(data_annotated)) {
    stop("❌ data_annotated must contain a Gene column")
  }
  
  if (!is.data.frame(sampleGroup)) {
    stop("❌ sampleGroup must be a data frame")
  }
  
  if (!"bioGroup" %in% colnames(sampleGroup)) {
    stop("❌ sampleGroup must contain a bioGroup column")
  }
  if (any(c("Global_RefNorm", "Local_RefNorm") %in% norm_types) && !"Context" %in% colnames(sampleGroup)) {
    stop("❌ sampleGroup must include Context column to use RefNorm")
  }
  
  # Validate norm_types
  valid_types <- c("noNorm",
                   "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
                   "Global_VSN", "Global_RefNorm", "Global_MNorm",
                   "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
                   "Local_VSN", "Local_RefNorm", "Local_MNorm")
  if (!all(norm_types %in% valid_types)) {
    stop("❌ norm_types contains invalid values; valid options: ", paste(valid_types, collapse = ", "))
  }
  
  # Ensure at least one global normalization
  if (!any(c("Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma", "Global_VSN", "Global_RefNorm", "Global_MNorm") %in% norm_types)) {
    stop("❌ Must include at least one global normalization method (Global_QNorm / Global_QNorm_Fix / Global_QNorm_Limma / Global_VSN / Global_RefNorm / Global_MNorm)")
  }
  
  cat("✓ Input validation passed\n")
  cat(sprintf("  - Data dimensions: %d rows × %d columns\n", nrow(data_annotated), ncol(data_annotated)))
  cat(sprintf("  - Sample count: %d\n", nrow(sampleGroup)))
  cat(sprintf("  - Normalization types: %s\n", paste(norm_types, collapse = ", ")))
  
  # 2. Identify data and annotation columns ####
  cat("\n[2] Identify data and annotation columns...\n")
  
  # Identify annotation columns (ending with "Localization")
  annotation_cols <- grep("Localization$", colnames(data_annotated), value = TRUE)
  
  if (length(annotation_cols) == 0) {
    stop("❌ No annotation columns found (should end with 'Localization')")
  }
  
  # Data columns: everything except Gene and annotation columns
  data_cols <- setdiff(colnames(data_annotated), c("Gene", annotation_cols))
  
  cat(sprintf("✓ Identification complete\n"))
  cat(sprintf("  - Data columns: %d\n", length(data_cols)))
  cat(sprintf("  - Annotation columns: %d (%s)\n", 
              length(annotation_cols), paste(annotation_cols, collapse = ", ")))
  
  # 3. Log2 transform ####
  cat("\n[3] Perform Log2 transform...\n")
  
  data_log2 <- data_annotated %>%
    mutate(across(all_of(data_cols), ~log2(.)))
  
  cat("✓ Log2 transform complete\n")
  
  # Save log2 data to CSV
  csv_file <- file.path(dir_config$output, "Module04_log2.csv")
  write.csv(data_log2, csv_file, row.names = FALSE)
  cat(sprintf("✓ Saved: %s\n", csv_file))
  
  # Assign colors for each bioGroup ####
  # Build mapping from samples to bioGroup
  sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
  
  # Verify data columns are present in sampleGroup
  valid_data_cols <- intersect(data_cols, names(sample_to_group))
  if (length(valid_data_cols) < length(data_cols)) {
    missing_cols <- setdiff(data_cols, valid_data_cols)
    warning("⚠ These data columns were not found in sampleGroup; using grey:", 
            paste(missing_cols, collapse = ", "))
  }
  
  # Map each data column to bioGroup
  col_groups <- sample_to_group[valid_data_cols]
  unique_groups <- unique(col_groups)
  n_groups <- length(unique_groups)
  
  cat(sprintf("  - Detected %d distinct bioGroups\n", n_groups))
  
  # Build color palette
  if (n_groups <= 12) {
    # Use a soft palette
    group_colors <- scales::hue_pal()(n_groups)
  } else {
    # Fallback rainbow palette
    group_colors <- rainbow(n_groups)
  }
  names(group_colors) <- unique_groups
  
  # Assign color to each data column
  box_colors <- character(length(data_cols))
  names(box_colors) <- data_cols
  for (col in data_cols) {
    if (col %in% names(col_groups)) {
      box_colors[col] <- group_colors[col_groups[col]]
    } else {
      box_colors[col] <- "grey80"  # grey for columns not found
    }
  }
  
  # Plot log2 boxplot colored by bioGroup
  pdf_file <- file.path(dir_config$output, "Module04_log2_boxplot.pdf")
  pdf(pdf_file, width = 10, height = 5)
  boxplot(data_log2[, data_cols], 
          cex.axis = 0.4, 
          las = 2, 
          main = "Log2 Transformed Data",
          col = box_colors[data_cols])
  # Add legend
  legend("topright", 
         legend = names(group_colors), 
         fill = group_colors, 
         cex = 0.6,
         title = "bioGroup")
  dev.off()
  cat(sprintf("✓ Saved: %s\n", pdf_file))
  
  # 4. Initialize result list ####
  standardized_data_list <- list()
  
  # noNorm: log2 transformed data only
  if ("noNorm" %in% norm_types) {
    standardized_data_list[["noNorm"]] <- data_log2
    cat("✓ Added: noNorm (log2 only)\n")
  }
  
  # 5. Global normalization ####
  if (any(c("Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma", "Global_VSN", "Global_RefNorm", "Global_MNorm") %in% norm_types)) {
    cat("\n[4] Run global normalization...\n")
    
    # Load preprocessCore
    if (!require(preprocessCore, quietly = TRUE)) {
      stop("❌ preprocessCore package required: BiocManager::install('preprocessCore')")
    }
    
    # Extract data matrix for normalization
    data_matrix <- as.matrix(data_log2[, data_cols])
    
    # 5.1 Global Quantile Normalization ####
    if ("Global_QNorm" %in% norm_types) {
      cat("  Running Global Quantile Normalization...\n")
      
      normalized_matrix <- normalize.quantiles(data_matrix)
      colnames(normalized_matrix) <- data_cols
      
      data_qnorm <- data_log2
      data_qnorm[, data_cols] <- as.data.frame(normalized_matrix)
      
      standardized_data_list[["Global_QNorm"]] <- data_qnorm
      
      # Save CSV
      csv_file <- file.path(dir_config$output, "Module04_Global_QNorm.csv")
      write.csv(data_qnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_QNorm complete; saved: %s\n", csv_file))
    }
    
    # 5.1b Global Quantile Normalization (NA-safe) ####
    if ("Global_QNorm_Fix" %in% norm_types) {
      cat("  Running Global Quantile Normalization (QNorm_Fix, mask NA then restore)...\n")
      
      normalized_matrix_fix <- quantile_normalize_preserve_na(data_matrix)
      colnames(normalized_matrix_fix) <- data_cols
      
      data_qnorm_fix <- data_log2
      data_qnorm_fix[, data_cols] <- as.data.frame(normalized_matrix_fix)
      
      standardized_data_list[["Global_QNorm_Fix"]] <- data_qnorm_fix
      
      # Save CSV
      csv_file <- file.path(dir_config$output, "Module04_Global_QNorm_Fix.csv")
      write.csv(data_qnorm_fix, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_QNorm_Fix complete; saved: %s\n", csv_file))
    }
    
    # 5.1c Global Quantile Normalization (limma, NA-safe) ####
    if ("Global_QNorm_Limma" %in% norm_types) {
      cat("  Running Global Quantile Normalization (limma)...\n")
      if (!require(limma, quietly = TRUE)) {
        stop("❌ limma package required: BiocManager::install('limma')")
      }
      normalized_matrix_limma <- limma::normalizeBetweenArrays(data_matrix, method = "quantile")
      colnames(normalized_matrix_limma) <- data_cols
      
      data_qnorm_limma <- data_log2
      data_qnorm_limma[, data_cols] <- as.data.frame(normalized_matrix_limma)
      
      standardized_data_list[["Global_QNorm_Limma"]] <- data_qnorm_limma
      
      csv_file <- file.path(dir_config$output, "Module04_Global_QNorm_Limma.csv")
      write.csv(data_qnorm_limma, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_QNorm_Limma complete; saved: %s\n", csv_file))
    }
    
    # 5.1d Global VSN ####
    if ("Global_VSN" %in% norm_types) {
      cat("  Running Global VSN...\n")
      vsn_input <- as.matrix(data_annotated[, data_cols])
      normalized_matrix_vsn <- vsn_transform_preserve_na(vsn_input)
      colnames(normalized_matrix_vsn) <- data_cols
      
      data_vsn <- data_log2
      data_vsn[, data_cols] <- as.data.frame(normalized_matrix_vsn)
      
      standardized_data_list[["Global_VSN"]] <- data_vsn
      
      csv_file <- file.path(dir_config$output, "Module04_Global_VSN.csv")
      write.csv(data_vsn, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_VSN complete; saved: %s\n", csv_file))
    }
    
    # 5.1e Global Reference Normalization (context-aware) ####
    if ("Global_RefNorm" %in% norm_types) {
      cat("  Running Global Reference Normalization (context-aware)...\n")
      
      sample_ctx_map <- setNames(sampleGroup$Context, sampleGroup$FinalName)
      ref_col <- choose_reference_sample(data_cols, data_matrix, sampleGroup, context_only_control = FALSE)
      if (is.null(ref_col)) {
        stop("❌ Unable to select reference sample (Global_RefNorm)")
      }
      cat(sprintf("    Reference sample: %s (Context=%s)\n", ref_col, sample_ctx_map[[ref_col]]))
      
      ref_vec <- data_matrix[, ref_col]
      refnorm_mat <- data_matrix
      for (cn in data_cols) {
        refnorm_mat[, cn] <- align_to_reference(refnorm_mat[, cn], ref_vec)
      }
      
      data_refnorm <- data_log2
      data_refnorm[, data_cols] <- as.data.frame(refnorm_mat)
      
      standardized_data_list[["Global_RefNorm"]] <- data_refnorm
      
      csv_file <- file.path(dir_config$output, "Module04_Global_RefNorm.csv")
      write.csv(data_refnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_RefNorm complete; saved: %s\n", csv_file))
    }
    
    # 5.2 Global Median Normalization ####
    if ("Global_MNorm" %in% norm_types) {
      cat("  Running Global Median Normalization...\n")
      
      # Compute median per sample
      sample_medians <- apply(data_matrix, 2, median, na.rm = TRUE)
      cat(sprintf("    Sample median range: %.2f - %.2f\n", 
                  min(sample_medians, na.rm = TRUE), 
                  max(sample_medians, na.rm = TRUE)))
      
      # Compute global median
      global_median <- median(sample_medians, na.rm = TRUE)
      cat(sprintf("    Global median: %.2f\n", global_median))
      
      # Normalize
      normalized_matrix <- sweep(data_matrix, 2, sample_medians / global_median, "/")
      
      data_mnorm <- data_log2
      data_mnorm[, data_cols] <- as.data.frame(normalized_matrix)
      
      standardized_data_list[["Global_MNorm"]] <- data_mnorm
      
      # Save CSV
      csv_file <- file.path(dir_config$output, "Module04_Global_MNorm.csv")
      write.csv(data_mnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_MNorm complete; saved: %s\n", csv_file))
    }
    
    # Plot global normalization boxplots
    pdf_file <- file.path(dir_config$output, "Module04_Global_Norm_boxplot.pdf")
    pdf(pdf_file, width = 10, height = 5)
    
    for (norm_name in c("Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma", "Global_VSN", "Global_RefNorm", "Global_MNorm")) {
      if (norm_name %in% names(standardized_data_list)) {
        boxplot(standardized_data_list[[norm_name]][, data_cols], 
                cex.axis = 0.4, 
                las = 2, 
                main = norm_name,
                col = box_colors[data_cols])
        # Add legend
        legend("topright", 
               legend = names(group_colors), 
               fill = group_colors, 
               cex = 0.6,
               title = "bioGroup")
      }
    }
    
    dev.off()
    cat(sprintf("✓ Saved: %s\n", pdf_file))
  }
  
  # 6. Local normalization ####
  if (any(c("Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma", "Local_VSN", "Local_RefNorm", "Local_MNorm") %in% norm_types)) {
    cat("\n[5] Run local normalization (by bioGroup)...\n")
    
    # Ensure preprocessCore is loaded
    if (!require(preprocessCore, quietly = TRUE)) {
      stop("❌ preprocessCore package required: BiocManager::install('preprocessCore')")
    }
    
    # Build sample to bioGroup mapping
    # sampleGroup$FinalName aligns to data_cols
    sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
    
    # Ensure all data columns have a mapped bioGroup
    missing_samples <- setdiff(data_cols, names(sample_to_group))
    if (length(missing_samples) > 0) {
      warning("⚠ These samples lack bioGroup in sampleGroup and will be skipped:",
              paste(missing_samples, collapse = ", "))
      data_cols <- setdiff(data_cols, missing_samples)
    }
    
    # Create group_vec
    group_vec <- sample_to_group[data_cols]
    
    cat(sprintf("  - Unique bioGroups: %d\n", length(unique(group_vec))))
    cat(sprintf("  - bioGroup distribution:\n"))
    for (grp in unique(group_vec)) {
      n_samples <- sum(group_vec == grp)
      cat(sprintf("    %s: %d samples\n", grp, n_samples))
    }
    
    # Extract data matrix
    data_matrix <- as.matrix(data_log2[, data_cols])
    
    # 6.1 Local Median Normalization ####
    if ("Local_MNorm" %in% norm_types) {
      cat("\n  Running Local Median Normalization...\n")
      
      # Initialize normalized matrix
      mnormalized_mat <- data_matrix
      
      # Normalize within each bioGroup
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    Processing %s (%d samples)...\n", grp, length(grp_cols)))
        
        # Extract sub-matrix for the group
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        
        # Median per column
        sample_medians <- apply(sub_mat, 2, median, na.rm = TRUE)
        
        # Group-level median
        group_global_median <- median(sample_medians, na.rm = TRUE)
        
        # Normalize
        normalized_sub_mat <- sweep(sub_mat, 2, sample_medians / group_global_median, "/")
        
        # Write back
        mnormalized_mat[, grp_cols] <- normalized_sub_mat
      }
      
      # Build final data frame
      data_local_mnorm <- data_log2
      data_local_mnorm[, data_cols] <- as.data.frame(mnormalized_mat)
      
      standardized_data_list[["Local_MNorm"]] <- data_local_mnorm
      
      # Save CSV
      csv_file <- file.path(dir_config$output, "Module04_Local_MNorm.csv")
      write.csv(data_local_mnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_MNorm complete; saved: %s\n", csv_file))
    }
    
    # 6.2 Local Quantile Normalization ####
    if ("Local_QNorm" %in% norm_types) {
      cat("\n  Running Local Quantile Normalization...\n")
      
      # Initialize normalized matrix
      qnormalized_mat <- data_matrix
      
      # Normalize within each bioGroup
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    Processing %s (%d samples)...\n", grp, length(grp_cols)))
        
        # Extract sub-matrix
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        
        # Quantile normalization
        normalized_sub_mat <- normalize.quantiles(as.matrix(sub_mat))
        colnames(normalized_sub_mat) <- grp_cols
        
        # Write back
        qnormalized_mat[, grp_cols] <- normalized_sub_mat
      }
      
      # Build final data frame
      data_local_qnorm <- data_log2
      data_local_qnorm[, data_cols] <- as.data.frame(qnormalized_mat)
      
      standardized_data_list[["Local_QNorm"]] <- data_local_qnorm
      
      # Save CSV
      csv_file <- file.path(dir_config$output, "Module04_Local_QNorm.csv")
      write.csv(data_local_qnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_QNorm complete; saved: %s\n", csv_file))
    }
    
    # 6.2b Local Quantile Normalization (NA-safe) ####
    if ("Local_QNorm_Fix" %in% norm_types) {
      cat("\n  Running Local Quantile Normalization (QNorm_Fix, mask NA then restore)...\n")
      
      qnormalized_fix_mat <- data_matrix
      
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    Processing %s (%d samples)...\n", grp, length(grp_cols)))
        
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        normalized_sub_mat <- quantile_normalize_preserve_na(sub_mat)
        colnames(normalized_sub_mat) <- grp_cols
        
        qnormalized_fix_mat[, grp_cols] <- normalized_sub_mat
      }
      
      data_local_qnorm_fix <- data_log2
      data_local_qnorm_fix[, data_cols] <- as.data.frame(qnormalized_fix_mat)
      
      standardized_data_list[["Local_QNorm_Fix"]] <- data_local_qnorm_fix
      
      csv_file <- file.path(dir_config$output, "Module04_Local_QNorm_Fix.csv")
      write.csv(data_local_qnorm_fix, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_QNorm_Fix complete; saved: %s\n", csv_file))
    }

    # 6.2c Local Quantile Normalization (limma) ####
    if ("Local_QNorm_Limma" %in% norm_types) {
      cat("\n  Running Local Quantile Normalization (limma)...\n")
      if (!require(limma, quietly = TRUE)) {
        stop("❌ limma package required: BiocManager::install('limma')")
      }
      
      qnormalized_limma_mat <- data_matrix
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    Processing %s (%d samples)...\n", grp, length(grp_cols)))
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        normalized_sub_mat <- limma::normalizeBetweenArrays(sub_mat, method = "quantile")
        colnames(normalized_sub_mat) <- grp_cols
        qnormalized_limma_mat[, grp_cols] <- normalized_sub_mat
      }
      
      data_local_qnorm_limma <- data_log2
      data_local_qnorm_limma[, data_cols] <- as.data.frame(qnormalized_limma_mat)
      
      standardized_data_list[["Local_QNorm_Limma"]] <- data_local_qnorm_limma
      
      csv_file <- file.path(dir_config$output, "Module04_Local_QNorm_Limma.csv")
      write.csv(data_local_qnorm_limma, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_QNorm_Limma complete; saved: %s\n", csv_file))
    }
    
    # 6.2d Local VSN ####
    if ("Local_VSN" %in% norm_types) {
      cat("\n  Running Local VSN...\n")
      
      vsn_mat <- data_matrix
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    Processing %s (%d samples)...\n", grp, length(grp_cols)))
        sub_mat <- as.matrix(data_annotated[, grp_cols, drop = FALSE])
        normalized_sub_mat <- vsn_transform_preserve_na(sub_mat)
        colnames(normalized_sub_mat) <- grp_cols
        vsn_mat[, grp_cols] <- normalized_sub_mat
      }
      
      data_local_vsn <- data_log2
      data_local_vsn[, data_cols] <- as.data.frame(vsn_mat)
      
      standardized_data_list[["Local_VSN"]] <- data_local_vsn
      
      csv_file <- file.path(dir_config$output, "Module04_Local_VSN.csv")
      write.csv(data_local_vsn, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_VSN complete; saved: %s\n", csv_file))
    }
    
    # 6.2e Local Reference Normalization (context-aware) ####
    if ("Local_RefNorm" %in% norm_types) {
      cat("\n  Running Local Reference Normalization (context-aware)...\n")
      
      refnorm_mat <- data_matrix
      sample_ctx_map <- setNames(sampleGroup$Context, sampleGroup$FinalName)
      
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        context_vec <- sample_ctx_map[grp_cols]
        only_control <- all(context_vec == "Control", na.rm = TRUE) && any(!is.na(context_vec))
        mode_label <- if (only_control) "Control mode (closest median)" else "Experiment first"
        cat(sprintf("    Processing %s (%d samples, %s)...\n", grp, length(grp_cols), mode_label))
        
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        ref_col <- choose_reference_sample(grp_cols, sub_mat, sampleGroup, context_only_control = only_control)
        if (is.null(ref_col)) {
          warning(sprintf("⚠ Group %s has no reference sample; skipping RefNorm", grp))
          next
        }
        cat(sprintf("      Reference sample: %s (Context=%s)\n", ref_col, sample_ctx_map[[ref_col]]))
        
        ref_vec <- sub_mat[, ref_col]
        for (cn in grp_cols) {
          refnorm_mat[, cn] <- align_to_reference(sub_mat[, cn], ref_vec)
        }
      }
      
      data_local_refnorm <- data_log2
      data_local_refnorm[, data_cols] <- as.data.frame(refnorm_mat)
      
      standardized_data_list[["Local_RefNorm"]] <- data_local_refnorm
      
      csv_file <- file.path(dir_config$output, "Module04_Local_RefNorm.csv")
      write.csv(data_local_refnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_RefNorm complete; saved: %s\n", csv_file))
    }
    
    # Plot local normalization boxplots
    pdf_file <- file.path(dir_config$output, "Module04_Local_Norm_boxplot.pdf")
    pdf(pdf_file, width = 10, height = 5)
    
    for (norm_name in c("Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma", "Local_VSN", "Local_RefNorm", "Local_MNorm")) {
      if (norm_name %in% names(standardized_data_list)) {
        boxplot(standardized_data_list[[norm_name]][, data_cols], 
                cex.axis = 0.4, 
                las = 2, 
                main = norm_name,
                col = box_colors[data_cols])
        # Add legend
        legend("topright", 
               legend = names(group_colors), 
               fill = group_colors, 
               cex = 0.6,
               title = "bioGroup")
      }
    }
    
    dev.off()
    cat(sprintf("✓ Saved: %s\n", pdf_file))
  }
  
  # 7. Combined comparison boxplot ####
  cat("\n[6] Generate combined comparison boxplot...\n")
  
  # Reorder according to requested sequence
  standardized_data_list <- standardized_data_list[intersect(norm_types, names(standardized_data_list))]
  
  pdf_file <- file.path(dir_config$output, "Module04_All_Norm_comparison_boxplot.pdf")
  pdf(pdf_file, width = 10, height = 5)
  
  for (norm_name in names(standardized_data_list)) {
    boxplot(standardized_data_list[[norm_name]][, data_cols], 
            cex.axis = 0.4, 
            las = 2, 
            main = norm_name,
            col = box_colors[data_cols])
    # Add legend
    legend("topright", 
           legend = names(group_colors), 
           fill = group_colors, 
           cex = 0.6,
           title = "bioGroup")
  }
  
  dev.off()
  cat(sprintf("✓ Saved: %s\n", pdf_file))
  
  # 8. Summary ####
  cat("\n=== Module 04 complete ===\n")
  cat(sprintf("✓ Generated %d normalization versions:\n", length(standardized_data_list)))
  for (norm_name in names(standardized_data_list)) {
    cat(sprintf("  - %s\n", norm_name))
  }
  
  cat("\nOutput files:\n")
  cat(sprintf("  - CSV files: Output/Module04_*.csv (%d files)\n", 
              length(standardized_data_list) + 1))  # +1 for log2
  cat(sprintf("  - PDF files: Output/Module04_*_boxplot.pdf (4 files)\n"))
  
  # Return result
  return(list(
    standardized_data_list = standardized_data_list,
    norm_types_used = names(standardized_data_list)
  ))
}
