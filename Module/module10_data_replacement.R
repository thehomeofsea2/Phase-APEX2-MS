# ============================================================================
# Module 10: Data Replacement (Module 5 global normalization + imputation)
# ============================================================================
# Aligns with CleanCode.R (2303-2361)
#
# Goals:
# - For Module 9 merged_data_A / merged_data_B (AfterROC_merged), perform targeted replacement:
#   keep all NA; replace all non-NA values with Module 5 Global_QNorm_Imputed;
#   if unavailable, use Global_MNorm_Imputed.
# - Generate sampling check reports (Check/ subdirectory) to verify replacement.
# - Export replacement workbooks (Method A and Method B).
#
# Input:
# - dir_config: includes output paths
# - sampleGroup: identify sample columns (FinalName)
# - merged_data_A / merged_data_B: Module 9 merged data (per version)
# - imputed_data_list: Module 5 outputs including Global_QNorm_Imputed / Global_MNorm_Imputed
# - selected_versions: versions to process (NULL = all from merged_data_*)
# - prefer_version: preferred imputed version (default "Global_QNorm_Imputed")
# - fallback_version: fallback imputed version (default "Global_MNorm_Imputed")
# - test_n: max genes for sampling check (default 20)
#
# Output:
# - Output/Module10_{method}_{version}_AfterFirstROC_Intersect.xlsx: replaced result
# - Output/Check/Module10_{method}_{version}_replacement_check.csv: sampling check table
# - Output/Check/Module10_{method}_{version}_replacement_summary.txt: replacement summary
#
# Return:
# list(
#   replaced_data_A = list(version -> dataframe),
#   replaced_data_B = list(version -> dataframe),
#   base_version_used = "<chosen imputed version>"
# )
# ============================================================================

module10_data_replacement <- function(dir_config,
                                      sampleGroup,
                                      merged_data_A,
                                      merged_data_B,
                                      imputed_data_list,
                                      selected_versions = NULL,
                                      prefer_version = "Global_QNorm_Imputed",
                                      fallback_version = "Global_MNorm_Imputed",
                                      test_n = 20) {
  cat("\n========================================\n")
  cat("Module 10: Data replacement (Module 5 global normalization + imputation)\n")
  cat("========================================\n")
  
  # Dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr package required")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("openxlsx package required")
  
  `%>%` <- dplyr::`%>%`
  
  # Create Check subdirectory
  check_dir <- file.path(dir_config$output, "Check")
  if (!dir.exists(check_dir)) {
    dir.create(check_dir, recursive = TRUE)
    cat(sprintf("✓ Created directory: %s\n", check_dir))
  }
  
  # Choose base imputed version
  base_version <- if (prefer_version %in% names(imputed_data_list)) {
    prefer_version
  } else if (fallback_version %in% names(imputed_data_list)) {
    fallback_version
  } else {
    stop(sprintf("%s or %s not found in imputed_data_list",
                 prefer_version, fallback_version))
  }
  cat(sprintf("✓ Base replacement version: %s\n", base_version))
  
  base_imputed <- imputed_data_list[[base_version]]
  if (!("Gene" %in% colnames(base_imputed))) {
    stop(sprintf("Version %s missing Gene column", base_version))
  }
  
  # Versions to process
  existing_versions <- unique(c(names(merged_data_A), names(merged_data_B)))
  if (is.null(selected_versions)) {
    selected_versions <- existing_versions
  } else {
    missing <- setdiff(selected_versions, existing_versions)
    if (length(missing) > 0) {
      stop(sprintf("Versions not present in Module 9 results: %s",
                   paste(missing, collapse = ", ")))
    }
  }
  cat(sprintf("✓ Versions to process: %s\n\n", paste(selected_versions, collapse = ", ")))
  
  # Helper: perform replacement on a merged dataframe
  replace_nonNA_cellwise <- function(df, ref_df, cols) {
    # Align by Gene and bring in reference values
    merged <- dplyr::left_join(df, ref_df, by = "Gene", suffix = c("", ".ref"))
    
    replaced_count <- 0L
    na_preserved <- 0L
    total_cells <- 0L
    
    for (col in cols) {
      ref_col <- paste0(col, ".ref")
      if (!(ref_col %in% colnames(merged))) next
      
      v_old <- merged[[col]]
      v_ref <- merged[[ref_col]]
      
      # Stats
      mask_non_na <- !is.na(v_old)
      total_cells <- total_cells + length(v_old)
      na_preserved <- na_preserved + sum(!mask_non_na, na.rm = TRUE)
      
      # Replacement: non-NA use ref; NA preserved
      v_new <- v_old
      v_new[mask_non_na & !is.na(v_ref)] <- v_ref[mask_non_na & !is.na(v_ref)]
      replaced_count <- replaced_count + sum(mask_non_na & !is.na(v_ref), na.rm = TRUE)
      merged[[col]] <- v_new
      
      # Drop helper column
      merged[[ref_col]] <- NULL
    }
    
    # Drop any remaining .ref columns
    merged <- merged[, !grepl("\\.ref$", colnames(merged)), drop = FALSE]
    
    return(list(
      data = merged,
      replaced_count = replaced_count,
      na_preserved = na_preserved,
      total_cells = total_cells
    ))
  }
  
  replaced_data_A <- list()
  replaced_data_B <- list()
  
  # Process each version
  for (ver in selected_versions) {
    cat(sprintf("════════════════════════════════════════\n"))
    cat(sprintf("Process version: %s\n", ver))
    cat(sprintf("════════════════════════════════════════\n"))
    
    # Prepare base table (only columns needed for this version)
    lfq_cols <- sampleGroup$FinalName
    base_cols <- c("Gene", lfq_cols)
    base_use <- base_imputed %>%
      dplyr::select(dplyr::any_of(base_cols)) %>%
      dplyr::distinct(Gene, .keep_all = TRUE)
    
    # --- Method A ---
    if (ver %in% names(merged_data_A)) {
    cat("\n[Method A]\n")
      df_A <- merged_data_A[[ver]]
      
      if (!is.data.frame(df_A) || nrow(df_A) == 0) {
        cat("  ⚠ Data empty; skipped\n")
        replaced_data_A[[ver]] <- df_A
      } else {
        # Identify columns to replace (LFQ)
        replace_cols <- intersect(lfq_cols, colnames(df_A))
        replace_cols <- intersect(replace_cols, colnames(base_use))
        
        if (length(replace_cols) == 0) {
          cat("  ⚠ No LFQ columns to replace; skipped\n")
          replaced_data_A[[ver]] <- df_A
        } else {
          cat(sprintf("  - Original data: %d proteins, %d columns\n", nrow(df_A), ncol(df_A)))
          cat(sprintf("  - Columns to replace: %d LFQ columns\n", length(replace_cols)))
          
          # Perform replacement
          result_A <- replace_nonNA_cellwise(df_A, base_use, replace_cols)
          replaced_data_A[[ver]] <- result_A$data
          
          cat(sprintf("  ✓ Replacement done: %d non-NA cells replaced\n", result_A$replaced_count))
          cat(sprintf("  ✓ NA preserved: %d/%d cells\n", result_A$na_preserved, result_A$total_cells))
          
          # Sampling check
          set.seed(123)
          show_cols <- head(replace_cols, 3)  # max 3 columns
          check_genes <- head(unique(df_A$Gene), test_n)
          
          if (length(check_genes) > 0 && length(show_cols) > 0) {
            # Original values
            orig_df <- df_A %>% 
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(orig_df) <- c("Gene", paste0("orig_", show_cols))
            
            # Replaced values
            repl_df <- result_A$data %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(repl_df) <- c("Gene", paste0("repl_", show_cols))
            
            # Reference values
            imp_df <- base_use %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(imp_df) <- c("Gene", paste0("imp_", show_cols))
            
            # Merge
            check_out <- Reduce(function(x, y) dplyr::left_join(x, y, by = "Gene"),
                               list(orig_df, repl_df, imp_df))
            
            # Long format
            long_list <- list()
            for (sc in show_cols) {
              long_list[[sc]] <- data.frame(
                Gene = check_out$Gene,
                Sample = sc,
                Original = check_out[[paste0("orig_", sc)]],
                Replaced = check_out[[paste0("repl_", sc)]],
                Imputed = check_out[[paste0("imp_", sc)]],
                stringsAsFactors = FALSE
              )
            }
            check_long <- do.call(rbind, long_list)
            
            # Validate: non-NA original -> Replaced equals Imputed
            check_long$CellMatch <- with(check_long, ifelse(
              !is.na(Original),
              abs(as.numeric(Replaced) - as.numeric(Imputed)) < 1e-9,
              NA
            ))
            
            # Save check CSV to Check/
            check_file <- file.path(check_dir, 
                                   sprintf("Module10_A_%s_replacement_check.csv", ver))
            utils::write.csv(check_long, check_file, row.names = FALSE)
            cat(sprintf("  ✓ Sampling check saved: Check/%s\n", basename(check_file)))
          }
          
          # Summary text
          summary_lines <- c(
            sprintf("Version: %s - Method A", ver),
            sprintf("Original data: %d proteins, %d columns", nrow(df_A), ncol(df_A)),
            sprintf("Columns replaced: %d LFQ columns", length(replace_cols)),
            sprintf("Replacement: %d non-NA cells replaced", result_A$replaced_count),
            sprintf("NA preserved: %d/%d cells", result_A$na_preserved, result_A$total_cells),
            sprintf("Base version: %s", base_version)
          )
          sum_file <- file.path(check_dir,
                               sprintf("Module10_A_%s_replacement_summary.txt", ver))
          writeLines(summary_lines, con = sum_file)
          cat(sprintf("  ✓ Replacement summary saved: Check/%s\n", basename(sum_file)))
          
          # Export replaced data
          wb <- openxlsx::createWorkbook()
          openxlsx::addWorksheet(wb, sheetName = "AfterROC_merged")
          openxlsx::writeData(wb, sheet = "AfterROC_merged", result_A$data)
          out_xlsx <- file.path(dir_config$output,
                               sprintf("Module10_A_%s_AfterFirstROC_Intersect.xlsx", ver))
          openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
          cat(sprintf("  ✓ Replaced data exported: %s\n\n", basename(out_xlsx)))
        }
      }
    } else {
      cat("\n[Method A]\n")
      cat("  ⚠ Note: Method A missing this version; skipped\n\n")
    }
    
  # --- Method B ---
  if (ver %in% names(merged_data_B)) {
    cat("[Method B]\n")
      df_B <- merged_data_B[[ver]]
      
      if (!is.data.frame(df_B) || nrow(df_B) == 0) {
        cat("  ⚠ Data empty; skipped\n")
        replaced_data_B[[ver]] <- df_B
      } else {
        # Identify columns to replace (LFQ)
        replace_cols <- intersect(lfq_cols, colnames(df_B))
        replace_cols <- intersect(replace_cols, colnames(base_use))
        
        if (length(replace_cols) == 0) {
          cat("  ⚠ No LFQ columns to replace; skipped\n")
          replaced_data_B[[ver]] <- df_B
        } else {
          cat(sprintf("  - Original data: %d proteins, %d columns\n", nrow(df_B), ncol(df_B)))
          cat(sprintf("  - Columns to replace: %d LFQ columns\n", length(replace_cols)))
          
          # Perform replacement
          result_B <- replace_nonNA_cellwise(df_B, base_use, replace_cols)
          replaced_data_B[[ver]] <- result_B$data
          
          cat(sprintf("  ✓ Replacement done: %d non-NA cells replaced\n", result_B$replaced_count))
          cat(sprintf("  ✓ NA preserved: %d/%d cells\n", result_B$na_preserved, result_B$total_cells))
          
          # Sampling check
          set.seed(123)
          show_cols <- head(replace_cols, 3)  # max 3 columns
          check_genes <- head(unique(df_B$Gene), test_n)
          
          if (length(check_genes) > 0 && length(show_cols) > 0) {
            # Original values
            orig_df <- df_B %>% 
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(orig_df) <- c("Gene", paste0("orig_", show_cols))
            
            # Replaced values
            repl_df <- result_B$data %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(repl_df) <- c("Gene", paste0("repl_", show_cols))
            
            # Reference values
            imp_df <- base_use %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(imp_df) <- c("Gene", paste0("imp_", show_cols))
            
            # Merge
            check_out <- Reduce(function(x, y) dplyr::left_join(x, y, by = "Gene"),
                               list(orig_df, repl_df, imp_df))
            
            # Long format
            long_list <- list()
            for (sc in show_cols) {
              long_list[[sc]] <- data.frame(
                Gene = check_out$Gene,
                Sample = sc,
                Original = check_out[[paste0("orig_", sc)]],
                Replaced = check_out[[paste0("repl_", sc)]],
                Imputed = check_out[[paste0("imp_", sc)]],
                stringsAsFactors = FALSE
              )
            }
            check_long <- do.call(rbind, long_list)
            
            # Validate: non-NA original -> Replaced equals Imputed
            check_long$CellMatch <- with(check_long, ifelse(
              !is.na(Original),
              abs(as.numeric(Replaced) - as.numeric(Imputed)) < 1e-9,
              NA
            ))
            
            # Save check CSV to Check/
            check_file <- file.path(check_dir, 
                                   sprintf("Module10_B_%s_replacement_check.csv", ver))
            utils::write.csv(check_long, check_file, row.names = FALSE)
            cat(sprintf("  ✓ Sampling check saved: Check/%s\n", basename(check_file)))
          }
          
          # Summary text
          summary_lines <- c(
            sprintf("Version: %s - Method B", ver),
            sprintf("Original data: %d proteins, %d columns", nrow(df_B), ncol(df_B)),
            sprintf("Columns replaced: %d LFQ columns", length(replace_cols)),
            sprintf("Replacement: %d non-NA cells replaced", result_B$replaced_count),
            sprintf("NA preserved: %d/%d cells", result_B$na_preserved, result_B$total_cells),
            sprintf("Base version: %s", base_version)
          )
          sum_file <- file.path(check_dir,
                               sprintf("Module10_B_%s_replacement_summary.txt", ver))
          writeLines(summary_lines, con = sum_file)
          cat(sprintf("  ✓ Replacement summary saved: Check/%s\n", basename(sum_file)))
          
          # Export replaced data
          wb <- openxlsx::createWorkbook()
          openxlsx::addWorksheet(wb, sheetName = "AfterROC_merged")
          openxlsx::writeData(wb, sheet = "AfterROC_merged", result_B$data)
          out_xlsx <- file.path(dir_config$output,
                               sprintf("Module10_B_%s_AfterFirstROC_Intersect.xlsx", ver))
          openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
          cat(sprintf("  ✓ Replaced data exported: %s\n\n", basename(out_xlsx)))
        }
      }
    } else {
      cat("[Method B]\n")
      cat("  ⚠ Note: Method B missing this version; skipped\n\n")
    }
  }
  
  cat("════════════════════════════════════════\n")
  cat("✓ Module 10 complete\n")
  cat("════════════════════════════════════════\n")
  
  return(list(
    replaced_data_A = replaced_data_A,
    replaced_data_B = replaced_data_B,
    base_version_used = base_version
  ))
}
