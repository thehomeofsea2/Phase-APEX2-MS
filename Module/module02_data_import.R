# ============================================================================
# Module 2: Data Import and Sample Table
# ============================================================================
# Features:
#   1. Read TSV data files
#   2. Clean column names automatically (remove common prefixes/suffixes)
#   3. Generate sampleGroup_template.csv (with placeholder hints)
#   4. Read user-provided sampleGroup.csv
#   5. Generate FinalName and rename data columns
# 
# Input: dir_config
# Output: data_raw.RData, sampleGroup.RData
# ============================================================================

#' Read data and generate template
#' 
#' @param dir_config Directory configuration
#' @param file_pattern TSV file pattern, default "_matrix.*\\.tsv$"
#' @return List containing data and template
module02_read_and_generate_template <- function(dir_config, file_pattern = "_matrix.*\\.tsv$") {
  
  cat("\n----------------------------------------\n")
  cat("Step 1: Read data files\n")
  cat("----------------------------------------\n")
  
  # Read TSV files
  setwd(dir_config$rawdata)
  myfiles <- list.files(pattern = file_pattern)
  
  if (length(myfiles) == 0) {
    stop(sprintf("✗ Error: no files matching '%s' found in %s", dir_config$rawdata, file_pattern))
  }
  
  cat(sprintf("Found %d file(s):\n", length(myfiles)))
  for (i in seq_along(myfiles)) {
    cat(sprintf("  [%d] %s\n", i, myfiles[i]))
  }
  
  # Read the first file (extend logic if multiple files are needed)
  detect_delim <- function(file_path, candidates = c("\t", ",", ";")) {
    lines_preview <- readr::read_lines(file_path, n_max = 20)
    if (length(lines_preview) == 0) return("\t")
    counts <- vapply(
      candidates,
      function(d) sum(lengths(regmatches(lines_preview, gregexpr(d, lines_preview, fixed = TRUE)))),
      numeric(1)
    )
    best_idx <- which.max(counts)
    if (all(counts == 0)) return("\t")
    candidates[best_idx]
  }
  
  target_file <- myfiles[1]
  delim <- detect_delim(target_file)
  cat(sprintf("Detected delimiter automatically: %s\n", ifelse(delim == "\t", "\\t", delim)))
  
  data_raw <- readr::read_delim(
    target_file,
    delim = delim,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  if (ncol(data_raw) == 1 && delim != ",") {
    cat("⚠ Only one column detected, retrying with comma delimiter...\n")
    data_raw <- readr::read_csv(target_file, show_col_types = FALSE, trim_ws = TRUE)
  }
  
  cat(sprintf("\n✓ Read file: %s\n", target_file))
  cat(sprintf("  Dimensions: %d rows × %d columns\n", nrow(data_raw), ncol(data_raw)))
  
  # --------------------------------------------------------------------------
  # Basic data cleaning
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Basic data cleaning\n")
  cat("----------------------------------------\n")
  
  # Standardize column names (trim leading/trailing spaces)
  colnames(data_raw) <- trimws(colnames(data_raw))
  if (colnames(data_raw)[1] != "Gene") {
    cat(sprintf("Renamed first column from '%s' to 'Gene'\n", colnames(data_raw)[1]))
    colnames(data_raw)[1] <- "Gene"
  }
  
  data_cols <- setdiff(colnames(data_raw), "Gene")
  
  # Convert data columns to numeric and drop NaN
  if (length(data_cols) > 0) {
    data_raw <- data_raw %>%
      mutate(across(all_of(data_cols), ~ {
        # Keep numeric values; trim/convert other types
        if (is.numeric(.)) return(.)
        suppressWarnings(as.numeric(trimws(as.character(.))))
      })) %>%
      mutate(across(all_of(data_cols), ~ na_if(., NaN)))
  }
  
  # Replace non-NA values < 1 with 1 (prevent negative log2)
  if (length(data_cols) > 0) {
    # Count how many values will be replaced
    n_replaced <- 0
    for (col in data_cols) {
      col_data <- data_raw[[col]]
      mask <- !is.na(col_data) & col_data < 1
      n_replaced <- n_replaced + sum(mask)
    }
    
    if (n_replaced > 0) {
      cat(sprintf("✓ Replaced %d non-NA values < 1 with 1 (avoid negative log2)\n", n_replaced))
      data_raw <- data_raw %>%
        mutate(across(all_of(data_cols), ~ {
          ifelse(!is.na(.) & . < 1, 1, .)
        }))
    } else {
      cat("✓ Data check: no replacement needed (all non-NA values >= 1)\n")
    }
  }
  
  # Clean Gene column (blanks, semicolons, duplicates)
  data_raw <- data_raw %>%
    mutate(Gene = trimws(as.character(Gene))) %>%
    filter(!is.na(Gene) & Gene != "") %>%
    filter(!grepl("[;；]", Gene))
  
  dup_removed <- nrow(data_raw) - nrow(data_raw %>% distinct(Gene, .keep_all = TRUE))
  data_raw <- data_raw %>% distinct(Gene, .keep_all = TRUE)
  if (dup_removed > 0) {
    cat(sprintf("✓ Removed %d rows with duplicate Gene\n", dup_removed))
  }
  
  # Remove rows where all data columns are NA
  if (length(data_cols) > 0 && nrow(data_raw) > 0) {
    all_na_mask <- apply(select(data_raw, all_of(data_cols)), 1, function(x) all(is.na(x)))
    removed_all_na <- sum(all_na_mask)
    if (removed_all_na > 0) {
      cat(sprintf("✓ Removed %d rows where all data columns are NA\n", removed_all_na))
      data_raw <- data_raw[!all_na_mask, , drop = FALSE]
    }
  }
  
  cat(sprintf("  Cleaned dimensions: %d rows × %d columns\n", nrow(data_raw), ncol(data_raw)))
  
  # --------------------------------------------------------------------------
  # Step 2: Clean column names automatically
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 2: Clean column names\n")
  cat("----------------------------------------\n")
  
  original_colnames <- colnames(data_raw)
  cat("Original column names:\n")
  print(head(original_colnames, 10))
  
  # First column defaults to Gene
  cleaned_colnames <- original_colnames
  cleaned_colnames[1] <- "Gene"
  
  # Detect and remove common prefix/suffix (excluding Gene)
  if (ncol(data_raw) > 1) {
    sample_cols <- original_colnames[-1]
    
    # Find common prefix
    common_prefix <- ""
    if (length(sample_cols) > 1) {
      # Use first and last strings to find the common prefix
      sorted_cols <- sort(sample_cols)
      first <- sorted_cols[1]
      last <- sorted_cols[length(sorted_cols)]
      
      for (i in 1:min(nchar(first), nchar(last))) {
        if (substr(first, 1, i) == substr(last, 1, i)) {
          common_prefix <- substr(first, 1, i)
        } else {
          break
        }
      }
    }
    
    # Remove common prefix if meaningful
    if (nchar(common_prefix) > 0 && all(startsWith(sample_cols, common_prefix))) {
      cleaned_sample_cols <- substr(sample_cols, nchar(common_prefix) + 1, nchar(sample_cols))
      # Ensure cleaned names still have content
      if (all(nchar(cleaned_sample_cols) > 0)) {
        cleaned_colnames[-1] <- cleaned_sample_cols
        cat(sprintf("\n✓ Removed common prefix: '%s'\n", common_prefix))
      }
    }
  }
  
  colnames(data_raw) <- cleaned_colnames
  cat("\nCleaned column names:\n")
  print(head(colnames(data_raw), 10))
  
  # --------------------------------------------------------------------------
  # Step 3: Generate sampleGroup template
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 3: Generate sampleGroup template\n")
  cat("----------------------------------------\n")
  
  # Get sample column names (excluding Gene)
  sample_names <- colnames(data_raw)[-1]
  n_samples <- length(sample_names)
  
  # Generate example values (ensure each option appears at least once)
  generate_example_values <- function(n, options) {
    n_opts <- length(options)
    if (n <= n_opts) {
      return(options[1:n])
    } else {
      # Spread options evenly
      repeats <- rep(ceiling(n / n_opts), n_opts)
      values <- rep(options, repeats)
      return(values[1:n])
    }
  }
  
  # Create template data frame
  template <- data.frame(
    OriginalName = sample_names,
    bioGroup = "",
    CatalyticGroup = generate_example_values(n_samples, c("Cata", "NoCat")),
    PLtype = generate_example_values(n_samples, c("Light", "H2O2", "PL")),
    Context = generate_example_values(n_samples, c("Experiment", "Control", "Spatial")),
    replicate = rep(1:3, length.out = n_samples),
    FirstROCgroup = generate_example_values(n_samples, c("A", "B", "C", "NA")),
    SecondROCgroup = generate_example_values(n_samples, c("A", "B", "C", "D", "E", "F")),
    Order = 1:n_samples,  # default order
    FinalName = "",
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("✓ Generated template: %d samples\n", nrow(template)))
  
  # Save template to working directory
  setwd(dir_config$root)
  write.csv(template, "Module02_sampleGroup_template.csv", row.names = FALSE)
  cat("\n✓ Template saved: Module02_sampleGroup_template.csv\n")
  cat("\n📝 Please fill in these columns:\n")
  cat("  - bioGroup: biological group name\n")
  cat("  - CatalyticGroup: Cata or NoCat\n")
  cat("  - PLtype: Light or H2O2 or PL\n")
  cat("  - Context: Experiment or Control or Spatial\n")
  cat("  - replicate: 1 or 2 or 3\n")
  cat("  - FirstROCgroup: A / B / C / A&B / NA\n")
  cat("  - SecondROCgroup: A / B / C / D / E / F\n")
  cat("  - Order: integer (1/2/3...) controls display order of data columns, 1 appears first\n")
  cat("  (FinalName is auto-generated; leave blank)\n")
  
  return(list(data = data_raw, template = template))
}


#' Read the completed sampleGroup_template and build sampleGroup
#' 
#' @param dir_config Directory configuration
#' @param data_raw Raw data
#' @return List containing data and sampleGroup
module02_process_samplegroup <- function(dir_config, data_raw) {
  
  cat("\n----------------------------------------\n")
  cat("Step 4: Read user-completed sampleGroup_template\n")
  cat("----------------------------------------\n")
  
  setwd(dir_config$root)
  
  # Check that the file exists
  template_file <- "Module02_sampleGroup_template.csv"
  if (!file.exists(template_file)) {
    stop(sprintf("✗ Error: file %s not found\nPlease run Step 1 to generate the template first", template_file))
  }
  
  # Read user-completed template
  sampleGroup <- read.csv(template_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ Loaded: %s\n", template_file))
  cat(sprintf("  Rows: %d\n", nrow(sampleGroup)))
  
  # Validate required columns
  required_cols <- c("OriginalName", "bioGroup", "CatalyticGroup", "PLtype", 
                     "Context", "replicate", "FirstROCgroup", "SecondROCgroup", "Order")
  missing_cols <- setdiff(required_cols, colnames(sampleGroup))
  if (length(missing_cols) > 0) {
    stop(sprintf("✗ Error: missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure values are filled
  empty_bioGroup <- sampleGroup$bioGroup == "" | is.na(sampleGroup$bioGroup)
  if (any(empty_bioGroup)) {
    stop(sprintf("✗ Error: bioGroup not filled on row(s) %s", 
                 paste(which(empty_bioGroup), collapse = ", ")))
  }
  
  # Validate Order column
  if (any(is.na(sampleGroup$Order))) {
    stop("✗ Error: Order column has empty values; please fill with integers")
  }
  if (!all(sampleGroup$Order == as.integer(sampleGroup$Order))) {
    stop("✗ Error: Order column must be integers")
  }
  
  # --------------------------------------------------------------------------
  # Step 5: Reorder sampleGroup and data columns by Order
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 5: Reorder data columns by Order\n")
  cat("----------------------------------------\n")
  
  # Sort sampleGroup by Order
  sampleGroup <- sampleGroup %>% arrange(Order)
  cat("✓ sampleGroup sorted by Order\n")
  
  # Get ordered column names (keep Gene first)
  ordered_sample_names <- sampleGroup$OriginalName
  
  # Verify all OriginalName entries exist in data_raw
  missing_samples <- setdiff(ordered_sample_names, colnames(data_raw)[-1])
  if (length(missing_samples) > 0) {
    stop(sprintf("✗ Error: samples in sampleGroup not found in data: %s", 
                 paste(missing_samples, collapse = ", ")))
  }
  
  # Reorder data_raw columns: Gene + ordered samples
  data_raw <- data_raw %>% select(Gene, all_of(ordered_sample_names))
  cat("✓ Data columns reordered by Order\n")
  cat("\nNew column order:\n")
  print(head(colnames(data_raw), 10))
  
  # --------------------------------------------------------------------------
  # Step 6: Generate FinalName
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 6: Generate FinalName\n")
  cat("----------------------------------------\n")
  
  # FinalName = bioGroup + _LFQ_ + replicate
  sampleGroup$FinalName <- paste0(sampleGroup$bioGroup, "_LFQ_", sampleGroup$replicate)
  
  cat("✓ FinalName created\n")
  cat("\nExample:\n")
  print(head(sampleGroup[, c("OriginalName", "bioGroup", "replicate", "FinalName", "Order")], 5))
  
  # --------------------------------------------------------------------------
  # Step 7: Rename data columns
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 7: Rename data columns\n")
  cat("----------------------------------------\n")
  
  # Create column name mapping
  name_mapping <- setNames(sampleGroup$FinalName, sampleGroup$OriginalName)
  
  # Rename (keep Gene column unchanged)
  current_colnames <- colnames(data_raw)
  new_colnames <- current_colnames
  for (i in 2:length(current_colnames)) {
    old_name <- current_colnames[i]
    if (old_name %in% names(name_mapping)) {
      new_colnames[i] <- name_mapping[old_name]
    }
  }
  
  colnames(data_raw) <- new_colnames
  cat("✓ Column renaming complete\n")
  cat("\nNew column names:\n")
  print(head(colnames(data_raw), 10))
  
  # --------------------------------------------------------------------------
  # Step 8: Save data
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 8: Save data\n")
  cat("----------------------------------------\n")
  
  # Save final sampleGroup to working directory (CSV)
  write.csv(sampleGroup, "Module02_sampleGroup.csv", row.names = FALSE)
  cat("✓ Saved: Module02_sampleGroup.csv (working directory, final version)\n")
  cat("  This file contains validated group info and generated FinalName\n")
  
  # Save CSV to Output directory
  output_file <- file.path(dir_config$output, "Module02_data_raw.csv")
  write.csv(data_raw, output_file, row.names = FALSE)
  cat(sprintf("✓ Saved: %s\n", output_file))
  
  cat("\n✓ Step 8 complete; returning data to main pipeline\n")
  
  return(list(data = data_raw, sampleGroup = sampleGroup))
}


#' Module 2 main function
#' 
#' @param dir_config Directory configuration
#' @param file_pattern TSV file pattern
#' @param auto_process Whether to auto process (FALSE generates template only)
module02_data_import <- function(dir_config, file_pattern = "_matrix.*\\.tsv$", auto_process = FALSE) {
  
  cat("\n========================================\n")
  cat("Module 2: Data Import and Sample Table\n")
  cat("========================================\n")
  
  # Steps 1-3: read data and generate template
  result <- module02_read_and_generate_template(dir_config, file_pattern)
  
  if (auto_process) {
    # If auto processing, try reading completed sampleGroup
    result <- module02_process_samplegroup(dir_config, result$data)
  } else {
    cat("\n⚠ Please complete the following steps:\n")
    cat("  1. Open sampleGroup_template.csv\n")
    cat("  2. Fill all required columns\n")
    cat("  3. Save as sampleGroup.csv\n")
    cat("  4. Continue running in main_pipeline.R\n")
  }
  
  return(invisible(result))
}

