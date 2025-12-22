# ============================================================================
# Module 11: Second differential analysis (based on Module 10 replacement results)
# ============================================================================
# Aligns with CleanCode.R (2377-2491) logic
#
# Goal:
# - For Module 10 generated replaced_data_A / replaced_data_B,
#   only keep samples with Context of Experiment or Spatial (remove Control)
# - Build comparison groups based on SecondROCgroup and PLtype:
#   * Exp vs Exp: compare all Experiment groups (PLtype is not restricted)
#   * Exp vs Spatial: only compare groups with the same PLtype
# - Use limma for differential analysis to obtain logFC and adj.P.Val
# - Produce FDR_combined_df_list_2nd (combined results of all comparisons)
#
# Input:
# - dir_config: contains paths such as output
# - sampleGroup: identifies sample grouping, Context, PLtype, and SecondROCgroup
# - replaced_data_A / replaced_data_B: data replacement outputs from Module 10
# - selected_versions: versions to process (NULL means all)
#
# Output:
# - FDR_combined_df_list_2nd: differential analysis results for all versions/methods
# - An Excel file per dataset (original topTable results)
# - A combined Excel file (all versions' logFC + adj.P.Val + annotation columns)
#
# ============================================================================

module11_diff_analysis2 <- function(
  dir_config,
  sampleGroup,
  replaced_data_A,
  replaced_data_B,
  selected_versions = NULL
) {
  
  cat("\n=== Module 11: Second differential analysis ===\n\n")
  
  # Load required packages
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("✗ Error: limma package is required\n  BiocManager::install('limma')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("✗ Error: dplyr package is required")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("✗ Error: openxlsx package is required")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("✗ Error: purrr package is required")
  }
  
  library(limma)
  library(dplyr)
  library(openxlsx)
  library(purrr)
  
  # ----------------------------------------------------------------------
  # 1. Prepare data list
  # ----------------------------------------------------------------------
  
  # Add suffix for method A and B data
  names(replaced_data_A) <- paste0(names(replaced_data_A), "_A")
  names(replaced_data_B) <- paste0(names(replaced_data_B), "_B")
  
  # Combine method A and B data
  all_data <- c(replaced_data_A, replaced_data_B)
  
  # Filter by selected_versions if provided
  if (!is.null(selected_versions)) {
    selected_names <- c(paste0(selected_versions, "_A"),
                       paste0(selected_versions, "_B"))
    all_data <- all_data[names(all_data) %in% selected_names]
  }
  
  if (length(all_data) == 0) {
    stop("✗ Error: no available data versions")
  }
  
  cat(sprintf("✓ Preparing to analyze %d data versions\n", length(all_data)))
  cat(sprintf("  Versions: %s\n", paste(names(all_data), collapse = ", ")))
  
  # ----------------------------------------------------------------------
  # 2. Extract sample grouping info from sampleGroup
  # ----------------------------------------------------------------------
  
  cat("\n--- Step 1: Extract sample grouping info ---\n")
  
  # Keep only samples with Context of Experiment or Spatial
  samples_exp_spatial <- sampleGroup %>%
    filter(Context %in% c("Experiment", "Spatial")) %>%
    arrange(Order)
  
  if (nrow(samples_exp_spatial) == 0) {
    stop("✗ Error: no Experiment or Spatial samples in sampleGroup")
  }
  
  cat(sprintf("✓ Found %d Experiment/Spatial samples\n", nrow(samples_exp_spatial)))
  
  # Extract unique bioGroup ordered by Order to keep consistent order
  bioGroup_info <- samples_exp_spatial %>%
    group_by(bioGroup) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(Order) %>%
    select(bioGroup, Context, PLtype, SecondROCgroup)
  
  bioGroups_selected <- bioGroup_info$bioGroup
  
  cat(sprintf("✓ Selected bioGroup (%d):\n", length(bioGroups_selected)))
  for (i in seq_along(bioGroups_selected)) {
    bg <- bioGroups_selected[i]
    info <- bioGroup_info %>% filter(bioGroup == bg)
    cat(sprintf("  %d. %s (Context: %s, PLtype: %s, SecondROCgroup: %s)\n",
                i, bg, info$Context, info$PLtype,
                ifelse(is.na(info$SecondROCgroup) | info$SecondROCgroup == "", 
                       "Empty", info$SecondROCgroup)))
  }
  
  # Build sample_info (for limma design matrix)
  sample_info <- data.frame(
    SampleName = samples_exp_spatial$FinalName,
    Group = samples_exp_spatial$bioGroup,
    stringsAsFactors = FALSE
  )
  
  # Set Group as factor following bioGroups_selected order
  sample_info$Group <- factor(sample_info$Group, levels = bioGroups_selected)
  
  cat(sprintf("✓ Built sample_info with %d samples\n", nrow(sample_info)))
  
  # ----------------------------------------------------------------------
  # 3. Build comparison matrix (based on Context and PLtype)
  # ----------------------------------------------------------------------
  
  cat("\n--- Step 2: Build comparison matrix ---\n")
  
  # Extract Experiment and Spatial bioGroups separately
  exp_groups <- bioGroup_info %>% filter(Context == "Experiment")
  spatial_groups <- bioGroup_info %>% filter(Context == "Spatial")
  
  cat(sprintf("✓ Experiment bioGroups (%d): %s\n",
              nrow(exp_groups),
              paste(exp_groups$bioGroup, collapse = ", ")))
  cat(sprintf("✓ Spatial bioGroups (%d): %s\n",
              nrow(spatial_groups),
              paste(spatial_groups$bioGroup, collapse = ", ")))
  
  # Build comparison list
  comparisons_list <- list()
  comparison_names <- c()
  
  # 1. Exp vs Exp comparisons (all Experiment groups, PLtype not restricted)
  if (nrow(exp_groups) >= 2) {
    for (i in 1:(nrow(exp_groups) - 1)) {
      for (j in (i + 1):nrow(exp_groups)) {
        exp1 <- exp_groups$bioGroup[i]
        exp2 <- exp_groups$bioGroup[j]
        comparisons_list[[length(comparisons_list) + 1]] <- c(exp1, exp2)
        # Simplify naming (remove _Light/_H2O2 suffix)
        name1 <- gsub("_Light$|_H2O2$", "", exp1)
        name2 <- gsub("_Light$|_H2O2$", "", exp2)
        comparison_names <- c(comparison_names, paste0(name1, "_vs_", name2))
      }
    }
  }
  
  # 2. Exp vs Spatial comparisons (require matching PLtype)
  for (i in 1:nrow(exp_groups)) {
    exp_bg <- exp_groups$bioGroup[i]
    exp_pltype <- exp_groups$PLtype[i]
    
    for (j in 1:nrow(spatial_groups)) {
      spatial_bg <- spatial_groups$bioGroup[j]
      spatial_pltype <- spatial_groups$PLtype[j]
      
      # Only compare when PLtype matches
      if (exp_pltype == spatial_pltype) {
        comparisons_list[[length(comparisons_list) + 1]] <- c(exp_bg, spatial_bg)
        # Simplify naming
        name1 <- gsub("_Light$|_H2O2$", "", exp_bg)
        name2 <- gsub("_Light$|_H2O2$", "", spatial_bg)
        comparison_names <- c(comparison_names, paste0(name1, "_vs_", name2))
      }
    }
  }
  
  names(comparisons_list) <- comparison_names
  
  cat(sprintf("\n✓ Built %d comparison groups:\n", length(comparisons_list)))
  for (i in seq_along(comparisons_list)) {
    cat(sprintf("  %d. %s (%s vs %s)\n", i, names(comparisons_list)[i],
                comparisons_list[[i]][1], comparisons_list[[i]][2]))
  }
  
  # ----------------------------------------------------------------------
  # 4. Run limma differential analysis for each data version
  # ----------------------------------------------------------------------
  
  cat("\n--- Step 3: Run differential analysis ---\n")
  
  FDR_combined_df_list_2nd <- list()
  
  for (j in seq_along(all_data)) {
    mydata <- all_data[[j]]
    myfilename <- names(all_data)[j]
    
    cat(sprintf("\n[%d/%d] Processing: %s\n", j, length(all_data), myfilename))
    cat(sprintf("------%s_start------\n", myfilename))
    
    # Extract data matrix (remove last 3 annotation columns)
    myrange <- ncol(mydata) - 3
    
    # Check data structure
    if (myrange < 2) {
      cat("  ⚠ Warning: insufficient data columns, skipping this version\n")
      next
    }
    
    # Extract expression matrix (columns 2 to myrange, column 1 is Gene)
    data_cols <- colnames(mydata)[2:myrange]
    
    # Check whether all required samples are present
    missing_cols <- setdiff(sample_info$SampleName, data_cols)
    if (length(missing_cols) > 0) {
      cat(sprintf("  ⚠ Warning: missing %d sample columns, skipping this version\n", length(missing_cols)))
      cat(sprintf("    Missing: %s\n", paste(head(missing_cols, 3), collapse = ", ")))
      next
    }
    
    # Extract expression matrix following sample_info order
    proteomics_data <- mydata %>% 
      select(all_of(sample_info$SampleName)) %>%
      as.matrix()
    
    rownames(proteomics_data) <- mydata$Gene
    
    cat(sprintf("  ✓ Extracted data matrix: %d genes x %d samples\n", 
                nrow(proteomics_data), ncol(proteomics_data)))
    
    # Create design matrix
    design <- model.matrix(~ 0 + Group, data = sample_info)
    colnames(design) <- levels(sample_info$Group)
    
    cat("  ✓ Created design matrix\n")
    
    # Fit linear model
    fit <- lmFit(proteomics_data, design)
    
    cat("  ✓ Fitted linear model\n")
    
    # Build contrast matrix expressions
    contrast_expressions <- character(length(comparisons_list))
    for (i in seq_along(comparisons_list)) {
      comp <- comparisons_list[[i]]
      comp_name <- names(comparisons_list)[i]
      # bioGroup names may contain special characters; wrap with backticks
      contrast_expressions[i] <- sprintf("%s = `%s` - `%s`", 
                                        comp_name, 
                                        comp[1], 
                                        comp[2])
    }
    
    # Dynamically build contrast matrix using makeContrasts
    contrast_matrix <- eval(parse(text = sprintf(
      "makeContrasts(%s, levels = design)",
      paste(contrast_expressions, collapse = ", ")
    )))
    
    cat("  ✓ Created contrast matrix\n")
    
    # Calculate contrasts and perform empirical Bayes adjustment
    fit_contrasts <- contrasts.fit(fit, contrast_matrix)
    fit_ebayes <- eBayes(fit_contrasts)
    
    cat("  ✓ Empirical Bayes adjustment completed\n")
    
    # Extract results for each comparison
    FDR_test_list <- list()
    Raw_FDR_test_list <- list()
    Comparision_Group <- colnames(contrast_matrix)
    
    for (i in seq_along(Comparision_Group)) {
      myGroup <- Comparision_Group[i]
      
      # Extract topTable results
      tem_dataframe <- topTable(fit_ebayes, coef = myGroup, n = Inf, adjust.method = "BH")
      tem_dataframe$Gene <- rownames(tem_dataframe)
      
      # Save raw results
      Raw_tem_dataframe <- tem_dataframe
      Raw_FDR_test_list[[paste0("RES_", myGroup)]] <- Raw_tem_dataframe
      
      # Keep only Gene, logFC, adj.P.Val
      tem_dataframe_simple <- tem_dataframe %>% 
        select(Gene, logFC, adj.P.Val)
      colnames(tem_dataframe_simple) <- c("Gene",
                                         paste0(myGroup, "_logFC"),
                                         paste0(myGroup, "_adj.P.Val"))
      
      FDR_test_list[[paste0("RES_", myGroup)]] <- tem_dataframe_simple
    }
    
    cat(sprintf("  ✓ Extracted differential analysis results for %d comparisons\n", length(FDR_test_list)))
    
    # Save raw topTable results (one sheet per comparison)
    wb <- createWorkbook()
    for (i in seq_along(Raw_FDR_test_list)) {
      sheet_name <- names(Raw_FDR_test_list)[i]
      # Excel sheet names are limited to 31 characters
      if (nchar(sheet_name) > 31) {
        sheet_name <- substr(sheet_name, 1, 31)
      }
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, Raw_FDR_test_list[[i]])
    }
    raw_file <- file.path(dir_config$output, 
                         paste0("Module11_Raw_FDR_test_list_", myfilename, ".xlsx"))
    saveWorkbook(wb, raw_file, overwrite = TRUE)
    cat(sprintf("  ✓ Exported raw topTable: %s\n", basename(raw_file)))
    
    # Merge logFC and adj.P.Val for all comparisons
    FDR_combined_df <- reduce(FDR_test_list, full_join, by = "Gene")
    
    # Add annotation columns (prefer *_Localization columns)
    preferred_annotation_cols <- c(
      "HaloMap_Localization",
      "GO_Localization",
      "MultiBait_Localization"
    )
    annotation_cols <- preferred_annotation_cols[
      preferred_annotation_cols %in% colnames(mydata)
    ]
    if (length(annotation_cols) == 0) {
      annotation_cols <- grep("_Localization$", colnames(mydata), value = TRUE)
    }
    if (length(annotation_cols) == 0) {
      stop(sprintf(
        "✗ Error: no *_Localization annotation columns found in dataset %s; please check input data.",
        myfilename
      ))
    }
    
    annotation_df <- mydata %>% select(Gene, all_of(annotation_cols))
    FDR_combined_df <- FDR_combined_df %>%
      left_join(annotation_df, by = "Gene")
    
    # Store in list
    FDR_combined_df_list_2nd[[myfilename]] <- FDR_combined_df
    
    cat(sprintf("  ✓ Merged results: %d genes x %d columns\n", 
                nrow(FDR_combined_df), ncol(FDR_combined_df)))
    cat(sprintf("------%s_finished------\n", myfilename))
  }
  
  # ----------------------------------------------------------------------
  # 5. Export combined differential analysis results
  # ----------------------------------------------------------------------
  
  cat("\n--- Step 4: Export combined results ---\n")
  
  if (length(FDR_combined_df_list_2nd) == 0) {
    cat("⚠ Warning: no results to export\n")
    return(list(
      FDR_combined_df_list_2nd = list(),
      sample_info = sample_info,
      comparisons_list = comparisons_list,
      bioGroups_selected = bioGroups_selected,
      bioGroup_info = bioGroup_info
    ))
  }
  
  wb <- createWorkbook()
  for (i in seq_along(FDR_combined_df_list_2nd)) {
    sheet_name <- names(FDR_combined_df_list_2nd)[i]
    # Excel sheet names are limited to 31 characters
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
    }
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, FDR_combined_df_list_2nd[[i]])
  }
  
  combined_file <- file.path(dir_config$output, 
                            "Module11_FC_FDR_2nd_combined.xlsx")
  saveWorkbook(wb, combined_file, overwrite = TRUE)
  
  cat(sprintf("✓ Exported combined differential analysis results:\n"))
  cat(sprintf("  %s\n", basename(combined_file)))
  cat(sprintf("  Contains %d data versions\n", length(FDR_combined_df_list_2nd)))
  
  # ----------------------------------------------------------------------
  # 6. Summary
  # ----------------------------------------------------------------------
  
  cat("\n=== Module 11 completed ===\n")
  cat(sprintf("✓ Analyzed %d data versions\n", length(FDR_combined_df_list_2nd)))
  cat(sprintf("✓ Built %d comparison groups\n", length(comparisons_list)))
  cat(sprintf("✓ Output file count: %d raw topTable + 1 combined result\n", 
              length(FDR_combined_df_list_2nd)))
  
  # Return results
  return(list(
    FDR_combined_df_list_2nd = FDR_combined_df_list_2nd,
    sample_info = sample_info,
    comparisons_list = comparisons_list,
    bioGroups_selected = bioGroups_selected,
    bioGroup_info = bioGroup_info
  ))
}
