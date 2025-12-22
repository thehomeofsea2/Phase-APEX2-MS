# Module 06: Heatmap System
# Purpose: Plot multiple heatmap types for imputed data, including correlation and localization-based heatmaps
# Author: CodeNorm Pipeline
# Date: 2024

#' Module 06: Heatmap System
#' 
#' @description
#' Plot heatmaps for imputed data, including:
#' 1. Correlation heatmap (with detailed customization)
#' 2. AllLocalization heatmap (all proteins)
#' 3. Heatmaps by Localization (SGs, Mitochondrion, Cytosol, Nuclear, Nuclear_Cytosol, Other)
#' 
#' @param dir_config Directory configuration list (from Module 1)
#' @param imputed_data_list Imputed data list (from Module 5)
#' @param sampleGroup Sample grouping information (from Module 2)
#' @param selected_versions Normalization version names to plot (vector), e.g., c("noNorm_Imputed", "Local_QNorm_Imputed")
#' @param heatmap_types Heatmap types to generate. Options: c("correlation", "all", "by_localization"). Default: c("all", "by_localization")
#' @param correlation_config Correlation heatmap config list (NULL uses defaults), including:
#'   - corr_min: minimum correlation (default 0.8)
#'   - corr_max: maximum correlation (default 1)
#'   - corr_center: center value (default 0.9)
#'   - col_low: low color (default "#4D97CD")
#'   - col_center: center color (default "white")
#'   - col_high: high color (default "#DB6968")
#'   - exclude_context: Context values to exclude (default NULL)
#'   - exclude_pltype: PLtype values to exclude (default NULL)
#'   - exclude_catalytic: CatalyticGroup values to exclude (default NULL)
#'   - include_biogroups: bioGroups to include explicitly (default NULL uses exclusions)
#' @param localization_columns Annotation column names for localization (auto-detected by default)
#' @param color_params Heatmap color parameters (for AllLocalization and by_localization heatmaps)
#' 
#' @return List containing:
#'   - heatmap_params: heatmap parameter info
#'   - annotation_info: annotation info
#' 
#' @export
module06_heatmap <- function(dir_config,
                             imputed_data_list,
                             sampleGroup,
                             selected_versions = NULL,
                             heatmap_types = c("all", "by_localization"),
                             correlation_config = NULL,
                             localization_columns = NULL,
                             color_params = NULL) {
  
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(scales)
  
  # 1. Input validation ####
  cat("\n[1] Validate input data...\n")
  
  if (!all(c("output") %in% names(dir_config))) {
    stop("❌ dir_config must include output path")
  }
  
  if (length(imputed_data_list) == 0) {
    stop("❌ imputed_data_list is empty")
  }
  
  required_cols <- c("FinalName", "bioGroup", "PLtype", "Context", "CatalyticGroup")
  if (!is.data.frame(sampleGroup) || !all(required_cols %in% names(sampleGroup))) {
    stop(sprintf("❌ sampleGroup must include: %s", paste(required_cols, collapse = ", ")))
  }
  
  # Auto-select versions if not provided
  if (is.null(selected_versions)) {
    if ("noNorm_Imputed" %in% names(imputed_data_list)) {
      selected_versions <- "noNorm_Imputed"
    } else {
      selected_versions <- names(imputed_data_list)[1]
    }
  }
  
  # Validate selected versions exist
  missing_versions <- setdiff(selected_versions, names(imputed_data_list))
  if (length(missing_versions) > 0) {
    stop(sprintf("❌ Versions not found: %s", paste(missing_versions, collapse = ", ")))
  }
  
  cat("✓ Input validation passed\n")
  cat(sprintf("  - Available datasets: %d\n", length(imputed_data_list)))
  cat(sprintf("  - Selected versions: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("  - Heatmap types: %s\n", paste(heatmap_types, collapse = ", ")))
  
  # 2. Set default correlation config ####
  if (is.null(correlation_config)) {
    correlation_config <- list(
      corr_min = 0.8,
      corr_max = 1,
      corr_center = 0.9,
      col_low = "#4D97CD",
      col_center = "white",
      col_high = "#DB6968",
      exclude_context = NULL,
      exclude_pltype = NULL,
      exclude_catalytic = NULL,
      include_biogroups = NULL
    )
  } else {
    # Fill missing defaults
    if (is.null(correlation_config$corr_min)) correlation_config$corr_min <- 0.8
    if (is.null(correlation_config$corr_max)) correlation_config$corr_max <- 1
    if (is.null(correlation_config$corr_center)) correlation_config$corr_center <- 0.9
    if (is.null(correlation_config$col_low)) correlation_config$col_low <- "#4D97CD"
    if (is.null(correlation_config$col_center)) correlation_config$col_center <- "white"
    if (is.null(correlation_config$col_high)) correlation_config$col_high <- "#DB6968"
  }
  
  # 3. generate_pheatmap_params helper ####
  generate_pheatmap_params <- function(mat, 
                                       custom_min = NULL, 
                                       custom_max = NULL, 
                                       custom_center = NULL, 
                                       col_low = "blue", 
                                       col_high = "red", 
                                       col_center = "white") {
    
    # Step 1: decide final data range
    data_min <- if (is.null(custom_min)) min(mat, na.rm = TRUE) else custom_min
    data_max <- if (is.null(custom_max)) max(mat, na.rm = TRUE) else custom_max
    
    # Step 2: decide color center
    center_val <- if (is.null(custom_center)) (data_min + data_max) / 2 else custom_center
    
    # Step 3: generate breaks
    n_breaks <- 100
    
    # If center is within range, build in two segments
    if (center_val > data_min && center_val < data_max) {
      # Left: data_min -> center_val
      breaks_left <- seq(data_min, center_val, length.out = n_breaks / 2)
      # Right: center_val -> data_max
      breaks_right <- seq(center_val, data_max, length.out = n_breaks / 2)
      # Merge and deduplicate
      my_breaks <- unique(c(breaks_left, breaks_right))
      
      # Corresponding color gradients
      colors_left <- colorRampPalette(c(col_low, col_center))(n_breaks / 2)
      colors_right <- colorRampPalette(c(col_center, col_high))(n_breaks / 2)
      my_colors <- c(colors_left, colors_right)
    } else {
      # Center outside range: direct span
      my_breaks <- seq(data_min, data_max, length.out = n_breaks)
      my_colors <- colorRampPalette(c(col_low, col_high))(n_breaks)
    }
    
    return(list(breaks = my_breaks, color = my_colors))
  }
  
  # 4. Filter samples for correlation heatmap ####
  if ("correlation" %in% heatmap_types) {
    cat("\n[2] Filter samples for correlation heatmap...\n")
    
    # If include_biogroups provided, use directly
    if (!is.null(correlation_config$include_biogroups)) {
      filtered_sampleGroup <- sampleGroup %>%
        filter(bioGroup %in% correlation_config$include_biogroups)
      
      cat(sprintf("  - Using specified bioGroup: %s\n", 
                  paste(correlation_config$include_biogroups, collapse = ", ")))
    } else {
      # Otherwise apply exclusion rules
      filtered_sampleGroup <- sampleGroup
      
      if (!is.null(correlation_config$exclude_context)) {
        filtered_sampleGroup <- filtered_sampleGroup %>%
          filter(!Context %in% correlation_config$exclude_context)
        cat(sprintf("  - Excluding Context: %s\n", 
                    paste(correlation_config$exclude_context, collapse = ", ")))
      }
      
      if (!is.null(correlation_config$exclude_pltype)) {
        filtered_sampleGroup <- filtered_sampleGroup %>%
          filter(!PLtype %in% correlation_config$exclude_pltype)
        cat(sprintf("  - Excluding PLtype: %s\n", 
                    paste(correlation_config$exclude_pltype, collapse = ", ")))
      }
      
      if (!is.null(correlation_config$exclude_catalytic)) {
        filtered_sampleGroup <- filtered_sampleGroup %>%
          filter(!CatalyticGroup %in% correlation_config$exclude_catalytic)
        cat(sprintf("  - Excluding CatalyticGroup: %s\n", 
                    paste(correlation_config$exclude_catalytic, collapse = ", ")))
      }
    }
    
    corr_samples <- filtered_sampleGroup$FinalName
    
    if (length(corr_samples) == 0) {
      warning("⚠ No samples left for correlation heatmap after filtering; skipping correlation heatmap")
      heatmap_types <- setdiff(heatmap_types, "correlation")
    } else {
      cat(sprintf("✓ Filtering complete, kept %d samples\n", length(corr_samples)))
      
      # 5. Prepare correlation heatmap annotations/colors ####
      cat("\n[3] Prepare correlation heatmap column annotations...\n")
      
      # Generate color maps for PLtype and Context
      unique_pltypes <- unique(filtered_sampleGroup$PLtype)
      unique_contexts <- unique(filtered_sampleGroup$Context)
      
      # PLtype colors: one per observed PLtype
      pltype_colors <- setNames(
        hue_pal()(length(unique_pltypes)),
        unique_pltypes
      )
      
      # Context colors: one per observed Context
      context_colors <- setNames(
        c("#CC0000", "#FFA07A")[1:length(unique_contexts)],
        unique_contexts
      )
      
      # Column annotation data frame (factor levels aligned with colors)
      Annotated_col_corr <- data.frame(
        PLtype = factor(filtered_sampleGroup$PLtype, levels = unique_pltypes),
        Context = factor(filtered_sampleGroup$Context, levels = unique_contexts),
        stringsAsFactors = FALSE
      )
      rownames(Annotated_col_corr) <- filtered_sampleGroup$FinalName
      
      # Column color rules
      # PLtype: Light=blues, H2O2=oranges, PL=greens
      # Context: Experiment=dark, Control=light
      # Spatial=black (highest priority)
      
      # Assign combined color per sample (PLtype × Context)
      sample_colors <- character(nrow(filtered_sampleGroup))
      names(sample_colors) <- filtered_sampleGroup$FinalName
      
      for (i in 1:nrow(filtered_sampleGroup)) {
        pltype <- filtered_sampleGroup$PLtype[i]
        context <- filtered_sampleGroup$Context[i]
        biogroup <- filtered_sampleGroup$bioGroup[i]
        
        # Spatial highest priority
        if (grepl("Spatial", biogroup, ignore.case = TRUE)) {
          sample_colors[i] <- "black"
        } else {
          # Otherwise map by PLtype and Context
          if (grepl("Light", pltype, ignore.case = TRUE)) {
            # Light - blues
            if (context == "Experiment") {
              sample_colors[i] <- "#0066CC"  # dark blue
            } else {
              sample_colors[i] <- "#87CEEB"  # light blue
            }
          } else if (grepl("H2O2", pltype, ignore.case = TRUE)) {
            # H2O2 - oranges
            if (context == "Experiment") {
              sample_colors[i] <- "#FF8C00"  # dark orange
            } else {
              sample_colors[i] <- "#FFD700"  # light orange/gold
            }
          } else if (grepl("PL", pltype, ignore.case = TRUE)) {
            # PL - greens
            if (context == "Experiment") {
              sample_colors[i] <- "#228B22"  # dark green
            } else {
              sample_colors[i] <- "#90EE90"  # light green
            }
          } else {
            # Other - greys
            if (context == "Experiment") {
              sample_colors[i] <- "#696969"  # dark grey
            } else {
              sample_colors[i] <- "#D3D3D3"  # light grey
            }
          }
        }
      }
      
      # Build mapping for Sample_Color
      unique_sample_colors <- unique(sample_colors)
      sample_color_mapping <- setNames(unique_sample_colors, unique_sample_colors)
      
      # Add per-sample color (for main heatmap column bar)
      Annotated_col_corr$Sample_Color <- factor(sample_colors[filtered_sampleGroup$FinalName], 
                                                 levels = unique_sample_colors)
      
      # Compose annotation color list
      anno_colors_corr <- list(
        PLtype = pltype_colors,
        Context = context_colors,
        Sample_Color = sample_color_mapping
      )
      
      cat("✓ Column annotations ready\n")
      cat(sprintf("  - PLtype categories: %d\n", length(unique_pltypes)))
      cat(sprintf("  - Context categories: %d\n", length(unique_contexts)))
    }
  }
  
  # 6. Process each selected version ####
  for (selected_version in selected_versions) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Processing version: %s\n", selected_version))
    cat(sprintf("========================================\n"))
    
    mydata <- imputed_data_list[[selected_version]]
    
    # Identify gene column
    gene_col <- if ("Gene" %in% names(mydata)) "Gene" else names(mydata)[1]
    
    # Identify annotation columns
    anno_cols <- grep("Localization|Annotation", names(mydata), value = TRUE)
    
    # Identify data columns
    data_cols <- setdiff(names(mydata), c(gene_col, anno_cols))
    data_cols <- data_cols[sapply(mydata[data_cols], is.numeric)]
    
    cat(sprintf("  - Gene column: %s\n", gene_col))
    cat(sprintf("  - Data columns: %d\n", length(data_cols)))
    cat(sprintf("  - Annotation columns: %d\n", length(anno_cols)))
    
    # 7. Draw correlation heatmap ####
    if ("correlation" %in% heatmap_types) {
      cat("\n[4] Plot correlation heatmap...\n")
      
      # Filter columns for correlation heatmap
      corr_data_cols <- intersect(data_cols, corr_samples)
      
      if (length(corr_data_cols) < 2) {
        warning("⚠ Fewer than 2 columns for correlation heatmap; skipping")
      } else {
        pdf_file <- file.path(dir_config$output, 
                             sprintf("Module06_Correlation_heatmap_%s.pdf", selected_version))
        pdf(pdf_file, width = 14, height = 10)
        
        # Prepare data matrix
        Pearson_matrix <- mydata[, corr_data_cols]
        rownames(Pearson_matrix) <- mydata[[gene_col]]
        
        # Compute Pearson correlation
        r <- cor(Pearson_matrix,
                 method = "pearson",
                 use = "complete.obs")
        
        # Prepare color params
        corr_params <- generate_pheatmap_params(
          r,
          custom_min = correlation_config$corr_min,
          custom_max = correlation_config$corr_max,
          custom_center = correlation_config$corr_center,
          col_low = correlation_config$col_low,
          col_center = correlation_config$col_center,
          col_high = correlation_config$col_high
        )
        
        # Filter column annotations
        Annotated_col_corr_filtered <- Annotated_col_corr[corr_data_cols, , drop = FALSE]
        
        # Plot heatmap
        pheatmap(r,
                 show_colnames = TRUE,
                 show_rownames = TRUE,
                 fontsize = 10,
                 color = corr_params$color,
                 breaks = corr_params$breaks,
                 border_color = NA,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 clustering_method = "ward.D2",
                 clustering_distance_cols = "correlation",
                 clustering_distance_rows = "correlation",
                 annotation_col = Annotated_col_corr_filtered,
                 annotation_row = Annotated_col_corr_filtered,
                 annotation_colors = anno_colors_corr,
                 main = paste("Correlation Heatmap -", selected_version))
        
        dev.off()
        cat(sprintf("✓ Saved: %s\n", pdf_file))
      }
    }
    
    # 8. Prepare annotations for AllLocalization and by_localization heatmaps ####
    if ("all" %in% heatmap_types || "by_localization" %in% heatmap_types) {
      cat("\n[5] Prepare AllLocalization heatmap data...\n")
      
      # Auto-set Localization columns
      if (is.null(localization_columns)) {
        localization_columns <- anno_cols
      }
      
      # Validate localization columns exist
      valid_loc_cols <- intersect(localization_columns, anno_cols)
      if (length(valid_loc_cols) == 0) {
        warning("⚠ No valid Localization columns found; will plot basic heatmaps only")
        localization_columns <- NULL
      } else {
        localization_columns <- valid_loc_cols
      }
      
      # Set default color params
      if (is.null(color_params)) {
        color_params <- list(
          custom_min = -4,
          custom_max = 4,
          custom_center = 0,
          col_low = "#4D97CD",
          col_center = "white",
          col_high = "#DB6968"
        )
      }
      
      # Compute color params
      custom_params <- generate_pheatmap_params(
        mydata[data_cols],
        custom_min = color_params$custom_min,
        custom_max = color_params$custom_max,
        custom_center = color_params$custom_center,
        col_low = color_params$col_low,
        col_center = color_params$col_center,
        col_high = color_params$col_high
      )
      
      # Remove NA rows
      data_clean <- mydata %>% na.omit()
      cat(sprintf("  - Original rows: %d, after removing NA: %d\n", nrow(mydata), nrow(data_clean)))
      
      if (nrow(data_clean) == 0) {
        warning("⚠ No data after removing NA; skipping AllLocalization heatmap")
      } else {
        # Prepare row annotations
        if (!is.null(localization_columns)) {
          Annotated_row <- data.frame(row.names = data_clean[[gene_col]])
          
          for (col in localization_columns) {
            Annotated_row[[col]] <- factor(
              data_clean[[col]],
              levels = c("SGs", "Nuclear", "Mitochondrion", "Cytosol", "Nuclear_Cytosol", "Other")
            )
          }
          
          # Row annotation colors
          anno_colors_list <- list()
          loc_colors <- c(
            SGs = "red",
            Nuclear = "blue",
            Mitochondrion = "black",
            Cytosol = "green",
            Nuclear_Cytosol = "yellow",
            Other = "white"
          )
          
          for (col in localization_columns) {
            anno_colors_list[[col]] <- loc_colors
          }
        } else {
          Annotated_row <- NULL
          anno_colors_list <- list()
        }
        
        # Prepare column annotations
        sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
        valid_samples <- intersect(data_cols, names(sample_to_group))
        
        if (length(valid_samples) > 0) {
          Annotated_col <- data.frame(
            Sample_Group = factor(sample_to_group[valid_samples]),
            row.names = valid_samples
          )
          
          # Column annotation colors
          unique_groups <- unique(sampleGroup$bioGroup)
          n_groups <- length(unique_groups)
          
          if (n_groups <= 12) {
            group_colors <- hue_pal()(n_groups)
          } else {
            group_colors <- rainbow(n_groups)
          }
          names(group_colors) <- unique_groups
          
          anno_colors_list$Sample_Group <- group_colors
        } else {
          Annotated_col <- NULL
        }
        
        # 9. Plot AllLocalization heatmap ####
        if ("all" %in% heatmap_types) {
          cat("\n[6] Plot AllLocalization heatmap...\n")
          
          pdf_file <- file.path(dir_config$output, 
                               sprintf("Module06_AllLocalization_heatmap_%s.pdf", selected_version))
          pdf(pdf_file, width = 10, height = 10)
          
          Dataforfigure <- data_clean[, data_cols]
          rownames(Dataforfigure) <- data_clean[[gene_col]]
          
          pheatmap(Dataforfigure,
                   cluster_cols = TRUE,
                   cluster_rows = TRUE,
                   clustering_method = "ward.D2",
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "correlation",
                   color = custom_params$color,
                   breaks = custom_params$breaks,
                   fontsize_row = 0.5,
                   scale = "row",
                   border_color = NA,
                   annotation_row = Annotated_row,
                   annotation_col = Annotated_col,
                   annotation_colors = anno_colors_list,
                   main = paste("Localization Heatmap -", selected_version))
          
          dev.off()
          cat(sprintf("✓ Saved: %s\n", pdf_file))
        }
        
        # 10. Plot heatmaps by Localization ####
        if ("by_localization" %in% heatmap_types && !is.null(localization_columns)) {
          cat("\n[7] Plot heatmaps by Localization...\n")
          
          loc_categories <- c("SGs", "Mitochondrion", "Cytosol", "Nuclear", "Nuclear_Cytosol", "Other")
          
          for (loc_col in localization_columns) {
            cat(sprintf("  - Processing %s...\n", loc_col))
            
            pdf_file <- file.path(dir_config$output, 
                                 sprintf("Module06_Different_%s_heatmap_%s.pdf", loc_col, selected_version))
            pdf(pdf_file, width = 10, height = 10)
            
            for (loc_cat in loc_categories) {
              mydata_filtered <- data_clean %>% filter(.data[[loc_col]] == loc_cat)
              
              if (nrow(mydata_filtered) > 0) {
                Dataforfigure <- mydata_filtered[, data_cols]
                rownames(Dataforfigure) <- mydata_filtered[[gene_col]]
                
                if (!is.null(Annotated_row)) {
                  Annotated_row_filtered <- Annotated_row[mydata_filtered[[gene_col]], , drop = FALSE]
                } else {
                  Annotated_row_filtered <- NULL
                }
                
                tryCatch({
                  pheatmap(Dataforfigure,
                           cluster_cols = TRUE,
                           cluster_rows = TRUE,
                           clustering_method = "ward.D2",
                           clustering_distance_cols = "correlation",
                           clustering_distance_rows = "correlation",
                           color = custom_params$color,
                           breaks = custom_params$breaks,
                           fontsize_row = 0.5,
                           scale = "row",
                           border_color = NA,
                           annotation_row = Annotated_row_filtered,
                           annotation_col = Annotated_col,
                           annotation_colors = anno_colors_list,
                           main = paste(loc_col, "Heatmap -", loc_cat, "-", selected_version))
                }, error = function(e) {
                  cat(sprintf("    ⚠ Skipped %s (%s): %s\n", loc_col, loc_cat, e$message))
                })
              }
            }
            
            dev.off()
            cat(sprintf("  ✓ Saved: %s\n", pdf_file))
          }
        }
      }
    }
  }
  
  # 11. Return result ####
  cat("\n[8] Heatmap generation complete\n")
  
  heatmap_info <- list(
    heatmap_params = list(
      selected_versions = selected_versions,
      heatmap_types = heatmap_types,
      localization_columns = localization_columns,
      correlation_config = correlation_config,
      color_params = color_params
    ),
    annotation_info = list(
      anno_cols = anno_cols,
      data_cols = data_cols,
      gene_col = gene_col
    )
  )
  
  cat("✓ Module 6 complete\n")
  
  return(invisible(heatmap_info))
}
