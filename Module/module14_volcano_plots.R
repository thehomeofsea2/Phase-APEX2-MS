## ============================================================================
## Module 14: Volcano plots based on Expr_FDR_df_list_2nd_SubSGs
## Aligns with CleanCode.R (2747-3322) logic but configurable
## ============================================================================

module14_default_category_ranges <- function() {
  list(
    Exp_vs_Exp = list(xlim = c(-2, 6), ylim = c(0, 13)),
    Exp_vs_Spatial = list(xlim = c(-2.5, 6), ylim = c(0, 13))
  )
}

module14_build_comparison_sets_from_metadata <- function(
    comparison_metadata,
    category_ranges = module14_default_category_ranges()
) {
  if (is.null(comparison_metadata) || nrow(comparison_metadata) == 0) {
    return(list())
  }

  metadata <- comparison_metadata
  metadata$category <- ifelse(
    is.na(metadata$category) | metadata$category == "",
    "Uncategorized",
    metadata$category
  )

  split_meta <- split(metadata, metadata$category)

  sets <- lapply(names(split_meta), function(cat) {
    subset_meta <- split_meta[[cat]]
    comps <- unique(subset_meta$comparison)
    comps <- comps[!is.na(comps) & comps != ""]
    if (length(comps) == 0) {
      return(NULL)
    }

    range_cfg <- category_ranges[[cat]]
    xlim <- if (!is.null(range_cfg) && !is.null(range_cfg$xlim)) range_cfg$xlim else NULL
    ylim <- if (!is.null(range_cfg) && !is.null(range_cfg$ylim)) range_cfg$ylim else NULL

    list(
      name = cat,
      fc_cols = paste0(comps, "_logFC"),
      fdr_cols = paste0(comps, "_adj.P.Val"),
      xlim = xlim,
      ylim = ylim
    )
  })

  Filter(Negate(is.null), sets)
}

module14_default_comparison_sets <- function(comparison_metadata = NULL) {
  if (!is.null(comparison_metadata) && nrow(comparison_metadata) > 0) {
    dynamic_sets <- module14_build_comparison_sets_from_metadata(comparison_metadata)
    if (length(dynamic_sets) > 0) {
      return(dynamic_sets)
    }
  }

  list(
    list(
      name = "Exp_vs_Exp",
      fc_cols = c(
        "K69A1B3_vs_K69C3_logFC",
        "K69A1B3_vs_C3_logFC",
        "K69C3_vs_C3_logFC"
      ),
      fdr_cols = c(
        "K69A1B3_vs_K69C3_adj.P.Val",
        "K69A1B3_vs_C3_adj.P.Val",
        "K69C3_vs_C3_adj.P.Val"
      ),
      xlim = c(-2, 6),
      ylim = c(0, 13)
    ),
    list(
      name = "Exp_vs_Spatial",
      fc_cols = c(
        "K69A1B3_vs_K20_logFC",
        "K69C3_vs_K20_logFC",
        "C3_vs_K73_logFC"
      ),
      fdr_cols = c(
        "K69A1B3_vs_K20_adj.P.Val",
        "K69C3_vs_K20_adj.P.Val",
        "C3_vs_K73_adj.P.Val"
      ),
      xlim = c(-2.5, 6),
      ylim = c(0, 13)
    )
  )
}

module14_collect_annotation_cols <- function(df_list) {
  unique(unlist(lapply(df_list, function(df) {
    grep("_Localization$", colnames(df), value = TRUE)
  })))
}

module14_collect_annotation_values <- function(df_list, column) {
  values <- unique(unlist(lapply(df_list, function(df) {
    if (column %in% colnames(df)) df[[column]] else NULL
  })))
  values[!is.na(values)]
}

module14_collect_comparisons <- function(df_list) {
  for (df in df_list) {
    logfc_cols <- grep("_logFC$", colnames(df), value = TRUE)
    fdr_cols <- grep("_adj\\.P\\.Val$", colnames(df), value = TRUE)
    candidates <- intersect(
      gsub("_logFC$", "", logfc_cols),
      gsub("_adj\\.P\\.Val$", "", fdr_cols)
    )
    if (length(candidates) > 0) {
      return(candidates)
    }
  }
  character(0)
}

module14_build_comparison_metadata <- function(
    comparisons_list,
    bioGroup_info
) {
  if (is.null(comparisons_list) || length(comparisons_list) == 0 ||
      is.null(bioGroup_info) || nrow(bioGroup_info) == 0) {
    return(NULL)
  }

  ctx_lookup <- bioGroup_info %>%
    select(bioGroup, Context) %>%
    distinct()

  rows <- lapply(seq_along(comparisons_list), function(idx) {
    name <- names(comparisons_list)[idx]
    groups <- comparisons_list[[idx]]
    g1 <- groups[1]
    g2 <- groups[2]
    ctx1 <- ctx_lookup$Context[match(g1, ctx_lookup$bioGroup)]
    ctx2 <- ctx_lookup$Context[match(g2, ctx_lookup$bioGroup)]
    ctx1 <- ifelse(is.na(ctx1), NA_character_, ctx1)
    ctx2 <- ifelse(is.na(ctx2), NA_character_, ctx2)

    category <- dplyr::case_when(
      identical(ctx1, "Experiment") && identical(ctx2, "Experiment") ~ "Exp_vs_Exp",
      identical(ctx1, "Experiment") && identical(ctx2, "Spatial") ~ "Exp_vs_Spatial",
      identical(ctx1, "Spatial") && identical(ctx2, "Experiment") ~ "Exp_vs_Spatial",
      TRUE ~ "Other"
    )

    data.frame(
      comparison = name,
      group1 = g1,
      group2 = g2,
      group1_context = ctx1,
      group2_context = ctx2,
      category = category,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

module14_print_available_options <- function(
    expr_fdr_df_list_2nd_subsgs,
    comparison_sets = NULL,
    comparison_metadata = NULL,
    max_preview = 15
) {
  if (is.null(comparison_sets)) {
    comparison_sets <- module14_default_comparison_sets(comparison_metadata)
  }
  if (length(expr_fdr_df_list_2nd_subsgs) == 0) {
    cat("[Module14] ⚠ Expr_FDR_df_list_2nd_SubSGs is empty; cannot show available options\n")
    return(invisible(FALSE))
  }

  versions <- names(expr_fdr_df_list_2nd_subsgs)
  cat("\n[Module14] Available dataset versions:\n")
  cat(sprintf("  %s\n", paste(versions, collapse = ", ")))

  annotation_cols <- module14_collect_annotation_cols(expr_fdr_df_list_2nd_subsgs)
  cat("\n[Module14] Available annotation columns and values:\n")
  if (length(annotation_cols) == 0) {
    cat("  (No *_Localization columns detected)\n")
  } else {
    for (col in annotation_cols) {
      values <- module14_collect_annotation_values(expr_fdr_df_list_2nd_subsgs, col)
      preview <- head(sort(values), max_preview)
      cat(sprintf("  - %s: %s", col, paste(preview, collapse = ", ")))
      extra <- length(values) - length(preview)
      if (extra > 0) cat(sprintf(" ... (+%d)", extra))
      cat("\n")
    }
  }

  comparisons <- module14_collect_comparisons(expr_fdr_df_list_2nd_subsgs)
  cat("\n[Module14] Available comparisons (paired logFC/FDR columns):\n")
  if (length(comparisons) == 0) {
    cat("  (No matching comparison columns found)\n")
  } else {
    cat(sprintf("  %s\n", paste(comparisons, collapse = ", ")))
  }

  if (!is.null(comparison_metadata) && nrow(comparison_metadata) > 0) {
    cat("\n[Module14] Comparison category counts:\n")
    for (cat_name in unique(comparison_metadata$category)) {
      subset <- comparison_metadata %>% filter(category == cat_name)
      cat(sprintf("  - %s (%d): %s\n",
                  cat_name,
                  nrow(subset),
                  paste(subset$comparison, collapse = ", ")))
    }
  }

  cat("\n[Module14] Default comparison groups:\n")
  if (length(comparison_sets) == 0) {
    cat("  (comparison_sets not defined)\n")
  } else {
    for (set in comparison_sets) {
      cat(sprintf("  - %s\n", set$name))
      cat(sprintf("      logFC: %s\n", paste(set$fc_cols, collapse = ", ")))
      cat(sprintf("      FDR : %s\n", paste(set$fdr_cols, collapse = ", ")))
    }
  }

  invisible(TRUE)
}

module14_try_load_module13 <- function(
    workspace_file = "Module13_workspace.RData"
) {
  if (exists("Expr_FDR_df_list_2nd_SubSGs", inherits = TRUE)) {
    return(FALSE)
  }
  if (!file.exists(workspace_file)) {
    cat(sprintf("[Module14] Note: %s not found; run Module 13 to generate this file\n",
                workspace_file))
    return(FALSE)
  }
  load(workspace_file, envir = .GlobalEnv)
  cat(sprintf("[Module14] ✓ Loaded %s; SubSG data available\n", workspace_file))
  TRUE
}

module14_autoprint_on_source <- function() {
  module14_try_load_module13()
  if (exists("Expr_FDR_df_list_2nd_SubSGs", inherits = TRUE)) {
    try(
      module14_print_available_options(
        expr_fdr_df_list_2nd_subsgs = get(
          "Expr_FDR_df_list_2nd_SubSGs",
          inherits = TRUE
        ),
        comparison_sets = module14_default_comparison_sets()
      ),
      silent = FALSE
    )
  } else {
    cat("[Module14] Note: Expr_FDR_df_list_2nd_SubSGs missing in the current environment; run Module 13 first\n")
  }
}

module14_volcano_plots <- function(
    dir_config,
    expr_fdr_df_list_2nd_subsgs,
    versions = NULL,
    annotation_column = "GO_Localization",
    annotation_color_map = c(
      "SGs" = "red",
      "Mitochondrion" = "black",
      "Nuclear" = "blue"
    ),
    use_threshold_colors = FALSE,
    logfc_threshold = 0.5,
    fdr_threshold = 0.05,
    threshold_color_above = "red",
    threshold_color_below = "grey70",
    label_mode = c("with", "without"),
    label_use_ggrepel = TRUE,
    label_use_segments = TRUE,
    label_segment_color = "grey60",
    label_annotations = c("SGs", "Mitochondrion"),
    label_size = 3,
    label_max_overlaps = 200,
    label_force = 1,
    label_box_padding = 0.25,
    comparison_sets = NULL,
    xlim_override = c(-3, 4),
    ylim_override = c(0, 7.5),
    pdf_width = 6,
    pdf_height = 5,
    output_prefix = "Module14_Step20",
    comparison_metadata = NULL,
    comparison_categories = NULL
) {

  cat("\n=== Module 14: Volcano plots (based on Expr_FDR_df_list_2nd_SubSGs) ===\n")

  required_pkgs <- c("dplyr", "ggplot2", "ggrepel", "scales", "rlang", "openxlsx")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("✗ Error: package %s is required", pkg))
    }
  }

  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(rlang)
  library(openxlsx)

  if (length(expr_fdr_df_list_2nd_subsgs) == 0) {
    stop("✗ Error: Expr_FDR_df_list_2nd_SubSGs is empty; cannot draw volcano plots")
  }

  comparison_sets_overview <- if (is.null(comparison_sets)) {
    module14_default_comparison_sets(comparison_metadata)
  } else {
    comparison_sets
  }
  if (length(comparison_sets_overview) == 0) {
    stop("✗ Error: comparison_sets could not be built; check inputs or provide a custom configuration")
  }

  module14_print_available_options(
    expr_fdr_df_list_2nd_subsgs = expr_fdr_df_list_2nd_subsgs,
    comparison_sets = comparison_sets_overview,
    comparison_metadata = comparison_metadata
  )

  available_versions <- names(expr_fdr_df_list_2nd_subsgs)
  cat(sprintf("Available datasets (total %d): %s\n",
              length(available_versions),
              paste(available_versions, collapse = ", ")))

  if (is.null(versions)) {
    versions <- available_versions
  }
  versions <- versions[versions %in% available_versions]
  if (length(versions) == 0) {
    stop("✗ Error: No matching versions available for plotting; check the versions parameter")
  }

  reference_df <- expr_fdr_df_list_2nd_subsgs[[versions[1]]]
  annotation_cols <- module14_collect_annotation_cols(expr_fdr_df_list_2nd_subsgs)

  comparison_candidates <- module14_collect_comparisons(list(reference_df))

  if (length(comparison_candidates) == 0) {
    stop("✗ Error: No paired *_logFC and *_adj.P.Val columns detected; cannot draw volcano plots")
  }

  if (!(annotation_column %in% colnames(reference_df))) {
    stop(sprintf("✗ Error: Annotation column %s does not exist; update the annotation_column parameter", annotation_column))
  }

  resolve_color_map <- function(color_map, column_name) {
    if (is.list(color_map) && !is.null(names(color_map))) {
      if (column_name %in% names(color_map)) {
        return(color_map[[column_name]])
      }
      warning(sprintf("⚠ Note: color_map does not include %s; attempting to merge all mappings", column_name))
      color_map <- unlist(color_map, use.names = TRUE)
    }
    if (is.null(names(color_map))) {
      stop("✗ Error: annotation_color_map must be named, e.g., c('SGs'='red')")
    }
    color_map
  }

  annotation_color_map <- resolve_color_map(annotation_color_map, annotation_column)

  allowed_categories <- comparison_categories
  if (!is.null(comparison_metadata) && is.null(allowed_categories)) {
    allowed_categories <- unique(comparison_metadata$category)
  }
  if (!is.null(allowed_categories) && is.null(comparison_metadata)) {
    warning("⚠ Note: comparison_categories supplied but comparison_metadata is missing; ignoring category filter")
    allowed_categories <- NULL
  }
  if (!is.null(allowed_categories)) {
    cat(sprintf("\n[Module14] Only plotting the following comparison categories: %s\n",
                paste(allowed_categories, collapse = ", ")))
  }

  default_xlim <- c(-3, 4)
  default_ylim <- c(0, 7.5)

  normalize_limit <- function(limit_vec, fallback_vec) {
    if (is.null(limit_vec)) {
      return(NULL)
    }
    if (length(limit_vec) != 2 || any(!is.finite(limit_vec))) {
      warning("⚠ Note: Axis range parameters are invalid; using default values")
      return(sort(fallback_vec))
    }
    sort(limit_vec)
  }

  xlim_override_norm <- normalize_limit(xlim_override, default_xlim)
  ylim_override_norm <- normalize_limit(ylim_override, default_ylim)

  apply_axis_defaults <- function(set) {
    if (!is.null(xlim_override_norm)) {
      set$xlim <- xlim_override_norm
    } else if (is.null(set$xlim)) {
      set$xlim <- default_xlim
    }
    if (!is.null(ylim_override_norm)) {
      set$ylim <- ylim_override_norm
    } else if (is.null(set$ylim)) {
      set$ylim <- default_ylim
    }
    set
  }

  comparison_sets <- lapply(comparison_sets_overview, apply_axis_defaults)
  if (!is.list(comparison_sets) || length(comparison_sets) == 0) {
    stop("✗ Error: comparison_sets must be a list containing at least one element")
  }
  if (!is.null(allowed_categories) && !is.null(comparison_metadata)) {
    get_set_names <- function(sets) {
      if (length(sets) == 0) return(character(0))
      vapply(
        sets,
        function(x) {
          if (!is.null(x$name)) x$name else "<unnamed>"
        },
        character(1)
      )
    }

    has_allowed_category <- function(set) {
      if (is.null(set$fc_cols)) {
        return(FALSE)
      }
      fc_names <- gsub("_logFC$", "", set$fc_cols)
      matched_idx <- match(fc_names, comparison_metadata$comparison)
      matched_categories <- comparison_metadata$category[matched_idx]
      any(!is.na(matched_categories) & matched_categories %in% allowed_categories)
    }

    filtered_sets <- Filter(has_allowed_category, comparison_sets)

    if (length(filtered_sets) == 0) {
      stop("✗ Error: comparison_categories did not match any comparison combinations; check the configuration")
    }

    skipped_names <- setdiff(get_set_names(comparison_sets), get_set_names(filtered_sets))
    if (length(skipped_names) > 0) {
      cat(sprintf(
        "[Module14] Skipping comparison combinations per comparison_categories: %s\n",
        paste(skipped_names, collapse = ", ")
      ))
    }

    comparison_sets <- filtered_sets
  }

  parse_label_mode <- function(modes) {
    allowed <- c("with", "without")
    modes <- unique(modes)
    matched <- modes[modes %in% allowed]
    if (length(matched) == 0) {
      stop("✗ Error: label_mode only supports 'with' or 'without'")
    }
    matched
  }
  label_mode <- parse_label_mode(label_mode)

  ensure_dir <- function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  ensure_dir(dir_config$output)

  pdf_records <- list()
  volcano_export_data <- list()

  build_plot_data <- function(df, fc_col, fdr_col) {
    df %>%
      filter(
        !is.na(.data[[fc_col]]),
        !is.na(.data[[fdr_col]])
      ) %>%
      mutate(
        log10FDR = -log10(abs(.data[[fdr_col]])),
        pass_threshold = abs(.data[[fc_col]]) >= logfc_threshold &
          abs(.data[[fdr_col]]) <= fdr_threshold
      )
  }

  for (mode in label_mode) {
    label_suffix <- ifelse(mode == "with", "withLabels", "noLabels")

    for (set in comparison_sets) {
      required_fields <- c("name", "fc_cols", "fdr_cols")
      missing_fields <- required_fields[!required_fields %in% names(set)]
      if (length(missing_fields) > 0) {
        warning(sprintf("⚠ Warning: Comparison configuration is missing fields: %s; skipping this configuration",
                        paste(missing_fields, collapse = ", ")))
        next
      }
      if (length(set$fc_cols) != length(set$fdr_cols)) {
        warning(sprintf("⚠ Warning: logFC/FDR column counts mismatch for comparison %s; skipping", set$name))
        next
      }

      pdf_file <- file.path(
        dir_config$output,
        sprintf(
          "%s_%s_%s_%s.pdf",
          output_prefix,
          gsub("[^A-Za-z0-9]+", "_", annotation_column),
          gsub("[^A-Za-z0-9]+", "_", set$name),
          label_suffix
        )
      )

      grDevices::pdf(pdf_file, width = pdf_width, height = pdf_height)
      tryCatch({
        for (dataset_name in versions) {
          df <- expr_fdr_df_list_2nd_subsgs[[dataset_name]]

          missing_cols <- setdiff(
            c(set$fc_cols, set$fdr_cols, annotation_column, "Gene"),
            colnames(df)
          )
          if (length(missing_cols) > 0) {
            warning(sprintf(
              "⚠ Warning: Dataset %s is missing columns %s; skipping this dataset",
              dataset_name,
              paste(missing_cols, collapse = ", ")
            ))
            next
          }

          for (idx in seq_along(set$fc_cols)) {
            fc_col <- set$fc_cols[idx]
            fdr_col <- set$fdr_cols[idx]

            comparison_name <- gsub("_logFC$", "", fc_col)
            comparison_category <- NULL
            if (!is.null(comparison_metadata)) {
              matched_row <- comparison_metadata %>%
                filter(comparison == comparison_name) %>%
                slice(1)
              if (nrow(matched_row) > 0) {
                comparison_category <- matched_row$category
              } else {
                comparison_category <- "Unknown"
              }
              if (!is.null(allowed_categories) &&
                  !(comparison_category %in% allowed_categories)) {
                next
              }
            }

            subset_df <- build_plot_data(df, fc_col, fdr_col)
            if (nrow(subset_df) == 0) {
              warning(sprintf(
                "⚠ Warning: Dataset %s lacks valid (logFC, FDR) data for comparison %s",
                dataset_name,
                fc_col
              ))
              next
            }

            subset_df <- subset_df %>%
              mutate(annotation_value = .data[[annotation_column]])

            if (use_threshold_colors) {
              subset_df <- subset_df %>%
                mutate(
                  point_color = ifelse(
                    pass_threshold,
                    threshold_color_above,
                    threshold_color_below
                  ),
                  label_color = point_color
                )
            } else {
              subset_df <- subset_df %>%
                mutate(
                  point_color = annotation_color_map[annotation_value],
                  label_color = point_color
                )
              subset_df$point_color[is.na(subset_df$point_color)] <- "grey70"
              subset_df$label_color[is.na(subset_df$label_color)] <- "grey70"
            }

            subset_df <- subset_df %>%
              mutate(
                is_significant = abs(.data[[fc_col]]) >= logfc_threshold &
                  abs(.data[[fdr_col]]) <= fdr_threshold,
                point_alpha = ifelse(pass_threshold, 1, 0.5),
                label_text = dplyr::case_when(
                  mode == "with" &
                    annotation_value %in% label_annotations &
                    is_significant ~ Gene,
                  TRUE ~ NA_character_
                ),
                label_color = ifelse(is.na(label_text), NA_character_, label_color)
              )

            category_label <- if (is.null(comparison_category)) "All" else comparison_category

            export_df <- subset_df %>%
              transmute(
                Gene = Gene,
                !!rlang::sym(fc_col) := .data[[fc_col]],
                !!rlang::sym(fdr_col) := .data[[fdr_col]],
                log10FDR = log10FDR,
                annotation_value = annotation_value,
                pass_threshold = pass_threshold,
                is_significant = is_significant,
                comparison_name = comparison_name,
                comparison_group = set$name,
                comparison_category = category_label,
                dataset = dataset_name,
                annotation_column = annotation_column
              )

            export_key <- sprintf("%s|%s|%s|%s", dataset_name, set$name, comparison_name, category_label)
            volcano_export_data[[export_key]] <- export_df

            p <- ggplot(
              subset_df,
              aes(x = .data[[fc_col]], y = log10FDR)
            ) +
              geom_point(
                aes(color = point_color, alpha = point_alpha),
                size = 0.7,
                show.legend = FALSE
              ) +
              geom_vline(
                xintercept = c(logfc_threshold, -logfc_threshold),
                linetype = "dashed",
                color = "black"
              ) +
              geom_hline(
                yintercept = -log10(fdr_threshold),
                linetype = "dashed",
                color = "black"
              ) +
              scale_color_identity() +
              scale_alpha_identity() +
              labs(
                title = sprintf(
                  "Volcano Plot (%s | %s)\n%s\n%s",
                  annotation_column,
                  category_label,
                  dataset_name,
                  fc_col
                ),
                x = "log2 Fold Change",
                y = "-log10(FDR)"
              ) +
              theme_minimal()

            if (!is.null(set$xlim)) {
              p <- p + xlim(set$xlim)
            }
            if (!is.null(set$ylim)) {
              p <- p + ylim(set$ylim)
            }

            if (mode == "with") {
              label_data <- subset_df %>% filter(!is.na(label_text))
              if (nrow(label_data) > 0) {
                segment_col <- if (label_use_segments) label_segment_color else NA
                if (label_use_ggrepel) {
                  p <- p +
                    geom_text_repel(
                      data = label_data,
                      aes(label = label_text, color = label_color),
                      size = label_size,
                      segment.color = segment_col,
                      show.legend = FALSE,
                      max.overlaps = label_max_overlaps,
                      force = label_force,
                      box.padding = label_box_padding
                    )
                } else {
                  p <- p +
                    geom_text(
                      data = label_data,
                      aes(label = label_text, color = label_color),
                      size = label_size,
                      show.legend = FALSE
                    )
                }
              }
            }

            print(p)
          }
        }
      },
      finally = {
        grDevices::dev.off()
      })

      pdf_records[[label_suffix]][[set$name]] <- pdf_file
      cat(sprintf("  ✓ Generated volcano plot: %s\n", basename(pdf_file)))
    }
  }

  volcano_source_file <- NULL
  if (length(volcano_export_data) > 0) {
    volcano_source_file <- file.path(
      dir_config$output,
      sprintf("%s_VolcanoSource.xlsx", output_prefix)
    )
    wb <- createWorkbook()
    existing_names <- character()
    for (key in names(volcano_export_data)) {
      sanitized_name <- gsub("[^A-Za-z0-9]+", "_", key)
      base_name <- substr(sanitized_name, 1, 31)
      sheet_name <- base_name
      if (sheet_name %in% existing_names) {
        suffix <- 1
        repeat {
          suffix_label <- paste0("_", suffix)
          max_base_len <- max(31 - nchar(suffix_label), 0)
          truncated_base <- if (max_base_len > 0) substr(sanitized_name, 1, max_base_len) else ""
          candidate <- paste0(truncated_base, suffix_label)
          if (!(candidate %in% existing_names)) {
            sheet_name <- candidate
            break
          }
          suffix <- suffix + 1
        }
      }
      existing_names <- c(existing_names, sheet_name)
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, x = volcano_export_data[[key]])
    }
    saveWorkbook(wb, file = volcano_source_file, overwrite = TRUE)
    cat(sprintf("  ✓ Exported volcano plot source data: %s\n", basename(volcano_source_file)))
  } else {
    cat("  ⚠ Note: No volcano plot source data exported (possibly no valid comparisons)\n")
  }

  cat("=== Module 14 completed ===\n")

  return(list(
    versions = versions,
    annotation_column = annotation_column,
    comparison_sets = vapply(comparison_sets, function(x) x$name, character(1)),
    comparison_categories_used = if (is.null(allowed_categories)) "All" else allowed_categories,
    pdf_files = pdf_records,
    source_data_file = volcano_source_file,
    summary = list(
      annotation_columns = annotation_cols,
      comparisons = comparison_candidates
    )
  ))
}

module14_autoprint_on_source()


