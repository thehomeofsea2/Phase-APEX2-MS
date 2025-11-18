## ============================================================================
## Module 14: Volcano plots based on Expr_FDR_df_list_2nd_SubSGs
## Aligns with CleanCode.R (2747-3322) logic but configurable
## ============================================================================

module14_default_comparison_sets <- function() {
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
    comparison_sets = module14_default_comparison_sets(),
    comparison_metadata = NULL,
    max_preview = 15
) {
  if (length(expr_fdr_df_list_2nd_subsgs) == 0) {
    cat("[Module14] ⚠ Expr_FDR_df_list_2nd_SubSGs 为空，无法展示可选项\n")
    return(invisible(FALSE))
  }

  versions <- names(expr_fdr_df_list_2nd_subsgs)
  cat("\n[Module14] 可用数据集版本：\n")
  cat(sprintf("  %s\n", paste(versions, collapse = ", ")))

  annotation_cols <- module14_collect_annotation_cols(expr_fdr_df_list_2nd_subsgs)
  cat("\n[Module14] 可选注释列及元素：\n")
  if (length(annotation_cols) == 0) {
    cat("  (未检测到 *_Localization 列)\n")
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
  cat("\n[Module14] 可用 Comparisons（logFC/FDR 成对列）：\n")
  if (length(comparisons) == 0) {
    cat("  (未找到匹配的比较列)\n")
  } else {
    cat(sprintf("  %s\n", paste(comparisons, collapse = ", ")))
  }

  if (!is.null(comparison_metadata) && nrow(comparison_metadata) > 0) {
    cat("\n[Module14] Comparison 分类统计：\n")
    for (cat_name in unique(comparison_metadata$category)) {
      subset <- comparison_metadata %>% filter(category == cat_name)
      cat(sprintf("  - %s (%d): %s\n",
                  cat_name,
                  nrow(subset),
                  paste(subset$comparison, collapse = ", ")))
    }
  }

  cat("\n[Module14] 默认 Comparison groups：\n")
  if (length(comparison_sets) == 0) {
    cat("  (未定义 comparison_sets)\n")
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
    cat(sprintf("[Module14] 提示：未找到 %s，可先运行 Module 13 以生成该文件\n",
                workspace_file))
    return(FALSE)
  }
  load(workspace_file, envir = .GlobalEnv)
  cat(sprintf("[Module14] ✓ 已加载 %s，检测 SubSG 数据可用\n", workspace_file))
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
    cat("[Module14] 提示：当前环境缺少 Expr_FDR_df_list_2nd_SubSGs，需先运行 Module 13\n")
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
    label_annotations = c("SGs", "Mitochondrion"),
    label_size = 3,
    label_max_overlaps = 15,
    comparison_sets = NULL,
    pdf_width = 6,
    pdf_height = 5,
    output_prefix = "Module14_Step20",
    comparison_metadata = NULL,
    comparison_categories = NULL
) {

  cat("\n=== Module 14: 火山图 (基于 Expr_FDR_df_list_2nd_SubSGs) ===\n")

  required_pkgs <- c("dplyr", "ggplot2", "ggrepel", "scales", "rlang", "openxlsx")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("✗ 错误：需要安装 %s 包", pkg))
    }
  }

  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(rlang)
  library(openxlsx)

  if (length(expr_fdr_df_list_2nd_subsgs) == 0) {
    stop("✗ 错误：Expr_FDR_df_list_2nd_SubSGs 为空，无法绘制火山图")
  }

  comparison_sets_overview <- if (is.null(comparison_sets)) {
    module14_default_comparison_sets()
  } else {
    comparison_sets
  }

  module14_print_available_options(
    expr_fdr_df_list_2nd_subsgs = expr_fdr_df_list_2nd_subsgs,
    comparison_sets = comparison_sets_overview,
    comparison_metadata = comparison_metadata
  )

  available_versions <- names(expr_fdr_df_list_2nd_subsgs)
  cat(sprintf("可用数据集（共 %d 个）：%s\n",
              length(available_versions),
              paste(available_versions, collapse = ", ")))

  if (is.null(versions)) {
    versions <- available_versions
  }
  versions <- versions[versions %in% available_versions]
  if (length(versions) == 0) {
    stop("✗ 错误：没有匹配的版本可供绘图，请检查 versions 参数")
  }

  reference_df <- expr_fdr_df_list_2nd_subsgs[[versions[1]]]
  annotation_cols <- module14_collect_annotation_cols(expr_fdr_df_list_2nd_subsgs)

  comparison_candidates <- module14_collect_comparisons(list(reference_df))

  if (length(comparison_candidates) == 0) {
    stop("✗ 错误：未检测到 *_logFC 与 *_adj.P.Val 成对的列，无法绘制火山图")
  }

  if (!(annotation_column %in% colnames(reference_df))) {
    stop(sprintf("✗ 错误：注释列 %s 不存在，请修改参数 annotation_column", annotation_column))
  }

  resolve_color_map <- function(color_map, column_name) {
    if (is.list(color_map) && !is.null(names(color_map))) {
      if (column_name %in% names(color_map)) {
        return(color_map[[column_name]])
      }
      warning(sprintf("⚠ 提示：color_map 未包含 %s，尝试合并所有映射", column_name))
      color_map <- unlist(color_map, use.names = TRUE)
    }
    if (is.null(names(color_map))) {
      stop("✗ 错误：annotation_color_map 需要具备名称，格式如 c('SGs'='red')")
    }
    color_map
  }

  annotation_color_map <- resolve_color_map(annotation_color_map, annotation_column)

  allowed_categories <- comparison_categories
  if (!is.null(comparison_metadata) && is.null(allowed_categories)) {
    allowed_categories <- unique(comparison_metadata$category)
  }
  if (!is.null(allowed_categories) && is.null(comparison_metadata)) {
    warning("⚠ 提示：提供了 comparison_categories 但缺少 comparison_metadata，忽略分类过滤")
    allowed_categories <- NULL
  }
  if (!is.null(allowed_categories)) {
    cat(sprintf("\n[Module14] 只绘制以下 comparison 类别：%s\n",
                paste(allowed_categories, collapse = ", ")))
  }

  comparison_sets <- comparison_sets_overview
  if (!is.list(comparison_sets) || length(comparison_sets) == 0) {
    stop("✗ 错误：comparison_sets 需要是包含至少一个元素的列表")
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
      stop("✗ 错误：comparison_categories 未匹配到任何比较组合，请检查配置")
    }

    skipped_names <- setdiff(get_set_names(comparison_sets), get_set_names(filtered_sets))
    if (length(skipped_names) > 0) {
      cat(sprintf(
        "[Module14] 根据 comparison_categories 跳过比较组合：%s\n",
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
      stop("✗ 错误：label_mode 只支持 'with' 或 'without'")
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
        warning(sprintf("⚠ 警告：比较配置缺少字段：%s，跳过该配置",
                        paste(missing_fields, collapse = ", ")))
        next
      }
      if (length(set$fc_cols) != length(set$fdr_cols)) {
        warning(sprintf("⚠ 警告：比较 %s 的 logFC/FDR 列数量不匹配，跳过", set$name))
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
              "⚠ 警告：数据集 %s 缺少列 %s，跳过该数据集",
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
                "⚠ 警告：数据集 %s 在比较 %s 缺少有效的 (logFC, FDR) 数据",
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
                aes(color = point_color),
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
                p <- p +
                  geom_text_repel(
                    data = label_data,
                    aes(label = label_text, color = label_color),
                    size = label_size,
                    segment.color = NA,
                    show.legend = FALSE,
                    max.overlaps = label_max_overlaps
                  )
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
      cat(sprintf("  ✓ 已生成火山图：%s\n", basename(pdf_file)))
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
    cat(sprintf("  ✓ 导出火山图源数据：%s\n", basename(volcano_source_file)))
  } else {
    cat("  ⚠ 提示：未生成任何火山图源数据导出（可能没有有效比较）\n")
  }

  cat("=== Module 14 完成 ===\n")

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


