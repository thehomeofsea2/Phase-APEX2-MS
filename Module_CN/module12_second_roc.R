## ============================================================================
## Module 12: 第二次 ROC 分析
## 对齐 CleanCode.R (2494-2702) 逻辑
## ============================================================================

module12_second_roc <- function(
    dir_config,
    FDR_combined_df_list_2nd,
    expr_data_list,
    comparison_cols = NULL,
    comparisons_list = NULL,
    annotation_column = "GO_Localization",
    tp_label = "SGs",
    fp_label = "Cytosol",
    desired_order = NULL,
    min_tp = 0.3,
    sample_columns = NULL,
    youden_ylim = c(0, 0.7)
) {
  cat("\n=== Module 12: 第二次ROC分析 ===\n\n")

  required_pkgs <- c("pROC", "openxlsx", "dplyr", "purrr", "stringr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("✗ 错误：需要安装 %s 包", pkg))
    }
  }

  library(pROC)
  library(openxlsx)
  library(dplyr)
  library(purrr)
  library(stringr)

  if (length(FDR_combined_df_list_2nd) == 0) {
    stop("✗ 错误：FDR_combined_df_list_2nd 为空，无法进行ROC分析")
  }
  if (length(expr_data_list) == 0) {
    stop("✗ 错误：expr_data_list 为空，无法合并表达矩阵")
  }

  # 名称对齐（默认按现有顺序，如果提供desired_order则按顺序筛选）
  common_names <- intersect(names(FDR_combined_df_list_2nd), names(expr_data_list))
  if (length(common_names) == 0) {
    stop("✗ 错误：FDR结果与表达矩阵的名称没有重叠，请检查输入列表命名")
  }

  if (!is.null(desired_order)) {
    desired_order <- desired_order[desired_order %in% common_names]
  }
  if (is.null(desired_order) || length(desired_order) == 0) {
    desired_order <- common_names
  }

  FDR_combined_df_list_2nd <- FDR_combined_df_list_2nd[desired_order]
  expr_data_list <- expr_data_list[desired_order]

  output_dir <- dir_config$output
  if (is.null(output_dir) || !dir.exists(output_dir)) {
    stop("✗ 错误：dir_config$output 未设置或目录不存在")
  }

  hardcoded_defaults <- c(
    "K69A1B3_vs_K69C3",
    "K69A1B3_vs_C3",
    "K69A1B3_vs_K20",
    "K69C3_vs_C3",
    "K69C3_vs_K20",
    "C3_vs_K73"
  )

  final_comparison_cols <- comparison_cols
  comparison_source <- NULL

  if (!is.null(final_comparison_cols) && length(final_comparison_cols) > 0) {
    comparison_source <- "自定义参数"
  }

  if ((is.null(final_comparison_cols) || length(final_comparison_cols) == 0) &&
      !is.null(comparisons_list) && length(comparisons_list) > 0) {
    comp_names <- names(comparisons_list)
    if (is.null(comp_names) || any(is.na(comp_names) | comp_names == "")) {
      comp_names <- vapply(
        seq_along(comparisons_list),
        function(idx) {
          comp <- comparisons_list[[idx]]
          if (!is.null(comp) && length(comp) >= 2) {
            paste0(comp[1], "_vs_", comp[2])
          } else {
            paste0("Comparison_", idx)
          }
        },
        character(1)
      )
    }
    final_comparison_cols <- comp_names
    comparison_source <- "Module11"
  }

  if (is.null(final_comparison_cols) || length(final_comparison_cols) == 0) {
    first_df <- FDR_combined_df_list_2nd[[desired_order[1]]]
    logfc_cols <- grep("_logFC$", colnames(first_df), value = TRUE)
    inferred <- sub("_logFC$", "", logfc_cols)
    if (length(inferred) > 0) {
      final_comparison_cols <- inferred
      comparison_source <- "FDR列自动识别"
    }
  }

  if (is.null(final_comparison_cols) || length(final_comparison_cols) == 0) {
    final_comparison_cols <- hardcoded_defaults
    comparison_source <- "默认列表"
  }

  if (length(final_comparison_cols) == 0) {
    stop("✗ 错误：未找到可用于ROC的比较列")
  }

  if (is.null(youden_ylim) || length(youden_ylim) != 2 || any(!is.finite(youden_ylim))) {
    warning("⚠ youden_ylim 无效，已重置为 c(0, 0.7)")
    youden_ylim <- c(0, 0.7)
  } else {
    youden_ylim <- sort(youden_ylim)
  }

  cat(sprintf("✓ 使用 %d 个比较组（来源：%s）\n", 
              length(final_comparison_cols), comparison_source))
  cat(sprintf("  比较列表：%s\n", paste(final_comparison_cols, collapse = ", ")))

  all_desired_thresholds_2nd <- list()
  Expr_FDR_df_list_2nd <- list()

  for (dataset_name in desired_order) {
    cat(sprintf("[Module12] 处理数据集：%s\n", dataset_name))
    fdr_df <- FDR_combined_df_list_2nd[[dataset_name]]
    expr_df <- expr_data_list[[dataset_name]]

    if (!("Gene" %in% colnames(fdr_df))) {
      stop(sprintf("✗ 错误：%s 缺少 Gene 列", dataset_name))
    }
    if (!("Gene" %in% colnames(expr_df))) {
      stop(sprintf("✗ 错误：表达矩阵 %s 缺少 Gene 列", dataset_name))
    }
    if (!(annotation_column %in% colnames(fdr_df))) {
      stop(sprintf("✗ 错误：%s 缺少注释列 %s", dataset_name, annotation_column))
    }

    roc_subset <- fdr_df %>%
      filter(.data[[annotation_column]] %in% c(tp_label, fp_label))

    if (nrow(roc_subset) == 0 ||
        length(intersect(unique(roc_subset[[annotation_column]]), c(tp_label, fp_label))) < 2) {
      warning(sprintf(
        "⚠ 警告：%s 中 %s 缺少 %s 或 %s（样本不足），跳过ROC计算",
        dataset_name, annotation_column, tp_label, fp_label
      ))
      next
    }

    roc_list <- list()
    roc_dfs <- list()
    roc_plot_titles <- character()

    for (comp in final_comparison_cols) {
      logfc_col <- paste0(comp, "_logFC")
      if (!(logfc_col %in% colnames(roc_subset))) {
        warning(sprintf("⚠ 警告：%s 缺少列 %s，跳过该比较", dataset_name, logfc_col))
        next
      }

      roc_obj <- roc(
        response = roc_subset[[annotation_column]],
        predictor = roc_subset[[logfc_col]],
        levels = c(fp_label, tp_label),
        direction = "<"
      )

      roc_expr <- as.character(roc_obj$call)[3]
      roc_key <- comp
      roc_list[[roc_key]] <- roc_obj
      roc_plot_titles <- c(roc_plot_titles, comp)

      roc_data <- data.frame(
        Threshold = roc_obj$thresholds,
        TP = roc_obj$sensitivities,
        FP = 1 - roc_obj$specificities,
        TP_FP = roc_obj$sensitivities - (1 - roc_obj$specificities)
      )
      colnames(roc_data) <- c(
        paste0(roc_expr, "_Threshold"),
        paste0(roc_expr, "_TP"),
        paste0(roc_expr, "_FP"),
        paste0(roc_expr, "TP_FP")
      )
      roc_dfs[[roc_key]] <- roc_data
    }

    if (length(roc_list) == 0) {
      warning(sprintf("⚠ 警告：%s 没有可用的比较，跳过该数据集", dataset_name))
      next
    }

    file_prefix <- "Module12_Step16"

    # 导出ROC数据
    roc_wb <- createWorkbook()
    for (roc_name in names(roc_dfs)) {
      addWorksheet(roc_wb, sheetName = substr(roc_name, 1, 31))
      writeData(roc_wb, sheet = substr(roc_name, 1, 31), roc_dfs[[roc_name]])
    }
    roc_file <- file.path(
      output_dir,
      sprintf("%s_%s_%s_%s_ROCdata.xlsx", file_prefix, dataset_name, fp_label, tp_label)
    )
    saveWorkbook(roc_wb, roc_file, overwrite = TRUE)
    cat(sprintf("  ✓ 导出ROC数据：%s\n", basename(roc_file)))

    # 绘制ROC曲线
    roc_pdf <- file.path(output_dir, sprintf("%s_1_%s_GO_CytosolROC.pdf", file_prefix, dataset_name))
    roc_pdf_open <- FALSE
    tryCatch({
      grDevices::pdf(roc_pdf, width = 12, height = 6)
      roc_pdf_open <- TRUE
      par(mfrow = c(2, 4))
      for (i in seq_along(roc_list)) {
        roc_obj <- roc_list[[i]]
        title_txt <- if (length(roc_plot_titles) >= i) roc_plot_titles[i] else names(roc_list)[i]
        plot(
          roc_obj,
          main = title_txt,
          xlab = "FPR",
          ylab = "TPR",
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
          cex.main = 1.5,
          cex.lab = 1.5,
          cex.axis = 1,
          col = "#DB6968",
          legacy.axes = TRUE
        )
      }
    }, finally = {
      if (roc_pdf_open) {
        grDevices::dev.off()
      }
    })
    cat(sprintf("  ✓ 导出ROC曲线PDF：%s\n", basename(roc_pdf)))

    # 绘制Youden折线图
    youden_pdf <- file.path(output_dir, sprintf("%s_1_%s_Youden_Index_Plots_BaseR.pdf", file_prefix, dataset_name))
    youden_pdf_open <- FALSE
    tryCatch({
      grDevices::pdf(youden_pdf, width = 12, height = 6)
      youden_pdf_open <- TRUE
      par(mfrow = c(2, 4))
      for (roc_name in names(roc_dfs)) {
        df <- roc_dfs[[roc_name]]
        tp_fp_col <- grep("TP_FP$", colnames(df), value = TRUE)[1]
        thres_col <- grep("Threshold$", colnames(df), value = TRUE)[1]
        if (is.na(tp_fp_col) || is.na(thres_col)) next
        optimal_idx <- which.max(df[[tp_fp_col]])
        optimal_point <- df[optimal_idx, , drop = FALSE]
        plot(
          x = df[[thres_col]],
          y = df[[tp_fp_col]],
          type = "l",
          col = "#DB6968",
          lwd = 2,
          main = roc_name,
          xlab = "Log2 FC",
          ylab = "TPR-FPR",
          cex.main = 1.5,
          cex.lab = 1.5,
          cex.axis = 1,
          ylim = youden_ylim
        )
        abline(v = optimal_point[[thres_col]], col = "#8B96AD", lty = 2)
        points(optimal_point[[thres_col]], optimal_point[[tp_fp_col]], col = "black", pch = 18, cex = 2)
        text(
          x = optimal_point[[thres_col]],
          y = optimal_point[[tp_fp_col]],
          labels = sprintf("(%.2f, %.2f)", optimal_point[[thres_col]], optimal_point[[tp_fp_col]]),
          pos = 4,
          col = "black",
          cex = 1.2
        )
      }
    }, finally = {
      if (youden_pdf_open) {
        grDevices::dev.off()
      }
    })
    cat(sprintf("  ✓ 导出Youden折线PDF：%s\n", basename(youden_pdf)))

    # 计算阈值
    threshold_vec <- numeric(length(roc_dfs))
    names(threshold_vec) <- names(roc_dfs)
    for (idx in seq_along(roc_dfs)) {
      df <- roc_dfs[[idx]]
      tp_fp_col <- grep("TP_FP$", colnames(df), value = TRUE)[1]
      thres_col <- grep("Threshold$", colnames(df), value = TRUE)[1]
      tp_col <- grep("_TP$", colnames(df), value = TRUE)[1]
      if (is.na(tp_fp_col) || is.na(thres_col) || is.na(tp_col)) {
        threshold_vec[idx] <- NA_real_
        next
      }
      max_tp_fp <- max(df[[tp_fp_col]], na.rm = TRUE)
      candidate_rows <- df[abs(df[[tp_fp_col]] - max_tp_fp) < .Machine$double.eps^0.5, , drop = FALSE]
      valid_rows <- candidate_rows[candidate_rows[[tp_col]] >= min_tp, , drop = FALSE]
      if (nrow(valid_rows) == 0) {
        threshold_vec[idx] <- 0
      } else {
        threshold_vec[idx] <- min(valid_rows[[thres_col]], na.rm = TRUE)
      }
    }
    all_desired_thresholds_2nd[[dataset_name]] <- threshold_vec

    # 构建 Expr_FDR_df_list_2nd
    annotation_cols <- grep("_Localization$", colnames(fdr_df), value = TRUE)
    if (length(annotation_cols) == 0 && ncol(fdr_df) >= 3) {
      annotation_cols <- tail(colnames(fdr_df), 3)
    }
    fdr_trimmed <- fdr_df %>%
      select(-all_of(annotation_cols))

    if (!is.null(sample_columns)) {
      sample_cols_use <- intersect(sample_columns, colnames(expr_df))
      loc_cols <- grep("_Localization$", colnames(expr_df), value = TRUE)
      expr_keep <- unique(c("Gene", sample_cols_use, loc_cols))
      expr_df <- expr_df %>% select(any_of(expr_keep))
    }

    merged_expr_fdr <- left_join(expr_df, fdr_trimmed, by = "Gene")
    Expr_FDR_df_list_2nd[[dataset_name]] <- merged_expr_fdr
  }

  if (length(all_desired_thresholds_2nd) == 0) {
    stop("✗ 错误：所有数据集都被跳过，未生成任何ROC结果")
  }

  # 导出阈值表
  threshold_wb <- createWorkbook()
  for (dataset_name in names(all_desired_thresholds_2nd)) {
    addWorksheet(threshold_wb, sheetName = substr(dataset_name, 1, 31))
    comp_names <- names(all_desired_thresholds_2nd[[dataset_name]])
    thr_vec <- as.numeric(all_desired_thresholds_2nd[[dataset_name]])
    thr_df <- data.frame(
      Comparison = comp_names,
      LogFC_Column = paste0(comp_names, "_logFC"),
      Annotation_Column = annotation_column,
      TP_Label = tp_label,
      FP_Label = fp_label,
      Threshold = thr_vec,
      stringsAsFactors = FALSE
    )
    writeData(threshold_wb, sheet = substr(dataset_name, 1, 31), thr_df)
  }
  threshold_file <- file.path(output_dir, "Module12_Step16_all_desired_thresholds_2nd.xlsx")
  saveWorkbook(threshold_wb, threshold_file, overwrite = TRUE)
  cat(sprintf("\n✓ 导出阈值工作簿：%s\n", basename(threshold_file)))

  # 导出 Expr_FDR_df_list_2nd
  expr_fdr_wb <- createWorkbook()
  for (dataset_name in names(Expr_FDR_df_list_2nd)) {
    addWorksheet(expr_fdr_wb, sheetName = substr(dataset_name, 1, 31))
    writeData(expr_fdr_wb, sheet = substr(dataset_name, 1, 31), Expr_FDR_df_list_2nd[[dataset_name]])
  }
  expr_fdr_file <- file.path(output_dir, "Module12_Step16_Expr_FDR_df_list_2nd.xlsx")
  saveWorkbook(expr_fdr_wb, expr_fdr_file, overwrite = TRUE)
  cat(sprintf("✓ 导出表达+FDR工作簿：%s\n\n", basename(expr_fdr_file)))

  cat("=== Module 12 完成 ===\n")
  cat(sprintf("✓ 处理数据集数量：%d\n", length(Expr_FDR_df_list_2nd)))

  return(list(
    FDR_combined_df_list_2nd = FDR_combined_df_list_2nd,
    Expr_FDR_df_list_2nd = Expr_FDR_df_list_2nd,
    all_desired_thresholds_2nd = all_desired_thresholds_2nd
  ))
}

