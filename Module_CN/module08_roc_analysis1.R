# ============================================================================
# Module 8: 第一次ROC分析
# ============================================================================
# 功能：
# 1. SubMito定位转化：将注释中的Mitochondrion替换为MitoCarta3的亚定位（MIM, Matrix, MOM, IMS）
# 2. ROC分析：使用pROC包进行ROC曲线分析
# 3. 输出ROC曲线图、Youden Index图、ROC数据表、最佳阈值表
#
# 输入：
# - Module07_workspace.RData（包含diff_results1, comparisons_used等）
# - MitoCarta3.0注释数据（用于SubMito转化）
#
# 输出：
# - Module08_workspace.RData（工作目录）
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

#' 应用SubMito定位转化
#'
#' @param data 包含注释列的数据框
#' @param mitocarta_anno MitoCarta注释数据
#' @param annotation_columns 需要转化的注释列名向量
#' @param enable_submito 是否启用SubMito转化（默认TRUE）
#' @return 转化后的数据框
apply_submito_transformation <- function(data, 
                                          mitocarta_anno, 
                                          annotation_columns = NULL,
                                          enable_submito = TRUE) {
  
  if (!enable_submito) {
    cat("  - SubMito转化已禁用，跳过\n")
    return(data)
  }
  
  # 自动检测注释列（如果未提供）
  if (is.null(annotation_columns)) {
    # 查找包含"Localization"的列
    annotation_columns <- grep("Localization", colnames(data), value = TRUE)
    if (length(annotation_columns) == 0) {
      cat("  - 警告：未找到注释列，跳过SubMito转化\n")
      return(data)
    }
  }
  
  cat(sprintf("  - 检测到 %d 个注释列需要转化\n", length(annotation_columns)))
  
  # 准备线粒体亚定位查找表
  mito_lookup <- mitocarta_anno %>%
    select(Gene = Symbol, mito_subClass = MitoCarta3.0_SubMitoLocalization) %>%
    filter(mito_subClass %in% c("MIM", "Matrix", "MOM", "IMS"))
  
  cat(sprintf("  - MitoCarta3亚定位类别: %s\n", 
              paste(unique(mito_lookup$mito_subClass), collapse = ", ")))
  
  # 使用left_join添加亚定位信息
  data_transformed <- data %>%
    left_join(mito_lookup, by = "Gene")
  
  # 对每个注释列进行转化
  for (col in annotation_columns) {
    if (!(col %in% colnames(data_transformed))) {
      cat(sprintf("  - 警告：列 %s 不存在，跳过\n", col))
      next
    }
    
    # 只替换值为"Mitochondrion"的行，保留其他注释（如Cytosol, Nuclear, SGs, Other等）
    # 这样不会破坏前面模块中设置的注释优先级
    data_transformed[[col]] <- ifelse(
      data_transformed[[col]] == "Mitochondrion" & !is.na(data_transformed$mito_subClass),
      data_transformed$mito_subClass,
      data_transformed[[col]]
    )
    
    cat(sprintf("  - 已转化列: %s（仅替换Mitochondrion -> SubMito）\n", col))
  }
  
  # 移除辅助列
  data_transformed <- data_transformed %>%
    select(-mito_subClass)
  
  cat("  ✓ SubMito转化完成\n")
  return(data_transformed)
}


#' 执行ROC分析
#'
#' @param data 包含logFC和注释列的数据框
#' @param logfc_columns logFC列名向量
#' @param annotation_column 用于ROC分析的注释列名
#' @param tp_label True Positive标签（如"SGs"）
#' @param fp_label False Positive标签（如"Matrix"）
#' @param direction ROC方向（默认"<"）
#' @return 包含ROC对象列表和数据框列表的list
perform_roc_analysis <- function(data,
                                  logfc_columns,
                                  annotation_column,
                                  tp_label,
                                  fp_label,
                                  direction = "<") {
  
  cat(sprintf("  - 执行ROC分析: %s vs %s\n", tp_label, fp_label))
  cat(sprintf("  - 注释列: %s\n", annotation_column))
  cat(sprintf("  - LogFC列数: %d\n", length(logfc_columns)))
  
  # 过滤数据：只保留TP和FP
  data_for_roc <- data %>%
    filter(!!sym(annotation_column) %in% c(tp_label, fp_label))
  
  n_tp <- sum(data_for_roc[[annotation_column]] == tp_label)
  n_fp <- sum(data_for_roc[[annotation_column]] == fp_label)
  
  cat(sprintf("  - %s 数量: %d\n", tp_label, n_tp))
  cat(sprintf("  - %s 数量: %d\n", fp_label, n_fp))
  
  if (n_tp < 5 || n_fp < 5) {
    cat("  - 警告：样本量不足（TP或FP < 5），跳过ROC分析\n")
    return(list(roc_objects = list(), roc_dfs = list()))
  }
  
  # 对每个logFC列执行ROC分析
  roc_objects <- list()
  roc_dfs <- list()
  
  for (i in seq_along(logfc_columns)) {
    logfc_col <- logfc_columns[i]
    
    if (!(logfc_col %in% colnames(data_for_roc))) {
      cat(sprintf("  - 警告：列 %s 不存在，跳过\n", logfc_col))
      next
    }
    
    # 移除NA值
    data_clean <- data_for_roc %>%
      filter(!is.na(!!sym(logfc_col)))
    
    if (nrow(data_clean) < 10) {
      cat(sprintf("  - 警告：%s 有效数据不足，跳过\n", logfc_col))
      next
    }
    
    tryCatch({
      # 执行ROC分析
      roc_obj <- roc(
        response = data_clean[[annotation_column]],
        predictor = data_clean[[logfc_col]],
        levels = c(fp_label, tp_label),
        direction = direction,
        quiet = TRUE
      )
      
      roc_name <- paste0("roc", i)
      roc_objects[[roc_name]] <- roc_obj
      
      # 提取ROC数据
      roc_data <- data.frame(
        Threshold = roc_obj$thresholds,
        TP = roc_obj$sensitivities,
        FP = 1 - roc_obj$specificities,
        TP_FP = roc_obj$sensitivities - (1 - roc_obj$specificities)
      )
      
      # 重命名列（使用logFC列名）
      colnames(roc_data) <- c(
        paste0(logfc_col, "_Threshold"),
        paste0(logfc_col, "_TP"),
        paste0(logfc_col, "_FP"),
        paste0(logfc_col, "_TP_FP")
      )
      
      roc_dfs[[roc_name]] <- roc_data
      
      cat(sprintf("  ✓ 完成: %s (AUC = %.3f)\n", logfc_col, auc(roc_obj)))
      
    }, error = function(e) {
      cat(sprintf("  - 错误：%s - %s\n", logfc_col, e$message))
    })
  }
  
  return(list(roc_objects = roc_objects, roc_dfs = roc_dfs))
}


#' 计算最佳阈值
#'
#' @param roc_dfs ROC数据框列表
#' @param min_tp 最小TP阈值（默认0.3）
#' @return 阈值向量
calculate_optimal_thresholds <- function(roc_dfs, min_tp = 0.3) {
  
  threshold_vec <- numeric(length(roc_dfs))
  comparison_names <- character(length(roc_dfs))
  
  for (i in seq_along(roc_dfs)) {
    df <- roc_dfs[[i]]
    
    # 动态获取列名
    tp_fp_col <- grep("TP_FP", colnames(df), value = TRUE)[1]
    thres_col <- grep("Threshold", colnames(df), value = TRUE)[1]
    tp_col <- grep("_TP$", colnames(df), value = TRUE)[1]
    
    # 提取Comparison名称（移除后缀）
    comparison_names[i] <- sub("_logFC_TP_FP$", "", tp_fp_col)
    
    # 找出TP_FP最大的行
    max_tp_fp <- max(df[[tp_fp_col]], na.rm = TRUE)
    candidate_rows <- df %>%
      filter(!!sym(tp_fp_col) == max_tp_fp)
    
    # 检查这些行中TP是否都 < min_tp
    if (all(candidate_rows[[tp_col]] < min_tp)) {
      threshold_vec[i] <- 0  # 不满足要求，设为0
    } else {
      # 满足TP >= min_tp的，找出最小的Threshold
      valid_candidates <- candidate_rows %>%
        filter(!!sym(tp_col) >= min_tp)
      threshold_vec[i] <- min(valid_candidates[[thres_col]], na.rm = TRUE)
    }
  }
  
  names(threshold_vec) <- comparison_names
  return(threshold_vec)
}


#' 绘制ROC曲线
#'
#' @param roc_objects ROC对象列表
#' @param roc_dfs ROC数据框列表（用于获取列名）
#' @param output_file 输出PDF文件路径
#' @param main_title 主标题
plot_roc_curves <- function(roc_objects, roc_dfs, output_file, main_title = "") {
  
  n_plots <- length(roc_objects)
  if (n_plots == 0) {
    cat("  - 警告：没有可绘制的ROC曲线\n")
    return(invisible(NULL))
  }
  
  # 计算布局
  n_cols <- min(4, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  pdf(output_file, width = 12, height = 6 * n_rows / 2)
  par(mfrow = c(n_rows, n_cols))
  
  for (i in seq_along(roc_objects)) {
    roc_obj <- roc_objects[[i]]
    
    # 从列名提取Comparison名称（和Youden图用同样的方法）
    roc_name <- names(roc_objects)[i]  # 例如 "roc1"
    df <- roc_dfs[[roc_name]]
    
    # 获取TP_FP列名，然后移除后缀得到Comparison名称
    tp_fp_col <- grep("_TP_FP$", colnames(df), value = TRUE)[1]
    # 移除"_logFC_TP_FP"后缀得到Comparison名称
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
    
    # 添加标题（pROC的plot不支持main参数，需要单独添加）
    title(main = predictor_name, cex.main = 1.5)
  }
  
  dev.off()
  cat(sprintf("  ✓ 已保存: %s\n", output_file))
}


#' 绘制Youden Index图
#'
#' @param roc_dfs ROC数据框列表
#' @param output_file 输出PDF文件路径
plot_youden_index <- function(roc_dfs, output_file) {
  
  n_plots <- length(roc_dfs)
  if (n_plots == 0) {
    cat("  - 警告：没有可绘制的Youden Index图\n")
    return(invisible(NULL))
  }
  
  # 计算布局
  n_cols <- min(4, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  pdf(output_file, width = 12, height = 6 * n_rows / 2)
  par(mfrow = c(n_rows, n_cols))
  
  for (i in seq_along(roc_dfs)) {
    df <- roc_dfs[[i]]
    
    # 动态获取列名
    tp_fp_col <- grep("_TP_FP$", colnames(df), value = TRUE)
    thres_col <- grep("_Threshold$", colnames(df), value = TRUE)
    
    # 提取干净的标题（移除"_logFC_TP_FP"后缀得到Comparison名称）
    main_title <- sub("_logFC_TP_FP$", "", tp_fp_col)
    
    # 找出最优点
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
    
    # 添加垂直线和最优点
    abline(v = optimal_point[[thres_col]], col = "#8B96AD", lty = 2)
    points(
      x = optimal_point[[thres_col]],
      y = optimal_point[[tp_fp_col]],
      col = "black",
      pch = 18,
      cex = 2
    )
    
    # 添加标注
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
  cat(sprintf("  ✓ 已保存: %s\n", output_file))
}


#' Module 8: 第一次ROC分析
#'
#' @param dir_config 目录配置
#' @param diff_results1 差异分析结果列表（来自Module 7）
#' @param annotation_references 注释参考数据（来自Module 3）
#' @param selected_versions 选择分析的数据版本（NULL表示所有版本）
#' @param roc_annotation_column 用于ROC分析的注释列（默认"GO_Localization"）
#' @param tp_label True Positive标签（默认"SGs"）
#' @param fp_label False Positive标签（默认"Matrix"）
#' @param enable_submito 是否启用SubMito转化（默认TRUE）
#' @param submito_annotation_columns SubMito转化的注释列（NULL表示自动检测）
#' @param min_tp 最小TP阈值（默认0.3）
#' @return 包含ROC分析结果的list
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
  cat("Module 8: 第一次ROC分析\n")
  cat("========================================\n")
  
  # 验证输入
  if (!all(c("reference", "output") %in% names(dir_config))) {
    stop("错误：dir_config 必须包含 'reference' 和 'output' 路径")
  }
  
  if (length(diff_results1) == 0) {
    stop("错误：diff_results1 为空，请先完成 Module 7")
  }
  
  # 选择要分析的版本
  if (is.null(selected_versions)) {
    selected_versions <- names(diff_results1)
  } else {
    # 验证版本名称
    missing_versions <- setdiff(selected_versions, names(diff_results1))
    if (length(missing_versions) > 0) {
      stop(sprintf("错误：未找到版本: %s", paste(missing_versions, collapse = ", ")))
    }
  }
  
  cat(sprintf("分析版本: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("ROC注释列: %s\n", roc_annotation_column))
  cat(sprintf("TP标签: %s\n", tp_label))
  cat(sprintf("FP标签: %s\n", fp_label))
  cat(sprintf("SubMito转化: %s\n", ifelse(enable_submito, "启用", "禁用")))
  
  # 获取MitoCarta注释
  mitocarta_anno <- annotation_references$MitoCarta
  if (is.null(mitocarta_anno) && enable_submito) {
    cat("  - 警告：未找到MitoCarta注释，禁用SubMito转化\n")
    enable_submito <- FALSE
  }
  
  # 存储结果
  roc_results <- list()
  all_thresholds <- list()
  data_with_submito_list <- list()
  expr_fdr_df_list <- list()
  
  # 对每个版本执行ROC分析
  for (version in selected_versions) {
    cat(sprintf("\n处理版本: %s\n", version))
    cat("----------------------------------------\n")
    
    # 获取差异分析数据（combined包含Gene + logFC列 + adj.P.Val列 + 注释列）
    diff_data <- diff_results1[[version]]$combined
    
    if (is.null(diff_data) || nrow(diff_data) == 0) {
      cat("  - 警告：数据为空，跳过\n")
      next
    }
    
    # 应用SubMito转化
    cat("步骤 1: SubMito定位转化\n")
    data_transformed <- apply_submito_transformation(
      data = diff_data,
      mitocarta_anno = mitocarta_anno,
      annotation_columns = submito_annotation_columns,
      enable_submito = enable_submito
    )
    
    # 保存转化后的数据
    data_with_submito_list[[version]] <- data_transformed
    output_csv <- file.path(dir_config$output, 
                            paste0("Module08_", version, "_data_with_SubMito.csv"))
    write.csv(data_transformed, output_csv, row.names = FALSE)
    cat(sprintf("  ✓ 已保存: %s\n", basename(output_csv)))
    
    # 生成 Expr_FDR_df（表达矩阵 + logFC/FDR + 注释），供后续背景扣除使用
    expr_matrix <- diff_results1[[version]]$expr_matrix
    if (!is.null(expr_matrix) && "Gene" %in% colnames(expr_matrix)) {
      expr_df <- expr_matrix %>%
        select(Gene, everything()) %>%
        distinct(Gene, .keep_all = TRUE)
      # 构建右表（只保留原始注释列与FC/FDR列，不含SubMito注释）
      anno_cols <- grep("_Localization$", colnames(diff_data), value = TRUE)
      fdr_cols <- grep("(_logFC$|_adj\\.P\\.Val$)", colnames(diff_data), value = TRUE)
      right_core <- diff_data %>% select(Gene, all_of(anno_cols), all_of(fdr_cols))
      # 最终表：Gene + 表达数据 + 原始注释 + FC/FDR（不含SubMito注释）
      expr_fdr_df <- expr_df %>% left_join(right_core, by = "Gene")
      expr_fdr_df_list[[version]] <- expr_fdr_df
      cat("  ✓ 已构建: Expr_FDR_df（表达矩阵 + 原始注释 + FC/FDR，无SubMito注释）\n")
    } else {
      cat("  - 警告：未找到expr_matrix或缺少Gene列，无法构建Expr_FDR_df\n")
    }
    
    # 检查注释列是否存在
    if (!(roc_annotation_column %in% colnames(data_transformed))) {
      cat(sprintf("  - 警告：注释列 %s 不存在，跳过ROC分析\n", roc_annotation_column))
      next
    }
    
    # 获取logFC列
    logfc_columns <- grep("_logFC$", colnames(data_transformed), value = TRUE)
    if (length(logfc_columns) == 0) {
      cat("  - 警告：未找到logFC列，跳过ROC分析\n")
      next
    }
    
    # 执行ROC分析
    cat("\n步骤 2: ROC分析\n")
    roc_analysis <- perform_roc_analysis(
      data = data_transformed,
      logfc_columns = logfc_columns,
      annotation_column = roc_annotation_column,
      tp_label = tp_label,
      fp_label = fp_label,
      direction = "<"
    )
    
    if (length(roc_analysis$roc_objects) == 0) {
      cat("  - 警告：ROC分析失败，跳过\n")
      next
    }
    
    # 计算最佳阈值
    cat("\n步骤 3: 计算最佳阈值\n")
    thresholds <- calculate_optimal_thresholds(roc_analysis$roc_dfs, min_tp = min_tp)
    all_thresholds[[version]] <- thresholds
    cat(sprintf("  - 阈值: %s\n", paste(round(thresholds, 3), collapse = ", ")))
    
    # 保存ROC数据
    cat("\n步骤 4: 保存ROC数据\n")
    output_xlsx <- file.path(dir_config$output, 
                             paste0("Module08_", version, "_ROC_data.xlsx"))
    wb <- createWorkbook()
    for (i in seq_along(roc_analysis$roc_dfs)) {
      sheet_name <- names(roc_analysis$roc_dfs)[i]
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, roc_analysis$roc_dfs[[i]])
    }
    saveWorkbook(wb, output_xlsx, overwrite = TRUE)
    cat(sprintf("  ✓ 已保存: %s\n", basename(output_xlsx)))
    
    # 绘制ROC曲线
    cat("\n步骤 5: 绘制ROC曲线\n")
    output_roc_pdf <- file.path(dir_config$output, 
                                 paste0("Module08_", version, "_ROC_curves.pdf"))
    plot_roc_curves(roc_analysis$roc_objects, roc_analysis$roc_dfs, output_roc_pdf, version)
    
    # 绘制Youden Index图
    cat("\n步骤 6: 绘制Youden Index图\n")
    output_youden_pdf <- file.path(dir_config$output, 
                                    paste0("Module08_", version, "_Youden_Index.pdf"))
    plot_youden_index(roc_analysis$roc_dfs, output_youden_pdf)
    
    # 保存结果
    roc_results[[version]] <- list(
      roc_objects = roc_analysis$roc_objects,
      roc_dfs = roc_analysis$roc_dfs,
      thresholds = thresholds,
      data_transformed = data_transformed
    )
    
    cat(sprintf("\n✓ 版本 %s 完成\n", version))
  }
  
  # 保存所有阈值
  if (length(all_thresholds) > 0) {
    cat("\n保存所有最佳阈值...\n")
    output_thresholds_xlsx <- file.path(dir_config$output, 
                                         "Module08_all_thresholds.xlsx")
    wb <- createWorkbook()
    for (version in names(all_thresholds)) {
      sheet_name <- substr(version, 1, 31)  # Excel限制
      addWorksheet(wb, sheetName = sheet_name)
      
      # 获取该版本的ROC数据框以提取Comparison名称
      version_roc_dfs <- roc_results[[version]]$roc_dfs
      comparison_names <- sapply(version_roc_dfs, function(df) {
        # 从列名中提取Comparison名称（和Youden图用同样的方法）
        tp_fp_col <- grep("_TP_FP$", colnames(df), value = TRUE)[1]
        # 移除"_logFC_TP_FP"后缀得到Comparison名称
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
    cat(sprintf("  ✓ 已保存: %s\n", basename(output_thresholds_xlsx)))
  }
  
  # 统一保存 Expr_FDR_df_list（每个版本一个sheet）
  if (length(expr_fdr_df_list) > 0) {
    cat("\n保存Expr_FDR_df_list...\n")
    expr_fdr_xlsx <- file.path(dir_config$output, "Module08_Expr_FDR_df_list.xlsx")
    wb_expr <- createWorkbook()
    for (version in names(expr_fdr_df_list)) {
      sheet_name <- substr(version, 1, 31)
      addWorksheet(wb_expr, sheetName = sheet_name)
      writeData(wb_expr, sheet = sheet_name, expr_fdr_df_list[[version]])
    }
    saveWorkbook(wb_expr, expr_fdr_xlsx, overwrite = TRUE)
    cat(sprintf("  ✓ 已保存: %s\n", basename(expr_fdr_xlsx)))
  }
  
  cat("\n========================================\n")
  cat("Module 8 完成\n")
  cat("========================================\n")
  
  return(list(
    roc_results = roc_results,
    all_thresholds = all_thresholds,
    data_with_submito = data_with_submito_list,
    expr_fdr_df_list = expr_fdr_df_list
  ))
}

