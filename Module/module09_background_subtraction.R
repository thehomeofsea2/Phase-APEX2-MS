#' Module 09: 背景扣除
#'
#' 根据Module 8的ROC阈值对数据进行过滤，自动生成Method A和Method B两种结果：
#' - Method A：对指定comparison不使用FDR过滤
#' - Method B：对所有comparison使用FDR过滤
#' - Context为"Experiment"的bioGroup：使用阈值过滤
#' - Context为"Spatial"的bioGroup：不使用阈值过滤，只选择列和去除NA
#' - Context为"Control"的bioGroup：不操作
#' - 生成每个bioGroup的统计汇总（count, mean, sum等）
#' - 绘制堆积条形图（数量比例和丰度比例）
#'
#' @param dir_config 目录配置列表（包含output等）
#' @param sampleGroup 样本分组信息（必须包含Context列）
#' @param diff_results 差异分析结果列表（来自Module 7）
#' @param roc_thresholds ROC阈值列表（来自Module 8）
#' @param comparison_info 比较信息列表（来自Module 7，包含comparison到bioGroup的映射）
#' @param data_with_submito Module 8输出的SubMito转化数据列表（可选）
#' @param selected_versions 需要分析的数据版本（NULL表示所有版本）
#' @param fdr_threshold FDR阈值（默认0.05）
#' @param no_fdr_comparisons Method A中不使用FDR过滤的comparison向量（默认NULL，即按Method B处理）
#' @param use_roc_threshold 是否使用ROC threshold（默认TRUE）
#' @param fixed_fc_threshold 如果不使用ROC threshold，使用的固定FC阈值（默认NULL）
#' @param min_valid_lfq 最小有效LFQ值个数（默认2，用于催化组过滤）
#' @param annotation_column 用于分组和作图的注释列（默认"GO_Localization"）
#' @param plot_all_annotations 是否显示所有注释（TRUE）还是只显示TP（FALSE）
#' @param tp_label TP标签（默认"SGs"）
#' @param tp_color TP颜色（默认"#DB6968"）
#'
#' @return 包含filtered_data_A, filtered_data_B, merged_data_A, merged_data_B的列表

module09_background_subtraction <- function(dir_config,
                                             sampleGroup,
                                             diff_results,
                                             roc_thresholds,
                                             comparison_info,
                                             data_with_submito = NULL,
                                             expr_fdr_df_list = NULL,
                                             selected_versions = NULL,
                                             fdr_threshold = 0.05,
                                             no_fdr_comparisons = NULL,
                                             use_roc_threshold = TRUE,
                                             fixed_fc_threshold = NULL,
                                             min_valid_lfq = 2,
                                             annotation_column = "GO_Localization",
                                             plot_all_annotations = FALSE,
                                             tp_label = "SGs",
                                             tp_color = "#DB6968") {
  
  cat("\n=== Module 09: 背景扣除 ===\n")
  
  # 加载必要的包
  if (!require("dplyr")) stop("需要安装dplyr包")
  if (!require("tidyr")) stop("需要安装tidyr包")
  if (!require("purrr")) stop("需要安装purrr包")
  if (!require("ggplot2")) stop("需要安装ggplot2包")
  if (!require("openxlsx")) stop("需要安装openxlsx包")
  
  # 参数验证
  if (!use_roc_threshold && is.null(fixed_fc_threshold)) {
    stop("如果不使用ROC threshold，必须提供 fixed_fc_threshold")
  }
  
  if (!"Context" %in% colnames(sampleGroup)) {
    stop("sampleGroup必须包含Context列")
  }
  
  # 确定要处理的版本
  if (is.null(selected_versions)) {
    selected_versions <- names(diff_results)
  }
  
  # 初始化结果列表
  filtered_data_A <- list()  # Method A结果
  filtered_data_B <- list()  # Method B结果
  merged_data_A <- list()    # Method A合并数据
  merged_data_B <- list()    # Method B合并数据
  
  # 对每个版本进行处理
  for (version in selected_versions) {
    cat(sprintf("\n处理版本: %s\n", version))
    cat(sprintf("  FDR阈值: %.2f\n", fdr_threshold))
    cat(sprintf("  最小有效LFQ: %d\n", min_valid_lfq))
    
    if (!version %in% names(diff_results)) {
      cat(sprintf("  - 警告：版本 %s 不存在于diff_results中，跳过\n", version))
      next
    }
    
    if (!version %in% names(roc_thresholds)) {
      cat(sprintf("  - 警告：版本 %s 不存在于roc_thresholds中，跳过\n", version))
      next
    }
    
    # 获取数据和阈值
    thresholds <- roc_thresholds[[version]]
    
    # 数据优先级：Expr_FDR_df（Module 8） > SubMito转化（Module 8） > Module 7 combined
    if (!is.null(expr_fdr_df_list) && version %in% names(expr_fdr_df_list)) {
      cat("  - 使用Module 8的Expr_FDR_df（表达矩阵 + logFC/FDR + 注释）\n")
      combined_data <- expr_fdr_df_list[[version]]
    } else if (!is.null(data_with_submito) && version %in% names(data_with_submito)) {
      cat("  - 使用Module 8的SubMito转化数据\n")
      combined_data <- data_with_submito[[version]]
    } else {
      cat("  - 提示：未找到Module 8的SubMito数据，使用Module 7的combined结果\n")
      combined_data <- diff_results[[version]]$combined
    }
    
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      cat("  - 警告：数据为空，跳过\n")
      next
    }
    
    # 若未使用Expr_FDR_df，则尝试附加ExprMatrix中的LFQ列
    if (is.null(expr_fdr_df_list) || !(version %in% names(expr_fdr_df_list))) {
      expr_matrix <- diff_results[[version]]$expr_matrix
      if (!is.null(expr_matrix) && "Gene" %in% colnames(expr_matrix)) {
        expr_lfq_cols <- intersect(sampleGroup$FinalName, colnames(expr_matrix))
        expr_df <- expr_matrix %>%
          select(Gene, all_of(expr_lfq_cols)) %>%
          distinct(Gene, .keep_all = TRUE)
        
        if (ncol(expr_df) > 1) {
          combined_data <- combined_data %>%
            left_join(expr_df, by = "Gene")
          cat(sprintf("  - 已附加 %d 个LFQ列\n", ncol(expr_df) - 1))
        } else {
          cat("  - 警告：ExprMatrix未提供匹配的LFQ列\n")
        }
      } else {
        cat("  - 警告：未找到ExprMatrix，无法附加LFQ列\n")
      }
    }
    
    # 识别所有bioGroup及其Context
    biogroups <- unique(sampleGroup$bioGroup)
    cat(sprintf("  发现 %d 个bioGroup\n", length(biogroups)))
    
    # 按Context分组bioGroups（Experiment组在前，Spatial组在后）
    exp_groups <- c()
    spatial_groups <- c()
    control_groups <- c()
    
    for (bg in biogroups) {
      bg_context <- unique(sampleGroup$Context[sampleGroup$bioGroup == bg])
      if (length(bg_context) > 1) {
        cat(sprintf("  - 警告：bioGroup %s 有多个Context，使用第一个: %s\n", bg, bg_context[1]))
        bg_context <- bg_context[1]
      }
      
      if (bg_context == "Experiment") {
        exp_groups <- c(exp_groups, bg)
      } else if (bg_context == "Spatial") {
        spatial_groups <- c(spatial_groups, bg)
      } else if (bg_context == "Control") {
        control_groups <- c(control_groups, bg)
      }
    }
    
    cat(sprintf("  Experiment组: %d 个\n", length(exp_groups)))
    cat(sprintf("  Spatial组: %d 个\n", length(spatial_groups)))
    cat(sprintf("  Control组: %d 个\n", length(control_groups)))
    
    # ========== Method A：部分comparison不用FDR ==========
    cat("\n  === Method A ===\n")
    if (is.null(no_fdr_comparisons) || length(no_fdr_comparisons) == 0) {
      cat("  - 提示：no_fdr_comparisons未指定，Method A将等同于Method B\n")
    } else {
      cat(sprintf("  - 不使用FDR过滤的comparison: %s\n", paste(no_fdr_comparisons, collapse=", ")))
    }
    
    result_A <- process_one_method(
      combined_data, thresholds, sampleGroup, comparison_info,
      exp_groups, spatial_groups, control_groups,
      fdr_threshold, no_fdr_comparisons, use_roc_threshold, 
      fixed_fc_threshold, min_valid_lfq, annotation_column,
      method_name = "A"
    )
    
    # ========== Method B：所有comparison都用FDR ==========
    cat("\n  === Method B ===\n")
    cat("  - 所有comparison都使用FDR过滤\n")
    
    result_B <- process_one_method(
      combined_data, thresholds, sampleGroup, comparison_info,
      exp_groups, spatial_groups, control_groups,
      fdr_threshold, NULL,  # Method B不指定no_fdr_comparisons
      use_roc_threshold, fixed_fc_threshold, min_valid_lfq, annotation_column,
      method_name = "B"
    )
    
    # 保存结果
    filtered_data_A[[version]] <- result_A$filtered_list
    filtered_data_B[[version]] <- result_B$filtered_list
    merged_data_A[[version]] <- result_A$merged_data
    merged_data_B[[version]] <- result_B$merged_data
    
    # 输出Excel文件
    output_results(result_A$filtered_list, version, "A", dir_config$output)
    output_results(result_B$filtered_list, version, "B", dir_config$output)
    # 额外输出：与源码一致的 summarise 与 merged 工作簿，并按组顺序排序
    output_summarise_workbook(result_A$filtered_list, version, "A", dir_config$output, exp_groups, spatial_groups)
    output_summarise_workbook(result_B$filtered_list, version, "B", dir_config$output, exp_groups, spatial_groups)
    output_merged_workbook(result_A$merged_data, version, "A", dir_config$output)
    output_merged_workbook(result_B$merged_data, version, "B", dir_config$output)
    
    # 绘制堆积条形图
    plot_stacked_bar(result_A$filtered_list, result_A$threshold_info,
                     version, "A", fdr_threshold, annotation_column,
                     plot_all_annotations, tp_label, tp_color, dir_config$output)
    
    plot_stacked_bar(result_B$filtered_list, result_B$threshold_info,
                     version, "B", fdr_threshold, annotation_column,
                     plot_all_annotations, tp_label, tp_color, dir_config$output)
  }
  
  cat("\n✓ Module 09 完成\n")
  
  return(list(
    filtered_data_A = filtered_data_A,
    filtered_data_B = filtered_data_B,
    merged_data_A = merged_data_A,
    merged_data_B = merged_data_B
  ))
}


#' 处理单个方法（A或B）
#'
#' @return 包含filtered_list, merged_data, threshold_info的列表

process_one_method <- function(combined_data, thresholds, sampleGroup, comparison_info,
                                exp_groups, spatial_groups, control_groups,
                                fdr_threshold, no_fdr_comparisons,
                                use_roc_threshold, fixed_fc_threshold, min_valid_lfq,
                                annotation_column, method_name) {
  
  filtered_list <- list()  # 用于保存所有过滤后的数据和汇总
  data_list <- list()      # 用于合并的数据列表（只包含过滤后的数据）
  threshold_info <- list() # 保存每个bioGroup使用的阈值信息
  
  # 按顺序处理：先Exp组，再Spatial组
  all_groups <- c(exp_groups, spatial_groups)
  
  for (bg in all_groups) {
    cat(sprintf("    处理 bioGroup: %s\n", bg))
    
    # 获取该bioGroup的Context信息
    bg_context <- unique(sampleGroup$Context[sampleGroup$bioGroup == bg])[1]
    
    # 获取该bioGroup对应的LFQ列名
    lfq_cols <- sampleGroup$FinalName[sampleGroup$bioGroup == bg]
    lfq_cols <- lfq_cols[lfq_cols %in% colnames(combined_data)]
    
    if (length(lfq_cols) == 0) {
      cat(sprintf("      - 警告：未找到对应的LFQ列，跳过\n"))
      next
    }
    
    # 如果是Spatial组，不使用阈值过滤
    if (bg_context == "Spatial") {
      cat(sprintf("      - Spatial组，不使用阈值过滤\n"))
      
      # 选择列：Gene, LFQ列，所有注释列
      annotation_cols <- grep("_Localization$", colnames(combined_data), value = TRUE)
      select_cols <- c("Gene", lfq_cols, annotation_cols)
      select_cols <- select_cols[select_cols %in% colnames(combined_data)]
      
      filtered_data <- combined_data %>%
        select(all_of(select_cols)) %>%
        na.omit()
      
      threshold_info[[bg]] <- "No threshold (Spatial)"
      
    } else if (bg_context == "Experiment") {
      # Experiment组，使用阈值过滤
      cat(sprintf("      - Experiment组，应用阈值过滤\n"))
      
      # 找出该bioGroup涉及的comparisons
      bg_comparisons_idx <- vapply(comparison_info, function(comp) {
        comp$exp_group == bg
      }, logical(1))
      candidate_info <- comparison_info[bg_comparisons_idx]
      
      if (length(candidate_info) == 0) {
        cat(sprintf("      - 警告：未找到对应的comparison，跳过\n"))
        next
      }
      
      bg_comparisons <- vapply(candidate_info, function(comp) {
        if (!is.null(comp$comparison_name) && comp$comparison_name != "") {
          comp$comparison_name
        } else if (!is.null(names(comp))) {
          names(comp)[1]
        } else {
          NA_character_
        }
      }, character(1))
      bg_comparisons <- bg_comparisons[!is.na(bg_comparisons)]
      
      if (length(bg_comparisons) == 0) {
        cat(sprintf("      - 警告：未找到对应的comparison，跳过\n"))
        next
      }
      
      cat(sprintf("      - 找到 %d 个comparison\n", length(bg_comparisons)))
      
      # 构建过滤条件
      filtered_data <- combined_data
      threshold_text <- c()
      
      for (comp in bg_comparisons) {
        logfc_col <- paste0(comp, "_logFC")
        fdr_col <- paste0(comp, "_adj.P.Val")
        
        if (!logfc_col %in% colnames(filtered_data) || !fdr_col %in% colnames(filtered_data)) {
          cat(sprintf("        - 警告：comparison %s 的列不存在，跳过\n", comp))
          next
        }
        
        # 确定使用的阈值
        if (use_roc_threshold) {
          # 使用ROC threshold
          threshold <- thresholds[comp]
          if (is.na(threshold)) {
            cat(sprintf("        - 警告：未找到comparison %s 的ROC阈值，跳过\n", comp))
            next
          }
        } else {
          # 使用固定FC阈值
          threshold <- fixed_fc_threshold
        }
        
        # 确定FDR阈值
        if (!is.null(no_fdr_comparisons) && comp %in% no_fdr_comparisons) {
          fdr_cutoff <- 1  # 不使用FDR过滤
          fdr_text <- "no FDR filter"
        } else {
          fdr_cutoff <- fdr_threshold
          fdr_text <- sprintf("FDR<%.2f", fdr_cutoff)
        }
        
        # 应用过滤
        filtered_data <- filtered_data %>%
          filter(!!sym(logfc_col) > threshold,
                 !!sym(fdr_col) < fdr_cutoff)
        
        threshold_text <- c(threshold_text, 
                            sprintf("%s: logFC>%.2f, %s", comp, threshold, fdr_text))
      }
      
      # 选择列：Gene, LFQ列，logFC/FDR列，所有注释列
      annotation_cols <- grep("_Localization$", colnames(combined_data), value = TRUE)
      logfc_cols <- paste0(bg_comparisons, "_logFC")
      fdr_cols <- paste0(bg_comparisons, "_adj.P.Val")
      select_cols <- c("Gene", lfq_cols, logfc_cols, fdr_cols, annotation_cols)
      select_cols <- select_cols[select_cols %in% colnames(filtered_data)]
      
      filtered_data <- filtered_data %>%
        select(all_of(select_cols))
      
      # 催化组有效值过滤（LFQ列中至少有min_valid_lfq个非NA值）
      if (length(lfq_cols) >= min_valid_lfq) {
        filtered_data <- filtered_data %>%
          filter(rowSums(!is.na(select(., all_of(lfq_cols)))) >= min_valid_lfq)
      }
      
      threshold_info[[bg]] <- paste(threshold_text, collapse = "; ")
      
    } else {
      # Control组，不操作
      cat(sprintf("      - Control组，跳过\n"))
      next
    }
    
    cat(sprintf("      - 过滤后保留 %d 个蛋白\n", nrow(filtered_data)))
    
    # 检查最小值（验证过滤是否成功）
    if (bg_context == "Experiment") {
      for (comp in bg_comparisons) {
        logfc_col <- paste0(comp, "_logFC")
        if (logfc_col %in% colnames(filtered_data)) {
          min_logfc <- min(filtered_data[[logfc_col]], na.rm = TRUE)
          cat(sprintf("      - %s 最小logFC: %.2f\n", comp, min_logfc))
        }
      }
    }
    
    # 保存过滤后的数据（先放数据）
    filtered_list[[bg]] <- filtered_data
    data_list[[bg]] <- filtered_data
    
    # 生成统计汇总（按annotation_column分组）
    if (annotation_column %in% colnames(filtered_data)) {
      summarise_data <- filtered_data %>%
        group_by(!!sym(annotation_column)) %>%
        summarise(
          count = n(),
          across(
            all_of(lfq_cols),
            list(
              mean = ~mean(.x, na.rm = TRUE),
              sum = ~sum(2^.x, na.rm = TRUE)  # 转回线性空间再求和
            ),
            .names = "{.col}_{.fn}"
          ),
          .groups = "drop"
        )
      
      # 计算百分比
      summarise_data <- summarise_data %>%
        mutate(Percent = count / sum(count))
      
      # 计算丰度百分比（每个LFQ列的sum百分比，然后取平均）
      sum_cols <- grep("_sum$", colnames(summarise_data), value = TRUE)
      if (length(sum_cols) > 0) {
        for (col in sum_cols) {
          percent_col <- paste0(col, "_percent")
          summarise_data[[percent_col]] <- summarise_data[[col]] / sum(summarise_data[[col]], na.rm = TRUE)
        }
        
        # 计算MeanSum（所有sum_percent的平均值）
        percent_cols <- grep("_sum_percent$", colnames(summarise_data), value = TRUE)
        if (length(percent_cols) > 0) {
          summarise_data <- summarise_data %>%
            mutate(MeanSum = rowMeans(select(., all_of(percent_cols)), na.rm = TRUE))
        }
      }
      
      # 保存统计汇总（后放汇总）
      summarise_name <- paste0(bg, "_Summarise")
      filtered_list[[summarise_name]] <- summarise_data
      cat(sprintf("      - 生成统计汇总：%d 个分组\n", nrow(summarise_data)))
    }
  }
  
  # 合并所有过滤后的数据（AfterROC_merged逻辑）
  merged_data <- NULL
  if (length(data_list) > 0) {
    # 提取所有注释列（从原始combined_data）
    annotation_cols <- grep("_Localization$", colnames(combined_data), value = TRUE)
    anno_data <- combined_data %>% 
      select(Gene, all_of(annotation_cols)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    # full_join合并所有过滤后的数据
    merged_data <- reduce(data_list, full_join, by = "Gene")
    
    # 只保留Gene和LFQ列
    lfq_cols_all <- c()
    for (bg in names(data_list)) {
      bg_lfq <- sampleGroup$FinalName[sampleGroup$bioGroup == bg]
      lfq_cols_all <- c(lfq_cols_all, bg_lfq)
    }
    lfq_cols_all <- lfq_cols_all[lfq_cols_all %in% colnames(merged_data)]
    
    merged_data <- merged_data %>%
      select(Gene, all_of(lfq_cols_all))
    
    # 添加注释列
    merged_data <- left_join(merged_data, anno_data, by = "Gene")
    
    cat(sprintf("    - 合并数据：%d 个蛋白, %d 列\n", nrow(merged_data), ncol(merged_data)))
  }
  
  return(list(
    filtered_list = filtered_list,
    merged_data = merged_data,
    threshold_info = threshold_info
  ))
}


#' 输出Excel文件
#'
#' @param filtered_list 过滤后的数据列表（包含数据和汇总）
#' @param version 数据版本
#' @param method 方法名称（A或B）
#' @param output_dir 输出目录

output_results <- function(filtered_list, version, method, output_dir) {
  if (length(filtered_list) == 0) {
    cat(sprintf("    - Method %s：没有数据输出\n", method))
    return(invisible(NULL))
  }
  
  output_xlsx <- file.path(output_dir, 
                           sprintf("Module09_%s_%s_FilteredData.xlsx", method, version))
  wb <- createWorkbook()
  
  # 按照源码顺序：先Exp组数据和汇总，再Spatial组数据和汇总
  for (name in names(filtered_list)) {
    sheet_name <- substr(name, 1, 31)  # Excel限制
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, filtered_list[[name]])
  }
  
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  cat(sprintf("    ✓ Method %s: 已保存 %s\n", method, basename(output_xlsx)))
}


#' 额外输出：按Exp→Spatial、数据→汇总顺序保存Summarise工作簿
#'
#' @param filtered_list 过滤后的数据与汇总列表
#' @param version 版本
#' @param method 方法（A/B）
#' @param output_dir 输出目录
#' @param exp_groups 实验组bioGroup向量
#' @param spatial_groups 空间组bioGroup向量
output_summarise_workbook <- function(filtered_list, version, method, output_dir, exp_groups, spatial_groups) {
  if (length(filtered_list) == 0) {
    return(invisible(NULL))
  }
  ordered_names <- c()
  # 1) 先所有“数据”sheet：Exp → Spatial
  for (bg in exp_groups) {
    if (bg %in% names(filtered_list)) ordered_names <- c(ordered_names, bg)
  }
  for (bg in spatial_groups) {
    if (bg %in% names(filtered_list)) ordered_names <- c(ordered_names, bg)
  }
  # 2) 再所有“Summarise”sheet：Exp → Spatial
  for (bg in exp_groups) {
    sum_name <- paste0(bg, "_Summarise")
    if (sum_name %in% names(filtered_list)) ordered_names <- c(ordered_names, sum_name)
  }
  for (bg in spatial_groups) {
    sum_name <- paste0(bg, "_Summarise")
    if (sum_name %in% names(filtered_list)) ordered_names <- c(ordered_names, sum_name)
  }
  # 回退：若有遗漏，追加其余项
  remaining <- setdiff(names(filtered_list), ordered_names)
  ordered_names <- c(ordered_names, remaining)
  
  output_xlsx <- file.path(output_dir,
                           sprintf("Module09_%s_%s_AfterROC_data_summarise.xlsx", method, version))
  wb <- createWorkbook()
  for (name in ordered_names) {
    sheet_name <- substr(name, 1, 31)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, filtered_list[[name]])
  }
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  cat(sprintf("    ✓ Method %s: 已保存 %s\n", method, basename(output_xlsx)))
}


#' 额外输出：合并后的merged数据工作簿
#'
#' @param merged_data 合并后的数据框
#' @param version 版本
#' @param method 方法（A/B）
#' @param output_dir 输出目录
output_merged_workbook <- function(merged_data, version, method, output_dir) {
  if (is.null(merged_data) || nrow(merged_data) == 0) {
    return(invisible(NULL))
  }
  output_xlsx <- file.path(output_dir,
                           sprintf("Module09_%s_%s_AfterROC_data_merged.xlsx", method, version))
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = substr(version, 1, 31))
  writeData(wb, sheet = substr(version, 1, 31), merged_data)
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  cat(sprintf("    ✓ Method %s: 已保存 %s\n", method, basename(output_xlsx)))
}


#' 绘制堆积条形图
#'
#' @param filtered_list 过滤后的数据列表（包含汇总数据）
#' @param threshold_info 阈值信息列表
#' @param version 数据版本
#' @param method 方法名称
#' @param fdr_threshold FDR阈值
#' @param annotation_column 注释列名
#' @param plot_all_annotations 是否显示所有注释
#' @param tp_label TP标签
#' @param tp_color TP颜色
#' @param output_dir 输出目录

plot_stacked_bar <- function(filtered_list, threshold_info, version, method,
                             fdr_threshold, annotation_column,
                             plot_all_annotations, tp_label, tp_color, output_dir) {
  
  # 提取汇总数据（名称包含"_Summarise"的）
  summarise_names <- grep("_Summarise$", names(filtered_list), value = TRUE)
  
  if (length(summarise_names) == 0) {
    cat(sprintf("    - Method %s：没有汇总数据用于作图\n", method))
    return(invisible(NULL))
  }
  
  # 准备作图数据
  plot_data_list <- list()
  
  for (name in summarise_names) {
    bg <- sub("_Summarise$", "", name)
    data <- filtered_list[[name]]
    
    if (!annotation_column %in% colnames(data)) {
      next
    }
    
    # 添加bioGroup列
    data <- data %>%
      mutate(bioGroup = bg,
             bioGroup_factor = factor(bg, levels = sub("_Summarise$", "", summarise_names)))
    
    plot_data_list[[bg]] <- data
  }
  
  if (length(plot_data_list) == 0) {
    cat(sprintf("    - Method %s：没有有效数据用于作图\n", method))
    return(invisible(NULL))
  }
  
  # 合并数据
  merged_data <- bind_rows(plot_data_list)
  
  # 准备阈值文本（用于图上方显示）
  threshold_text <- sapply(names(threshold_info), function(bg) {
    paste0(bg, ": ", threshold_info[[bg]])
  })
  threshold_subtitle <- paste(threshold_text, collapse = "\n")
  
  # 提取TP数据（用于标注）
  tp_data <- merged_data %>%
    filter(!!sym(annotation_column) == tp_label)
  
  # 设置颜色
  if (plot_all_annotations) {
    # 显示所有注释，TP为红色，其他使用调色板
    all_annotations <- unique(merged_data[[annotation_column]])
    fill_levels <- c(tp_label, setdiff(all_annotations, tp_label))
    color_values <- setNames(
      c(tp_color, scales::hue_pal()(length(all_annotations) - 1)),
      fill_levels
    )
    
    # 调整因子顺序，使TP在底部
    merged_data[[annotation_column]] <- factor(
      merged_data[[annotation_column]],
      levels = fill_levels
    )
  } else {
    # 只显示TP，其他为灰色
    all_annotations <- unique(merged_data[[annotation_column]])
    fill_levels <- c(tp_label, setdiff(all_annotations, tp_label))
    color_values <- setNames(
      c(tp_color, rep("grey", length(all_annotations) - 1)),
      fill_levels
    )
    
    # 调整因子顺序，使TP在底部
    merged_data[[annotation_column]] <- factor(
      merged_data[[annotation_column]],
      levels = fill_levels
    )
  }
  
  # 明确堆叠顺序：TP在底部（order越小越靠下）
  merged_data <- merged_data %>%
    mutate(.stack_order = ifelse(!!sym(annotation_column) == tp_label, 0L, 1L))
  
  # 绘制数量比例图
  output_file <- file.path(output_dir, 
                           sprintf("Module09_%s_%s_CountPercent_StackBar.pdf", 
                                   method, version))
  pdf(output_file, width = 10, height = 7)
  
  p <- ggplot(merged_data, aes(x = bioGroup_factor, y = Percent, fill = !!sym(annotation_column), order = .stack_order)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_text(data = tp_data,
              aes(x = bioGroup_factor, y = 1,
                  label = sprintf("%s:\n%d\n(%.1f%%)", 
                                  tp_label, count, Percent * 100)),
              inherit.aes = FALSE,
              vjust = -0.25,
              size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(title = sprintf("%s - Count Percent (Method %s, FDR<%.2f)", 
                         version, method, fdr_threshold),
         subtitle = threshold_subtitle,
         y = "Percent",
         x = "bioGroup") +
    scale_fill_manual(values = color_values, breaks = fill_levels, limits = fill_levels) +
    guides(fill = guide_legend(title = sprintf("Localization\n(%s)", annotation_column))) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0, size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  print(p)
  dev.off()
  cat(sprintf("    ✓ Method %s: 已保存 %s\n", method, basename(output_file)))
  
  # 绘制丰度比例图（如果有MeanSum列）
  if ("MeanSum" %in% colnames(merged_data)) {
    output_file <- file.path(output_dir, 
                             sprintf("Module09_%s_%s_AbundancePercent_StackBar.pdf", 
                                     method, version))
    pdf(output_file, width = 10, height = 7)
    
    tp_data_abundance <- merged_data %>%
      filter(!!sym(annotation_column) == tp_label)
    
    p <- ggplot(merged_data, aes(x = bioGroup_factor, y = MeanSum, fill = !!sym(annotation_column), order = .stack_order)) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(data = tp_data_abundance,
                aes(x = bioGroup_factor, y = 1,
                    label = sprintf("%s:\n%d\n(%.1f%%)", 
                                    tp_label, count, MeanSum * 100)),
                inherit.aes = FALSE,
                vjust = -0.25,
                size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(title = sprintf("%s - Abundance Percent (Method %s, FDR<%.2f)", 
                           version, method, fdr_threshold),
           subtitle = threshold_subtitle,
           y = "Percent Abundance",
           x = "bioGroup") +
      scale_fill_manual(values = color_values, breaks = fill_levels, limits = fill_levels) +
      guides(fill = guide_legend(title = sprintf("Localization\n(%s)", annotation_column))) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    print(p)
    dev.off()
    cat(sprintf("    ✓ Method %s: 已保存 %s\n", method, basename(output_file)))
  }
  
  return(invisible(NULL))
}
