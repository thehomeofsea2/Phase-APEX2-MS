# ============================================================================
# Module 7: 第一次差异分析
# ============================================================================
# 功能：
#   1. 基于sampleGroup$FirstROCgroup自动构建对比组
#   2. 使用limma包进行差异表达分析
#   3. 使用eBayes经验贝叶斯调整
#   4. 提取logFC和adj.P.Val（FDR校正）
#   5. 输出合并的差异分析结果和表达矩阵
#
# FirstROCgroup逻辑：
#   - 包含相同元素的组为同一组（如A、A/B、A&B都属于A组）
#   - 在同一组内，Context为Experiment的bioGroup vs Context为Control的bioGroup
#   - 例如：
#     * FirstROCgroup=A, Context=Experiment：K69A1B3_Light
#     * FirstROCgroup=A, Context=Control：K69A1B3_noLight, A1B3_Light
#     * FirstROCgroup=A&B, Context=Control：K69_Light（同时属于A组和B组）
#     * 对比：K69A1B3_Light vs K69A1B3_noLight
#             K69A1B3_Light vs K69_Light
#             K69A1B3_Light vs A1B3_Light
#
# 输入：
#   - dir_config: 目录配置（必须包含output）
#   - imputed_data_list: Module05输出的填补后数据列表
#   - sampleGroup: 样本分组信息（必须包含FinalName, bioGroup, FirstROCgroup, Context）
#   - selected_versions: 需要分析的数据版本（默认NULL表示所有版本）
#
# 输出：
#   - Module07_workspace.RData: 包含Module01-06所有数据 + diff_results1
#   - CSV文件（Output目录）：
#     - Module07_<version>_DiffAnalysis.csv - 合并的差异分析结果
#     - Module07_<version>_ExprMatrix.csv - 表达矩阵（Gene × Samples）
#     - Module07_<version>_LogFC_Wide.csv - LogFC宽表（Gene × Comparisons）
#   - XLSX文件（Output目录）：
#     - Module07_<version>_DiffAnalysis_Full.xlsx - 完整的topTable结果（多sheet）
#
# 依赖包：limma, dplyr, tidyr, openxlsx

module07_diff_analysis1 <- function(dir_config,
                                     imputed_data_list,
                                     sampleGroup,
                                     selected_versions = NULL) {
  
  # 加载必需的包
  if (!require("limma", quietly = TRUE)) {
    stop("需要安装limma包: BiocManager::install('limma')")
  }
  if (!require("dplyr", quietly = TRUE)) {
    stop("需要安装dplyr包: install.packages('dplyr')")
  }
  if (!require("tidyr", quietly = TRUE)) {
    stop("需要安装tidyr包: install.packages('tidyr')")
  }
  if (!require("openxlsx", quietly = TRUE)) {
    stop("需要安装openxlsx包: install.packages('openxlsx')")
  }
  
  cat("\n========================================\n")
  cat("Module 7: 第一次差异分析\n")
  cat("========================================\n")
  
  # 验证输入
  if (!("output" %in% names(dir_config))) {
    stop("dir_config必须包含'output'路径")
  }
  
  if (is.null(imputed_data_list) || length(imputed_data_list) == 0) {
    stop("imputed_data_list不能为空")
  }
  
  # 验证sampleGroup必需列
  required_cols <- c("FinalName", "bioGroup", "FirstROCgroup", "Context")
  missing_cols <- setdiff(required_cols, colnames(sampleGroup))
  if (length(missing_cols) > 0) {
    stop(paste("sampleGroup缺少必需列:", paste(missing_cols, collapse = ", ")))
  }
  
  # 如果未指定版本，使用所有版本
  if (is.null(selected_versions)) {
    selected_versions <- names(imputed_data_list)
  }
  
  # 验证选择的版本是否存在
  missing_versions <- setdiff(selected_versions, names(imputed_data_list))
  if (length(missing_versions) > 0) {
    stop(paste("以下版本不存在:", paste(missing_versions, collapse = ", ")))
  }
  
  cat(sprintf("  - 分析版本: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("  - 样本数: %d\n", nrow(sampleGroup)))
  
  # ========================================
  # 步骤1：解析FirstROCgroup，构建对比组
  # ========================================
  cat("\n步骤1：解析FirstROCgroup并构建对比组...\n")
  
  # 过滤掉FirstROCgroup为NA或空的样本
  sampleGroup_filtered <- sampleGroup %>%
    filter(!is.na(FirstROCgroup) & FirstROCgroup != "")
  
  cat(sprintf("  - 参与第一次差异分析的样本数: %d\n", nrow(sampleGroup_filtered)))
  
  # 解析FirstROCgroup，提取所有独立元素
  # 例如："A&B" 或 "A/B" -> c("A", "B")
  parse_roc_group <- function(roc_str) {
    if (is.na(roc_str) || roc_str == "") return(character(0))
    # 分割符可能是 &, /, 或其他
    elements <- unlist(strsplit(roc_str, "[&/,;|]+"))
    # 去除空格
    elements <- trimws(elements)
    # 去除空元素
    elements <- elements[elements != ""]
    return(unique(elements))
  }
  
  # 为每个样本提取FirstROCgroup元素
  sampleGroup_filtered$ROCgroups <- lapply(sampleGroup_filtered$FirstROCgroup, parse_roc_group)
  
  # 获取所有唯一的ROC组元素
  all_roc_elements <- unique(unlist(sampleGroup_filtered$ROCgroups))
  cat(sprintf("  - 检测到的ROC组: %s\n", paste(all_roc_elements, collapse = ", ")))
  
  # 为每个ROC组构建对比
  comparisons <- list()
  comparison_info <- list()
  
  for (roc_element in all_roc_elements) {
    cat(sprintf("\n  处理ROC组: %s\n", roc_element))
    
    # 找出属于该ROC组的所有bioGroup
    belongs_to_group <- sapply(sampleGroup_filtered$ROCgroups, function(x) roc_element %in% x)
    group_samples <- sampleGroup_filtered[belongs_to_group, ]
    
    cat(sprintf("    - 该组内样本数: %d\n", nrow(group_samples)))
    cat(sprintf("    - 该组内bioGroup: %s\n", paste(unique(group_samples$bioGroup), collapse = ", ")))
    
    # 找出Experiment组（Context为Experiment）
    exp_groups <- group_samples %>%
      filter(Context == "Experiment") %>%
      pull(bioGroup) %>%
      unique()
    
    # 找出Control组（Context为Control）
    ctrl_groups <- group_samples %>%
      filter(Context == "Control") %>%
      pull(bioGroup) %>%
      unique()
    
    cat(sprintf("    - Experiment组: %s\n", 
                ifelse(length(exp_groups) > 0, paste(exp_groups, collapse = ", "), "无")))
    cat(sprintf("    - Control组: %s\n", 
                ifelse(length(ctrl_groups) > 0, paste(ctrl_groups, collapse = ", "), "无")))
    
    # 构建对比：每个Experiment组 vs 每个Control组
    if (length(exp_groups) > 0 && length(ctrl_groups) > 0) {
      for (exp_group in exp_groups) {
        for (ctrl_group in ctrl_groups) {
          comp <- c(exp_group, ctrl_group)
          comparisons[[length(comparisons) + 1]] <- comp
          
          # 保存对比信息
          comparison_info[[length(comparison_info) + 1]] <- list(
            roc_group = roc_element,
            exp_group = exp_group,
            ctrl_group = ctrl_group,
            comparison_name = paste0(make.names(exp_group), "_vs_", make.names(ctrl_group))
          )
          
          cat(sprintf("    - 添加对比: %s vs %s\n", exp_group, ctrl_group))
        }
      }
    } else {
      cat(sprintf("    - 警告：ROC组 %s 缺少Experiment或Control组，跳过\n", roc_element))
    }
  }
  
  if (length(comparisons) == 0) {
    stop("未能构建任何对比组，请检查sampleGroup的FirstROCgroup和Context列")
  }
  
  cat(sprintf("\n✓ 总共构建了 %d 个对比组\n", length(comparisons)))
  
  # 存储所有版本的差异分析结果
  diff_results1 <- list()
  
  # ========================================
  # 步骤2：对每个版本进行差异分析
  # ========================================
  for (version in selected_versions) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("分析版本: %s\n", version))
    cat(sprintf("========================================\n"))
    
    # 获取数据
    data_imputed <- imputed_data_list[[version]]
    
    # 识别数据列和注释列
    data_cols <- intersect(sampleGroup$FinalName, colnames(data_imputed))
    anno_cols <- setdiff(colnames(data_imputed), data_cols)
    
    if (length(data_cols) == 0) {
      warning(sprintf("版本 %s 没有找到匹配的数据列，跳过", version))
      next
    }
    
    cat(sprintf("  - 数据列数: %d\n", length(data_cols)))
    cat(sprintf("  - 注释列数: %d\n", length(anno_cols)))
    
    # 提取表达矩阵（只保留数据列）
    expr_matrix <- as.matrix(data_imputed[, data_cols, drop = FALSE])
    
    # 设置行名为Gene（假设有Gene列）
    if ("Gene" %in% anno_cols) {
      rownames(expr_matrix) <- data_imputed$Gene
      gene_col <- data_imputed$Gene
    } else {
      warning("数据中没有Gene列，使用行号作为行名")
      rownames(expr_matrix) <- paste0("Protein_", 1:nrow(expr_matrix))
      gene_col <- rownames(expr_matrix)
    }
    
    # 输出表达矩阵
    expr_matrix_df <- data.frame(Gene = gene_col, expr_matrix, check.names = FALSE)
    expr_matrix_file <- file.path(dir_config$output, 
                                   paste0("Module07_", version, "_ExprMatrix.csv"))
    write.csv(expr_matrix_df, expr_matrix_file, row.names = FALSE)
    cat(sprintf("✓ 已保存表达矩阵: %s\n", basename(expr_matrix_file)))
    
    # 确保sampleGroup的顺序与expr_matrix列顺序一致
    sampleGroup_ordered <- sampleGroup[match(data_cols, sampleGroup$FinalName), ]
    
    # 构建设计矩阵
    cat("  - 构建设计矩阵...\n")
    Group <- factor(sampleGroup_ordered$bioGroup)
    design <- model.matrix(~ 0 + Group)
    colnames(design) <- levels(Group)
    
    # 拟合线性模型
    cat("  - 拟合线性模型...\n")
    fit <- limma::lmFit(expr_matrix, design)
    
    # 构建对比矩阵
    cat("  - 构建对比矩阵...\n")
    
    # 获取设计矩阵中实际存在的bioGroup
    available_groups <- colnames(design)
    cat(sprintf("  - 设计矩阵中的bioGroup: %s\n", paste(available_groups, collapse = ", ")))
    
    # 过滤掉不在设计矩阵中的对比组
    valid_comparisons <- list()
    contrast_strings <- character()
    contrast_names <- character()
    
    for (i in seq_along(comparisons)) {
      comp <- comparisons[[i]]
      group1 <- comp[1]
      group2 <- comp[2]
      
      # 检查两个组是否都在设计矩阵中
      group1_safe <- make.names(group1)
      group2_safe <- make.names(group2)
      
      if (!(group1_safe %in% available_groups)) {
        cat(sprintf("  - 警告：跳过对比 %s vs %s（%s不在数据中）\n", group1, group2, group1))
        next
      }
      if (!(group2_safe %in% available_groups)) {
        cat(sprintf("  - 警告：跳过对比 %s vs %s（%s不在数据中）\n", group1, group2, group2))
        next
      }
      
      # 创建对比名称（GroupA_vs_GroupB）
      contrast_name <- paste0(group1_safe, "_vs_", group2_safe)
      contrast_names <- c(contrast_names, contrast_name)
      
      # 创建对比字符串（GroupA - GroupB）
      contrast_strings <- c(contrast_strings, paste0(group1_safe, " - ", group2_safe))
      
      # 保存有效对比
      valid_comparisons[[length(valid_comparisons) + 1]] <- comp
    }
    
    if (length(contrast_names) == 0) {
      warning(sprintf("版本 %s 没有有效的对比组，跳过", version))
      next
    }
    
    names(contrast_strings) <- contrast_names
    
    # 构建对比矩阵
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_strings,
      levels = design
    )
    
    cat(sprintf("  - 有效对比组数: %d\n", length(contrast_names)))
    for (name in contrast_names) {
      cat(sprintf("    - %s\n", name))
    }
    
    # 计算对比并进行经验贝叶斯调整
    cat("  - 经验贝叶斯调整...\n")
    fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
    fit_ebayes <- limma::eBayes(fit_contrasts)
    
    # 提取结果
    cat("  - 提取差异分析结果...\n")
    
    # 获取对比矩阵的实际列名（limma生成的名称）
    actual_coef_names <- colnames(contrast_matrix)
    cat(sprintf("  - 对比矩阵列名: %s\n", paste(actual_coef_names, collapse = ", ")))
    
    # 存储每个对比的结果
    FDR_test_list <- list()
    Raw_FDR_test_list <- list()
    LogFC_list <- list()  # 用于构建LogFC宽表
    successful_comparisons_list <- list()  # 保存成功的对比
    
    for (i in seq_along(contrast_names)) {
      contrast_name <- contrast_names[i]
      actual_coef_name <- actual_coef_names[i]  # 使用对比矩阵的实际列名
      
      # 提取topTable结果（使用实际的系数名称）
      tryCatch({
        top_table <- limma::topTable(fit_ebayes, 
                                     coef = actual_coef_name,  # 使用实际的系数名称
                                     n = Inf, 
                                     adjust.method = "BH")
        
        # 添加Gene列
        top_table$Gene <- rownames(top_table)
        
        # 保存完整结果（使用友好的名称作为键）
        Raw_FDR_test_list[[contrast_name]] <- top_table
        
        # 提取logFC和adj.P.Val
        result_df <- top_table %>%
          select(Gene, logFC, adj.P.Val) %>%
          rename(!!paste0(contrast_name, "_logFC") := logFC,
                 !!paste0(contrast_name, "_adj.P.Val") := adj.P.Val)
        
        FDR_test_list[[contrast_name]] <- result_df
        
        # 提取logFC用于宽表
        logfc_df <- top_table %>%
          select(Gene, logFC) %>%
          rename(!!contrast_name := logFC)
        
        LogFC_list[[contrast_name]] <- logfc_df
        
        # 保存成功的对比
        successful_comparisons_list[[contrast_name]] <- valid_comparisons[[i]]
        
        cat(sprintf("  ✓ 完成对比: %s\n", contrast_name))
        
      }, error = function(e) {
        cat(sprintf("  - 错误：对比 %s 失败 - %s\n", contrast_name, e$message))
      })
    }
    
    # 检查是否有成功的对比
    if (length(FDR_test_list) == 0) {
      warning(sprintf("版本 %s 没有任何成功的对比，跳过", version))
      next
    }
    
    # 合并所有对比的结果
    cat(sprintf("  - 合并对比结果（成功 %d/%d 个对比）...\n", 
                length(FDR_test_list), length(contrast_names)))
    FDR_combined_df <- Reduce(function(x, y) full_join(x, y, by = "Gene"), FDR_test_list)
    
    # 添加注释列（如果存在）
    if (length(anno_cols) > 0) {
      anno_df <- data_imputed %>% select(all_of(anno_cols))
      FDR_combined_df <- FDR_combined_df %>%
        left_join(anno_df, by = "Gene")
    }
    
    # 保存合并结果为CSV
    csv_file <- file.path(dir_config$output, 
                          paste0("Module07_", version, "_DiffAnalysis.csv"))
    write.csv(FDR_combined_df, csv_file, row.names = FALSE)
    cat(sprintf("✓ 已保存合并结果: %s\n", basename(csv_file)))
    
    # 合并LogFC宽表
    if (length(LogFC_list) > 0) {
      LogFC_wide <- Reduce(function(x, y) full_join(x, y, by = "Gene"), LogFC_list)
    } else {
      LogFC_wide <- data.frame(Gene = gene_col)
    }
    logfc_wide_file <- file.path(dir_config$output, 
                                  paste0("Module07_", version, "_LogFC_Wide.csv"))
    write.csv(LogFC_wide, logfc_wide_file, row.names = FALSE)
    cat(sprintf("✓ 已保存LogFC宽表: %s\n", basename(logfc_wide_file)))
    
    # 保存完整结果为XLSX（多sheet）
    xlsx_file <- file.path(dir_config$output, 
                           paste0("Module07_", version, "_DiffAnalysis_Full.xlsx"))
    wb <- openxlsx::createWorkbook()
    
    # 添加合并结果sheet
    openxlsx::addWorksheet(wb, "Combined")
    openxlsx::writeData(wb, "Combined", FDR_combined_df)
    
    # 添加LogFC宽表sheet
    openxlsx::addWorksheet(wb, "LogFC_Wide")
    openxlsx::writeData(wb, "LogFC_Wide", LogFC_wide)
    
    # 添加表达矩阵sheet
    openxlsx::addWorksheet(wb, "ExprMatrix")
    openxlsx::writeData(wb, "ExprMatrix", expr_matrix_df)
    
    # 添加每个对比的完整结果sheet
    for (contrast_name in names(Raw_FDR_test_list)) {
      sheet_name <- substr(contrast_name, 1, 31)  # Excel sheet名称限制31字符
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, Raw_FDR_test_list[[contrast_name]])
    }
    
    openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
    cat(sprintf("✓ 已保存完整结果: %s\n", basename(xlsx_file)))
    
    # 存储结果
    diff_results1[[version]] <- list(
      combined = FDR_combined_df,
      logfc_wide = LogFC_wide,
      expr_matrix = expr_matrix_df,
      raw_results = Raw_FDR_test_list,
      comparisons = successful_comparisons_list,  # 只保存成功的对比
      comparison_info = comparison_info,
      contrast_names = names(Raw_FDR_test_list)  # 只保存成功的对比名称
    )
  }
  
  cat("\n✓ Module 7 完成\n")
  cat(sprintf("  - 分析版本数: %d\n", length(diff_results1)))
  cat(sprintf("  - 总对比组数: %d\n", length(comparisons)))
  
  # 返回结果
  return(list(
    diff_results1 = diff_results1,
    comparisons_used = comparisons,  # 返回所有构建的对比（可能某些版本中不全部有效）
    comparison_info = comparison_info
  ))
}
