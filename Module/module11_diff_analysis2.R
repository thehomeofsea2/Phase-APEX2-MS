# ============================================================================
# Module 11: 第二次差异分析（基于Module 10的数据替换结果）
# ============================================================================
# 对齐 CleanCode.R (2377-2491) 逻辑
#
# 目标：
# - 对 Module 10 生成的 replaced_data_A / replaced_data_B
#   只选择 Context 为 Experiment 或 Spatial 的样本（去除 Control）
# - 根据 SecondROCgroup 和 PLtype 构建比较组：
#   * Exp vs Exp：所有 Experiment 组之间的比较（不限PLtype）
#   * Exp vs Spatial：只比较 PLtype 相同的组
# - 使用 limma 进行差异分析，得到 logFC 和 adj.P.Val
# - 生成 FDR_combined_df_list_2nd（合并所有比较结果）
#
# 输入：
# - dir_config：包含 output 等目录路径
# - sampleGroup：用于识别样本分组、Context、PLtype 和 SecondROCgroup
# - replaced_data_A / replaced_data_B：Module 10 输出的数据替换结果
# - selected_versions：需要处理的版本（NULL 表示全部）
#
# 输出：
# - FDR_combined_df_list_2nd：包含所有版本和方法的差异分析结果
# - 每个数据集对应一个 Excel 文件（包含原始 topTable 结果）
# - 合并的 Excel 文件（包含所有版本的 logFC + adj.P.Val + 注释列）
#
# ============================================================================

module11_diff_analysis2 <- function(
  dir_config,
  sampleGroup,
  replaced_data_A,
  replaced_data_B,
  selected_versions = NULL
) {
  
  cat("\n=== Module 11: 第二次差异分析 ===\n\n")
  
  # 加载必需包
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("✗ 错误：需要安装 limma 包\n  BiocManager::install('limma')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("✗ 错误：需要安装 dplyr 包")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("✗ 错误：需要安装 openxlsx 包")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("✗ 错误：需要安装 purrr 包")
  }
  
  library(limma)
  library(dplyr)
  library(openxlsx)
  library(purrr)
  
  # ----------------------------------------------------------------------
  # 1. 准备数据列表
  # ----------------------------------------------------------------------
  
  # 给 A 和 B 方法的数据添加后缀
  names(replaced_data_A) <- paste0(names(replaced_data_A), "_A")
  names(replaced_data_B) <- paste0(names(replaced_data_B), "_B")
  
  # 合并 A 和 B 方法的数据
  all_data <- c(replaced_data_A, replaced_data_B)
  
  # 如果指定了 selected_versions，则筛选
  if (!is.null(selected_versions)) {
    selected_names <- c(paste0(selected_versions, "_A"),
                       paste0(selected_versions, "_B"))
    all_data <- all_data[names(all_data) %in% selected_names]
  }
  
  if (length(all_data) == 0) {
    stop("✗ 错误：没有可用的数据版本")
  }
  
  cat(sprintf("✓ 准备分析 %d 个数据版本\n", length(all_data)))
  cat(sprintf("  版本：%s\n", paste(names(all_data), collapse = ", ")))
  
  # ----------------------------------------------------------------------
  # 2. 从 sampleGroup 提取样本分组信息
  # ----------------------------------------------------------------------
  
  cat("\n--- 步骤 1: 提取样本分组信息 ---\n")
  
  # 只保留 Context 为 Experiment 或 Spatial 的样本
  samples_exp_spatial <- sampleGroup %>%
    filter(Context %in% c("Experiment", "Spatial")) %>%
    arrange(Order)
  
  if (nrow(samples_exp_spatial) == 0) {
    stop("✗ 错误：sampleGroup 中没有 Experiment 或 Spatial 样本")
  }
  
  cat(sprintf("✓ 找到 %d 个 Experiment/Spatial 样本\n", nrow(samples_exp_spatial)))
  
  # 提取唯一的 bioGroup（按 Order 排序，保证顺序一致）
  bioGroup_info <- samples_exp_spatial %>%
    group_by(bioGroup) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(Order) %>%
    select(bioGroup, Context, PLtype, SecondROCgroup)
  
  bioGroups_selected <- bioGroup_info$bioGroup
  
  cat(sprintf("✓ 选择的 bioGroup (%d):\n", length(bioGroups_selected)))
  for (i in seq_along(bioGroups_selected)) {
    bg <- bioGroups_selected[i]
    info <- bioGroup_info %>% filter(bioGroup == bg)
    cat(sprintf("  %d. %s (Context: %s, PLtype: %s, SecondROCgroup: %s)\n", 
                i, bg, info$Context, info$PLtype,
                ifelse(is.na(info$SecondROCgroup) | info$SecondROCgroup == "", 
                       "空", info$SecondROCgroup)))
  }
  
  # 构建 sample_info（用于 limma 设计矩阵）
  sample_info <- data.frame(
    SampleName = samples_exp_spatial$FinalName,
    Group = samples_exp_spatial$bioGroup,
    stringsAsFactors = FALSE
  )
  
  # 设置 Group 为因子，按 bioGroups_selected 的顺序
  sample_info$Group <- factor(sample_info$Group, levels = bioGroups_selected)
  
  cat(sprintf("✓ 构建 sample_info，共 %d 个样本\n", nrow(sample_info)))
  
  # ----------------------------------------------------------------------
  # 3. 构建比较矩阵（基于 Context 和 PLtype）
  # ----------------------------------------------------------------------
  
  cat("\n--- 步骤 2: 构建比较矩阵 ---\n")
  
  # 分别提取 Experiment 和 Spatial 的 bioGroup
  exp_groups <- bioGroup_info %>% filter(Context == "Experiment")
  spatial_groups <- bioGroup_info %>% filter(Context == "Spatial")
  
  cat(sprintf("✓ Experiment bioGroups (%d): %s\n", 
              nrow(exp_groups), 
              paste(exp_groups$bioGroup, collapse = ", ")))
  cat(sprintf("✓ Spatial bioGroups (%d): %s\n", 
              nrow(spatial_groups), 
              paste(spatial_groups$bioGroup, collapse = ", ")))
  
  # 构建比较列表
  comparisons_list <- list()
  comparison_names <- c()
  
  # 1. Exp vs Exp 的比较（所有Experiment组之间，不限PLtype）
  if (nrow(exp_groups) >= 2) {
    for (i in 1:(nrow(exp_groups) - 1)) {
      for (j in (i + 1):nrow(exp_groups)) {
        exp1 <- exp_groups$bioGroup[i]
        exp2 <- exp_groups$bioGroup[j]
        comparisons_list[[length(comparisons_list) + 1]] <- c(exp1, exp2)
        # 简化命名（去掉_Light/_H2O2后缀）
        name1 <- gsub("_Light$|_H2O2$", "", exp1)
        name2 <- gsub("_Light$|_H2O2$", "", exp2)
        comparison_names <- c(comparison_names, paste0(name1, "_vs_", name2))
      }
    }
  }
  
  # 2. Exp vs Spatial 的比较（必须PLtype相同）
  for (i in 1:nrow(exp_groups)) {
    exp_bg <- exp_groups$bioGroup[i]
    exp_pltype <- exp_groups$PLtype[i]
    
    for (j in 1:nrow(spatial_groups)) {
      spatial_bg <- spatial_groups$bioGroup[j]
      spatial_pltype <- spatial_groups$PLtype[j]
      
      # 只有PLtype相同才比较
      if (exp_pltype == spatial_pltype) {
        comparisons_list[[length(comparisons_list) + 1]] <- c(exp_bg, spatial_bg)
        # 简化命名
        name1 <- gsub("_Light$|_H2O2$", "", exp_bg)
        name2 <- gsub("_Light$|_H2O2$", "", spatial_bg)
        comparison_names <- c(comparison_names, paste0(name1, "_vs_", name2))
      }
    }
  }
  
  names(comparisons_list) <- comparison_names
  
  cat(sprintf("\n✓ 构建 %d 个比较组：\n", length(comparisons_list)))
  for (i in seq_along(comparisons_list)) {
    cat(sprintf("  %d. %s (%s vs %s)\n", i, names(comparisons_list)[i],
                comparisons_list[[i]][1], comparisons_list[[i]][2]))
  }
  
  # ----------------------------------------------------------------------
  # 4. 对每个数据版本进行 limma 差异分析
  # ----------------------------------------------------------------------
  
  cat("\n--- 步骤 3: 执行差异分析 ---\n")
  
  FDR_combined_df_list_2nd <- list()
  
  for (j in seq_along(all_data)) {
    mydata <- all_data[[j]]
    myfilename <- names(all_data)[j]
    
    cat(sprintf("\n[%d/%d] 处理：%s\n", j, length(all_data), myfilename))
    cat(sprintf("------%s_start------\n", myfilename))
    
    # 提取数据矩阵（去除最后3列注释列）
    myrange <- ncol(mydata) - 3
    
    # 检查数据结构
    if (myrange < 2) {
      cat("  ⚠ 警告：数据列数不足，跳过此版本\n")
      next
    }
    
    # 提取表达矩阵（从第2列到myrange列，第1列是Gene）
    data_cols <- colnames(mydata)[2:myrange]
    
    # 检查是否包含所有需要的样本
    missing_cols <- setdiff(sample_info$SampleName, data_cols)
    if (length(missing_cols) > 0) {
      cat(sprintf("  ⚠ 警告：缺少 %d 个样本列，跳过此版本\n", length(missing_cols)))
      cat(sprintf("    缺少：%s\n", paste(head(missing_cols, 3), collapse = ", ")))
      next
    }
    
    # 提取表达矩阵（按sample_info的顺序）
    proteomics_data <- mydata %>% 
      select(all_of(sample_info$SampleName)) %>%
      as.matrix()
    
    rownames(proteomics_data) <- mydata$Gene
    
    cat(sprintf("  ✓ 提取数据矩阵：%d genes × %d samples\n", 
                nrow(proteomics_data), ncol(proteomics_data)))
    
    # 创建设计矩阵
    design <- model.matrix(~ 0 + Group, data = sample_info)
    colnames(design) <- levels(sample_info$Group)
    
    cat("  ✓ 创建设计矩阵\n")
    
    # 拟合线性模型
    fit <- lmFit(proteomics_data, design)
    
    cat("  ✓ 拟合线性模型\n")
    
    # 构建对比矩阵表达式
    contrast_expressions <- character(length(comparisons_list))
    for (i in seq_along(comparisons_list)) {
      comp <- comparisons_list[[i]]
      comp_name <- names(comparisons_list)[i]
      # bioGroup 名称可能包含特殊字符，用反引号包围
      contrast_expressions[i] <- sprintf("%s = `%s` - `%s`", 
                                        comp_name, 
                                        comp[1], 
                                        comp[2])
    }
    
    # 使用 makeContrasts 动态构建对比矩阵
    contrast_matrix <- eval(parse(text = sprintf(
      "makeContrasts(%s, levels = design)",
      paste(contrast_expressions, collapse = ", ")
    )))
    
    cat("  ✓ 创建对比矩阵\n")
    
    # 计算对比并进行经验贝叶斯调整
    fit_contrasts <- contrasts.fit(fit, contrast_matrix)
    fit_ebayes <- eBayes(fit_contrasts)
    
    cat("  ✓ 经验贝叶斯调整\n")
    
    # 提取每个对比的结果
    FDR_test_list <- list()
    Raw_FDR_test_list <- list()
    Comparision_Group <- colnames(contrast_matrix)
    
    for (i in seq_along(Comparision_Group)) {
      myGroup <- Comparision_Group[i]
      
      # 提取 topTable 结果
      tem_dataframe <- topTable(fit_ebayes, coef = myGroup, n = Inf, adjust.method = "BH")
      tem_dataframe$Gene <- rownames(tem_dataframe)
      
      # 保存原始结果
      Raw_tem_dataframe <- tem_dataframe
      Raw_FDR_test_list[[paste0("RES_", myGroup)]] <- Raw_tem_dataframe
      
      # 只保留 Gene, logFC, adj.P.Val
      tem_dataframe_simple <- tem_dataframe %>% 
        select(Gene, logFC, adj.P.Val)
      colnames(tem_dataframe_simple) <- c("Gene",
                                         paste0(myGroup, "_logFC"),
                                         paste0(myGroup, "_adj.P.Val"))
      
      FDR_test_list[[paste0("RES_", myGroup)]] <- tem_dataframe_simple
    }
    
    cat(sprintf("  ✓ 提取 %d 个比较的差异分析结果\n", length(FDR_test_list)))
    
    # 保存原始 topTable 结果（每个比较一个 sheet）
    wb <- createWorkbook()
    for (i in seq_along(Raw_FDR_test_list)) {
      sheet_name <- names(Raw_FDR_test_list)[i]
      # Excel sheet 名称最多 31 个字符
      if (nchar(sheet_name) > 31) {
        sheet_name <- substr(sheet_name, 1, 31)
      }
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, Raw_FDR_test_list[[i]])
    }
    raw_file <- file.path(dir_config$output, 
                         paste0("Module11_Raw_FDR_test_list_", myfilename, ".xlsx"))
    saveWorkbook(wb, raw_file, overwrite = TRUE)
    cat(sprintf("  ✓ 导出原始 topTable：%s\n", basename(raw_file)))
    
    # 合并所有比较的 logFC 和 adj.P.Val
    FDR_combined_df <- reduce(FDR_test_list, full_join, by = "Gene")
    
    # 添加注释列（优先使用 *_Localization 命名的列）
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
        "✗ 错误：数据集 %s 中找不到任何 *_Localization 注释列，请检查输入数据。",
        myfilename
      ))
    }
    
    annotation_df <- mydata %>% select(Gene, all_of(annotation_cols))
    FDR_combined_df <- FDR_combined_df %>%
      left_join(annotation_df, by = "Gene")
    
    # 存储到列表
    FDR_combined_df_list_2nd[[myfilename]] <- FDR_combined_df
    
    cat(sprintf("  ✓ 合并结果：%d genes × %d columns\n", 
                nrow(FDR_combined_df), ncol(FDR_combined_df)))
    cat(sprintf("------%s_finished------\n", myfilename))
  }
  
  # ----------------------------------------------------------------------
  # 5. 导出合并的差异分析结果
  # ----------------------------------------------------------------------
  
  cat("\n--- 步骤 4: 导出合并结果 ---\n")
  
  if (length(FDR_combined_df_list_2nd) == 0) {
    cat("⚠ 警告：没有可导出的结果\n")
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
    # Excel sheet 名称最多 31 个字符
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
    }
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, FDR_combined_df_list_2nd[[i]])
  }
  
  combined_file <- file.path(dir_config$output, 
                            "Module11_FC_FDR_2nd_combined.xlsx")
  saveWorkbook(wb, combined_file, overwrite = TRUE)
  
  cat(sprintf("✓ 导出合并的差异分析结果：\n"))
  cat(sprintf("  %s\n", basename(combined_file)))
  cat(sprintf("  包含 %d 个数据版本\n", length(FDR_combined_df_list_2nd)))
  
  # ----------------------------------------------------------------------
  # 6. 汇总统计
  # ----------------------------------------------------------------------
  
  cat("\n=== Module 11 完成 ===\n")
  cat(sprintf("✓ 分析了 %d 个数据版本\n", length(FDR_combined_df_list_2nd)))
  cat(sprintf("✓ 构建了 %d 个比较组\n", length(comparisons_list)))
  cat(sprintf("✓ 输出文件数量：%d 个原始 topTable + 1 个合并结果\n", 
              length(FDR_combined_df_list_2nd)))
  
  # 返回结果
  return(list(
    FDR_combined_df_list_2nd = FDR_combined_df_list_2nd,
    sample_info = sample_info,
    comparisons_list = comparisons_list,
    bioGroups_selected = bioGroups_selected,
    bioGroup_info = bioGroup_info
  ))
}
