# ============================================================================
# Module 2: 数据读取与分组表
# ============================================================================
# 功能：
#   1. 读取TSV数据文件
#   2. 自动清理列名（去除共同前缀/后缀）
#   3. 生成sampleGroup_template.csv（含提示元素）
#   4. 读取用户填写的sampleGroup.csv
#   5. 生成FinalName并重命名数据列
# 
# 输入：dir_config
# 输出：data_raw.RData, sampleGroup.RData
# ============================================================================

#' 读取数据并生成模板
#' 
#' @param dir_config 目录配置
#' @param file_pattern TSV文件匹配模式，默认"_matrix.*\\.tsv$"
#' @return 包含data和template的列表
module02_read_and_generate_template <- function(dir_config, file_pattern = "_matrix.*\\.tsv$") {
  
  cat("\n----------------------------------------\n")
  cat("步骤1: 读取数据文件\n")
  cat("----------------------------------------\n")
  
  # 读取TSV文件
  setwd(dir_config$rawdata)
  myfiles <- list.files(pattern = file_pattern)
  
  if (length(myfiles) == 0) {
    stop(sprintf("✗ 错误：在 %s 中未找到匹配 '%s' 的文件", dir_config$rawdata, file_pattern))
  }
  
  cat(sprintf("找到 %d 个文件:\n", length(myfiles)))
  for (i in seq_along(myfiles)) {
    cat(sprintf("  [%d] %s\n", i, myfiles[i]))
  }
  
  # 读取第一个文件（如有多个文件，可扩展此逻辑）
  data_raw <- read_tsv(myfiles[1], show_col_types = FALSE)
  cat(sprintf("\n✓ 已读取文件: %s\n", myfiles[1]))
  cat(sprintf("  维度: %d 行 × %d 列\n", nrow(data_raw), ncol(data_raw)))
  
  # --------------------------------------------------------------------------
  # 步骤2: 自动清理列名
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤2: 自动清理列名\n")
  cat("----------------------------------------\n")
  
  original_colnames <- colnames(data_raw)
  cat("原始列名:\n")
  print(head(original_colnames, 10))
  
  # 第一列默认为Gene
  cleaned_colnames <- original_colnames
  cleaned_colnames[1] <- "Gene"
  
  # 检测并去除共同前缀/后缀（除Gene列外）
  if (ncol(data_raw) > 1) {
    sample_cols <- original_colnames[-1]
    
    # 查找共同前缀
    common_prefix <- ""
    if (length(sample_cols) > 1) {
      # 使用第一个和最后一个字符串找共同前缀
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
    
    # 去除共同前缀（如果存在且有意义）
    if (nchar(common_prefix) > 0 && all(startsWith(sample_cols, common_prefix))) {
      cleaned_sample_cols <- substr(sample_cols, nchar(common_prefix) + 1, nchar(sample_cols))
      # 检查清理后是否还有有效内容
      if (all(nchar(cleaned_sample_cols) > 0)) {
        cleaned_colnames[-1] <- cleaned_sample_cols
        cat(sprintf("\n✓ 去除共同前缀: '%s'\n", common_prefix))
      }
    }
  }
  
  colnames(data_raw) <- cleaned_colnames
  cat("\n清理后列名:\n")
  print(head(colnames(data_raw), 10))
  
  # --------------------------------------------------------------------------
  # 步骤3: 生成sampleGroup模板
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤3: 生成sampleGroup模板\n")
  cat("----------------------------------------\n")
  
  # 获取样本列名（除Gene外）
  sample_names <- colnames(data_raw)[-1]
  n_samples <- length(sample_names)
  
  # 生成智能示例值（确保所有选项至少出现一次）
  generate_example_values <- function(n, options) {
    n_opts <- length(options)
    if (n <= n_opts) {
      return(options[1:n])
    } else {
      # 均匀分配
      repeats <- rep(ceiling(n / n_opts), n_opts)
      values <- rep(options, repeats)
      return(values[1:n])
    }
  }
  
  # 创建模板数据框
  template <- data.frame(
    OriginalName = sample_names,
    bioGroup = "",
    CatalyticGroup = generate_example_values(n_samples, c("Cata", "NoCat")),
    PLtype = generate_example_values(n_samples, c("Light", "H2O2", "PL")),
    Context = generate_example_values(n_samples, c("Experiment", "Control", "Spatial")),
    replicate = rep(1:3, length.out = n_samples),
    FirstROCgroup = generate_example_values(n_samples, c("A", "B", "C", "NA")),
    SecondROCgroup = generate_example_values(n_samples, c("A", "B", "C", "D", "E", "F")),
    Order = 1:n_samples,  # 默认顺序
    FinalName = "",
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("✓ 生成模板: %d 个样本\n", nrow(template)))
  
  # 保存模板到工作目录
  setwd(dir_config$root)
  write.csv(template, "Module02_sampleGroup_template.csv", row.names = FALSE)
  cat("\n✓ 模板已保存: Module02_sampleGroup_template.csv\n")
  cat("\n📝 请填写以下列:\n")
  cat("  - bioGroup: 生物学分组名称\n")
  cat("  - CatalyticGroup: Cata 或 NoCat\n")
  cat("  - PLtype: Light 或 H2O2 或 PL\n")
  cat("  - Context: Experiment 或 Control 或 Spatial\n")
  cat("  - replicate: 1 或 2 或 3\n")
  cat("  - FirstROCgroup: A / B / C / A&B / NA\n")
  cat("  - SecondROCgroup: A / B / C / D / E / F\n")
  cat("  - Order: 整数（1/2/3...），控制数据列展示顺序，1最先出现\n")
  cat("  (FinalName 将自动生成，请留空)\n")
  
  return(list(data = data_raw, template = template))
}


#' 读取用户填写的sampleGroup_template并生成sampleGroup
#' 
#' @param dir_config 目录配置
#' @param data_raw 原始数据
#' @return 包含data和sampleGroup的列表
module02_process_samplegroup <- function(dir_config, data_raw) {
  
  cat("\n----------------------------------------\n")
  cat("步骤4: 读取用户填写的sampleGroup_template\n")
  cat("----------------------------------------\n")
  
  setwd(dir_config$root)
  
  # 检查文件是否存在
  template_file <- "Module02_sampleGroup_template.csv"
  if (!file.exists(template_file)) {
    stop(sprintf("✗ 错误：未找到文件 %s\n请先运行步骤1生成模板", template_file))
  }
  
  # 读取用户填写的模板文件
  sampleGroup <- read.csv(template_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ 已读取: %s\n", template_file))
  cat(sprintf("  行数: %d\n", nrow(sampleGroup)))
  
  # 验证必需列
  required_cols <- c("OriginalName", "bioGroup", "CatalyticGroup", "PLtype", 
                     "Context", "replicate", "FirstROCgroup", "SecondROCgroup", "Order")
  missing_cols <- setdiff(required_cols, colnames(sampleGroup))
  if (length(missing_cols) > 0) {
    stop(sprintf("✗ 错误：缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # 验证是否填写
  empty_bioGroup <- sampleGroup$bioGroup == "" | is.na(sampleGroup$bioGroup)
  if (any(empty_bioGroup)) {
    stop(sprintf("✗ 错误：第 %s 行的 bioGroup 未填写", 
                 paste(which(empty_bioGroup), collapse = ", ")))
  }
  
  # 验证Order列
  if (any(is.na(sampleGroup$Order))) {
    stop("✗ 错误：Order列存在空值，请填写整数")
  }
  if (!all(sampleGroup$Order == as.integer(sampleGroup$Order))) {
    stop("✗ 错误：Order列必须为整数")
  }
  
  # --------------------------------------------------------------------------
  # 步骤5: 根据Order排序sampleGroup和重排数据列
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤5: 根据Order重排数据列\n")
  cat("----------------------------------------\n")
  
  # 按Order排序sampleGroup
  sampleGroup <- sampleGroup %>% arrange(Order)
  cat("✓ sampleGroup已按Order排序\n")
  
  # 获取排序后的列名顺序（保持Gene在第一列）
  ordered_sample_names <- sampleGroup$OriginalName
  
  # 检查是否所有OriginalName都在data_raw中
  missing_samples <- setdiff(ordered_sample_names, colnames(data_raw)[-1])
  if (length(missing_samples) > 0) {
    stop(sprintf("✗ 错误：sampleGroup中的样本在数据中不存在: %s", 
                 paste(missing_samples, collapse = ", ")))
  }
  
  # 重排data_raw的列：Gene列 + 按Order排序的样本列
  data_raw <- data_raw %>% select(Gene, all_of(ordered_sample_names))
  cat("✓ 数据列已按Order重新排列\n")
  cat("\n新的列顺序:\n")
  print(head(colnames(data_raw), 10))
  
  # --------------------------------------------------------------------------
  # 步骤6: 生成FinalName
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤6: 生成FinalName\n")
  cat("----------------------------------------\n")
  
  # FinalName = bioGroup + _LFQ_ + replicate
  sampleGroup$FinalName <- paste0(sampleGroup$bioGroup, "_LFQ_", sampleGroup$replicate)
  
  cat("✓ FinalName 生成完成\n")
  cat("\n示例:\n")
  print(head(sampleGroup[, c("OriginalName", "bioGroup", "replicate", "FinalName", "Order")], 5))
  
  # --------------------------------------------------------------------------
  # 步骤7: 重命名数据列
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤7: 重命名数据列\n")
  cat("----------------------------------------\n")
  
  # 创建列名映射
  name_mapping <- setNames(sampleGroup$FinalName, sampleGroup$OriginalName)
  
  # 重命名（保持Gene列不变）
  current_colnames <- colnames(data_raw)
  new_colnames <- current_colnames
  for (i in 2:length(current_colnames)) {
    old_name <- current_colnames[i]
    if (old_name %in% names(name_mapping)) {
      new_colnames[i] <- name_mapping[old_name]
    }
  }
  
  colnames(data_raw) <- new_colnames
  cat("✓ 列名重命名完成\n")
  cat("\n新列名:\n")
  print(head(colnames(data_raw), 10))
  
  # --------------------------------------------------------------------------
  # 步骤8: 保存数据
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤8: 保存数据\n")
  cat("----------------------------------------\n")
  
  # 保存最终的sampleGroup到工作目录（CSV格式）
  write.csv(sampleGroup, "Module02_sampleGroup.csv", row.names = FALSE)
  cat("✓ 已保存: Module02_sampleGroup.csv (工作目录，最终版本)\n")
  cat("  该文件包含验证后的分组信息和生成的FinalName\n")
  
  # 保存CSV到Output目录
  output_file <- file.path(dir_config$output, "Module02_data_raw.csv")
  write.csv(data_raw, output_file, row.names = FALSE)
  cat(sprintf("✓ 已保存: %s\n", output_file))
  
  cat("\n✓ 步骤8完成，返回数据到主程序\n")
  
  return(list(data = data_raw, sampleGroup = sampleGroup))
}


#' Module 2 主函数
#' 
#' @param dir_config 目录配置
#' @param file_pattern TSV文件匹配模式
#' @param auto_process 是否自动处理（FALSE则只生成模板）
module02_data_import <- function(dir_config, file_pattern = "_matrix.*\\.tsv$", auto_process = FALSE) {
  
  cat("\n========================================\n")
  cat("Module 2: 数据读取与分组表\n")
  cat("========================================\n")
  
  # 步骤1-3: 读取数据并生成模板
  result <- module02_read_and_generate_template(dir_config, file_pattern)
  
  if (auto_process) {
    # 如果自动处理，尝试读取已填写的sampleGroup
    result <- module02_process_samplegroup(dir_config, result$data)
  } else {
    cat("\n⚠ 请完成以下步骤:\n")
    cat("  1. 打开 sampleGroup_template.csv\n")
    cat("  2. 填写所有必需列\n")
    cat("  3. 另存为 sampleGroup.csv\n")
    cat("  4. 在 main_pipeline.R 中继续运行\n")
  }
  
  return(invisible(result))
}

