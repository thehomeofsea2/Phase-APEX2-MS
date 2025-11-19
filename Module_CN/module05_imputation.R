# Module 05: 缺失值填补
# 功能：对标准化后的数据进行缺失值填补（主要针对NoCat组）
# 作者：CodeNorm Pipeline
# 日期：2024

#' Module 05: 缺失值填补
#' 
#' @description
#' 对标准化后的数据进行缺失值填补，使用Perseus方法
#' 主要针对CatalyticGroup为NoCat的组，Cat组可选择是否填补
#' 
#' @param dir_config 目录配置列表（来自Module 1）
#' @param standardized_data_list 标准化数据列表（来自Module 4）
#' @param sampleGroup 样本分组信息（来自Module 2）
#' @param impute_cat_mean 是否对Cat组的n_valid=2情况用平均值填补（默认FALSE）
#' @param random_seed 随机数种子，用于保证结果可重复（默认123）
#' 
#' @return 列表，包含：
#'   - imputed_data_list: 包含所有填补后版本的列表
#'   - imputation_params: 填补参数信息
#' 
#' @details
#' 填补策略：
#' 1. Perseus参数计算：imputation_mean = col_mean - 1.8 * col_sd
#'                     imputation_sd = 0.3 * col_sd
#' 2. NoCat组填补规则：
#'    - n_valid == 2: 用组内平均值填补
#'    - n_valid < 2: 用Perseus参数随机抽样填补
#' 3. Cat组：默认不填补，用户可选择n_valid==2时用平均值填补
#' 4. 输出填补前后对比boxplot
#' 
#' @export
module05_imputation <- function(dir_config, 
                                standardized_data_list,
                                sampleGroup,
                                impute_cat_mean = FALSE,
                                random_seed = 123) {
  
  # 0. 加载必需包 ####
  require(tidyverse)
  require(dplyr)
  require(tidyr)
  
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("Module 05: 缺失值填补\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # 1. 输入验证 ####
  cat("\n[1] 验证输入数据...\n")
  
  if (!all(c("reference", "output") %in% names(dir_config))) {
    stop("❌ dir_config必须包含reference和output路径")
  }
  
  if (!is.list(standardized_data_list) || length(standardized_data_list) == 0) {
    stop("❌ standardized_data_list必须是非空列表")
  }
  
  required_cols <- c("FinalName", "bioGroup", "CatalyticGroup")
  if (!all(required_cols %in% names(sampleGroup))) {
    stop("❌ sampleGroup必须包含：", paste(required_cols, collapse = ", "))
  }
  
  cat("✓ 输入验证通过\n")
  cat(sprintf("  - 待处理数据集: %d 个\n", length(standardized_data_list)))
  cat(sprintf("  - Cat组填补策略: %s\n", ifelse(impute_cat_mean, "n_valid=2时填补", "不填补")))
  cat(sprintf("  - 随机数种子: %d\n", random_seed))
  
  # 2. 准备分组信息 ####
  cat("\n[2] 准备分组信息...\n")
  
  # 从第一个数据集中识别数据列和注释列
  first_data <- standardized_data_list[[1]]
  data_cols <- sampleGroup$FinalName
  anno_cols <- setdiff(names(first_data), c("Gene", data_cols))
  
  cat(sprintf("✓ 识别到 %d 个数据列\n", length(data_cols)))
  cat(sprintf("✓ 识别到 %d 个注释列: %s\n", 
              length(anno_cols), 
              paste(anno_cols, collapse = ", ")))
  
  # 识别NoCat和Cat组
  cat_groups <- sampleGroup %>%
    filter(CatalyticGroup == "Cat") %>%
    pull(bioGroup) %>%
    unique()
  
  nocat_groups <- sampleGroup %>%
    filter(CatalyticGroup == "NoCat") %>%
    pull(bioGroup) %>%
    unique()
  
  cat(sprintf("✓ Cat组 (%d个): %s\n", 
              length(cat_groups), 
              paste(cat_groups, collapse = ", ")))
  cat(sprintf("✓ NoCat组 (%d个): %s\n", 
              length(nocat_groups), 
              paste(nocat_groups, collapse = ", ")))
  
  # 3. 定义填补函数 ####
  impute_missing_values <- function(data, data_name, nocat_groups, impute_cat) {
    
    cat(sprintf("\n  处理: %s\n", data_name))
    
    # 计算每列的填补参数（Perseus方法）
    cat("    - 计算Perseus填补参数...\n")
    
    perseus_params <- data %>%
      select(all_of(data_cols)) %>%
      pivot_longer(everything(), names_to = "sample", values_to = "intensity") %>%
      filter(!is.na(intensity)) %>%
      group_by(sample) %>%
      summarise(
        col_mean = mean(intensity),
        col_sd = sd(intensity),
        .groups = "drop"
      ) %>%
      mutate(
        imputation_mean = col_mean - 1.8 * col_sd,
        imputation_sd = 0.3 * col_sd
      )
    
    # 为每个样本匹配bioGroup
    sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
    
    # 执行填补
    cat("    - 执行缺失值填补...\n")
    
    # 设置随机种子保证可重复性
    set.seed(random_seed)
    
    data_imputed <- data %>%
      pivot_longer(
        cols = all_of(data_cols),
        names_to = "sample",
        values_to = "intensity",
        values_drop_na = FALSE
      ) %>%
      mutate(group = sample_to_group[sample]) %>%
      group_by(Gene, group) %>%
      mutate(
        n_valid = sum(!is.na(intensity)),
        mean_valid_within_group = mean(intensity, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      left_join(perseus_params, by = "sample") %>%
      rowwise() %>%
      mutate(
        intensity_imputed = case_when(
          # NoCat组：n_valid=2用平均值，<2用Perseus
          group %in% nocat_groups & n_valid == 2 & is.na(intensity) ~ mean_valid_within_group,
          group %in% nocat_groups & n_valid < 2 & is.na(intensity) ~ rnorm(1, mean = imputation_mean, sd = imputation_sd),
          # Cat组：根据用户选择
          impute_cat & !(group %in% nocat_groups) & n_valid == 2 & is.na(intensity) ~ mean_valid_within_group,
          # 其他情况保持原值
          TRUE ~ intensity
        )
      ) %>%
      ungroup() %>%
      select(Gene, sample, intensity_imputed, all_of(anno_cols)) %>%
      distinct() %>%
      pivot_wider(
        names_from = sample,
        values_from = intensity_imputed
      ) %>%
      select(Gene, all_of(data_cols), all_of(anno_cols))
    
    # 统计填补情况
    n_imputed <- sum(is.na(data[, data_cols])) - sum(is.na(data_imputed[, data_cols]))
    cat(sprintf("    ✓ 填补了 %d 个缺失值\n", n_imputed))
    
    return(data_imputed)
  }
  
  # 4. 对所有数据集执行填补 ####
  cat("\n[3] 执行缺失值填补...\n")
  
  imputed_data_list <- list()
  
  for (i in seq_along(standardized_data_list)) {
    data_name <- names(standardized_data_list)[i]
    data <- standardized_data_list[[i]]
    
    # 执行填补
    data_imputed <- impute_missing_values(data, data_name, nocat_groups, impute_cat_mean)
    
    # 保存到列表
    imputed_name <- paste0(data_name, "_Imputed")
    imputed_data_list[[imputed_name]] <- data_imputed
    
    # 保存CSV
    csv_file <- file.path(dir_config$output, paste0("Module05_", imputed_name, ".csv"))
    write.csv(data_imputed, csv_file, row.names = FALSE)
    cat(sprintf("    ✓ 已保存: %s\n", basename(csv_file)))
  }
  
  cat(sprintf("\n✓ 共处理 %d 个数据集\n", length(imputed_data_list)))
  
  # 5. 生成对比boxplot ####
  cat("\n[4] 生成填补前后对比boxplot...\n")
  
  # 为不同bioGroup分配颜色（与Module 4一致）
  sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
  col_groups <- sample_to_group[data_cols]
  unique_groups <- unique(col_groups)
  n_groups <- length(unique_groups)
  
  if (n_groups <= 12) {
    group_colors <- scales::hue_pal()(n_groups)
  } else {
    group_colors <- rainbow(n_groups)
  }
  
  names(group_colors) <- unique_groups
  box_colors <- group_colors[col_groups]
  
  pdf_file <- file.path(dir_config$output, "Module05_Imputation_comparison_boxplot.pdf")
  pdf(pdf_file, width = 12, height = 6)
  
  for (i in seq_along(standardized_data_list)) {
    original_name <- names(standardized_data_list)[i]
    imputed_name <- names(imputed_data_list)[i]
    
    original_data <- standardized_data_list[[i]]
    imputed_data <- imputed_data_list[[i]]
    
    # 设置双图布局
    par(mfrow = c(1, 2))
    
    # 填补前
    boxplot(original_data[, data_cols], 
            log = "y", 
            cex.axis = 0.4, 
            las = 2, 
            main = paste0(original_name, "\n(Before Imputation)"),
            col = box_colors,
            border = "black")
    
    # 填补后
    boxplot(imputed_data[, data_cols], 
            log = "y", 
            cex.axis = 0.4, 
            las = 2, 
            main = paste0(imputed_name, "\n(After Imputation)"),
            col = box_colors,
            border = "black")
    
    par(mfrow = c(1, 1))
  }
  
  dev.off()
  cat(sprintf("✓ 已保存: %s\n", basename(pdf_file)))
  
  # 6. 验证注释列完整性 ####
  cat("\n[5] 验证注释列完整性...\n")
  
  for (i in seq_along(imputed_data_list)) {
    original_data <- standardized_data_list[[i]]
    imputed_data <- imputed_data_list[[i]]
    
    for (anno_col in anno_cols) {
      if (!all(original_data[[anno_col]] == imputed_data[[anno_col]], na.rm = TRUE)) {
        warning(sprintf("⚠ %s的%s列可能有变化", 
                       names(imputed_data_list)[i], anno_col))
      }
    }
  }
  
  cat("✓ 注释列完整性验证通过\n")
  
  # 7. 生成填补参数摘要 ####
  cat("\n[6] 生成填补参数摘要...\n")
  
  imputation_params <- list(
    nocat_groups = nocat_groups,
    cat_groups = cat_groups,
    impute_cat_mean = impute_cat_mean,
    random_seed = random_seed,
    data_cols = data_cols,
    anno_cols = anno_cols
  )
  
  # 8. 返回结果 ####
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("Module 05 完成\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("\n生成的数据集：\n")
  for (name in names(imputed_data_list)) {
    cat(sprintf("  - %s\n", name))
  }
  cat("\n")
  
  return(list(
    imputed_data_list = imputed_data_list,
    imputation_params = imputation_params
  ))
}

