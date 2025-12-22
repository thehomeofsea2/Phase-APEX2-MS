# Module 04: 数据标准化
# 功能：Log2转化、全局标准化（QNorm/QNorm_Fix/QNorm_Limma/VSN/RefNorm/MNorm）、局部标准化（Local QNorm/QNorm_Fix/QNorm_Limma/VSN/RefNorm/MNorm）
# 作者：CodeNorm Pipeline
# 日期：2024

# 分位数标准化的NA安全版本：先以每列中位数填补，再标准化，最后恢复NA
quantile_normalize_preserve_na <- function(mat) {
  na_mask <- is.na(mat)
  col_medians <- apply(mat, 2, function(x) {
    med <- median(x, na.rm = TRUE)
    if (is.na(med)) 0 else med
  })
  
  filled_mat <- mat
  for (i in seq_len(ncol(filled_mat))) {
    filled_mat[is.na(filled_mat[, i]), i] <- col_medians[i]
  }
  
  normalized <- normalize.quantiles(filled_mat)
  colnames(normalized) <- colnames(mat)
  rownames(normalized) <- rownames(mat)
  
  normalized[na_mask] <- NA
  normalized
}

# 将目标列对齐到参考分布（RefNorm辅助函数）
align_to_reference <- function(vec, ref_vec) {
  if (all(is.na(vec))) return(vec)
  ref_sorted <- sort(ref_vec, na.last = NA)
  if (length(ref_sorted) == 0) return(vec)
  ref_probs <- seq(0, 1, length.out = length(ref_sorted))
  
  idx <- which(!is.na(vec))
  if (length(idx) == 0) return(vec)
  
  vec_non_na <- vec[idx]
  vec_order <- order(vec_non_na)
  
  target_vals <- approx(ref_probs, ref_sorted, xout = seq(0, 1, length.out = length(idx)), rule = 2)$y
  
  new_vals <- vec_non_na
  new_vals[vec_order] <- target_vals
  
  vec[idx] <- new_vals
  vec
}

# 选择参考样本（RefNorm）
choose_reference_sample <- function(cols, data_matrix, sampleGroup, context_only_control = FALSE) {
  ctx_map <- setNames(sampleGroup$Context, sampleGroup$FinalName)
  contexts <- ctx_map[cols]
  valid_counts <- apply(data_matrix[, cols, drop = FALSE], 2, function(x) sum(!is.na(x)))
  
  pick_by_max <- function(pool_cols) {
    if (length(pool_cols) == 0) return(NULL)
    pool_counts <- valid_counts[pool_cols]
    pool_cols[which.max(pool_counts)]
  }
  
  if (context_only_control) {
    med_count <- median(valid_counts, na.rm = TRUE)
    if (is.na(med_count)) med_count <- 0
    diffs <- abs(valid_counts - med_count)
    ref_col <- cols[which.min(diffs)]
    return(ref_col)
  }
  
  exp_cols <- cols[contexts == "Experiment"]
  ref_col <- pick_by_max(exp_cols)
  if (is.null(ref_col)) {
    ref_col <- pick_by_max(cols)
  }
  ref_col
}

# VSN转换，保持全NA行位置（避免警告“only NA elements”）
vsn_transform_preserve_na <- function(mat) {
  if (!require(vsn, quietly = TRUE)) {
    stop("❌ 需要安装vsn包：BiocManager::install('vsn')")
  }
  row_all_na <- apply(mat, 1, function(x) all(is.na(x)))
  mat_use <- mat[!row_all_na, , drop = FALSE]
  if (nrow(mat_use) == 0) return(mat)
  min_val <- suppressWarnings(min(mat_use, na.rm = TRUE))
  if (is.finite(min_val) && min_val <= 0) {
    shift <- abs(min_val) + 1e-3
    mat_use <- mat_use + shift
  }
  fit_vsn <- vsn::vsn2(mat_use)
  transformed <- predict(fit_vsn, newdata = mat_use)
  out <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  out[!row_all_na, ] <- transformed
  colnames(out) <- colnames(mat)
  rownames(out) <- rownames(mat)
  out
}

#' Module 04: 数据标准化
#' 
#' @description
#' 对数据进行Log2转化和多种标准化处理，包括全局标准化和局部标准化
#' 
#' @param dir_config 目录配置列表（来自Module 1）
#' @param data_annotated 带注释的数据框（来自Module 3）
#' @param sampleGroup 样本分组信息（来自Module 2）
#' @param norm_types 字符向量，指定需要的标准化类型
#'   可选值："noNorm", "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
#'            "Global_VSN", "Global_RefNorm", "Global_MNorm",
#'            "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
#'            "Local_VSN", "Local_RefNorm", "Local_MNorm"
#'   默认：c("noNorm", "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
#'            "Global_VSN", "Global_RefNorm", "Global_MNorm",
#'            "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
#'            "Local_VSN", "Local_RefNorm", "Local_MNorm")
#'   注意：必须至少包含一个全局标准化（Global_QNorm / Global_QNorm_Fix / Global_QNorm_Limma / Global_VSN / Global_RefNorm / Global_MNorm）
#' 
#' @return 列表，包含：
#'   - standardized_data_list: 包含所有标准化版本的列表
#'   - norm_types_used: 实际使用的标准化类型
#' 
#' @details
#' - Log2转化：对数值列进行log2转换
#' - Global标准化：对所有样本统一标准化
#' - Local标准化：在同一bioGroup内分别标准化，然后合并
#' - 输出boxplot展示不同标准化方法的效果
#' 
#' @export
module04_standardization <- function(dir_config, 
                                     data_annotated, 
                                     sampleGroup,
                                     norm_types = c("noNorm",
                                                   "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
                                                   "Global_VSN", "Global_RefNorm", "Global_MNorm",
                                                   "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
                                                   "Local_VSN", "Local_RefNorm", "Local_MNorm")) {
  
  cat("\n=== Module 04: 数据标准化 ===\n")
  
  # 1. 输入验证 ####
  cat("\n[1] 验证输入数据...\n")
  
  if (!all(c("reference", "output") %in% names(dir_config))) {
    stop("❌ dir_config必须包含reference和output路径")
  }
  
  if (!is.data.frame(data_annotated)) {
    stop("❌ data_annotated必须是数据框")
  }
  
  if (!"Gene" %in% colnames(data_annotated)) {
    stop("❌ data_annotated必须包含Gene列")
  }
  
  if (!is.data.frame(sampleGroup)) {
    stop("❌ sampleGroup必须是数据框")
  }
  
  if (!"bioGroup" %in% colnames(sampleGroup)) {
    stop("❌ sampleGroup必须包含bioGroup列")
  }
  if (any(c("Global_RefNorm", "Local_RefNorm") %in% norm_types) && !"Context" %in% colnames(sampleGroup)) {
    stop("❌ 使用RefNorm需要sampleGroup包含Context列")
  }
  
  # 验证norm_types
  valid_types <- c("noNorm",
                   "Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma",
                   "Global_VSN", "Global_RefNorm", "Global_MNorm",
                   "Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma",
                   "Local_VSN", "Local_RefNorm", "Local_MNorm")
  if (!all(norm_types %in% valid_types)) {
    stop("❌ norm_types包含无效值，有效值为：", paste(valid_types, collapse = ", "))
  }
  
  # 检查是否至少有一个全局标准化
  if (!any(c("Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma", "Global_VSN", "Global_RefNorm", "Global_MNorm") %in% norm_types)) {
    stop("❌ 必须至少包含一个全局标准化方法（Global_QNorm / Global_QNorm_Fix / Global_QNorm_Limma / Global_VSN / Global_RefNorm / Global_MNorm）")
  }
  
  cat("✓ 输入验证通过\n")
  cat(sprintf("  - 数据维度: %d行 × %d列\n", nrow(data_annotated), ncol(data_annotated)))
  cat(sprintf("  - 样本数量: %d\n", nrow(sampleGroup)))
  cat(sprintf("  - 标准化类型: %s\n", paste(norm_types, collapse = ", ")))
  
  # 2. 识别数据列和注释列 ####
  cat("\n[2] 识别数据列和注释列...\n")
  
  # 识别注释列（包含"Localization"的列）
  annotation_cols <- grep("Localization$", colnames(data_annotated), value = TRUE)
  
  if (length(annotation_cols) == 0) {
    stop("❌ 未找到注释列（应包含'Localization'后缀）")
  }
  
  # 数据列：除了Gene列和注释列之外的所有列
  data_cols <- setdiff(colnames(data_annotated), c("Gene", annotation_cols))
  
  cat(sprintf("✓ 识别完成\n"))
  cat(sprintf("  - 数据列数量: %d\n", length(data_cols)))
  cat(sprintf("  - 注释列数量: %d (%s)\n", 
              length(annotation_cols), paste(annotation_cols, collapse = ", ")))
  
  # 3. Log2转化 ####
  cat("\n[3] 执行Log2转化...\n")
  
  data_log2 <- data_annotated %>%
    mutate(across(all_of(data_cols), ~log2(.)))
  
  cat("✓ Log2转化完成\n")
  
  # 保存log2数据到CSV
  csv_file <- file.path(dir_config$output, "Module04_log2.csv")
  write.csv(data_log2, csv_file, row.names = FALSE)
  cat(sprintf("✓ 已保存: %s\n", csv_file))
  
  # 为不同bioGroup分配颜色 ####
  # 构建样本到bioGroup的映射
  sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
  
  # 验证数据列是否在sampleGroup中
  valid_data_cols <- intersect(data_cols, names(sample_to_group))
  if (length(valid_data_cols) < length(data_cols)) {
    missing_cols <- setdiff(data_cols, valid_data_cols)
    warning("⚠ 以下数据列未在sampleGroup中找到，将使用灰色：", 
            paste(missing_cols, collapse = ", "))
  }
  
  # 为每个数据列分配bioGroup
  col_groups <- sample_to_group[valid_data_cols]
  unique_groups <- unique(col_groups)
  n_groups <- length(unique_groups)
  
  cat(sprintf("  - 识别到 %d 个不同的bioGroup\n", n_groups))
  
  # 生成颜色调色板（使用rainbow或其他调色板）
  if (n_groups <= 12) {
    # 使用Set3调色板（柔和的颜色）
    group_colors <- scales::hue_pal()(n_groups)
  } else {
    # 使用rainbow调色板
    group_colors <- rainbow(n_groups)
  }
  names(group_colors) <- unique_groups
  
  # 为每个数据列分配颜色
  box_colors <- character(length(data_cols))
  names(box_colors) <- data_cols
  for (col in data_cols) {
    if (col %in% names(col_groups)) {
      box_colors[col] <- group_colors[col_groups[col]]
    } else {
      box_colors[col] <- "grey80"  # 未找到的列使用灰色
    }
  }
  
  # 绘制log2 boxplot（按bioGroup着色）
  pdf_file <- file.path(dir_config$output, "Module04_log2_boxplot.pdf")
  pdf(pdf_file, width = 10, height = 5)
  boxplot(data_log2[, data_cols], 
          cex.axis = 0.4, 
          las = 2, 
          main = "Log2 Transformed Data",
          col = box_colors[data_cols])
  # 添加图例
  legend("topright", 
         legend = names(group_colors), 
         fill = group_colors, 
         cex = 0.6,
         title = "bioGroup")
  dev.off()
  cat(sprintf("✓ 已保存: %s\n", pdf_file))
  
  # 4. 初始化结果列表 ####
  standardized_data_list <- list()
  
  # noNorm: log2转化后的数据
  if ("noNorm" %in% norm_types) {
    standardized_data_list[["noNorm"]] <- data_log2
    cat("✓ 已添加: noNorm (仅Log2转化)\n")
  }
  
  # 5. 全局标准化 ####
  if (any(c("Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma", "Global_VSN", "Global_RefNorm", "Global_MNorm") %in% norm_types)) {
    cat("\n[4] 执行全局标准化...\n")
    
    # 加载preprocessCore包
    if (!require(preprocessCore, quietly = TRUE)) {
      stop("❌ 需要安装preprocessCore包：BiocManager::install('preprocessCore')")
    }
    
    # 提取数据矩阵（用于标准化）
    data_matrix <- as.matrix(data_log2[, data_cols])
    
    # 5.1 Global Quantile Normalization ####
    if ("Global_QNorm" %in% norm_types) {
      cat("  执行Global Quantile Normalization...\n")
      
      normalized_matrix <- normalize.quantiles(data_matrix)
      colnames(normalized_matrix) <- data_cols
      
      data_qnorm <- data_log2
      data_qnorm[, data_cols] <- as.data.frame(normalized_matrix)
      
      standardized_data_list[["Global_QNorm"]] <- data_qnorm
      
      # 保存CSV
      csv_file <- file.path(dir_config$output, "Module04_Global_QNorm.csv")
      write.csv(data_qnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_QNorm完成，已保存: %s\n", csv_file))
    }
    
    # 5.1b Global Quantile Normalization (NA安全版本) ####
    if ("Global_QNorm_Fix" %in% norm_types) {
      cat("  执行Global Quantile Normalization（QNorm_Fix，伪装NA后恢复）...\n")
      
      normalized_matrix_fix <- quantile_normalize_preserve_na(data_matrix)
      colnames(normalized_matrix_fix) <- data_cols
      
      data_qnorm_fix <- data_log2
      data_qnorm_fix[, data_cols] <- as.data.frame(normalized_matrix_fix)
      
      standardized_data_list[["Global_QNorm_Fix"]] <- data_qnorm_fix
      
      # 保存CSV
      csv_file <- file.path(dir_config$output, "Module04_Global_QNorm_Fix.csv")
      write.csv(data_qnorm_fix, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_QNorm_Fix完成，已保存: %s\n", csv_file))
    }
    
    # 5.1c Global Quantile Normalization (limma, NA安全) ####
    if ("Global_QNorm_Limma" %in% norm_types) {
      cat("  执行Global Quantile Normalization（limma）...\n")
      if (!require(limma, quietly = TRUE)) {
        stop("❌ 需要安装limma包：BiocManager::install('limma')")
      }
      normalized_matrix_limma <- limma::normalizeBetweenArrays(data_matrix, method = "quantile")
      colnames(normalized_matrix_limma) <- data_cols
      
      data_qnorm_limma <- data_log2
      data_qnorm_limma[, data_cols] <- as.data.frame(normalized_matrix_limma)
      
      standardized_data_list[["Global_QNorm_Limma"]] <- data_qnorm_limma
      
      csv_file <- file.path(dir_config$output, "Module04_Global_QNorm_Limma.csv")
      write.csv(data_qnorm_limma, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_QNorm_Limma完成，已保存: %s\n", csv_file))
    }
    
    # 5.1d Global VSN ####
    if ("Global_VSN" %in% norm_types) {
      cat("  执行Global VSN...\n")
      vsn_input <- as.matrix(data_annotated[, data_cols])
      normalized_matrix_vsn <- vsn_transform_preserve_na(vsn_input)
      colnames(normalized_matrix_vsn) <- data_cols
      
      data_vsn <- data_log2
      data_vsn[, data_cols] <- as.data.frame(normalized_matrix_vsn)
      
      standardized_data_list[["Global_VSN"]] <- data_vsn
      
      csv_file <- file.path(dir_config$output, "Module04_Global_VSN.csv")
      write.csv(data_vsn, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_VSN完成，已保存: %s\n", csv_file))
    }
    
    # 5.1e Global Reference Normalization（上下文感知） ####
    if ("Global_RefNorm" %in% norm_types) {
      cat("  执行Global Reference Normalization（上下文感知）...\n")
      
      sample_ctx_map <- setNames(sampleGroup$Context, sampleGroup$FinalName)
      ref_col <- choose_reference_sample(data_cols, data_matrix, sampleGroup, context_only_control = FALSE)
      if (is.null(ref_col)) {
        stop("❌ 无法选择参考样本（Global_RefNorm）")
      }
      cat(sprintf("    参考样本: %s (Context=%s)\n", ref_col, sample_ctx_map[[ref_col]]))
      
      ref_vec <- data_matrix[, ref_col]
      refnorm_mat <- data_matrix
      for (cn in data_cols) {
        refnorm_mat[, cn] <- align_to_reference(refnorm_mat[, cn], ref_vec)
      }
      
      data_refnorm <- data_log2
      data_refnorm[, data_cols] <- as.data.frame(refnorm_mat)
      
      standardized_data_list[["Global_RefNorm"]] <- data_refnorm
      
      csv_file <- file.path(dir_config$output, "Module04_Global_RefNorm.csv")
      write.csv(data_refnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_RefNorm完成，已保存: %s\n", csv_file))
    }
    
    # 5.2 Global Median Normalization ####
    if ("Global_MNorm" %in% norm_types) {
      cat("  执行Global Median Normalization...\n")
      
      # 计算每个样本的中位数
      sample_medians <- apply(data_matrix, 2, median, na.rm = TRUE)
      cat(sprintf("    样本中位数范围: %.2f - %.2f\n", 
                  min(sample_medians, na.rm = TRUE), 
                  max(sample_medians, na.rm = TRUE)))
      
      # 计算全局中位数
      global_median <- median(sample_medians, na.rm = TRUE)
      cat(sprintf("    全局中位数: %.2f\n", global_median))
      
      # 标准化
      normalized_matrix <- sweep(data_matrix, 2, sample_medians / global_median, "/")
      
      data_mnorm <- data_log2
      data_mnorm[, data_cols] <- as.data.frame(normalized_matrix)
      
      standardized_data_list[["Global_MNorm"]] <- data_mnorm
      
      # 保存CSV
      csv_file <- file.path(dir_config$output, "Module04_Global_MNorm.csv")
      write.csv(data_mnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Global_MNorm完成，已保存: %s\n", csv_file))
    }
    
    # 绘制全局标准化boxplot
    pdf_file <- file.path(dir_config$output, "Module04_Global_Norm_boxplot.pdf")
    pdf(pdf_file, width = 10, height = 5)
    
    for (norm_name in c("Global_QNorm", "Global_QNorm_Fix", "Global_QNorm_Limma", "Global_VSN", "Global_RefNorm", "Global_MNorm")) {
      if (norm_name %in% names(standardized_data_list)) {
        boxplot(standardized_data_list[[norm_name]][, data_cols], 
                cex.axis = 0.4, 
                las = 2, 
                main = norm_name,
                col = box_colors[data_cols])
        # 添加图例
        legend("topright", 
               legend = names(group_colors), 
               fill = group_colors, 
               cex = 0.6,
               title = "bioGroup")
      }
    }
    
    dev.off()
    cat(sprintf("✓ 已保存: %s\n", pdf_file))
  }
  
  # 6. 局部标准化 ####
  if (any(c("Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma", "Local_VSN", "Local_RefNorm", "Local_MNorm") %in% norm_types)) {
    cat("\n[5] 执行局部标准化（按bioGroup分组）...\n")
    
    # 确保preprocessCore已加载
    if (!require(preprocessCore, quietly = TRUE)) {
      stop("❌ 需要安装preprocessCore包：BiocManager::install('preprocessCore')")
    }
    
    # 构建样本到bioGroup的映射
    # sampleGroup$FinalName对应data_cols
    sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
    
    # 验证所有数据列都有对应的bioGroup
    missing_samples <- setdiff(data_cols, names(sample_to_group))
    if (length(missing_samples) > 0) {
      warning("⚠ 以下样本未在sampleGroup中找到对应的bioGroup，将被跳过：",
              paste(missing_samples, collapse = ", "))
      data_cols <- setdiff(data_cols, missing_samples)
    }
    
    # 创建group_vec
    group_vec <- sample_to_group[data_cols]
    
    cat(sprintf("  - 唯一bioGroup数量: %d\n", length(unique(group_vec))))
    cat(sprintf("  - bioGroup分布:\n"))
    for (grp in unique(group_vec)) {
      n_samples <- sum(group_vec == grp)
      cat(sprintf("    %s: %d个样本\n", grp, n_samples))
    }
    
    # 提取数据矩阵
    data_matrix <- as.matrix(data_log2[, data_cols])
    
    # 6.1 Local Median Normalization ####
    if ("Local_MNorm" %in% norm_types) {
      cat("\n  执行Local Median Normalization...\n")
      
      # 初始化标准化矩阵
      mnormalized_mat <- data_matrix
      
      # 对每个bioGroup分别标准化
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    处理 %s (%d个样本)...\n", grp, length(grp_cols)))
        
        # 提取该组的子矩阵
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        
        # 计算每列的中位数
        sample_medians <- apply(sub_mat, 2, median, na.rm = TRUE)
        
        # 计算该组的全局中位数
        group_global_median <- median(sample_medians, na.rm = TRUE)
        
        # 标准化
        normalized_sub_mat <- sweep(sub_mat, 2, sample_medians / group_global_median, "/")
        
        # 写回结果矩阵
        mnormalized_mat[, grp_cols] <- normalized_sub_mat
      }
      
      # 构建最终数据框
      data_local_mnorm <- data_log2
      data_local_mnorm[, data_cols] <- as.data.frame(mnormalized_mat)
      
      standardized_data_list[["Local_MNorm"]] <- data_local_mnorm
      
      # 保存CSV
      csv_file <- file.path(dir_config$output, "Module04_Local_MNorm.csv")
      write.csv(data_local_mnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_MNorm完成，已保存: %s\n", csv_file))
    }
    
    # 6.2 Local Quantile Normalization ####
    if ("Local_QNorm" %in% norm_types) {
      cat("\n  执行Local Quantile Normalization...\n")
      
      # 初始化标准化矩阵
      qnormalized_mat <- data_matrix
      
      # 对每个bioGroup分别标准化
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    处理 %s (%d个样本)...\n", grp, length(grp_cols)))
        
        # 提取子矩阵
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        
        # Quantile normalization
        normalized_sub_mat <- normalize.quantiles(as.matrix(sub_mat))
        colnames(normalized_sub_mat) <- grp_cols
        
        # 写回标准化矩阵
        qnormalized_mat[, grp_cols] <- normalized_sub_mat
      }
      
      # 构建最终数据框
      data_local_qnorm <- data_log2
      data_local_qnorm[, data_cols] <- as.data.frame(qnormalized_mat)
      
      standardized_data_list[["Local_QNorm"]] <- data_local_qnorm
      
      # 保存CSV
      csv_file <- file.path(dir_config$output, "Module04_Local_QNorm.csv")
      write.csv(data_local_qnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_QNorm完成，已保存: %s\n", csv_file))
    }
    
    # 6.2b Local Quantile Normalization (NA安全版本) ####
    if ("Local_QNorm_Fix" %in% norm_types) {
      cat("\n  执行Local Quantile Normalization（QNorm_Fix，伪装NA后恢复）...\n")
      
      qnormalized_fix_mat <- data_matrix
      
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    处理 %s (%d个样本)...\n", grp, length(grp_cols)))
        
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        normalized_sub_mat <- quantile_normalize_preserve_na(sub_mat)
        colnames(normalized_sub_mat) <- grp_cols
        
        qnormalized_fix_mat[, grp_cols] <- normalized_sub_mat
      }
      
      data_local_qnorm_fix <- data_log2
      data_local_qnorm_fix[, data_cols] <- as.data.frame(qnormalized_fix_mat)
      
      standardized_data_list[["Local_QNorm_Fix"]] <- data_local_qnorm_fix
      
      csv_file <- file.path(dir_config$output, "Module04_Local_QNorm_Fix.csv")
      write.csv(data_local_qnorm_fix, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_QNorm_Fix完成，已保存: %s\n", csv_file))
    }

    # 6.2c Local Quantile Normalization (limma) ####
    if ("Local_QNorm_Limma" %in% norm_types) {
      cat("\n  执行Local Quantile Normalization（limma）...\n")
      if (!require(limma, quietly = TRUE)) {
        stop("❌ 需要安装limma包：BiocManager::install('limma')")
      }
      
      qnormalized_limma_mat <- data_matrix
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    处理 %s (%d个样本)...\n", grp, length(grp_cols)))
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        normalized_sub_mat <- limma::normalizeBetweenArrays(sub_mat, method = "quantile")
        colnames(normalized_sub_mat) <- grp_cols
        qnormalized_limma_mat[, grp_cols] <- normalized_sub_mat
      }
      
      data_local_qnorm_limma <- data_log2
      data_local_qnorm_limma[, data_cols] <- as.data.frame(qnormalized_limma_mat)
      
      standardized_data_list[["Local_QNorm_Limma"]] <- data_local_qnorm_limma
      
      csv_file <- file.path(dir_config$output, "Module04_Local_QNorm_Limma.csv")
      write.csv(data_local_qnorm_limma, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_QNorm_Limma完成，已保存: %s\n", csv_file))
    }
    
    # 6.2d Local VSN ####
    if ("Local_VSN" %in% norm_types) {
      cat("\n  执行Local VSN...\n")
      
      vsn_mat <- data_matrix
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        cat(sprintf("    处理 %s (%d个样本)...\n", grp, length(grp_cols)))
        sub_mat <- as.matrix(data_annotated[, grp_cols, drop = FALSE])
        normalized_sub_mat <- vsn_transform_preserve_na(sub_mat)
        colnames(normalized_sub_mat) <- grp_cols
        vsn_mat[, grp_cols] <- normalized_sub_mat
      }
      
      data_local_vsn <- data_log2
      data_local_vsn[, data_cols] <- as.data.frame(vsn_mat)
      
      standardized_data_list[["Local_VSN"]] <- data_local_vsn
      
      csv_file <- file.path(dir_config$output, "Module04_Local_VSN.csv")
      write.csv(data_local_vsn, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_VSN完成，已保存: %s\n", csv_file))
    }
    
    # 6.2e Local Reference Normalization（上下文感知） ####
    if ("Local_RefNorm" %in% norm_types) {
      cat("\n  执行Local Reference Normalization（上下文感知）...\n")
      
      refnorm_mat <- data_matrix
      sample_ctx_map <- setNames(sampleGroup$Context, sampleGroup$FinalName)
      
      for (grp in unique(group_vec)) {
        grp_cols <- names(group_vec[group_vec == grp])
        context_vec <- sample_ctx_map[grp_cols]
        only_control <- all(context_vec == "Control", na.rm = TRUE) && any(!is.na(context_vec))
        mode_label <- if (only_control) "Control模式（中位数就近）" else "Experiment优先"
        cat(sprintf("    处理 %s (%d个样本, %s)...\n", grp, length(grp_cols), mode_label))
        
        sub_mat <- data_matrix[, grp_cols, drop = FALSE]
        ref_col <- choose_reference_sample(grp_cols, sub_mat, sampleGroup, context_only_control = only_control)
        if (is.null(ref_col)) {
          warning(sprintf("⚠ 组 %s 未找到参考样本，跳过RefNorm", grp))
          next
        }
        cat(sprintf("      参考样本: %s (Context=%s)\n", ref_col, sample_ctx_map[[ref_col]]))
        
        ref_vec <- sub_mat[, ref_col]
        for (cn in grp_cols) {
          refnorm_mat[, cn] <- align_to_reference(sub_mat[, cn], ref_vec)
        }
      }
      
      data_local_refnorm <- data_log2
      data_local_refnorm[, data_cols] <- as.data.frame(refnorm_mat)
      
      standardized_data_list[["Local_RefNorm"]] <- data_local_refnorm
      
      csv_file <- file.path(dir_config$output, "Module04_Local_RefNorm.csv")
      write.csv(data_local_refnorm, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Local_RefNorm完成，已保存: %s\n", csv_file))
    }
    
    # 绘制局部标准化boxplot
    pdf_file <- file.path(dir_config$output, "Module04_Local_Norm_boxplot.pdf")
    pdf(pdf_file, width = 10, height = 5)
    
    for (norm_name in c("Local_QNorm", "Local_QNorm_Fix", "Local_QNorm_Limma", "Local_VSN", "Local_RefNorm", "Local_MNorm")) {
      if (norm_name %in% names(standardized_data_list)) {
        boxplot(standardized_data_list[[norm_name]][, data_cols], 
                cex.axis = 0.4, 
                las = 2, 
                main = norm_name,
                col = box_colors[data_cols])
        # 添加图例
        legend("topright", 
               legend = names(group_colors), 
               fill = group_colors, 
               cex = 0.6,
               title = "bioGroup")
      }
    }
    
    dev.off()
    cat(sprintf("✓ 已保存: %s\n", pdf_file))
  }
  
  # 7. 综合对比boxplot ####
  cat("\n[6] 生成综合对比boxplot...\n")
  
  # 按照用户指定的顺序重新排列
  standardized_data_list <- standardized_data_list[intersect(norm_types, names(standardized_data_list))]
  
  pdf_file <- file.path(dir_config$output, "Module04_All_Norm_comparison_boxplot.pdf")
  pdf(pdf_file, width = 10, height = 5)
  
  for (norm_name in names(standardized_data_list)) {
    boxplot(standardized_data_list[[norm_name]][, data_cols], 
            cex.axis = 0.4, 
            las = 2, 
            main = norm_name,
            col = box_colors[data_cols])
    # 添加图例
    legend("topright", 
           legend = names(group_colors), 
           fill = group_colors, 
           cex = 0.6,
           title = "bioGroup")
  }
  
  dev.off()
  cat(sprintf("✓ 已保存: %s\n", pdf_file))
  
  # 8. 总结输出 ####
  cat("\n=== Module 04 完成 ===\n")
  cat(sprintf("✓ 已生成 %d 种标准化版本:\n", length(standardized_data_list)))
  for (norm_name in names(standardized_data_list)) {
    cat(sprintf("  - %s\n", norm_name))
  }
  
  cat("\n输出文件:\n")
  cat(sprintf("  - CSV文件: Output/Module04_*.csv (%d个)\n", 
              length(standardized_data_list) + 1))  # +1 for log2
  cat(sprintf("  - PDF文件: Output/Module04_*_boxplot.pdf (4个)\n"))
  
  # 返回结果
  return(list(
    standardized_data_list = standardized_data_list,
    norm_types_used = names(standardized_data_list)
  ))
}

