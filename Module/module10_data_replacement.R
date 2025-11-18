# ============================================================================
# Module 10: 数据替换（基于Module 5的全局标准化+填补结果）
# ============================================================================
# 对齐 CleanCode.R (2303-2361) 逻辑
#
# 目标：
# - 对 Module 9 生成的 merged_data_A / merged_data_B（AfterROC_merged）
#   执行定点替换：保留所有 NA，所有非 NA 值用 Module 5 的
#   Global_QNorm_Imputed 对应值替换；若无该版本，则用 Global_MNorm_Imputed。
# - 生成抽样校验报告（存放在Check/子目录），验证替换是否成功。
# - 导出替换后的结果工作簿（Method A和Method B各一份）。
#
# 输入：
# - dir_config：包含 output 等目录路径
# - sampleGroup：用于识别样本列（FinalName）
# - merged_data_A / merged_data_B：Module 9 输出的合并数据（按 version 的列表）
# - imputed_data_list：Module 5 输出，包含 Global_QNorm_Imputed / Global_MNorm_Imputed 等版本
# - selected_versions：需要处理的版本（NULL 表示使用 merged_data_* 的全部版本）
# - prefer_version：优先使用的填补数据版本（默认 "Global_QNorm_Imputed"）
# - fallback_version：备选填补数据版本（默认 "Global_MNorm_Imputed"）
# - test_n：抽样校验的最大基因数（默认 20）
#
# 输出：
# - Output/Module10_{method}_{version}_AfterFirstROC_Intersect.xlsx：替换后的结果
# - Output/Check/Module10_{method}_{version}_replacement_check.csv：抽样校验表
# - Output/Check/Module10_{method}_{version}_replacement_summary.txt：替换计数汇总
#
# 返回：
# list(
#   replaced_data_A = list(version -> dataframe),
#   replaced_data_B = list(version -> dataframe),
#   base_version_used = "<选择的填补版本>"
# )
# ============================================================================

module10_data_replacement <- function(dir_config,
                                      sampleGroup,
                                      merged_data_A,
                                      merged_data_B,
                                      imputed_data_list,
                                      selected_versions = NULL,
                                      prefer_version = "Global_QNorm_Imputed",
                                      fallback_version = "Global_MNorm_Imputed",
                                      test_n = 20) {
  cat("\n========================================\n")
  cat("Module 10: 数据替换（使用Module 5全局标准化+填补结果）\n")
  cat("========================================\n")
  
  # 依赖
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("需要安装 dplyr 包")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("需要安装 openxlsx 包")
  
  `%>%` <- dplyr::`%>%`
  
  # 创建Check子目录
  check_dir <- file.path(dir_config$output, "Check")
  if (!dir.exists(check_dir)) {
    dir.create(check_dir, recursive = TRUE)
    cat(sprintf("✓ 创建目录：%s\n", check_dir))
  }
  
  # 选择基准填补版本
  base_version <- if (prefer_version %in% names(imputed_data_list)) {
    prefer_version
  } else if (fallback_version %in% names(imputed_data_list)) {
    fallback_version
  } else {
    stop(sprintf("在 imputed_data_list 中未找到 %s 或 %s",
                 prefer_version, fallback_version))
  }
  cat(sprintf("✓ 基准替换版本：%s\n", base_version))
  
  base_imputed <- imputed_data_list[[base_version]]
  if (!("Gene" %in% colnames(base_imputed))) {
    stop(sprintf("版本 %s 缺少 Gene 列", base_version))
  }
  
  # 需要处理的 version 集合
  existing_versions <- unique(c(names(merged_data_A), names(merged_data_B)))
  if (is.null(selected_versions)) {
    selected_versions <- existing_versions
  } else {
    missing <- setdiff(selected_versions, existing_versions)
    if (length(missing) > 0) {
      stop(sprintf("以下版本不存在于 Module 9 结果中：%s",
                   paste(missing, collapse = ", ")))
    }
  }
  cat(sprintf("✓ 处理版本：%s\n\n", paste(selected_versions, collapse = ", ")))
  
  # 工具：对一个 merged dataframe 执行替换
  replace_nonNA_cellwise <- function(df, ref_df, cols) {
    # 按Gene对齐，带入参考值
    merged <- dplyr::left_join(df, ref_df, by = "Gene", suffix = c("", ".ref"))
    
    replaced_count <- 0L
    na_preserved <- 0L
    total_cells <- 0L
    
    for (col in cols) {
      ref_col <- paste0(col, ".ref")
      if (!(ref_col %in% colnames(merged))) next
      
      v_old <- merged[[col]]
      v_ref <- merged[[ref_col]]
      
      # 统计
      mask_non_na <- !is.na(v_old)
      total_cells <- total_cells + length(v_old)
      na_preserved <- na_preserved + sum(!mask_non_na, na.rm = TRUE)
      
      # 替换逻辑：非NA用ref替换；NA保持
      v_new <- v_old
      v_new[mask_non_na & !is.na(v_ref)] <- v_ref[mask_non_na & !is.na(v_ref)]
      replaced_count <- replaced_count + sum(mask_non_na & !is.na(v_ref), na.rm = TRUE)
      merged[[col]] <- v_new
      
      # 去掉辅助列
      merged[[ref_col]] <- NULL
    }
    
    # 去掉所有残余.ref列
    merged <- merged[, !grepl("\\.ref$", colnames(merged)), drop = FALSE]
    
    return(list(
      data = merged,
      replaced_count = replaced_count,
      na_preserved = na_preserved,
      total_cells = total_cells
    ))
  }
  
  replaced_data_A <- list()
  replaced_data_B <- list()
  
  # 逐版本执行
  for (ver in selected_versions) {
    cat(sprintf("════════════════════════════════════════\n"))
    cat(sprintf("处理版本：%s\n", ver))
    cat(sprintf("════════════════════════════════════════\n"))
    
    # 准备基准表（只选择该版本需要的列）
    lfq_cols <- sampleGroup$FinalName
    base_cols <- c("Gene", lfq_cols)
    base_use <- base_imputed %>%
      dplyr::select(dplyr::any_of(base_cols)) %>%
      dplyr::distinct(Gene, .keep_all = TRUE)
    
    # --- Method A ---
    if (ver %in% names(merged_data_A)) {
      cat("\n【Method A】\n")
      df_A <- merged_data_A[[ver]]
      
      if (!is.data.frame(df_A) || nrow(df_A) == 0) {
        cat("  ⚠ 数据为空，跳过\n")
        replaced_data_A[[ver]] <- df_A
      } else {
        # 识别要替换的列（LFQ列）
        replace_cols <- intersect(lfq_cols, colnames(df_A))
        replace_cols <- intersect(replace_cols, colnames(base_use))
        
        if (length(replace_cols) == 0) {
          cat("  ⚠ 无可替换的LFQ列，跳过\n")
          replaced_data_A[[ver]] <- df_A
        } else {
          cat(sprintf("  - 原始数据：%d 个蛋白, %d 列\n", nrow(df_A), ncol(df_A)))
          cat(sprintf("  - 替换列数：%d 个LFQ列\n", length(replace_cols)))
          
          # 执行替换
          result_A <- replace_nonNA_cellwise(df_A, base_use, replace_cols)
          replaced_data_A[[ver]] <- result_A$data
          
          cat(sprintf("  ✓ 替换完成：%d 个非NA单元格被替换\n", result_A$replaced_count))
          cat(sprintf("  ✓ 保留NA：%d/%d 单元格\n", result_A$na_preserved, result_A$total_cells))
          
          # 抽样校验
          set.seed(123)
          show_cols <- head(replace_cols, 3)  # 最多3列
          check_genes <- head(unique(df_A$Gene), test_n)
          
          if (length(check_genes) > 0 && length(show_cols) > 0) {
            # 原始值
            orig_df <- df_A %>% 
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(orig_df) <- c("Gene", paste0("orig_", show_cols))
            
            # 替换后值
            repl_df <- result_A$data %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(repl_df) <- c("Gene", paste0("repl_", show_cols))
            
            # 参考值
            imp_df <- base_use %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(imp_df) <- c("Gene", paste0("imp_", show_cols))
            
            # 合并
            check_out <- Reduce(function(x, y) dplyr::left_join(x, y, by = "Gene"),
                               list(orig_df, repl_df, imp_df))
            
            # 转换为长表
            long_list <- list()
            for (sc in show_cols) {
              long_list[[sc]] <- data.frame(
                Gene = check_out$Gene,
                Sample = sc,
                Original = check_out[[paste0("orig_", sc)]],
                Replaced = check_out[[paste0("repl_", sc)]],
                Imputed = check_out[[paste0("imp_", sc)]],
                stringsAsFactors = FALSE
              )
            }
            check_long <- do.call(rbind, long_list)
            
            # 判断匹配（非NA原值的Replaced应等于Imputed）
            check_long$CellMatch <- with(check_long, ifelse(
              !is.na(Original),
              abs(as.numeric(Replaced) - as.numeric(Imputed)) < 1e-9,
              NA
            ))
            
            # 保存校验CSV到Check/子目录
            check_file <- file.path(check_dir, 
                                   sprintf("Module10_A_%s_replacement_check.csv", ver))
            utils::write.csv(check_long, check_file, row.names = FALSE)
            cat(sprintf("  ✓ 抽样校验导出：Check/%s\n", basename(check_file)))
          }
          
          # 汇总文本
          summary_lines <- c(
            sprintf("版本：%s - Method A", ver),
            sprintf("原始数据：%d 蛋白, %d 列", nrow(df_A), ncol(df_A)),
            sprintf("替换列数：%d LFQ列", length(replace_cols)),
            sprintf("替换结果：%d 个非NA单元格被替换", result_A$replaced_count),
            sprintf("NA保留：%d/%d 单元格", result_A$na_preserved, result_A$total_cells),
            sprintf("基准版本：%s", base_version)
          )
          sum_file <- file.path(check_dir,
                               sprintf("Module10_A_%s_replacement_summary.txt", ver))
          writeLines(summary_lines, con = sum_file)
          cat(sprintf("  ✓ 替换汇总导出：Check/%s\n", basename(sum_file)))
          
          # 导出替换后的数据
          wb <- openxlsx::createWorkbook()
          openxlsx::addWorksheet(wb, sheetName = "AfterROC_merged")
          openxlsx::writeData(wb, sheet = "AfterROC_merged", result_A$data)
          out_xlsx <- file.path(dir_config$output,
                               sprintf("Module10_A_%s_AfterFirstROC_Intersect.xlsx", ver))
          openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
          cat(sprintf("  ✓ 导出替换结果：%s\n\n", basename(out_xlsx)))
        }
      }
    } else {
      cat("\n【Method A】\n")
      cat("  ⚠ 提示：Method A 无该版本，跳过\n\n")
    }
    
    # --- Method B ---
    if (ver %in% names(merged_data_B)) {
      cat("【Method B】\n")
      df_B <- merged_data_B[[ver]]
      
      if (!is.data.frame(df_B) || nrow(df_B) == 0) {
        cat("  ⚠ 数据为空，跳过\n")
        replaced_data_B[[ver]] <- df_B
      } else {
        # 识别要替换的列（LFQ列）
        replace_cols <- intersect(lfq_cols, colnames(df_B))
        replace_cols <- intersect(replace_cols, colnames(base_use))
        
        if (length(replace_cols) == 0) {
          cat("  ⚠ 无可替换的LFQ列，跳过\n")
          replaced_data_B[[ver]] <- df_B
        } else {
          cat(sprintf("  - 原始数据：%d 个蛋白, %d 列\n", nrow(df_B), ncol(df_B)))
          cat(sprintf("  - 替换列数：%d 个LFQ列\n", length(replace_cols)))
          
          # 执行替换
          result_B <- replace_nonNA_cellwise(df_B, base_use, replace_cols)
          replaced_data_B[[ver]] <- result_B$data
          
          cat(sprintf("  ✓ 替换完成：%d 个非NA单元格被替换\n", result_B$replaced_count))
          cat(sprintf("  ✓ 保留NA：%d/%d 单元格\n", result_B$na_preserved, result_B$total_cells))
          
          # 抽样校验
          set.seed(123)
          show_cols <- head(replace_cols, 3)  # 最多3列
          check_genes <- head(unique(df_B$Gene), test_n)
          
          if (length(check_genes) > 0 && length(show_cols) > 0) {
            # 原始值
            orig_df <- df_B %>% 
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(orig_df) <- c("Gene", paste0("orig_", show_cols))
            
            # 替换后值
            repl_df <- result_B$data %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(repl_df) <- c("Gene", paste0("repl_", show_cols))
            
            # 参考值
            imp_df <- base_use %>%
              dplyr::filter(Gene %in% check_genes) %>%
              dplyr::select(Gene, dplyr::all_of(show_cols))
            names(imp_df) <- c("Gene", paste0("imp_", show_cols))
            
            # 合并
            check_out <- Reduce(function(x, y) dplyr::left_join(x, y, by = "Gene"),
                               list(orig_df, repl_df, imp_df))
            
            # 转换为长表
            long_list <- list()
            for (sc in show_cols) {
              long_list[[sc]] <- data.frame(
                Gene = check_out$Gene,
                Sample = sc,
                Original = check_out[[paste0("orig_", sc)]],
                Replaced = check_out[[paste0("repl_", sc)]],
                Imputed = check_out[[paste0("imp_", sc)]],
                stringsAsFactors = FALSE
              )
            }
            check_long <- do.call(rbind, long_list)
            
            # 判断匹配（非NA原值的Replaced应等于Imputed）
            check_long$CellMatch <- with(check_long, ifelse(
              !is.na(Original),
              abs(as.numeric(Replaced) - as.numeric(Imputed)) < 1e-9,
              NA
            ))
            
            # 保存校验CSV到Check/子目录
            check_file <- file.path(check_dir, 
                                   sprintf("Module10_B_%s_replacement_check.csv", ver))
            utils::write.csv(check_long, check_file, row.names = FALSE)
            cat(sprintf("  ✓ 抽样校验导出：Check/%s\n", basename(check_file)))
          }
          
          # 汇总文本
          summary_lines <- c(
            sprintf("版本：%s - Method B", ver),
            sprintf("原始数据：%d 蛋白, %d 列", nrow(df_B), ncol(df_B)),
            sprintf("替换列数：%d LFQ列", length(replace_cols)),
            sprintf("替换结果：%d 个非NA单元格被替换", result_B$replaced_count),
            sprintf("NA保留：%d/%d 单元格", result_B$na_preserved, result_B$total_cells),
            sprintf("基准版本：%s", base_version)
          )
          sum_file <- file.path(check_dir,
                               sprintf("Module10_B_%s_replacement_summary.txt", ver))
          writeLines(summary_lines, con = sum_file)
          cat(sprintf("  ✓ 替换汇总导出：Check/%s\n", basename(sum_file)))
          
          # 导出替换后的数据
          wb <- openxlsx::createWorkbook()
          openxlsx::addWorksheet(wb, sheetName = "AfterROC_merged")
          openxlsx::writeData(wb, sheet = "AfterROC_merged", result_B$data)
          out_xlsx <- file.path(dir_config$output,
                               sprintf("Module10_B_%s_AfterFirstROC_Intersect.xlsx", ver))
          openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
          cat(sprintf("  ✓ 导出替换结果：%s\n\n", basename(out_xlsx)))
        }
      }
    } else {
      cat("【Method B】\n")
      cat("  ⚠ 提示：Method B 无该版本，跳过\n\n")
    }
  }
  
  cat("════════════════════════════════════════\n")
  cat("✓ Module 10 完成\n")
  cat("════════════════════════════════════════\n")
  
  return(list(
    replaced_data_A = replaced_data_A,
    replaced_data_B = replaced_data_B,
    base_version_used = base_version
  ))
}
