# Module 06: 热图系统
# 功能：为填补后的数据绘制多种类型的热图，包括相关性热图和按定位分类的热图
# 作者：CodeNorm Pipeline
# 日期：2024

#' Module 06: 热图系统
#' 
#' @description
#' 为填补后的数据绘制热图，包括：
#' 1. 相关性热图（Correlation heatmap）- 支持详细定制
#' 2. AllLocalization热图（所有蛋白）
#' 3. 按Localization分类的热图（SGs, Mitochondrion, Cytosol, Nuclear, Nuclear_Cytosol, Other）
#' 
#' @param dir_config 目录配置列表（来自Module 1）
#' @param imputed_data_list 填补后的数据列表（来自Module 5）
#' @param sampleGroup 样本分组信息（来自Module 2）
#' @param selected_versions 选择用于绘制热图的标准化版本名称（支持向量）
#'   例如：c("noNorm_Imputed", "Local_QNorm_Imputed")
#' @param heatmap_types 需要绘制的热图类型向量
#'   可选值：c("correlation", "all", "by_localization")
#'   默认：c("all", "by_localization")
#' @param correlation_config 相关性热图配置列表（默认NULL使用默认配置）
#'   包含：
#'   - corr_min: 相关性最小值（默认0.8）
#'   - corr_max: 相关性最大值（默认1）
#'   - corr_center: 相关性中心值（默认0.9）
#'   - col_low: 低值颜色（默认"#4D97CD"）
#'   - col_center: 中心颜色（默认"white"）
#'   - col_high: 高值颜色（默认"#DB6968"）
#'   - exclude_context: 排除的Context向量（默认NULL）
#'   - exclude_pltype: 排除的PLtype向量（默认NULL）
#'   - exclude_catalytic: 排除的CatalyticGroup向量（默认NULL）
#'   - include_biogroups: 直接指定包含的bioGroup向量（默认NULL，使用排除规则）
#' @param localization_columns 用于分类的注释列名（默认自动检测）
#' @param color_params 热图颜色参数列表（用于AllLocalization和by_localization热图）
#' 
#' @return 列表，包含：
#'   - heatmap_params: 热图参数信息
#'   - annotation_info: 注释信息
#' 
#' @export
module06_heatmap <- function(dir_config,
                             imputed_data_list,
                             sampleGroup,
                             selected_versions = NULL,
                             heatmap_types = c("all", "by_localization"),
                             correlation_config = NULL,
                             localization_columns = NULL,
                             color_params = NULL) {
  
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(scales)
  
  # 1. 输入验证 ####
  cat("\n[1] 验证输入数据...\n")
  
  if (!all(c("output") %in% names(dir_config))) {
    stop("❌ dir_config必须包含output路径")
  }
  
  if (length(imputed_data_list) == 0) {
    stop("❌ imputed_data_list为空")
  }
  
  required_cols <- c("FinalName", "bioGroup", "PLtype", "Context", "CatalyticGroup")
  if (!is.data.frame(sampleGroup) || !all(required_cols %in% names(sampleGroup))) {
    stop(sprintf("❌ sampleGroup必须包含以下列: %s", paste(required_cols, collapse = ", ")))
  }
  
  # 自动选择版本
  if (is.null(selected_versions)) {
    if ("noNorm_Imputed" %in% names(imputed_data_list)) {
      selected_versions <- "noNorm_Imputed"
    } else {
      selected_versions <- names(imputed_data_list)[1]
    }
  }
  
  # 验证所有选择的版本存在
  missing_versions <- setdiff(selected_versions, names(imputed_data_list))
  if (length(missing_versions) > 0) {
    stop(sprintf("❌ 以下版本不存在: %s", paste(missing_versions, collapse = ", ")))
  }
  
  cat("✓ 输入验证通过\n")
  cat(sprintf("  - 可用数据集: %d 个\n", length(imputed_data_list)))
  cat(sprintf("  - 选择版本: %s\n", paste(selected_versions, collapse = ", ")))
  cat(sprintf("  - 热图类型: %s\n", paste(heatmap_types, collapse = ", ")))
  
  # 2. 设置默认相关性配置 ####
  if (is.null(correlation_config)) {
    correlation_config <- list(
      corr_min = 0.8,
      corr_max = 1,
      corr_center = 0.9,
      col_low = "#4D97CD",
      col_center = "white",
      col_high = "#DB6968",
      exclude_context = NULL,
      exclude_pltype = NULL,
      exclude_catalytic = NULL,
      include_biogroups = NULL
    )
  } else {
    # 补充缺失的默认值
    if (is.null(correlation_config$corr_min)) correlation_config$corr_min <- 0.8
    if (is.null(correlation_config$corr_max)) correlation_config$corr_max <- 1
    if (is.null(correlation_config$corr_center)) correlation_config$corr_center <- 0.9
    if (is.null(correlation_config$col_low)) correlation_config$col_low <- "#4D97CD"
    if (is.null(correlation_config$col_center)) correlation_config$col_center <- "white"
    if (is.null(correlation_config$col_high)) correlation_config$col_high <- "#DB6968"
  }
  
  # 3. generate_pheatmap_params函数（内部函数）####
  generate_pheatmap_params <- function(mat, 
                                       custom_min = NULL, 
                                       custom_max = NULL, 
                                       custom_center = NULL, 
                                       col_low = "blue", 
                                       col_high = "red", 
                                       col_center = "white") {
    
    # 步骤1: 决定最终的数据范围
    data_min <- if (is.null(custom_min)) min(mat, na.rm = TRUE) else custom_min
    data_max <- if (is.null(custom_max)) max(mat, na.rm = TRUE) else custom_max
    
    # 步骤2: 决定最终的颜色中点
    center_val <- if (is.null(custom_center)) (data_min + data_max) / 2 else custom_center
    
    # 步骤3: 生成breaks
    n_breaks <- 100
    
    # 如果center在范围内，分两段生成
    if (center_val > data_min && center_val < data_max) {
      # 左侧：data_min -> center_val
      breaks_left <- seq(data_min, center_val, length.out = n_breaks / 2)
      # 右侧：center_val -> data_max
      breaks_right <- seq(center_val, data_max, length.out = n_breaks / 2)
      # 合并并去重
      my_breaks <- unique(c(breaks_left, breaks_right))
      
      # 对应的颜色渐变
      colors_left <- colorRampPalette(c(col_low, col_center))(n_breaks / 2)
      colors_right <- colorRampPalette(c(col_center, col_high))(n_breaks / 2)
      my_colors <- c(colors_left, colors_right)
    } else {
      # center不在范围内，直接从data_min到data_max
      my_breaks <- seq(data_min, data_max, length.out = n_breaks)
      my_colors <- colorRampPalette(c(col_low, col_high))(n_breaks)
    }
    
    return(list(breaks = my_breaks, color = my_colors))
  }
  
  # 4. 筛选用于相关性热图的样本 ####
  if ("correlation" %in% heatmap_types) {
    cat("\n[2] 筛选相关性热图样本...\n")
    
    # 如果指定了include_biogroups，直接使用
    if (!is.null(correlation_config$include_biogroups)) {
      filtered_sampleGroup <- sampleGroup %>%
        filter(bioGroup %in% correlation_config$include_biogroups)
      
      cat(sprintf("  - 使用指定的bioGroup: %s\n", 
                  paste(correlation_config$include_biogroups, collapse = ", ")))
    } else {
      # 否则使用排除规则
      filtered_sampleGroup <- sampleGroup
      
      if (!is.null(correlation_config$exclude_context)) {
        filtered_sampleGroup <- filtered_sampleGroup %>%
          filter(!Context %in% correlation_config$exclude_context)
        cat(sprintf("  - 排除Context: %s\n", 
                    paste(correlation_config$exclude_context, collapse = ", ")))
      }
      
      if (!is.null(correlation_config$exclude_pltype)) {
        filtered_sampleGroup <- filtered_sampleGroup %>%
          filter(!PLtype %in% correlation_config$exclude_pltype)
        cat(sprintf("  - 排除PLtype: %s\n", 
                    paste(correlation_config$exclude_pltype, collapse = ", ")))
      }
      
      if (!is.null(correlation_config$exclude_catalytic)) {
        filtered_sampleGroup <- filtered_sampleGroup %>%
          filter(!CatalyticGroup %in% correlation_config$exclude_catalytic)
        cat(sprintf("  - 排除CatalyticGroup: %s\n", 
                    paste(correlation_config$exclude_catalytic, collapse = ", ")))
      }
    }
    
    corr_samples <- filtered_sampleGroup$FinalName
    
    if (length(corr_samples) == 0) {
      warning("⚠ 筛选后没有样本用于相关性热图，跳过相关性热图")
      heatmap_types <- setdiff(heatmap_types, "correlation")
    } else {
      cat(sprintf("✓ 筛选完成，保留 %d 个样本\n", length(corr_samples)))
      
      # 5. 准备相关性热图的列注释和颜色 ####
      cat("\n[3] 准备相关性热图列注释...\n")
      
      # 为PLtype和Context分别生成颜色映射
      unique_pltypes <- unique(filtered_sampleGroup$PLtype)
      unique_contexts <- unique(filtered_sampleGroup$Context)
      
      # PLtype颜色：为每个实际出现的PLtype生成颜色
      pltype_colors <- setNames(
        hue_pal()(length(unique_pltypes)),
        unique_pltypes
      )
      
      # Context颜色：根据实际出现的Context生成
      context_colors <- setNames(
        c("#CC0000", "#FFA07A")[1:length(unique_contexts)],
        unique_contexts
      )
      
      # 创建列注释数据框（确保factor levels与颜色定义匹配）
      Annotated_col_corr <- data.frame(
        PLtype = factor(filtered_sampleGroup$PLtype, levels = unique_pltypes),
        Context = factor(filtered_sampleGroup$Context, levels = unique_contexts),
        stringsAsFactors = FALSE
      )
      rownames(Annotated_col_corr) <- filtered_sampleGroup$FinalName
      
      # 定义列注释颜色规则
      # PLtype颜色：Light=蓝色系，H2O2=橘黄色系，PL=绿色系
      # Context颜色：Experiment=深色，Control=浅色
      # Spatial=黑色（优先级最高）
      
      # 为每个样本分配组合颜色（PLtype×Context）
      sample_colors <- character(nrow(filtered_sampleGroup))
      names(sample_colors) <- filtered_sampleGroup$FinalName
      
      for (i in 1:nrow(filtered_sampleGroup)) {
        pltype <- filtered_sampleGroup$PLtype[i]
        context <- filtered_sampleGroup$Context[i]
        biogroup <- filtered_sampleGroup$bioGroup[i]
        
        # Spatial优先级最高
        if (grepl("Spatial", biogroup, ignore.case = TRUE)) {
          sample_colors[i] <- "black"
        } else {
          # 根据PLtype和Context组合分配颜色
          if (grepl("Light", pltype, ignore.case = TRUE)) {
            # Light - 蓝色系
            if (context == "Experiment") {
              sample_colors[i] <- "#0066CC"  # 深蓝
            } else {
              sample_colors[i] <- "#87CEEB"  # 浅蓝（天蓝）
            }
          } else if (grepl("H2O2", pltype, ignore.case = TRUE)) {
            # H2O2 - 橘黄色系
            if (context == "Experiment") {
              sample_colors[i] <- "#FF8C00"  # 深橙
            } else {
              sample_colors[i] <- "#FFD700"  # 浅橙（金色）
            }
          } else if (grepl("PL", pltype, ignore.case = TRUE)) {
            # PL - 绿色系
            if (context == "Experiment") {
              sample_colors[i] <- "#228B22"  # 深绿
            } else {
              sample_colors[i] <- "#90EE90"  # 浅绿
            }
          } else {
            # 其他 - 灰色系
            if (context == "Experiment") {
              sample_colors[i] <- "#696969"  # 深灰
            } else {
              sample_colors[i] <- "#D3D3D3"  # 浅灰
            }
          }
        }
      }
      
      # 为Sample_Color创建唯一的颜色映射
      unique_sample_colors <- unique(sample_colors)
      sample_color_mapping <- setNames(unique_sample_colors, unique_sample_colors)
      
      # 添加样本级别的颜色（用于主热图的列颜色条）
      Annotated_col_corr$Sample_Color <- factor(sample_colors[filtered_sampleGroup$FinalName], 
                                                 levels = unique_sample_colors)
      
      # 组装注释颜色列表
      anno_colors_corr <- list(
        PLtype = pltype_colors,
        Context = context_colors,
        Sample_Color = sample_color_mapping
      )
      
      cat("✓ 列注释准备完成\n")
      cat(sprintf("  - PLtype类别: %d 个\n", length(unique_pltypes)))
      cat(sprintf("  - Context类别: %d 个\n", length(unique_contexts)))
    }
  }
  
  # 6. 处理每个选择的版本 ####
  for (selected_version in selected_versions) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("处理版本: %s\n", selected_version))
    cat(sprintf("========================================\n"))
    
    mydata <- imputed_data_list[[selected_version]]
    
    # 识别基因列
    gene_col <- if ("Gene" %in% names(mydata)) "Gene" else names(mydata)[1]
    
    # 识别注释列
    anno_cols <- grep("Localization|Annotation", names(mydata), value = TRUE)
    
    # 识别数据列
    data_cols <- setdiff(names(mydata), c(gene_col, anno_cols))
    data_cols <- data_cols[sapply(mydata[data_cols], is.numeric)]
    
    cat(sprintf("  - 基因列: %s\n", gene_col))
    cat(sprintf("  - 数据列: %d 个\n", length(data_cols)))
    cat(sprintf("  - 注释列: %d 个\n", length(anno_cols)))
    
    # 7. 绘制相关性热图 ####
    if ("correlation" %in% heatmap_types) {
      cat("\n[4] 绘制相关性热图...\n")
      
      # 筛选相关性热图的列
      corr_data_cols <- intersect(data_cols, corr_samples)
      
      if (length(corr_data_cols) < 2) {
        warning("⚠ 用于相关性热图的列少于2个，跳过")
      } else {
        pdf_file <- file.path(dir_config$output, 
                             sprintf("Module06_Correlation_heatmap_%s.pdf", selected_version))
        pdf(pdf_file, width = 14, height = 10)
        
        # 准备数据矩阵
        Pearson_matrix <- mydata[, corr_data_cols]
        rownames(Pearson_matrix) <- mydata[[gene_col]]
        
        # 计算Pearson相关性
        r <- cor(Pearson_matrix,
                 method = "pearson",
                 use = "complete.obs")
        
        # 准备颜色参数
        corr_params <- generate_pheatmap_params(
          r,
          custom_min = correlation_config$corr_min,
          custom_max = correlation_config$corr_max,
          custom_center = correlation_config$corr_center,
          col_low = correlation_config$col_low,
          col_center = correlation_config$col_center,
          col_high = correlation_config$col_high
        )
        
        # 筛选列注释
        Annotated_col_corr_filtered <- Annotated_col_corr[corr_data_cols, , drop = FALSE]
        
        # 绘制热图
        pheatmap(r,
                 show_colnames = TRUE,
                 show_rownames = TRUE,
                 fontsize = 10,
                 color = corr_params$color,
                 breaks = corr_params$breaks,
                 border_color = NA,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 clustering_method = "ward.D2",
                 clustering_distance_cols = "correlation",
                 clustering_distance_rows = "correlation",
                 annotation_col = Annotated_col_corr_filtered,
                 annotation_row = Annotated_col_corr_filtered,
                 annotation_colors = anno_colors_corr,
                 main = paste("Correlation Heatmap -", selected_version))
        
        dev.off()
        cat(sprintf("✓ 已保存: %s\n", pdf_file))
      }
    }
    
    # 8. 准备AllLocalization和by_localization热图的注释 ####
    if ("all" %in% heatmap_types || "by_localization" %in% heatmap_types) {
      cat("\n[5] 准备AllLocalization热图数据...\n")
      
      # 自动设置Localization列
      if (is.null(localization_columns)) {
        localization_columns <- anno_cols
      }
      
      # 验证Localization列存在
      valid_loc_cols <- intersect(localization_columns, anno_cols)
      if (length(valid_loc_cols) == 0) {
        warning("⚠ 未找到有效的Localization列，将只绘制基础热图")
        localization_columns <- NULL
      } else {
        localization_columns <- valid_loc_cols
      }
      
      # 设置默认颜色参数
      if (is.null(color_params)) {
        color_params <- list(
          custom_min = -4,
          custom_max = 4,
          custom_center = 0,
          col_low = "#4D97CD",
          col_center = "white",
          col_high = "#DB6968"
        )
      }
      
      # 计算颜色参数
      custom_params <- generate_pheatmap_params(
        mydata[data_cols],
        custom_min = color_params$custom_min,
        custom_max = color_params$custom_max,
        custom_center = color_params$custom_center,
        col_low = color_params$col_low,
        col_center = color_params$col_center,
        col_high = color_params$col_high
      )
      
      # 去除NA行
      data_clean <- mydata %>% na.omit()
      cat(sprintf("  - 原始行数: %d, 去除NA后: %d\n", nrow(mydata), nrow(data_clean)))
      
      if (nrow(data_clean) == 0) {
        warning("⚠ 去除NA后没有数据，跳过AllLocalization热图")
      } else {
        # 准备行注释
        if (!is.null(localization_columns)) {
          Annotated_row <- data.frame(row.names = data_clean[[gene_col]])
          
          for (col in localization_columns) {
            Annotated_row[[col]] <- factor(
              data_clean[[col]],
              levels = c("SGs", "Nuclear", "Mitochondrion", "Cytosol", "Nuclear_Cytosol", "Other")
            )
          }
          
          # 行注释颜色
          anno_colors_list <- list()
          loc_colors <- c(
            SGs = "red",
            Nuclear = "blue",
            Mitochondrion = "black",
            Cytosol = "green",
            Nuclear_Cytosol = "yellow",
            Other = "white"
          )
          
          for (col in localization_columns) {
            anno_colors_list[[col]] <- loc_colors
          }
        } else {
          Annotated_row <- NULL
          anno_colors_list <- list()
        }
        
        # 准备列注释
        sample_to_group <- setNames(sampleGroup$bioGroup, sampleGroup$FinalName)
        valid_samples <- intersect(data_cols, names(sample_to_group))
        
        if (length(valid_samples) > 0) {
          Annotated_col <- data.frame(
            Sample_Group = factor(sample_to_group[valid_samples]),
            row.names = valid_samples
          )
          
          # 列注释颜色
          unique_groups <- unique(sampleGroup$bioGroup)
          n_groups <- length(unique_groups)
          
          if (n_groups <= 12) {
            group_colors <- hue_pal()(n_groups)
          } else {
            group_colors <- rainbow(n_groups)
          }
          names(group_colors) <- unique_groups
          
          anno_colors_list$Sample_Group <- group_colors
        } else {
          Annotated_col <- NULL
        }
        
        # 9. 绘制AllLocalization热图 ####
        if ("all" %in% heatmap_types) {
          cat("\n[6] 绘制AllLocalization热图...\n")
          
          pdf_file <- file.path(dir_config$output, 
                               sprintf("Module06_AllLocalization_heatmap_%s.pdf", selected_version))
          pdf(pdf_file, width = 10, height = 10)
          
          Dataforfigure <- data_clean[, data_cols]
          rownames(Dataforfigure) <- data_clean[[gene_col]]
          
          pheatmap(Dataforfigure,
                   cluster_cols = TRUE,
                   cluster_rows = TRUE,
                   clustering_method = "ward.D2",
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "correlation",
                   color = custom_params$color,
                   breaks = custom_params$breaks,
                   fontsize_row = 0.5,
                   scale = "row",
                   border_color = NA,
                   annotation_row = Annotated_row,
                   annotation_col = Annotated_col,
                   annotation_colors = anno_colors_list,
                   main = paste("Localization Heatmap -", selected_version))
          
          dev.off()
          cat(sprintf("✓ 已保存: %s\n", pdf_file))
        }
        
        # 10. 按Localization分类绘制热图 ####
        if ("by_localization" %in% heatmap_types && !is.null(localization_columns)) {
          cat("\n[7] 按Localization分类绘制热图...\n")
          
          loc_categories <- c("SGs", "Mitochondrion", "Cytosol", "Nuclear", "Nuclear_Cytosol", "Other")
          
          for (loc_col in localization_columns) {
            cat(sprintf("  - 处理 %s...\n", loc_col))
            
            pdf_file <- file.path(dir_config$output, 
                                 sprintf("Module06_Different_%s_heatmap_%s.pdf", loc_col, selected_version))
            pdf(pdf_file, width = 10, height = 10)
            
            for (loc_cat in loc_categories) {
              mydata_filtered <- data_clean %>% filter(.data[[loc_col]] == loc_cat)
              
              if (nrow(mydata_filtered) > 0) {
                Dataforfigure <- mydata_filtered[, data_cols]
                rownames(Dataforfigure) <- mydata_filtered[[gene_col]]
                
                if (!is.null(Annotated_row)) {
                  Annotated_row_filtered <- Annotated_row[mydata_filtered[[gene_col]], , drop = FALSE]
                } else {
                  Annotated_row_filtered <- NULL
                }
                
                tryCatch({
                  pheatmap(Dataforfigure,
                           cluster_cols = TRUE,
                           cluster_rows = TRUE,
                           clustering_method = "ward.D2",
                           clustering_distance_cols = "correlation",
                           clustering_distance_rows = "correlation",
                           color = custom_params$color,
                           breaks = custom_params$breaks,
                           fontsize_row = 0.5,
                           scale = "row",
                           border_color = NA,
                           annotation_row = Annotated_row_filtered,
                           annotation_col = Annotated_col,
                           annotation_colors = anno_colors_list,
                           main = paste(loc_col, "Heatmap -", loc_cat, "-", selected_version))
                }, error = function(e) {
                  cat(sprintf("    ⚠ 跳过 %s（%s）: %s\n", loc_col, loc_cat, e$message))
                })
              }
            }
            
            dev.off()
            cat(sprintf("  ✓ 已保存: %s\n", pdf_file))
          }
        }
      }
    }
  }
  
  # 11. 返回结果 ####
  cat("\n[8] 热图绘制完成\n")
  
  heatmap_info <- list(
    heatmap_params = list(
      selected_versions = selected_versions,
      heatmap_types = heatmap_types,
      localization_columns = localization_columns,
      correlation_config = correlation_config,
      color_params = color_params
    ),
    annotation_info = list(
      anno_cols = anno_cols,
      data_cols = data_cols,
      gene_col = gene_col
    )
  )
  
  cat("✓ Module 6 完成\n")
  
  return(invisible(heatmap_info))
}
