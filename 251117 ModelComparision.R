#军科院 SGs 样品分析 
#1.加载包 #### R4.1
library(readr)
library(tidyverse)
library(readxl)
library(purrr)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(gridExtra)
library(pheatmap)
library(pROC)
library(purrr)
library(openxlsx)
library(preprocessCore)
library(limma)
library(colourpicker)
library(RColorBrewer)
library(writexl)


#2. 工作目录和数据加载 ####
setwd("D:/Phase_APEX2/Raw_ForFigures/BioMap")
#source("../202501017 3DAUC models base.R")
##替换为任意一个能够完成主流程至FC/Scatter plot的环境中
load("D:/Phase_APEX2/Raw_ForFigures/BioMap/ToStep21_all_data.RData") 


names(ForStep21)

names(ForStep21$noMBR_QNorm_New)

#修改ForStep21列表中两个数据框的列名
##为第2-16列在末尾数字前添加LFQ_前缀
# 创建列名修改函数
renameColumnsWithLFQ <- function(df) {
  # 获取当前列名
  col_names <- colnames(df)
  # 对第2-16列进行重命名
  for (i in 2:16) {
    if (i <= length(col_names)) {
      # 使用正则表达式在末尾数字前插入LFQ_
      col_names[i] <- gsub("([a-zA-Z_]+)([0-9]+)$", "\\1LFQ_\\2", col_names[i])
    }
  }
  colnames(df) <- col_names
  return(df)
}

# 应用函数到ForStep21的所有数据框
#for (i in seq_along(ForStep21)) {
#  ForStep21[[i]] <- renameColumnsWithLFQ(ForStep21[[i]])
#}

str(ForStep21)

## 分组结构
group_info <- list(
  CD3EAP = list(
    samples = c("CD3EAP_LFQ_1", "CD3EAP_LFQ_2"),
    logFC  = "CD3EAP_vs_NLS_logFC",
    logFC_FDR="CD3EAP_vs_NLS_adj.P.Val"
  ),
  NOP56 = list(
    samples = c("NOP56_LFQ_1", "NOP56_LFQ_2"),
    logFC  = "NOP56_vs_NLS_logFC",
    logFC_FDR="NOP56_vs_NLS_adj.P.Val"
  ),
  NPM1 = list(
    samples = c("NPM1_LFQ_1", "NPM1_LFQ_2"),
    logFC  = "NPM1_vs_NLS_logFC",
    logFC_FDR="NPM1_vs_NLS_adj.P.Val"
  ),
  POLR1E = list(
    samples = c("POLR1E_LFQ_1", "POLR1E_LFQ_2"),
    logFC  = "POLR1E_vs_NLS_logFC",
    logFC_FDR="POLR1E_vs_NLS_adj.P.Val"
  ),
  ZNF330 = list(
    samples = c("ZNF330_LFQ_1", "ZNF330_LFQ_2"),
    logFC  = "ZNF330_vs_NLS_logFC",
    logFC_FDR="ZNF330_vs_NLS_adj.P.Val"
  )
)

ForStep21=ForStep21[c("noMBR_QNorm_New")]

str(ForStep21)

for (i in seq_along(ForStep21)) {
  ForStep21[[i]]=ForStep21[[i]] %>% mutate(MultiBait_Localization=Sub_HPA_Localization)
}

str(ForStep21)

unique(ForStep21$noMBR_QNorm_New$Sub_HPA_Localization)
myTP_vector=c("Nucleolus","Other&Nucleolus","Nuclear&Nucleolus","Nuclear&Cytosol&Nucleolus","Cytosol&Nucleolus")
myTP_vector %in% unique(ForStep21$noMBR_QNorm_New$Sub_HPA_Localization)
ForStep19=ForStep21

#Step29 多模型复杂需求 （Gemini版本）#####

# ===================================================================
#
###    配置接口 - 用户可选设置 ####
#
# ===================================================================

#1. 组别选择配置 ####
## 可选组别：c("K69A1B3", "K69C3","C3")
## 设置为NULL表示使用所有组别，或者指定特定组别的向量
names(group_info)
selected_groups <- names(group_info)


#2. 模型选择配置 ####
## 可选模型类型：
## 1D模型: "1D_Abundance", "1D_logFC", "1D_PPI_Global", "1D_PPI_Local"
## 2D模型: "2D_Abundance_logFC", "2D_Abundance_PPI_Global", "2D_Abundance_PPI_Local", "2D_logFC_PPI_Global", "2D_logFC_PPI_Local"
## 3D模型: "3D_with_Global_PPI", "3D_with_Local_PPI"
## 传统方法: "logFC>0.5 & FDR<0.05"
## 设置为NULL表示运行所有模型，或者指定特定模型的向量
selected_models <- c(
  "logFC>0.5 & FDR<0.05",
  "1D_logFC",
  "1D_Abundance",
  "2D_Abundance_logFC")  # 可选: c("1D_Abundance", "2D_Abundance_logFC") 等


#3. 输出和可视化控制 ####
## 是否生成各类文件和图表
enable_cytoscape_export <- TRUE    # 是否导出Cytoscape文件
enable_visualization <- TRUE       # 是否生成可视化图表
enable_robustness_analysis <- TRUE # 是否进行稳健性分析
enable_group_comparison <- TRUE    # 是否进行组间比较

## 可视化子图控制
enable_plot_auc_summary <- TRUE                # 图1：AUC_summary分组森林图
enable_plot_auc_summary_vs_logFC <- TRUE       # 图1b：与1D_logFC比较的分组森林图
enable_plot_group_comparison <- TRUE           # 图2：跨邻近方法比较（柱状图+热图）
enable_plot_intersection_dumbbell <- TRUE      # 图3：样本交集公平性哑铃图
enable_plot_robustness_scatter <- TRUE         # 图4：性能-稳定性散点图
enable_plot_controlled_comparison <- TRUE      # 图5：控制变量比较增益图
enable_plot_localization_distribution <- FALSE  # 图6：候选蛋白定位分布堆叠柱状图
enable_plot_controlled_comparison_vs_logFC <- TRUE  # 图7：与1D_logFC对比的控制变量增益图

#4. 并行计算配置 ####
## 并行核心数设置（建议保持12，总24核心中使用一半）
parallel_cores <- 20

#4.1 控制变量比较图配置 ####
## 图5a-5c：指定要展示的模型
models_for_controlled_comparison <- c(
  "1D_Abundance",
  "1D_logFC", 
  "2D_Abundance_logFC"
  # 传统模型会自动作为参考，无需在此列出
)

## 图5a-5c：坐标轴范围设置
controlled_comparison_xlim <- c(-0.15, 0.5)  # Specificity/Sensitivity增益的x轴范围
controlled_comparison_ylim <- c(-0.15, 0.5)  # Sensitivity增益的y轴范围（用于5c图）

## 图5a-5c：Bootstrap重采样次数设置
n_bootstrap <- 500  # 可选: 200(快速测试), 500(常规), 1000(发表), 2000+(严格)
                     # 注意: 1000次 × 6组 × 3模型 = 约36,000次重采样，需要较长时间

#5. 验证和初始化配置 ####
cat("=== 配置验证开始 ===\n")

# 验证组别选择
available_groups <- names(group_info)
if (is.null(selected_groups)) {
  selected_groups <- available_groups
  cat("✓ 使用所有可用组别:", paste(selected_groups, collapse = ", "), "\n")
} else {
  invalid_groups <- setdiff(selected_groups, available_groups)
  if (length(invalid_groups) > 0) {
    stop("❌ 错误：无效的组别选择 - ", paste(invalid_groups, collapse = ", "), 
         "\n   可用组别：", paste(available_groups, collapse = ", "))
  }
  cat("✓ 已选择组别:", paste(selected_groups, collapse = ", "), "\n")
}

# 验证模型选择
available_models <- c(
  "1D_Abundance", "1D_logFC", "1D_PPI_Global", "1D_PPI_Local",
  "2D_Abundance_logFC", "2D_Abundance_PPI_Global", "2D_Abundance_PPI_Local", 
  "2D_logFC_PPI_Global", "2D_logFC_PPI_Local",
  "3D_with_Global_PPI", "3D_with_Local_PPI",
  "logFC>0.5 & FDR<0.05"
)
if (is.null(selected_models)) {
  selected_models <- setdiff(available_models, "logFC>0.5 & FDR<0.05")  # 默认排除传统方法
  cat("✓ 使用所有可用模型（不含传统方法）:", length(selected_models), "个\n")
} else {
  invalid_models <- setdiff(selected_models, available_models)
  if (length(invalid_models) > 0) {
    stop("❌ 错误：无效的模型选择 - ", paste(invalid_models, collapse = ", "), 
         "\n   可用模型：", paste(available_models, collapse = ", "))
  }
  cat("✓ 已选择模型:", length(selected_models), "个 -", paste(selected_models, collapse = ", "), "\n")
}

# 验证控制变量比较配置
if (!is.null(models_for_controlled_comparison)) {
  invalid_cc_models <- setdiff(models_for_controlled_comparison, available_models)
  if (length(invalid_cc_models) > 0) {
    warning("控制变量比较图配置中有无效模型: ", paste(invalid_cc_models, collapse = ", "))
    models_for_controlled_comparison <- intersect(models_for_controlled_comparison, available_models)
  }
  cat("✓ 控制变量比较图展示模型:", paste(models_for_controlled_comparison, collapse = ", "), "\n")
  cat("✓ 坐标轴范围 - X轴:", paste(controlled_comparison_xlim, collapse = " 至 "), "\n")
  cat("✓ 坐标轴范围 - Y轴:", paste(controlled_comparison_ylim, collapse = " 至 "), "\n")
  if (exists("n_bootstrap")) {
    cat("✓ Bootstrap重采样次数:", n_bootstrap, "次\n")
  }
}

# 显示其他配置
cat("✓ Cytoscape导出:", ifelse(enable_cytoscape_export, "启用", "禁用"), "\n")
cat("✓ 可视化生成:", ifelse(enable_visualization, "启用", "禁用"), "\n")
plot_toggle_status <- c(
  "图1 AUC分组森林图" = enable_plot_auc_summary,
  "图1b AUC vs 1D_logFC" = enable_plot_auc_summary_vs_logFC,
  "图2 组间比较" = enable_plot_group_comparison,
  "图3 样本交集哑铃图" = enable_plot_intersection_dumbbell,
  "图4 性能-稳定性散点图" = enable_plot_robustness_scatter,
  "图5 控制变量增益" = enable_plot_controlled_comparison,
  "图6 候选蛋白定位分布" = enable_plot_localization_distribution,
  "图7 控制变量增益 vs 1D_logFC" = enable_plot_controlled_comparison_vs_logFC
)
cat("  可视化子图控制:\n")
for (toggle_name in names(plot_toggle_status)) {
  cat("   -", toggle_name, ":", ifelse(plot_toggle_status[[toggle_name]], "启用", "禁用"), "\n")
}
cat("✓ 稳健性分析:", ifelse(enable_robustness_analysis, "启用", "禁用"), "\n")
cat("✓ 组间比较:", ifelse(enable_group_comparison, "启用", "禁用"), "\n")
cat("✓ 并行核心数:", parallel_cores, "\n")

cat("=== 配置验证完成 ===\n")
cat("=== 使用说明 ===\n")
cat("要修改配置，请编辑上述变量：\n")
cat("1. selected_groups: 选择特定组别，如 c(\"E7A2B4\", \"C2\")\n")
cat("2. selected_models: 选择特定模型，如 c(\"1D_Abundance\", \"2D_Abundance_logFC\")\n")
cat("3. enable_* 变量: 控制是否生成特定类型的输出文件和图表\n")
cat("4. parallel_cores: 调整并行计算的核心数\n")
cat("=== 开始分析 ===\n\n")

# ===================================================================
#
###    Phase 1: Data Calculation and Result Export #####
#
# ===================================================================

# -------------------------------------------------------------------
# 1. Setup and Configuration
# -------------------------------------------------------------------

# Load necessary libraries
# Ensure packages are installed: install.packages(c("dplyr", "tidyr", "pROC", "STRINGdb", "igraph", "future", "furrr", "progressr"))
library(dplyr)
library(tidyr)
library(pROC)
library(STRINGdb)
library(igraph)
library(future)
library(furrr)
library(progressr)

# Define a functional error handler as per your preference
my_error_handler <- function(e) {
  cat("错误发生:", conditionMessage(e), "\n")
}

# Setup parallel processing
# Using 20 cores as requested
plan(multisession, workers = 20)
cat(paste("并行计算已启动，将使用", nbrOfWorkers(), "个核心。\n"))

# --- Key Parameters ---
# This parameter can be adjusted for the local PPI calculation
local_ppi_window_size <- 10

# Create output directories
output_dir <- "Step29_ModelComparision"
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "Cytoscape_Files"), showWarnings = FALSE)

# -------------------------------------------------------------------
# 2. Input Data Preparation
# -------------------------------------------------------------------
# This section assumes 'ForStep19' and 'group_info' are loaded in your environment.

# -------------------------------------------------------------------
# 3. PPI Data Pre-processing (本地STRING数据库版本) ####
# -------------------------------------------------------------------
cat("开始预处理PPI数据（使用本地STRING文件）...\n")

# === 配置：本地STRING文件路径 ===
LOCAL_STRING_DIR <- file.path(getwd(), "Reference", "StringDb")
LOCAL_STRING_ALIAS_FILE <- file.path(LOCAL_STRING_DIR, "9606.protein.aliases.v12.0.txt")
LOCAL_STRING_LINKS_FILE <- file.path(LOCAL_STRING_DIR, "9606.protein.links.detailed.v12.0.txt")
SCORE_THRESHOLD <- 400  # 与原在线版本一致

# 验证文件路径
cat(paste("STRING文件目录:", LOCAL_STRING_DIR, "\n"))
cat(paste("别名文件:", LOCAL_STRING_ALIAS_FILE, "\n"))
cat(paste("互作文件:", LOCAL_STRING_LINKS_FILE, "\n"))

# === 步骤1：加载STRING别名映射表（只需加载一次）====
cat("正在加载STRING别名映射表...\n")
if (!file.exists(LOCAL_STRING_ALIAS_FILE)) {
  stop(paste("错误：找不到映射文件", LOCAL_STRING_ALIAS_FILE))
}

string_aliases <- read_tsv(LOCAL_STRING_ALIAS_FILE, 
                           col_names = c("STRING_id", "alias", "source"),
                           col_types = "ccc",
                           skip = 1,
                           show_col_types = FALSE)

# 过滤出基因名来源的别名（保留多个来源以提高映射率）
string_gene_map <- string_aliases %>%
  filter(source %in% c("BioMart_HUGO", "Ensembl_gene_name", "BLAST_UniProt_GN", 
                       "Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)")) %>%
  select(STRING_id, gene = alias) %>%
  distinct()

cat(paste("  ✅ 加载了", nrow(string_gene_map), "条基因映射记录\n"))

# === 步骤2：加载STRING互作网络（只需加载一次）====
cat("正在加载STRING互作网络...\n")
if (!file.exists(LOCAL_STRING_LINKS_FILE)) {
  stop(paste("错误：找不到互作文件", LOCAL_STRING_LINKS_FILE))
}

# STRING文件使用空格分隔，需要明确指定
string_links <- read.table(LOCAL_STRING_LINKS_FILE, 
                           header = TRUE,           # 第一行是列名
                           sep = " ",               # 空格分隔
                           stringsAsFactors = FALSE,
                           comment.char = "",       # 不处理注释
                           quote = "")              # 不处理引号

# 显示实际的列名以便调试
cat(paste("  文件实际列名:", paste(colnames(string_links)[1:5], "..."), "\n"))
cat(paste("  共", ncol(string_links), "列,", format(nrow(string_links), big.mark=","), "行\n"))

# 自动识别列名（STRING数据库可能使用不同的列名格式）
col_names <- colnames(string_links)

# 识别第一个蛋白列（通常是第1列）
protein1_col <- col_names[1]
# 识别第二个蛋白列（通常是第2列）  
protein2_col <- col_names[2]
# 识别combined_score列（可能叫combined_score或其他名称）
score_cols <- col_names[grepl("combined", col_names, ignore.case = TRUE)]
if (length(score_cols) == 0) {
  # 如果没有combined列，尝试找score列
  score_cols <- col_names[grepl("score", col_names, ignore.case = TRUE)]
  if (length(score_cols) == 0) {
    stop("错误：在文件中找不到分数列")
  }
}
score_col <- score_cols[1]

cat(paste("  使用列: [", protein1_col, "] - [", protein2_col, "] - [", score_col, "]\n"))

# 提取并重命名列
string_links <- string_links %>%
  select(protein1 = all_of(protein1_col), 
         protein2 = all_of(protein2_col), 
         combined_score = all_of(score_col)) %>%
  mutate(combined_score = as.numeric(combined_score))

# 应用分数阈值过滤
string_links_filtered <- string_links %>%
  filter(combined_score >= SCORE_THRESHOLD)

cat(paste("  ✅ 加载了", format(nrow(string_links_filtered), big.mark=","), 
          "条互作记录（分数 ≥", SCORE_THRESHOLD, "）\n"))

# === 步骤3：为每个数据源处理PPI ====
ppi_data <- list(dfs = ForStep19, interactions = list(), ppi_global = list())

for(df_name in names(ForStep19)) {
  cat(paste("正在处理数据源的PPI信息:", df_name, "\n"))
  
  tryCatch({
    # 获取蛋白列表
    protein_list <- unique(ForStep19[[df_name]]$Gene)
    cat(paste("  - 数据源包含", length(protein_list), "个蛋白\n"))
    
    # 映射基因名到STRING ID
    mapped_proteins <- string_gene_map %>%
      filter(gene %in% protein_list) %>%
      group_by(gene) %>%
      slice(1) %>%  # 每个基因只保留第一个STRING ID
      ungroup()
    
    mapping_rate <- round(nrow(mapped_proteins) / length(protein_list) * 100, 1)
    cat(paste("  - 成功映射", nrow(mapped_proteins), "/", length(protein_list), 
              "个蛋白 (", mapping_rate, "%)\n"))
    
    if (nrow(mapped_proteins) == 0) {
      cat("  ⚠️  警告：没有蛋白被成功映射到STRING ID，跳过此数据源\n")
      ppi_data$interactions[[df_name]] <- data.frame(
        protein1 = character(0), 
        protein2 = character(0), 
        combined_score = numeric(0)
      )
      ppi_data$ppi_global[[df_name]] <- data.frame(
        Gene = ForStep19[[df_name]]$Gene, 
        ppi_score_global = 0
      )
      next
    }
    
    # 筛选互作关系（只保留映射到的蛋白）
    valid_ids <- mapped_proteins$STRING_id
    interactions_filtered <- string_links_filtered %>%
      filter(protein1 %in% valid_ids & protein2 %in% valid_ids)
    
    # 将STRING ID转回基因名，并缩放分数到0-1
    ppi_data$interactions[[df_name]] <- interactions_filtered %>%
      mutate(combined_score = combined_score / 1000) %>%  # 缩放到0-1
      left_join(mapped_proteins %>% select(STRING_id, gene), 
                by = c("protein1" = "STRING_id"), 
                relationship = "many-to-many") %>%
      rename(protein1_gene = gene) %>%
      left_join(mapped_proteins %>% select(STRING_id, gene), 
                by = c("protein2" = "STRING_id"), 
                relationship = "many-to-many") %>%
      rename(protein2_gene = gene) %>%
      filter(!is.na(protein1_gene) & !is.na(protein2_gene)) %>%
      select(protein1 = protein1_gene, protein2 = protein2_gene, combined_score) %>%
      distinct()
    
    cat(paste("  - 获得", format(nrow(ppi_data$interactions[[df_name]]), big.mark=","), 
              "条互作关系\n"))
    
    # 计算全局PPI分数（每个蛋白的总互作强度）
    scores1 <- ppi_data$interactions[[df_name]] %>% 
      group_by(protein1) %>% 
      summarise(total_score = sum(combined_score)) %>% 
      rename(Gene = protein1)
    
    scores2 <- ppi_data$interactions[[df_name]] %>% 
      group_by(protein2) %>% 
      summarise(total_score = sum(combined_score)) %>% 
      rename(Gene = protein2)
    
    ppi_data$ppi_global[[df_name]] <- bind_rows(scores1, scores2) %>%
      group_by(Gene) %>%
      summarise(ppi_score_global = sum(total_score))
    
    cat(paste("  - 计算了", nrow(ppi_data$ppi_global[[df_name]]), 
              "个蛋白的全局PPI分数\n"))
    
  }, error = my_error_handler)
}

cat("✅ PPI数据预处理完成（本地版本）。\n\n")


# -------------------------------------------------------------------
# 4. Core Worker Function for Parallel Execution
# -------------------------------------------------------------------

process_group <- function(df_name, group, ppi_data) {
  
  # Load libraries within each worker for safety
  library(dplyr)
  library(pROC)
  
  # Define myTP_vector within function to avoid global variable dependency
  myTP_vector=c("Nucleolus","Other&Nucleolus","Nuclear&Nucleolus","Nuclear&Cytosol&Nucleolus","Cytosol&Nucleolus")
  
  # Define local_ppi_window_size within function
  local_ppi_window_size <- 10
  
  # --- Unpack data for the current task ---
  # Add safety checks for data availability
  if (is.null(ppi_data$dfs[[df_name]])) {
    warning(paste("No data found for df_name:", df_name))
    return(NULL)
  }
  
  df <- ppi_data$dfs[[df_name]]
  interactions_df <- ppi_data$interactions[[df_name]]
  ppi_feature_global <- ppi_data$ppi_global[[df_name]]
  
  # Handle missing interactions data gracefully
  if (is.null(interactions_df)) {
    warning(paste("No PPI interactions data for:", df_name))
    interactions_df <- data.frame(protein1 = character(0), protein2 = character(0), combined_score = numeric(0))
  }
  
  if (is.null(ppi_feature_global)) {
    warning(paste("No PPI global features for:", df_name))
    ppi_feature_global <- data.frame(Gene = character(0), ppi_score_global = numeric(0))
  }
  
  # --- Data Preparation ---
  # Check if group exists in group_info
  if (is.null(group_info[[group]])) {
    warning(paste("No group info found for group:", group))
    return(NULL)
  }
  
  samples <- group_info[[group]]$samples
  logFC_col <- group_info[[group]]$logFC
  fdr_col <- group_info[[group]]$logFC_FDR
  
  # Verify required columns exist
  if (is.null(samples) || is.null(logFC_col) || is.null(fdr_col)) {
    warning(paste("Missing required columns for group:", group))
    return(NULL)
  }
  
  # Check if required columns exist in dataframe
  missing_cols <- setdiff(c(samples, logFC_col, fdr_col, "MultiBait_Localization", "Gene"), names(df))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in dataframe for", df_name, "-", group, ":", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  df_processed <- df %>%
    mutate(mean_group = rowMeans(select(., all_of(samples)), na.rm = TRUE)) %>%
    mutate(mean_group_medianNorm = mean_group - median(mean_group, na.rm = TRUE)) %>%
    mutate(
      # Define ground truth based on MultiBait_Localization
      plot_group = if_else(
        MultiBait_Localization %in% myTP_vector,
        "TP", "Other"
      ),
      # Set factor levels to ensure consistent modeling direction
      plot_group = factor(plot_group, levels = c("Other", "TP"))
    )
  
  # Filter for proteins with valid data points for modeling
  # Corrected the !is.na call to be more robust
  plot_data <- df_processed %>% 
    filter(!is.na(mean_group_medianNorm) & !is.na(.data[[logFC_col]]) & !is.na(plot_group))
  
  # If no data remains after filtering, exit early
  if(nrow(plot_data) == 0) return(NULL)
  
  # --- Local PPI Feature Calculation ---
  sorted_genes <- plot_data %>% arrange(desc(mean_group_medianNorm)) %>% pull(Gene)
  
  local_ppi_scores <- sapply(seq_along(sorted_genes), function(i) {
    current_gene <- sorted_genes[i]
    start_index <- max(1, i - local_ppi_window_size)
    end_index <- min(length(sorted_genes), i + local_ppi_window_size)
    neighbors <- sorted_genes[start_index:end_index]
    
    # Sum scores of interactions with neighbors in the window
    local_score <- interactions_df %>% 
      filter((protein1 == current_gene & protein2 %in% neighbors) | (protein2 == current_gene & protein1 %in% neighbors)) %>% 
      pull(combined_score) %>% 
      sum()
    return(local_score)
  })
  
  ppi_feature_local <- data.frame(Gene = sorted_genes, ppi_score_local = local_ppi_scores)
  
  # --- Merge all features and perform log-transformation ---
  plot_data_full <- plot_data %>%
    left_join(ppi_feature_global, by = "Gene") %>% 
    left_join(ppi_feature_local, by = "Gene") %>%
    mutate(
      # Convert NA PPI scores to 0 (for proteins with no known interactions)
      ppi_score_global = ifelse(is.na(ppi_score_global), 0, ppi_score_global),
      ppi_score_local = ifelse(is.na(ppi_score_local), 0, ppi_score_local),
      # Log1p transform to handle skewness, as discussed
      log_ppi_score_global = log1p(ppi_score_global),
      log_ppi_score_local = log1p(ppi_score_local)
    )
  
  # --- Model Building, Evaluation, and Protein Scoring ---
  all_aucs <- list()
  model_predictions <- data.frame(protein_id = plot_data_full$Gene)
  all_proteins_scored_list <- list()
  auc_summary_df <- NULL
  roc_objects <- list()
  
  # Ensure there are at least two groups to compare for modeling
  if (n_distinct(plot_data_full$plot_group) < 2 || any(table(plot_data_full$plot_group) < 5)) {
    warning(paste("Skipping modeling for", df_name, "-", group, "due to insufficient data in one or both classes."))
    return(NULL)
  }
  
  tryCatch({
    # Define all available models
    all_models_available <- list(
      "1D_Abundance" = list(type="1D", predictor="mean_group_medianNorm"),
      "1D_logFC" = list(type="1D", predictor=logFC_col),
      "1D_PPI_Global" = list(type="1D", predictor="log_ppi_score_global"),
      "1D_PPI_Local" = list(type="1D", predictor="log_ppi_score_local"),
      "2D_Abundance_logFC" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col)),
      "2D_Abundance_PPI_Global" = list(type="GLM", formula="plot_group ~ mean_group_medianNorm + log_ppi_score_global"),
      "2D_Abundance_PPI_Local" = list(type="GLM", formula="plot_group ~ mean_group_medianNorm + log_ppi_score_local"),
      "2D_logFC_PPI_Global" = list(type="GLM", formula=paste("plot_group ~", logFC_col, "+ log_ppi_score_global")),
      "2D_logFC_PPI_Local" = list(type="GLM", formula=paste("plot_group ~", logFC_col, "+ log_ppi_score_local")),
      "3D_with_Global_PPI" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col, "+ log_ppi_score_global")),
      "3D_with_Local_PPI" = list(type="GLM", formula=paste("plot_group ~ mean_group_medianNorm +", logFC_col, "+ log_ppi_score_local"))
    )
    
    # Filter models based on user selection
    models_to_run <- all_models_available[intersect(names(all_models_available), selected_models)]
    
    # Use the pre-initialized roc_objects list

    for(model_name in names(models_to_run)) {
      model_info <- models_to_run[[model_name]]
      predictor_values <- NULL
      
      # Train model and get predictor values
      if(model_info$type == "1D") {
        model <- glm(as.formula(paste("plot_group ~", model_info$predictor)), data = plot_data_full, family = "binomial")
        predictor_values <- plot_data_full[[model_info$predictor]]
        model_predictions[[paste0("prob_", model_name)]] <- predict(model, type = "response")
      } else { # GLM for 2D/3D
        model <- glm(as.formula(model_info$formula), data = plot_data_full, family = "binomial")
        predictor_values <- predict(model, type = "response")
        model_predictions[[paste0("prob_", model_name)]] <- predictor_values
      }
      
      # Perform ROC analysis and store the roc object
      roc_objects[[model_name]] <- roc(response = plot_data_full$plot_group, predictor = predictor_values, 
                                       levels = c("Other", "TP"), quiet = TRUE, direction = "<")
      
      # Get best threshold using Youden index
      best_coords <- coords(roc_objects[[model_name]], "best", best.method = "youden", ret = "threshold")
      threshold_value <- best_coords$threshold
      
      # Score all proteins based on the model and create the detailed output table
      scored_proteins_df <- plot_data_full %>%
        mutate(
          model_score = if(model_info$type == "1D") .data[[model_info$predictor]] else predict(model, type = "response"),
          Model_best_threshold = threshold_value,
          PassThreshold = if_else(model_score > threshold_value, "Pass", "Fail"),
          isCandidate = if_else(plot_group == "Other" & PassThreshold == "Pass", "TP_candidate", NA_character_)
        ) %>%
        select(
          protein_id = Gene, 
          group = plot_group,
          contains("Localization"),
          abundance = mean_group_medianNorm, 
          logFC = all_of(logFC_col),
          ppi_score_global, ppi_score_local,
          log_ppi_score_global, log_ppi_score_local,
          model_score, 
          Model_best_threshold,
          PassThreshold, 
          isCandidate
        )
      all_proteins_scored_list[[model_name]] <- scored_proteins_df
    }

    # --- Add Traditional and Control Methods to ROC objects list ---
    # 只有在用户选择了传统方法时才计算
    if ("logFC>0.5 & FDR<0.05" %in% selected_models) {
      # Handle NA values in logFC and FDR columns gracefully
      traditional_predictor <- ifelse(
        !is.na(plot_data_full[[logFC_col]]) & !is.na(plot_data_full[[fdr_col]]) & 
        plot_data_full[[logFC_col]] > 0.5 & plot_data_full[[fdr_col]] < 0.05, 
        1, 0
      )
      
      # Only create ROC if there's variation in the predictor
      if (length(unique(traditional_predictor)) > 1) {
        roc_objects[["logFC>0.5 & FDR<0.05"]] <- roc(response = plot_data_full$plot_group, predictor = traditional_predictor, levels = c("Other", "TP"), quiet = TRUE, direction = "<")
      } else {
        warning(paste("Traditional method has no variation for", df_name, "-", group, ", skipping ROC creation"))
      }
    }
    
    set.seed(123)
    random_predictor <- runif(nrow(plot_data_full))
    roc_objects[["Control"]] <- roc(response = plot_data_full$plot_group, predictor = random_predictor, levels = c("Other", "TP"), quiet = TRUE, direction = "<")

    # --- Calculate AUC, CIs, and P-values vs. Traditional Method ---
    # Get AUC scores and CIs into a single data frame
    auc_summary_df <- purrr::map_df(roc_objects, ~ {
      ci_vals <- pROC::ci.auc(.x, method = "delong")
      tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
    }, .id = "ModelName")

    # Get p-values vs. the traditional model
    traditional_roc <- roc_objects[["logFC>0.5 & FDR<0.05"]]
    p_values_vs_trad <- purrr::map_dbl(roc_objects, ~ {
      # Don't compare the traditional model to itself or if traditional_roc is NULL
      if (identical(.x, traditional_roc) || is.null(traditional_roc)) return(NA)
      tryCatch({
        test_result <- roc.test(.x, traditional_roc, method = "delong")
        test_result$p.value
      }, error = function(e) { 
        warning(paste("ROC test failed:", e$message))
        NA 
      })
    })
    
    # Add p-values to the summary dataframe
    auc_summary_df <- auc_summary_df %>%
      mutate(p_value_vs_traditional = p_values_vs_trad)
    
    # --- Calculate P-values vs. 1D_logFC Model (for Figure 1b) ---
    logFC_roc <- roc_objects[["1D_logFC"]]
    p_values_vs_logFC <- purrr::map_dbl(roc_objects, ~ {
      # Don't compare 1D_logFC to itself or if logFC_roc is NULL
      if (identical(.x, logFC_roc) || is.null(logFC_roc)) return(NA)
      tryCatch({
        test_result <- roc.test(.x, logFC_roc, method = "delong")
        test_result$p.value
      }, error = function(e) { 
        warning(paste("ROC test vs 1D_logFC failed:", e$message))
        NA 
      })
    })
    
    # Add p-values vs 1D_logFC to the summary dataframe
    auc_summary_df <- auc_summary_df %>%
      mutate(p_value_vs_1D_logFC = p_values_vs_logFC)

    # --- New Analysis: Fixed Sensitivity/Specificity Comparison ---
    # 控制变量比较：固定Sensitivity比较Specificity，或固定Specificity比较Sensitivity
    controlled_comparison_results <- NULL
    
    if (!is.null(traditional_roc)) {
      # 获取传统方法的性能指标
      traditional_coords <- coords(traditional_roc, x = 1, input = "threshold", 
                                   ret = c("sensitivity", "specificity"))
      traditional_sens <- traditional_coords$sensitivity
      traditional_spec <- traditional_coords$specificity
      
      # 对每个连续模型进行控制变量比较
      controlled_comparison_list <- list()
      
      for (model_name in setdiff(names(roc_objects), c("logFC>0.5 & FDR<0.05", "Control"))) {
        model_roc <- roc_objects[[model_name]]
        
        tryCatch({
          # 1. 固定Sensitivity，比较Specificity
          # 在ROC曲线上找到最接近传统方法sensitivity的点
          model_coords_all <- coords(model_roc, x = "all", ret = c("threshold", "sensitivity", "specificity"))
          
          # 找到最接近传统方法sensitivity的点
          sens_diff <- abs(model_coords_all$sensitivity - traditional_sens)
          closest_sens_idx <- which.min(sens_diff)
          model_spec_at_fixed_sens <- model_coords_all$specificity[closest_sens_idx]
          actual_sens_used <- model_coords_all$sensitivity[closest_sens_idx]
          
          # 2. 固定Specificity，比较Sensitivity
          # 找到最接近传统方法specificity的点
          spec_diff <- abs(model_coords_all$specificity - traditional_spec)
          closest_spec_idx <- which.min(spec_diff)
          model_sens_at_fixed_spec <- model_coords_all$sensitivity[closest_spec_idx]
          actual_spec_used <- model_coords_all$specificity[closest_spec_idx]
          
          # 计算增益
          specificity_gain <- model_spec_at_fixed_sens - traditional_spec
          sensitivity_gain <- model_sens_at_fixed_spec - traditional_sens
          
          # 使用bootstrap方法计算置信区间和p值
          # Bootstrap for specificity comparison (fixed sensitivity)
          set.seed(123)
          n_boot <- if(exists("n_bootstrap")) n_bootstrap else 1000  # 使用全局配置的bootstrap次数，默认1000
          boot_spec_diff <- numeric(n_boot)
          
          # 获取原始数据
          response_vec <- as.numeric(model_roc$response) - 1  # Convert to 0/1
          predictor_traditional <- as.numeric(traditional_roc$predictor)
          predictor_model <- as.numeric(model_roc$predictor)
          
          for (b in 1:n_boot) {
            # Bootstrap sample
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_trad <- predictor_traditional[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            # Calculate traditional method specificity
            boot_trad_pred_binary <- ifelse(boot_pred_trad >= 1, 1, 0)
            boot_trad_spec <- sum(boot_response == 0 & boot_trad_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            # For model, find threshold that gives similar sensitivity to traditional
            boot_trad_sens <- sum(boot_response == 1 & boot_trad_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            # Find model threshold
            thresholds <- sort(unique(boot_pred_model))
            best_thresh <- thresholds[1]
            min_sens_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
              sens_diff_boot <- abs(boot_model_sens - boot_trad_sens)
              if (sens_diff_boot < min_sens_diff) {
                min_sens_diff <- sens_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            boot_spec_diff[b] <- boot_model_spec - boot_trad_spec
          }
          
          # Calculate CI and p-value for specificity gain
          spec_gain_ci_lower <- quantile(boot_spec_diff, 0.025, na.rm = TRUE)
          spec_gain_ci_upper <- quantile(boot_spec_diff, 0.975, na.rm = TRUE)
          # 使用保守的P值估计：(x+1)/(n+1)，避免P=0的歧义
          # 这是Agresti-Coull调整的简化版本
          n_valid_boot <- sum(!is.na(boot_spec_diff))
          n_unfavorable <- sum(boot_spec_diff <= 0, na.rm = TRUE)
          spec_gain_p_value <- (n_unfavorable + 1) / (n_valid_boot + 1)
          
          # Bootstrap for sensitivity comparison (fixed specificity)
          boot_sens_diff <- numeric(n_boot)
          
          for (b in 1:n_boot) {
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_trad <- predictor_traditional[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            # Calculate traditional method sensitivity
            boot_trad_pred_binary <- ifelse(boot_pred_trad >= 1, 1, 0)
            boot_trad_sens <- sum(boot_response == 1 & boot_trad_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            # For model, find threshold that gives similar specificity to traditional
            boot_trad_spec <- sum(boot_response == 0 & boot_trad_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            # Find model threshold
            thresholds <- sort(unique(boot_pred_model), decreasing = TRUE)
            best_thresh <- thresholds[1]
            min_spec_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
              spec_diff_boot <- abs(boot_model_spec - boot_trad_spec)
              if (spec_diff_boot < min_spec_diff) {
                min_spec_diff <- spec_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            boot_sens_diff[b] <- boot_model_sens - boot_trad_sens
          }
          
          # Calculate CI and p-value for sensitivity gain
          sens_gain_ci_lower <- quantile(boot_sens_diff, 0.025, na.rm = TRUE)
          sens_gain_ci_upper <- quantile(boot_sens_diff, 0.975, na.rm = TRUE)
          # 使用保守的P值估计：(x+1)/(n+1)，避免P=0的歧义
          n_valid_boot_sens <- sum(!is.na(boot_sens_diff))
          n_unfavorable_sens <- sum(boot_sens_diff <= 0, na.rm = TRUE)
          sens_gain_p_value <- (n_unfavorable_sens + 1) / (n_valid_boot_sens + 1)
          
          # 存储结果
          controlled_comparison_list[[model_name]] <- data.frame(
            ModelName = model_name,
            # Traditional method metrics
            Traditional_Sensitivity = traditional_sens,
            Traditional_Specificity = traditional_spec,
            # Fixed Sensitivity comparison
            Fixed_Sensitivity = actual_sens_used,
            Model_Specificity_at_Fixed_Sens = model_spec_at_fixed_sens,
            Specificity_Gain = specificity_gain,
            Specificity_Gain_CI_Lower = spec_gain_ci_lower,
            Specificity_Gain_CI_Upper = spec_gain_ci_upper,
            Specificity_Gain_P_Value = spec_gain_p_value,
            # Fixed Specificity comparison
            Fixed_Specificity = actual_spec_used,
            Model_Sensitivity_at_Fixed_Spec = model_sens_at_fixed_spec,
            Sensitivity_Gain = sensitivity_gain,
            Sensitivity_Gain_CI_Lower = sens_gain_ci_lower,
            Sensitivity_Gain_CI_Upper = sens_gain_ci_upper,
            Sensitivity_Gain_P_Value = sens_gain_p_value,
            stringsAsFactors = FALSE
          )
          
        }, error = function(e) {
          warning(paste("Controlled comparison failed for", model_name, ":", e$message))
        })
      }
      
      # 合并所有模型的控制变量比较结果
      if (length(controlled_comparison_list) > 0) {
        controlled_comparison_results <- bind_rows(controlled_comparison_list)
      }
    }
    
    # --- New Analysis: Comparison vs 1D_logFC Model (for Figure 7a-7c) ---
    # 与1D_logFC模型的控制变量比较
    controlled_comparison_vs_logFC <- NULL
    
    logFC_roc <- roc_objects[["1D_logFC"]]
    
    if (!is.null(logFC_roc)) {
      # 获取1D_logFC模型在Youden最佳阈值下的性能指标
      logFC_best_coords <- coords(logFC_roc, "best", best.method = "youden", 
                                   ret = c("threshold", "sensitivity", "specificity"))
      logFC_sens <- logFC_best_coords$sensitivity
      logFC_spec <- logFC_best_coords$specificity
      
      # 对每个其他连续模型进行控制变量比较
      controlled_comparison_vs_logFC_list <- list()
      
      for (model_name in setdiff(names(roc_objects), c("1D_logFC", "logFC>0.5 & FDR<0.05", "Control"))) {
        model_roc <- roc_objects[[model_name]]
        
        tryCatch({
          # 1. 固定Sensitivity，比较Specificity
          model_coords_all <- coords(model_roc, x = "all", ret = c("threshold", "sensitivity", "specificity"))
          
          sens_diff <- abs(model_coords_all$sensitivity - logFC_sens)
          closest_sens_idx <- which.min(sens_diff)
          model_spec_at_fixed_sens <- model_coords_all$specificity[closest_sens_idx]
          actual_sens_used <- model_coords_all$sensitivity[closest_sens_idx]
          
          # 2. 固定Specificity，比较Sensitivity
          spec_diff <- abs(model_coords_all$specificity - logFC_spec)
          closest_spec_idx <- which.min(spec_diff)
          model_sens_at_fixed_spec <- model_coords_all$sensitivity[closest_spec_idx]
          actual_spec_used <- model_coords_all$specificity[closest_spec_idx]
          
          # 计算增益
          specificity_gain <- model_spec_at_fixed_sens - logFC_spec
          sensitivity_gain <- model_sens_at_fixed_spec - logFC_sens
          
          # 使用bootstrap方法计算置信区间和p值
          set.seed(123)
          n_boot <- if(exists("n_bootstrap")) n_bootstrap else 1000
          boot_spec_diff <- numeric(n_boot)
          
          # 获取原始数据
          response_vec <- as.numeric(model_roc$response) - 1
          predictor_logFC <- as.numeric(logFC_roc$predictor)
          predictor_model <- as.numeric(model_roc$predictor)
          
          # Bootstrap for specificity comparison (fixed sensitivity)
          for (b in 1:n_boot) {
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_logFC <- predictor_logFC[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            # Calculate logFC model specificity at its best threshold
            logFC_best_thresh <- logFC_best_coords$threshold
            boot_logFC_pred_binary <- ifelse(boot_pred_logFC >= logFC_best_thresh, 1, 0)
            boot_logFC_spec <- sum(boot_response == 0 & boot_logFC_pred_binary == 0) / max(1, sum(boot_response == 0))
            boot_logFC_sens <- sum(boot_response == 1 & boot_logFC_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            # Find model threshold that gives similar sensitivity
            thresholds <- sort(unique(boot_pred_model))
            best_thresh <- thresholds[1]
            min_sens_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
              sens_diff_boot <- abs(boot_model_sens - boot_logFC_sens)
              if (sens_diff_boot < min_sens_diff) {
                min_sens_diff <- sens_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            boot_spec_diff[b] <- boot_model_spec - boot_logFC_spec
          }
          
          spec_gain_ci_lower <- quantile(boot_spec_diff, 0.025, na.rm = TRUE)
          spec_gain_ci_upper <- quantile(boot_spec_diff, 0.975, na.rm = TRUE)
          # 使用保守的P值估计：(x+1)/(n+1)
          n_valid_boot <- sum(!is.na(boot_spec_diff))
          n_unfavorable <- sum(boot_spec_diff <= 0, na.rm = TRUE)
          spec_gain_p_value <- (n_unfavorable + 1) / (n_valid_boot + 1)
          
          # Bootstrap for sensitivity comparison (fixed specificity)
          boot_sens_diff <- numeric(n_boot)
          
          for (b in 1:n_boot) {
            boot_idx <- sample(length(response_vec), replace = TRUE)
            boot_response <- response_vec[boot_idx]
            boot_pred_logFC <- predictor_logFC[boot_idx]
            boot_pred_model <- predictor_model[boot_idx]
            
            boot_logFC_pred_binary <- ifelse(boot_pred_logFC >= logFC_best_thresh, 1, 0)
            boot_logFC_sens <- sum(boot_response == 1 & boot_logFC_pred_binary == 1) / max(1, sum(boot_response == 1))
            boot_logFC_spec <- sum(boot_response == 0 & boot_logFC_pred_binary == 0) / max(1, sum(boot_response == 0))
            
            # Find model threshold that gives similar specificity
            thresholds <- sort(unique(boot_pred_model), decreasing = TRUE)
            best_thresh <- thresholds[1]
            min_spec_diff <- Inf
            
            for (thresh in thresholds) {
              boot_model_pred_binary <- ifelse(boot_pred_model >= thresh, 1, 0)
              boot_model_spec <- sum(boot_response == 0 & boot_model_pred_binary == 0) / max(1, sum(boot_response == 0))
              spec_diff_boot <- abs(boot_model_spec - boot_logFC_spec)
              if (spec_diff_boot < min_spec_diff) {
                min_spec_diff <- spec_diff_boot
                best_thresh <- thresh
              }
            }
            
            boot_model_pred_binary <- ifelse(boot_pred_model >= best_thresh, 1, 0)
            boot_model_sens <- sum(boot_response == 1 & boot_model_pred_binary == 1) / max(1, sum(boot_response == 1))
            
            boot_sens_diff[b] <- boot_model_sens - boot_logFC_sens
          }
          
          sens_gain_ci_lower <- quantile(boot_sens_diff, 0.025, na.rm = TRUE)
          sens_gain_ci_upper <- quantile(boot_sens_diff, 0.975, na.rm = TRUE)
          # 使用保守的P值估计：(x+1)/(n+1)
          n_valid_boot_sens <- sum(!is.na(boot_sens_diff))
          n_unfavorable_sens <- sum(boot_sens_diff <= 0, na.rm = TRUE)
          sens_gain_p_value <- (n_unfavorable_sens + 1) / (n_valid_boot_sens + 1)
          
          # 存储结果
          controlled_comparison_vs_logFC_list[[model_name]] <- data.frame(
            ModelName = model_name,
            # 1D_logFC model metrics
            Baseline_Sensitivity = logFC_sens,
            Baseline_Specificity = logFC_spec,
            # Fixed Sensitivity comparison
            Fixed_Sensitivity = actual_sens_used,
            Model_Specificity_at_Fixed_Sens = model_spec_at_fixed_sens,
            Specificity_Gain = specificity_gain,
            Specificity_Gain_CI_Lower = spec_gain_ci_lower,
            Specificity_Gain_CI_Upper = spec_gain_ci_upper,
            Specificity_Gain_P_Value = spec_gain_p_value,
            # Fixed Specificity comparison
            Fixed_Specificity = actual_spec_used,
            Model_Sensitivity_at_Fixed_Spec = model_sens_at_fixed_spec,
            Sensitivity_Gain = sensitivity_gain,
            Sensitivity_Gain_CI_Lower = sens_gain_ci_lower,
            Sensitivity_Gain_CI_Upper = sens_gain_ci_upper,
            Sensitivity_Gain_P_Value = sens_gain_p_value,
            stringsAsFactors = FALSE
          )
          
        }, error = function(e) {
          warning(paste("Controlled comparison vs 1D_logFC failed for", model_name, ":", e$message))
        })
      }
      
      # 合并所有模型相对于1D_logFC的比较结果
      if (length(controlled_comparison_vs_logFC_list) > 0) {
        controlled_comparison_vs_logFC <- bind_rows(controlled_comparison_vs_logFC_list)
      }
    }

    # Add traditional method to the scored list (for consistency in Phase 3) only if ROC was created
    # Note: Threshold is binary (1), so candidates are simply those that pass
    if ("logFC>0.5 & FDR<0.05" %in% names(roc_objects) && "logFC>0.5 & FDR<0.05" %in% selected_models) {
      all_proteins_scored_list[["logFC>0.5 & FDR<0.05"]] <- plot_data_full %>%
          mutate(
            model_score = traditional_predictor,
            Model_best_threshold = 1,
            PassThreshold = if_else(model_score >= 1, "Pass", "Fail"),
            isCandidate = if_else(plot_group == "Other" & PassThreshold == "Pass", "TP_candidate", NA_character_)
          ) %>%
          select(
            protein_id = Gene, group = plot_group, contains("Localization"),
            abundance = mean_group_medianNorm, logFC = all_of(logFC_col),
            ppi_score_global, ppi_score_local, log_ppi_score_global, log_ppi_score_local,
            model_score, Model_best_threshold, PassThreshold, isCandidate
          )
    }

  }, error = function(e) { 
    warning(paste("Model evaluation failed for", df_name, "-", group, ":", e$message))
    # Return empty objects to prevent downstream errors
    auc_summary_df <<- data.frame()
    roc_objects <<- list()
    all_proteins_scored_list <<- list()
  })
  
  # --- Package all results for return ---
  return(list(
    df_name = df_name, 
    group_name = group,
    auc_summary_df = if(is.null(auc_summary_df)) data.frame() else auc_summary_df, # AUC, CIs, and p-value vs traditional
    roc_objects = roc_objects, # Return full ROC objects for cross-group comparison
    node_attributes = plot_data_full,
    model_predictions = model_predictions, 
    candidates_full_list = all_proteins_scored_list,
    interactions = interactions_df,
    controlled_comparison = if(is.null(controlled_comparison_results)) data.frame() else controlled_comparison_results, # Controlled variable comparison vs traditional
    controlled_comparison_vs_logFC = if(is.null(controlled_comparison_vs_logFC)) data.frame() else controlled_comparison_vs_logFC # Controlled variable comparison vs 1D_logFC
  ))
}


# -------------------------------------------------------------------
# 5. Execute Parallel Analysis
# -------------------------------------------------------------------
cat("开始执行所有组合的并行分析...\n")

# Enable progress bar tracking
handlers(global = TRUE)
handlers("progress")

# Create a grid of tasks based on selected groups
tasks <- expand_grid(df_name = names(ForStep19), group = selected_groups)

with_progress({
  # Setup the progressor
  p <- progressr::progressor(steps = nrow(tasks))
  
  # Use future_pmap to iterate over all combinations in parallel
  all_results <- future_pmap(list(
    df_name = tasks$df_name,
    group = tasks$group
  ), function(df_name, group) {
    # Run the main worker function
    result <- tryCatch({
      process_group(df_name, group, ppi_data)
    }, error = function(e) {
      warning(paste("Task failed for", df_name, "-", group, ":", e$message))
      return(NULL)
    })
    
    # Update progress only if progressor is still active
    tryCatch({
      if (!is.null(p)) {
        p(message = paste("已完成", df_name, "-", group))
      }
    }, error = function(e) {
      # Silently ignore progress update errors - this is expected when progressor completes
    })
    
    return(result)
  }, .options = furrr_options(seed = TRUE)) # Use seed for reproducibility
  
  # 为all_results添加有意义的名称
  names(all_results) <- paste(tasks$df_name, tasks$group, sep = "_")
})
 

cat("\n所有并行计算任务已完成。\n")


# -------------------------------------------------------------------
# 6. Consolidate Results and Export Files
# -------------------------------------------------------------------
cat("正在汇总结果并导出文件...\n")

# --- Consolidate AUC Summaries and Apply Multiple-Test Correction ---
# Extract the auc_summary_df from each task result
all_auc_summaries <- purrr::map_df(all_results, ~ {
  if (is.null(.x) || !is.list(.x) || is.null(.x$auc_summary_df)) {
    return(NULL)
  }
  .x$auc_summary_df
}, .id = "Analysis_Group")

# Apply BH correction to p-values within each analysis group
final_auc_summary <- all_auc_summaries %>% 
  group_by(Analysis_Group) %>% 
  mutate(
    q_value_vs_traditional = p.adjust(p_value_vs_traditional, method = "BH")
  ) %>% 
  ungroup() %>%
  select(Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, p_value_vs_traditional, q_value_vs_traditional, everything())

# --- Export Main AUC Summary ---
writexl::write_xlsx(list(AUC_Summary = final_auc_summary), path = file.path(output_dir, "AUC_summary.xlsx"))

# --- Consolidate and Export Candidate Lists and Cytoscape Files ---
final_candidates_list <- list()

# Get the task list again for matching failed tasks
tasks_check <- expand_grid(df_name = names(ForStep19), group = names(group_info))

for(i in seq_along(all_results)) {
  result <- all_results[[i]]
  
  # Check for failed tasks
  if (is.null(result) || !is.list(result) || is.null(result$candidates_full_list)) {
    failed_task <- tasks_check[i, ]
    warning(paste("任务", failed_task$df_name, "-", failed_task$group, "的结果不完整，部分导出文件可能缺失。"))
    next
  }
  
  df_name <- result$df_name
  group <- result$group_name

  # Consolidate scored protein list
  candidates_by_model <- result$candidates_full_list
  if (length(candidates_by_model) > 0) {
    for (model_name in names(candidates_by_model)) {
      candidates_df <- candidates_by_model[[model_name]]
      if (nrow(candidates_df) > 0) {
        candidates_df$DataAnalysisMethod <- df_name
        candidates_df$ExperimentGroup <- group
        candidates_df$Model <- model_name
        final_candidates_list[[paste(df_name, group, model_name)]] <- candidates_df
      }
    }
  }
  
  # Export files for Cytoscape (如果启用)
  if (enable_cytoscape_export) {
    cytoscape_path_prefix <- file.path(output_dir, "Cytoscape_Files", paste0(df_name, "_", group))
    dir.create(dirname(cytoscape_path_prefix), recursive = TRUE, showWarnings = FALSE)
    logFC_col <- group_info[[group]]$logFC
    
    # Node attributes file
    node_attributes_for_export <- result$node_attributes %>% 
      left_join(result$model_predictions, by = c("Gene" = "protein_id")) %>% 
      select(
        protein_id = Gene, 
        group = plot_group,
        contains("Localization"),
        logFC = all_of(logFC_col),
        everything()
      )
    write.csv(node_attributes_for_export, 
              paste0(cytoscape_path_prefix, "_node_attributes.csv"), 
              row.names = FALSE,quote=FALSE)
    
    # Edge attributes file
    write.csv(result$interactions, paste0(cytoscape_path_prefix, "_edge_attributes.csv"), 
              row.names = FALSE,quote=FALSE)
  }
}
#install.packages("writexl")
library("writexl")

# Combine and export candidate list
final_candidates_df <- bind_rows(final_candidates_list)
writexl::write_xlsx(list(Candidates = final_candidates_df), path = file.path(output_dir, "detailed_candidates_by_model.xlsx"))

cat("\n所有结果文件已成功导出到目录:", output_dir, "\n")

# --- Consolidate and Export Controlled Comparison Results ---
cat("\n正在汇总控制变量比较结果（固定Sensitivity/Specificity）...\n")

all_controlled_comparisons <- purrr::map_df(all_results, ~ {
  if (is.null(.x) || !is.list(.x) || is.null(.x$controlled_comparison)) {
    return(NULL)
  }
  .x$controlled_comparison
}, .id = "Analysis_Group")

if (nrow(all_controlled_comparisons) > 0) {
  # 添加Df_Name列
  group_to_df_map_controlled <- purrr::map_df(all_results, ~tibble::tibble(Df_Name = .x$df_name), .id = "Analysis_Group") %>% dplyr::distinct()
  
  final_controlled_comparison <- all_controlled_comparisons %>%
    dplyr::left_join(group_to_df_map_controlled, by = "Analysis_Group") %>%
    select(Df_Name, Analysis_Group, everything()) %>%
    # 应用BH多重检验校正
    group_by(Analysis_Group) %>%
    mutate(
      Specificity_Gain_Q_Value = p.adjust(Specificity_Gain_P_Value, method = "BH"),
      Sensitivity_Gain_Q_Value = p.adjust(Sensitivity_Gain_P_Value, method = "BH"),
      # 添加P值理论下限说明（基于bootstrap次数）
      P_Value_Lower_Limit = 1 / (n_bootstrap + 1),
      # 添加显著性标记
      Specificity_Significant = if_else(Specificity_Gain_Q_Value < 0.05, "Yes", "No"),
      Sensitivity_Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "Yes", "No")
    ) %>%
    ungroup()
  
  # 导出到Excel
  writexl::write_xlsx(
    list(Controlled_Comparison = final_controlled_comparison), 
    path = file.path(output_dir, "Controlled_Variable_Comparison.xlsx")
  )
  
  cat("控制变量比较结果已导出到 Controlled_Variable_Comparison.xlsx\n")
  cat("  - 固定Sensitivity比较Specificity：", sum(!is.na(final_controlled_comparison$Specificity_Gain)), "个比较\n")
  cat("  - 固定Specificity比较Sensitivity：", sum(!is.na(final_controlled_comparison$Sensitivity_Gain)), "个比较\n")
} else {
  cat("警告: 未能从分析结果中提取任何控制变量比较信息。\n")
}

# --- Consolidate and Export Controlled Comparison vs 1D_logFC Results ---
cat("\n正在汇总与1D_logFC模型的控制变量比较结果...\n")

all_controlled_comparisons_vs_logFC <- purrr::map_df(all_results, ~ {
  if (is.null(.x) || !is.list(.x) || is.null(.x$controlled_comparison_vs_logFC)) {
    return(NULL)
  }
  .x$controlled_comparison_vs_logFC
}, .id = "Analysis_Group")

if (nrow(all_controlled_comparisons_vs_logFC) > 0) {
  # 添加Df_Name列
  group_to_df_map_logFC <- purrr::map_df(all_results, ~tibble::tibble(Df_Name = .x$df_name), .id = "Analysis_Group") %>% dplyr::distinct()
  
  final_controlled_comparison_vs_logFC <- all_controlled_comparisons_vs_logFC %>%
    dplyr::left_join(group_to_df_map_logFC, by = "Analysis_Group") %>%
    select(Df_Name, Analysis_Group, everything()) %>%
    # 应用BH多重检验校正
    group_by(Analysis_Group) %>%
    mutate(
      Specificity_Gain_Q_Value = p.adjust(Specificity_Gain_P_Value, method = "BH"),
      Sensitivity_Gain_Q_Value = p.adjust(Sensitivity_Gain_P_Value, method = "BH"),
      # 添加P值理论下限说明（基于bootstrap次数）
      P_Value_Lower_Limit = 1 / (n_bootstrap + 1),
      # 添加显著性标记
      Specificity_Significant = if_else(Specificity_Gain_Q_Value < 0.05, "Yes", "No"),
      Sensitivity_Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "Yes", "No")
    ) %>%
    ungroup()
  
  # 导出到Excel
  writexl::write_xlsx(
    list(Controlled_Comparison_vs_1D_logFC = final_controlled_comparison_vs_logFC), 
    path = file.path(output_dir, "Controlled_Variable_Comparison_vs_1D_logFC.xlsx")
  )
  
  cat("与1D_logFC模型的控制变量比较结果已导出到 Controlled_Variable_Comparison_vs_1D_logFC.xlsx\n")
  cat("  - 固定Sensitivity比较Specificity：", sum(!is.na(final_controlled_comparison_vs_logFC$Specificity_Gain)), "个比较\n")
  cat("  - 固定Specificity比较Sensitivity：", sum(!is.na(final_controlled_comparison_vs_logFC$Sensitivity_Gain)), "个比较\n")
} else {
  cat("警告: 未能从分析结果中提取任何与1D_logFC的控制变量比较信息。\n")
}


# ===================================================================
#
###    Phase 1.1: Data Export and Robustness Checks #####
#
# ===================================================================

cat("\n\n--- 开始执行 Phase 1.1: 数据导出与稳健性检验 ---\n")

# --- 0. Consolidate, Correct, and Export All Model AUCs ---
cat("\n--- 正在汇总、校正并导出所有模型的AUC结果... ---\n")

# Extract the auc_summary dataframe from each result, adding the group name
all_auc_summaries_raw <- purrr::map_dfr(all_results, ~ .x$auc_summary_df, .id = "Analysis_Group")

# Create a mapping from Analysis_Group (list names) to Df_Name (field in list)
# This is crucial for linking results back to their original dataset (e.g., APEX, BioID)
group_to_df_map <- purrr::map_df(all_results, ~tibble::tibble(Df_Name = .x$df_name), .id = "Analysis_Group") %>% dplyr::distinct()

# Check if the summary data frame is not empty
if (nrow(all_auc_summaries_raw) > 0) {
  
  # 检测是否为单次比较情况
  total_comparisons <- length(selected_groups) * length(setdiff(selected_models, "Control"))
  single_comparison_flag <- (total_comparisons == 1)
  
  if (single_comparison_flag) {
    cat("✓ 检测到单次比较情况 - 只有1个组别使用1个模型，p值无需多重检验校正\n")
  } else {
    cat("✓ 检测到多次比较情况 - 共", total_comparisons, "次比较，将应用BH多重检验校正\n")
  }
  
  # Apply BH correction to p-values within each analysis group to get q-values
  # This is the definitive, final summary table for AUCs.
  final_auc_summary <- all_auc_summaries_raw %>% 
    group_by(Analysis_Group) %>% 
    mutate(
      q_value_vs_traditional = if (single_comparison_flag) {
        # 单次比较时，q值等于p值
        p_value_vs_traditional
      } else {
        # 多次比较时，应用BH校正
        p.adjust(p_value_vs_traditional, method = "BH")
      }
    ) %>%
    ungroup() %>% 
    # Add the Df_Name column by joining the map
    dplyr::left_join(group_to_df_map, by = "Analysis_Group") %>%
    # Reorder columns for clarity, bringing Df_Name to the front
    select(Df_Name, Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, p_value_vs_traditional, q_value_vs_traditional, everything()) %>%
    # 添加标记列表示是否为单次比较
    mutate(single_comparison = single_comparison_flag)

  # Export the final, corrected summary to its own Excel file
  auc_summary_path <- file.path(output_dir, "AUC_summary.xlsx")
  writexl::write_xlsx(list("AUC_Summary" = final_auc_summary), path = auc_summary_path)
  cat("所有模型的AUC汇总结果（已进行多重检验校正）已导出到 AUC_summary.xlsx\n")
  
  # --- Export AUC Summary with comparison to 1D_logFC (for Figure 1b) ---
  cat("\n--- 正在创建与1D_logFC比较的AUC汇总结果... ---\n")
  
  # 检查是否有p_value_vs_1D_logFC列
  if ("p_value_vs_1D_logFC" %in% names(final_auc_summary)) {
    
    # 创建针对1D_logFC比较的汇总表
    final_auc_summary_b <- final_auc_summary %>%
      # 只保留与1D_logFC比较相关的列
      select(Df_Name, Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, p_value_vs_1D_logFC) %>%
      # 过滤掉1D_logFC自己和Control
      filter(!ModelName %in% c("1D_logFC", "Control")) %>%
      # 应用BH校正到vs 1D_logFC的p值
      group_by(Analysis_Group) %>%
      mutate(
        q_value_vs_1D_logFC = p.adjust(p_value_vs_1D_logFC, method = "BH"),
        # 添加显著性标记
        is_significant_vs_1D_logFC = if_else(
          !is.na(q_value_vs_1D_logFC) & q_value_vs_1D_logFC < 0.05,
          "Yes",
          "No"
        )
      ) %>%
      ungroup() %>%
      # 重新排列列顺序
      select(Df_Name, Analysis_Group, ModelName, AUC, AUC_Lower_CI, AUC_Upper_CI, 
             p_value_vs_1D_logFC, q_value_vs_1D_logFC, is_significant_vs_1D_logFC, everything())
    
    # 导出AUC_Summary_b
    auc_summary_b_path <- file.path(output_dir, "AUC_Summary_b.xlsx")
    writexl::write_xlsx(list("AUC_Summary_vs_1D_logFC" = final_auc_summary_b), path = auc_summary_b_path)
    cat("与1D_logFC模型比较的AUC汇总结果已导出到 AUC_Summary_b.xlsx\n")
    
  } else {
    cat("警告: 未找到 p_value_vs_1D_logFC 列，跳过 AUC_Summary_b.xlsx 的生成。\n")
  }
  
} else {
  cat("警告: 未能从分析结果中提取任何AUC摘要信息，跳过 AUC_summary.xlsx 的生成。\n")
}


# --- 1. Type 2 Robustness (Model Stability) Analysis ---
if (enable_robustness_analysis) {
  cat("\n执行第二类稳健性分析（模型稳定性）...\n")

  # Calculate Mean and StdDev of AUCs for each group across all models
  # This now uses the `final_auc_summary` table which contains q-values
  group_robustness_summary <- final_auc_summary %>%
    # Exclude the control and traditional methods from robustness calculation
    filter(!ModelName %in% c("Control", "logFC>0.5 & FDR<0.05")) %>%
    group_by(Df_Name, Analysis_Group) %>%
    summarise(
      Mean_AUC = mean(AUC, na.rm = TRUE),
      AUC_SD = sd(AUC, na.rm = TRUE),
      Num_Models_Tested = n(),
      .groups = 'drop'  # Avoid nested grouping issues
    ) %>%
    arrange(desc(Mean_AUC), AUC_SD) # Sort by best and most stable

  # Export the robustness summary
  writexl::write_xlsx(
    list(Group_Robustness_Summary = group_robustness_summary),
    path = file.path(output_dir, "Group_Robustness_Summary.xlsx")
  )

  cat("\n模型稳定性结果已导出到 Group_Robustness_Summary.xlsx\n")
} else {
  cat("\n稳健性分析已禁用，跳过相关计算。\n")
}


# --- New Section: Cross-Group (Inter-Method) Comparison ---
if (enable_group_comparison) {
  cat("\n开始执行跨邻近方法比较分析...\n")

# Define the baseline model for this comparison
baseline_model_for_comparison <- "2D_Abundance_logFC"

# Group results by the original dataframe name (e.g., "ForStep19_APEX", "ForStep19_BioID2")
# Ensure Phase 1.2 data copies exist when running this block alone
if (!exists("all_results_p1_2")) {
  all_results_p1_2 <- rlang::duplicate(all_results, shallow = FALSE)
}
if (!exists("group_info_p1_2") && exists("group_info")) {
  group_info_p1_2 <- rlang::duplicate(group_info, shallow = FALSE)
}

results_by_df <- all_results_p1_2 %>% 
  purrr::keep(~!is.null(.x)) %>% # Filter out failed runs
  split(., purrr::map_chr(., "df_name"))

# List to store the comparison results for each dataframe
all_group_comparisons <- list()

for (df_name in names(results_by_df)) {
  
  # Extract the ROC objects for the baseline model from each group
  roc_objects_for_comparison <- results_by_df[[df_name]] %>% 
    purrr::map(~.x$roc_objects[[baseline_model_for_comparison]]) %>% 
    # Name the list elements by their group name for clarity
    setNames(purrr::map_chr(results_by_df[[df_name]], "group_name")) %>% 
    # Remove any NULL entries if a model failed for a specific group
    purrr::compact()

  if (length(roc_objects_for_comparison) < 2) {
    warning(paste("Skipping group comparison for", df_name, "as it has fewer than 2 valid groups for the model", baseline_model_for_comparison))
    next
  }
  
  # --- Perform pairwise DeLong's test ---
  group_names <- names(roc_objects_for_comparison)
  p_value_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                           dimnames = list(group_names, group_names))
  
  group_pairs <- combn(group_names, 2, simplify = FALSE)
  
  p_values <- purrr::map_dbl(group_pairs, ~ {
    tryCatch({
      test_result <- roc.test(roc_objects_for_comparison[[.x[1]]], roc_objects_for_comparison[[.x[2]]], method = "delong")
      test_result$p.value
    }, error = function(e) { NA })
  })
  
  for (i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    p_value_matrix[pair[1], pair[2]] <- p_values[i]
    p_value_matrix[pair[2], pair[1]] <- p_values[i]
  }
  diag(p_value_matrix) <- 1.0

  # --- Consolidate results into a summary table ---
  # Get AUC and CIs for each group
  auc_ci_summary <- purrr::map_df(roc_objects_for_comparison, ~{
    ci_vals <- pROC::ci.auc(.x, method = "delong")
    tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
  }, .id = "GroupName")
  
  # Convert p-value matrix to a long format and apply BH correction
  p_values_long_temp <- as.data.frame(p_value_matrix) %>% 
    tibble::rownames_to_column("Group1") %>% 
    tidyr::pivot_longer(-Group1, names_to = "Group2", values_to = "p_value") %>% 
    filter(Group1 < Group2) # Avoid duplicates and self-comparisons
  
  p_values_long <- p_values_long_temp %>%
    mutate(
      q_value = if (nrow(p_values_long_temp) == 1) {
        # 单次组间比较时，q值等于p值
        p_value
      } else {
        # 多次组间比较时，应用BH校正
        p.adjust(p_value, method = "BH")
      },
      single_group_comparison = (nrow(p_values_long_temp) == 1)
    )

  # Store results
  all_group_comparisons[[df_name]] <- list(
    auc_ci_summary = auc_ci_summary,
    pairwise_tests = p_values_long
  )
}

# --- Export Cross-Group Comparison Results ---
if (length(all_group_comparisons) > 0) {
  # Create a single Excel file with sheets for each type of result
  export_list_group_comp <- list()
  
  # Add AUC/CI summaries for each dataset
  auc_sheets <- purrr::map(all_group_comparisons, "auc_ci_summary")
  names(auc_sheets) <- paste("C_AUC", names(auc_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, auc_sheets)
  
  # Add pairwise test results for each dataset
  pairwise_sheets <- purrr::map(all_group_comparisons, "pairwise_tests")
  names(pairwise_sheets) <- paste("C_PVal", names(pairwise_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, pairwise_sheets)
  
  writexl::write_xlsx(export_list_group_comp, path = file.path(output_dir, "Group_Comparison_Summary.xlsx"))
  cat("\n跨邻近方法比较结果已导出到 Group_Comparison_Summary.xlsx\n")
} else {
  cat("\n未能生成任何跨邻近方法比较结果。\n")
}

} else {
  cat("\n组间比较分析已禁用，跳过相关计算。\n")
}

# --- New Section: Type 1 Robustness (Sample Intersection Fairness) Analysis ---
if (enable_robustness_analysis) {
  cat("\n开始执行第一类稳健性分析（样本交集公平性）...\n")

# Define the baseline model again for clarity
baseline_model_for_intersection <- "2D_Abundance_logFC"

# List to store the intersection comparison results
all_intersection_comparisons <- list()

for (df_name in names(results_by_df)) {
  
  current_groups_results <- results_by_df[[df_name]]
  
  # --- Find intersecting proteins for the current df_name ---
  list_of_protein_ids <- current_groups_results %>% 
    purrr::map(~.x$node_attributes$Gene)
  
  intersecting_proteins <- Reduce(intersect, list_of_protein_ids)
  
  cat(paste("\n对于", df_name, "，找到", length(intersecting_proteins), "个所有组共有的蛋白用于交集分析。\n"))
  
  if (length(intersecting_proteins) < 20) { # Heuristic check
      warning(paste("Skipping intersection analysis for", df_name, "due to very small protein intersection (", length(intersecting_proteins), ")"))
      next
  }

  # --- Re-run model on intersected data for each group ---
  roc_objects_intersect <- list()
  
  for(result in current_groups_results) {
    group_name <- result$group_name
    
    # Filter the original full data to the intersection
    intersected_data <- result$node_attributes %>% 
      filter(Gene %in% intersecting_proteins)
    
    # Check if there are enough samples in both classes
    if (n_distinct(intersected_data$plot_group) < 2 || any(table(intersected_data$plot_group) < 5)) {
      warning(paste("Skipping group", group_name, "in intersection analysis for", df_name, "due to insufficient data post-intersection."))
      next
    }
    
    # Re-fit the baseline model on the intersected data
    logFC_col <- group_info[[group_name]]$logFC
    formula_str <- paste("plot_group ~ mean_group_medianNorm +", logFC_col)
    
    tryCatch({
      model_intersect <- glm(as.formula(formula_str), data = intersected_data, family = "binomial")
      predictor_intersect <- predict(model_intersect, type = "response")
      
      # Calculate and store the new ROC object
      roc_objects_intersect[[group_name]] <- roc(
        response = intersected_data$plot_group, 
        predictor = predictor_intersect, 
        levels = c("Other", "TP"), 
        quiet = TRUE, 
        direction = "<"
      )
    }, error = function(e) {
      warning(paste("Model re-fitting failed for group", group_name, "in intersection analysis for", df_name, ":", e$message))
    })
  }
  
  if (length(roc_objects_intersect) < 2) {
    warning(paste("Skipping intersection comparison for", df_name, "as it has fewer than 2 valid groups after re-fitting."))
    next
  }
  
  # --- Perform pairwise DeLong's test on intersection ROCs ---
  group_names_intersect <- names(roc_objects_intersect)
  p_matrix_intersect <- matrix(NA, nrow = length(group_names_intersect), ncol = length(group_names_intersect),
                               dimnames = list(group_names_intersect, group_names_intersect))
  
  group_pairs_intersect <- combn(group_names_intersect, 2, simplify = FALSE)
  
  p_values_intersect <- purrr::map_dbl(group_pairs_intersect, ~ {
    tryCatch({
      roc.test(roc_objects_intersect[[.x[1]]], roc_objects_intersect[[.x[2]]], method = "delong")$p.value
    }, error = function(e) { NA })
  })
  
  for (i in seq_along(group_pairs_intersect)) {
    pair <- group_pairs_intersect[[i]]
    p_matrix_intersect[pair[1], pair[2]] <- p_values_intersect[i]
    p_matrix_intersect[pair[2], pair[1]] <- p_values_intersect[i]
  }
  diag(p_matrix_intersect) <- 1.0
  
  # --- Consolidate and store results ---
  auc_ci_intersect <- purrr::map_df(roc_objects_intersect, ~{
    ci_vals <- pROC::ci.auc(.x, method = "delong")
    tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
  }, .id = "GroupName")
  
  p_values_long_intersect <- as.data.frame(p_matrix_intersect) %>% 
    tibble::rownames_to_column("Group1") %>% 
    tidyr::pivot_longer(-Group1, names_to = "Group2", values_to = "p_value") %>% 
    filter(Group1 < Group2) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
    
  all_intersection_comparisons[[df_name]] <- list(
    auc_ci_summary = auc_ci_intersect,
    pairwise_tests = p_values_long_intersect
  )
}

  # --- Export Intersection Analysis Results ---
if (length(all_intersection_comparisons) > 0) {
  export_list_intersect <- list()
  
  auc_sheets_intersect <- purrr::map(all_intersection_comparisons, "auc_ci_summary")
  names(auc_sheets_intersect) <- paste("I_AUC", names(auc_sheets_intersect), sep = "_")
  export_list_intersect <- c(export_list_intersect, auc_sheets_intersect)
  
  pairwise_sheets_intersect <- purrr::map(all_intersection_comparisons, "pairwise_tests")
  names(pairwise_sheets_intersect) <- paste("I_PVal", names(pairwise_sheets_intersect), sep = "_")
  export_list_intersect <- c(export_list_intersect, pairwise_sheets_intersect)
  
  writexl::write_xlsx(export_list_intersect, path = file.path(output_dir, "Group_Comparison_Intersection_Summary.xlsx"))
  cat("\n样本交集公平性分析结果已导出到 Group_Comparison_Intersection_Summary.xlsx\n")
} else {
  cat("\n未能生成任何样本交集公平性分析结果。\n")
}

} else {
  cat("\n稳健性分析已禁用，跳过样本交集公平性分析。\n")
}

# ===================================================================
#
###    Phase 1.2: Visualization (Isolated Data) #####
#
# ===================================================================

if (enable_visualization) {
  cat("\n\n--- 开始执行 Phase 1.2: 结果可视化 ---\n")

  # To prevent the analysis in this section from accidentally modifying the original
  # `all_results` and `group_info` objects, which are needed for Phase 2
  # and Phase 3, we create deep copies of them specifically for this section.
  all_results_p1_2 <- rlang::duplicate(all_results, shallow = FALSE)
  group_info_p1_2 <- rlang::duplicate(group_info, shallow = FALSE)


# Create a directory for plots
plots_dir <- file.path(output_dir, "Plots")
dir.create(plots_dir, showWarnings = FALSE)

# --- 1. Visualization for AUC_summary.xlsx: Faceted Forest Plot ---
if (enable_plot_auc_summary) {
  cat("\n正在为 AUC_summary.xlsx 创建分组森林图...\n")

  # Check if the final summary data frame exists and is not empty
  if (exists("final_auc_summary") && nrow(final_auc_summary) > 0) {

    # Prepare data for plotting directly from the final_auc_summary dataframe
    plot_data_auc <- final_auc_summary %>%
      # Exclude the random control as it clutters the plot
      filter(ModelName != "Control") %>%
      # Create a more descriptive label for faceting
      mutate(Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n")) %>%
      # For each group, find the AUC of the traditional method to use as a reference line
      group_by(Analysis_Group) %>%
      mutate(
        traditional_auc = tryCatch(AUC[ModelName == "logFC>0.5 & FDR<0.05"][1], error = function(e) NA),
        # The q-value is already calculated, so we just create the significance flag
        is_significant = if_else(
          !is.na(q_value_vs_traditional) & q_value_vs_traditional < 0.05, 
          "Significant (q < 0.05)", 
          "Not Significant"
        )
      ) %>%
      ungroup() %>%
      # Make ModelName a factor to control plotting order
      mutate(ModelName = factor(ModelName, levels = rev(c("logFC>0.5 & FDR<0.05", setdiff(unique(ModelName), "logFC>0.5 & FDR<0.05")))))

    # Generate the plot
    forest_plot <- ggplot(plot_data_auc, aes(y = ModelName, x = AUC, xmin = AUC_Lower_CI, xmax = AUC_Upper_CI)) +
      # Add reference line for the traditional method's AUC
      geom_vline(aes(xintercept = traditional_auc), linetype = "dashed", color = "grey50", na.rm = TRUE) +
      # Add points and error bars, colored by significance
      geom_point(aes(color = is_significant), size = 2.5) +
      geom_errorbarh(aes(color = is_significant), height = 0.2) +
      # Facet by the new combined label
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      # Customize colors and labels
      scale_color_manual(values = c("Not Significant" = "black", "Significant (q < 0.05)" = "#d62728")) +
      labs(
        title = "Model Performance Comparison within Each Analysis Group",
        subtitle = "AUCs compared to the traditional method (dashed line). Red models are significantly better (q < 0.05).",
        x = "Area Under Curve (AUC) with 95% Confidence Interval",
        y = "Model",
        color = "Significance vs. Traditional"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      # Set x-axis limits for better visualization
      coord_cartesian(xlim = c(0.4, 1.0))

    # Save the plot
    ggsave(
      filename = file.path(plots_dir, "1_AUC_Faceted_Forest_Plot.pdf"),
      plot = forest_plot,
      width = 16,
      height = 12,
      bg = "white"
    )

    cat("分组森林图已保存到:", file.path(plots_dir, "1_AUC_Faceted_Forest_Plot.pdf"), "\n")

  } else {
    cat("警告: 未找到 AUC_summary.xlsx，跳过分组森林图的生成。\n")
  }
} else {
  cat("\nAUC_summary分组森林图生成已禁用，跳过该图表。\n")
}


# --- 1b. Visualization for AUC_Summary_b.xlsx: Faceted Forest Plot vs 1D_logFC ---
if (enable_plot_auc_summary_vs_logFC) {
  cat("\n正在为 AUC_Summary_b.xlsx 创建分组森林图（vs 1D_logFC）...\n")

  # Check if the AUC_Summary_b data exists
  auc_summary_b_path <- file.path(output_dir, "AUC_Summary_b.xlsx")
  if (file.exists(auc_summary_b_path)) {
    auc_summary_b_data <- readxl::read_excel(auc_summary_b_path)
    
    if (nrow(auc_summary_b_data) > 0) {
      
      # Prepare data for plotting
      # 首先创建1D_logFC的AUC参考表
      logFC_auc_reference <- final_auc_summary %>%
        filter(ModelName == "1D_logFC") %>%
        select(Analysis_Group, logFC_auc = AUC)
      
      plot_data_auc_b <- auc_summary_b_data %>%
        # Create a more descriptive label for faceting
        mutate(Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n")) %>%
        # Join with 1D_logFC AUC as reference line
        left_join(logFC_auc_reference, by = "Analysis_Group") %>%
        # Create significance flag based on q_value_vs_1D_logFC
        mutate(
          is_significant = if_else(
            !is.na(q_value_vs_1D_logFC) & q_value_vs_1D_logFC < 0.05,
            "Significant (q < 0.05)",
            "Not Significant"
          )
        ) %>%
        # Make ModelName a factor to control plotting order
        mutate(ModelName = factor(ModelName, levels = rev(c("logFC>0.5 & FDR<0.05", "1D_Abundance", "2D_Abundance_logFC", setdiff(unique(ModelName), c("logFC>0.5 & FDR<0.05", "1D_Abundance", "2D_Abundance_logFC"))))))
      
      # Generate the forest plot
      forest_plot_b <- ggplot(plot_data_auc_b, aes(y = ModelName, x = AUC, xmin = AUC_Lower_CI, xmax = AUC_Upper_CI)) +
        # Add reference line for 1D_logFC's AUC
        geom_vline(aes(xintercept = logFC_auc), linetype = "dashed", color = "grey50", na.rm = TRUE) +
        # Add points and error bars, colored by significance
        geom_point(aes(color = is_significant), size = 2.5) +
        geom_errorbarh(aes(color = is_significant), height = 0.2) +
        # Facet by the combined label
        facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
        # Customize colors and labels
        scale_color_manual(values = c("Not Significant" = "black", "Significant (q < 0.05)" = "#d62728")) +
        labs(
          title = "Model Performance Comparison vs 1D_logFC within Each Analysis Group",
          subtitle = "AUCs compared to the 1D_logFC model (dashed line). Red models are significantly better (q < 0.05).",
          x = "Area Under Curve (AUC) with 95% Confidence Interval",
          y = "Model",
          color = "Significance vs. 1D_logFC"
        ) +
        theme_bw(base_size = 12) +
        theme(
          legend.position = "bottom",
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 7, face = "bold"),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        ) +
        # Set x-axis limits for better visualization
        coord_cartesian(xlim = c(0.4, 1.0))
      
      # Save the plot
      ggsave(
        filename = file.path(plots_dir, "1b_AUC_Faceted_Forest_Plot_vs_1D_logFC.pdf"),
        plot = forest_plot_b,
        width = 16,
        height = 12,
        bg = "white"
      )
      
      cat("与1D_logFC比较的分组森林图已保存到:", file.path(plots_dir, "1b_AUC_Faceted_Forest_Plot_vs_1D_logFC.pdf"), "\n")
      
    } else {
      cat("警告: AUC_Summary_b.xlsx 为空，跳过图1b的生成。\n")
    }
  } else {
    cat("警告: 未找到 AUC_Summary_b.xlsx，跳过图1b的生成。\n")
  }
} else {
  cat("\n与1D_logFC比较的分组森林图生成已禁用，跳过该图表。\n")
}


# --- 2. Visualization for Group_Comparison_Summary.xlsx: Bar Plots and Heatmaps ---
group_comp_path <- file.path(output_dir, "Group_Comparison_Summary.xlsx")
if (enable_plot_group_comparison) {
  cat("\n正在为 Group_Comparison_Summary.xlsx 创建可视化图表...\n")

  if (file.exists(group_comp_path)) {
    
    # Get all sheet names from the Excel file
    sheet_names <- readxl::excel_sheets(group_comp_path)
    
    # Identify the base names for each analysis (e.g., "ForStep19_APEX")
    base_names <- unique(sub("^(C_AUC_|C_PVal_)", "", sheet_names))

    for (base_name in base_names) {
      cat(paste("  - 正在处理组比较:", base_name, "\n"))
      
      auc_sheet_name <- paste0("C_AUC_", base_name)
      pval_sheet_name <- paste0("C_PVal_", base_name)
      
      # --- a) Ordered Bar Plot with CIs ---
      if (auc_sheet_name %in% sheet_names) {
        auc_data <- readxl::read_excel(group_comp_path, sheet = auc_sheet_name)
        
        if (nrow(auc_data) > 0) {
          # 按照selected_groups的顺序排列，确保只包含实际存在的组别
          groups_from_data <- unique(auc_data$GroupName)
          ordered_groups <- intersect(selected_groups, groups_from_data)
          missing_groups <- setdiff(groups_from_data, selected_groups)
          final_group_order <- c(ordered_groups, missing_groups)
          
          # 将GroupName转换为有序因子
          auc_data <- auc_data %>%
            mutate(GroupName = factor(GroupName, levels = final_group_order))
          
          bar_plot <- ggplot(auc_data, aes(x = GroupName, y = AUC)) +
            geom_col(fill = "#0072B2", width = 0.7) +
            geom_errorbar(aes(ymin = AUC_Lower_CI, ymax = AUC_Upper_CI), width = 0.2, color = "gray30") +
            labs(
              title = paste("Cross-Group AUC Comparison for:", base_name),
              subtitle = paste("Using baseline model '2D_Abundance_logFC' - Groups ordered by selection"),
              x = "Proximity Labeling Method (Group)",
              y = "AUC (with 95% CI)"
            ) +
            theme_minimal(base_size = 12) +
            theme(
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5, size = 3.5, color = "black") + 
            ylim(NA, max(auc_data$AUC_Upper_CI, na.rm = TRUE) * 1.05) # Adjust y-axis limit

          ggsave(
            filename = file.path(plots_dir, paste0("2a_Group_Comparison_Bar_Plot_", base_name, ".pdf")),
            plot = bar_plot,
            width = 10, height = 8, bg = "white"
          )
        }
      }
      
      # --- b) Pairwise Comparison Heatmap ---
      if (pval_sheet_name %in% sheet_names) {
        pval_data <- readxl::read_excel(group_comp_path, sheet = pval_sheet_name)
        
        if (nrow(pval_data) > 0) {
          # 检测是否为单次比较
          is_single_comparison <- nrow(pval_data) == 1 || 
                                 (exists("single_group_comparison") && any(pval_data$single_group_comparison == TRUE, na.rm = TRUE))
          
          # 确定使用p值还是q值，以及相应的标签
          value_column <- if (is_single_comparison) "p_value" else "q_value"
          value_label <- if (is_single_comparison) "p-value" else "q-value (BH-adj.)"
          subtitle_text <- if (is_single_comparison) {
            "Single comparison - using p-values (no multiple testing correction needed)"
          } else {
            "Multiple comparisons - using BH-adjusted q-values"
          }
          
          # 使用selected_groups的顺序，确保只包含实际存在的组别
          all_groups_from_data <- union(pval_data$Group1, pval_data$Group2)
          ordered_groups <- intersect(selected_groups, all_groups_from_data)
          # 如果有数据中存在但不在selected_groups中的组别，追加到末尾
          missing_groups <- setdiff(all_groups_from_data, selected_groups)
          all_groups <- c(ordered_groups, missing_groups)
          
          full_pval_data <- pval_data %>% 
            # Add the reverse pairs for a symmetric matrix
            bind_rows(rename(., Group1 = Group2, Group2 = Group1)) %>% 
            # Add diagonal comparisons
            bind_rows(tibble(Group1 = all_groups, Group2 = all_groups, 
                            p_value = NA, q_value = NA, single_group_comparison = NA)) %>%
            # 将组别转换为有序因子，按照selected_groups的顺序
            mutate(
              Group1 = factor(Group1, levels = all_groups),
              Group2 = factor(Group2, levels = all_groups)
            )

          # 动态选择填充的值
          heatmap_plot <- ggplot(full_pval_data, aes(x = Group1, y = Group2, fill = .data[[value_column]])) +
            geom_tile(color = "white", size = 0.5) +
            # Use a color scale that highlights significant values
            scale_fill_gradientn(
              colors = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF"),
              na.value = "grey80",
              limits = c(0, 1),
              name = value_label
            ) +
            geom_text(aes(label = ifelse(is.na(.data[[value_column]]), "", sprintf("%.3f", .data[[value_column]]))), 
                     size = 3, color = "white") +
            coord_fixed() + # Ensure tiles are square
            # 设置轴的顺序，取消默认的字母排序
            scale_x_discrete(limits = all_groups) +
            scale_y_discrete(limits = rev(all_groups)) + # y轴反转以匹配矩阵显示习惯
            labs(
              title = paste("Pairwise DeLong Test", value_label, "for:", base_name),
              subtitle = subtitle_text,
              x = "", y = ""
            ) +
            theme_minimal(base_size = 12) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right"
            )

          ggsave(
            filename = file.path(plots_dir, paste0("2b_Group_Comparison_Heatmap_", base_name, ".pdf")),
            plot = heatmap_plot,
            width = 10, height = 9, bg = "white"
          )
        }
      }
    }
    cat("跨邻近方法比较的可视化图表已保存。\n")
  } else {
    cat("警告: 未找到 Group_Comparison_Summary.xlsx，跳过相关图表的生成。\n")
  }
} else {
  cat("\n跨邻近方法比较图表生成已禁用，跳过该图表。\n")
}


# --- 3. Visualization for Intersection Analysis: Dumbbell Plot ---
intersection_comp_path <- file.path(output_dir, "Group_Comparison_Intersection_Summary.xlsx")
if (enable_plot_intersection_dumbbell) {
  cat("\n正在为样本交集公平性分析创建哑铃图...\n")

  # We need both files to create the comparison plot
  if (file.exists(group_comp_path) && file.exists(intersection_comp_path)) {

    # Get sheet names from both files
    union_sheets <- readxl::excel_sheets(group_comp_path)
    intersection_sheets <- readxl::excel_sheets(intersection_comp_path)
    
    # Identify the base names for which both analyses were successful
    union_base_names <- unique(sub("^C_AUC_", "", union_sheets[grep("^C_AUC_", union_sheets)]))
    intersection_base_names <- unique(sub("^I_AUC_", "", intersection_sheets[grep("^I_AUC_", intersection_sheets)]))
    common_base_names <- intersect(union_base_names, intersection_base_names)

    for (base_name in common_base_names) {
      cat(paste("  - 正在为以下对象创建哑铃图:", base_name, "\n"))
      
      # Load data from both analyses
      auc_union <- readxl::read_excel(group_comp_path, sheet = paste0("C_AUC_", base_name)) %>%
        mutate(AnalysisType = "Union (All Proteins)")
      
      auc_intersection <- readxl::read_excel(intersection_comp_path, sheet = paste0("I_AUC_", base_name)) %>%
        mutate(AnalysisType = "Intersection (Common Proteins)")
        
      # Combine the data
      combined_auc_data <- bind_rows(auc_union, auc_intersection)

      if (nrow(combined_auc_data) > 0) {
        # 按照selected_groups的顺序排列，确保只包含实际存在的组别
        groups_from_data <- unique(combined_auc_data$GroupName)
        ordered_groups <- intersect(selected_groups, groups_from_data)
        missing_groups <- setdiff(groups_from_data, selected_groups)
        final_group_order <- c(ordered_groups, missing_groups)
        
        # 将GroupName转换为有序因子
        combined_auc_data <- combined_auc_data %>%
          mutate(GroupName = factor(GroupName, levels = final_group_order))
        
        dumbbell_plot <- ggplot(combined_auc_data, aes(x = AUC, y = GroupName)) +
          # The line connecting the two points
          geom_line(aes(group = GroupName), color = "gray", linewidth = 1.5) +
          # The points for each analysis type
          geom_point(aes(color = AnalysisType), size = 4) +
          # Customize colors
          scale_color_manual(values = c("Union (All Proteins)" = "#1f77b4", "Intersection (Common Proteins)" = "#ff7f0e")) +
          labs(
            title = paste("Robustness to Protein Set for:", base_name),
            subtitle = "Comparing AUC on all available proteins (Union) vs. only common proteins (Intersection) - Groups ordered by selection",
            x = "AUC Score",
            y = "Proximity Labeling Method (Group)",
            color = "Analysis Set"
          ) +
          theme_bw(base_size = 12) +
          theme(
            legend.position = "top",
            panel.grid.major.y = element_line(linetype = "dashed")
          )

        ggsave(
          filename = file.path(plots_dir, paste0("3_Intersection_Dumbbell_Plot_", base_name, ".pdf")),
          plot = dumbbell_plot,
          width = 11, height = 8, bg = "white"
        )
      }
    }
    cat("样本交集公平性分析的哑铃图已保存。\n")
  } else {
    cat("警告: 未找到 Group_Comparison_Summary.xlsx 或 Group_Comparison_Intersection_Summary.xlsx，跳过哑铃图的生成。\n")
  }
} else {
  cat("\n样本交集公平性哑铃图生成已禁用，跳过该图表。\n")
}


# --- 4. Visualization for Model Choice Robustness: Performance vs. Stability Plot ---
robustness_path <- file.path(output_dir, "Group_Robustness_Summary.xlsx")
if (enable_plot_robustness_scatter) {
  cat("\n正在为模型选择稳健性创建性能-稳定性散点图...\n")

  if (file.exists(robustness_path)) {
    robustness_data <- readxl::read_excel(robustness_path)

    if (nrow(robustness_data) > 0) {
      # Create a new column for combined labels
      robustness_data <- robustness_data %>%
        mutate(Combined_Label = paste(Df_Name, Analysis_Group, sep = " - "))

      perf_stability_plot <- ggplot(robustness_data, aes(x = AUC_SD, y = Mean_AUC)) +
        geom_point(aes(color = Combined_Label), size = 5, alpha = 0.8) +
        # Add text labels with a bit of offset
        ggrepel::geom_text_repel(aes(label = Combined_Label), box.padding = 0.5, max.overlaps = Inf) +
        labs(
          title = "Performance vs. Stability of Proximity Labeling Methods",
          subtitle = "Across all defined 1D, 2D, and 3D models",
          x = "AUC Standard Deviation (Stability, lower is better)",
          y = "Mean AUC (Performance, higher is better)",
          caption = "Each point represents a method's average performance and variability across different models."
        ) +
        scale_color_viridis_d() + # Use a nice color scale
        theme_bw(base_size = 12) +
        theme(legend.position = "none") # Hide legend as labels are on the plot

      ggsave(
        filename = file.path(plots_dir, "4_Performance_vs_Stability_Scatter_Plot.pdf"),
        plot = perf_stability_plot,
        width = 12, height = 10, bg = "white"
      )
      cat("性能-稳定性散点图已保存。\n")
    } else {
      cat("警告: Group_Robustness_Summary.xlsx 为空，跳过散点图的生成。\n")
    }
  } else {
    cat("警告: 未找到 Group_Robustness_Summary.xlsx，跳过散点图的生成。\n")
  }
} else {
  cat("\n模型选择稳健性散点图生成已禁用，跳过该图表。\n")
}



# --- 5. Visualization for Controlled Variable Comparison ---
controlled_comp_path <- file.path(output_dir, "Controlled_Variable_Comparison.xlsx")
if (enable_plot_controlled_comparison) {
  cat("\n正在为控制变量比较创建增益图...\n")

if (file.exists(controlled_comp_path)) {
  controlled_comp_data <- readxl::read_excel(controlled_comp_path)
  
  if (nrow(controlled_comp_data) > 0) {
    
    # --- 5a. Fixed Sensitivity - Specificity Gain Plot ---
    cat("  - 正在创建固定Sensitivity下的Specificity增益图...\n")
    
    # 筛选指定的模型
    spec_gain_data <- controlled_comp_data %>%
      filter(ModelName %in% models_for_controlled_comparison) %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Specificity_Gain_Q_Value < 0.05, "显著 (q<0.05)", "不显著"),
        ModelName = factor(ModelName, levels = rev(models_for_controlled_comparison))
      )
    
    spec_gain_plot <- ggplot(spec_gain_data, aes(y = ModelName, x = Specificity_Gain, 
                                                   xmin = Specificity_Gain_CI_Lower, 
                                                   xmax = Specificity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("不显著" = "black", "显著 (q<0.05)" = "#d62728")) +
      labs(
        title = "控制变量比较：固定Sensitivity下的Specificity增益",
        subtitle = "在传统方法的Sensitivity水平下，连续模型的Specificity增益（正值表示改进）",
        x = "Specificity增益（相对于传统方法）",
        y = "模型",
        color = "统计显著性",
        caption = "误差棒表示95% Bootstrap置信区间"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "5a_Controlled_Comparison_Specificity_Gain.pdf"),
      plot = spec_gain_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 5b. Fixed Specificity - Sensitivity Gain Plot ---
    cat("  - 正在创建固定Specificity下的Sensitivity增益图...\n")
    
    # 筛选指定的模型
    sens_gain_data <- controlled_comp_data %>%
      filter(ModelName %in% models_for_controlled_comparison) %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "显著 (q<0.05)", "不显著"),
        ModelName = factor(ModelName, levels = rev(models_for_controlled_comparison))
      )
    
    sens_gain_plot <- ggplot(sens_gain_data, aes(y = ModelName, x = Sensitivity_Gain, 
                                                   xmin = Sensitivity_Gain_CI_Lower, 
                                                   xmax = Sensitivity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("不显著" = "black", "显著 (q<0.05)" = "#d62728")) +
      labs(
        title = "控制变量比较：固定Specificity下的Sensitivity增益",
        subtitle = "在传统方法的Specificity水平下，连续模型的Sensitivity增益（正值表示改进）",
        x = "Sensitivity增益（相对于传统方法）",
        y = "模型",
        color = "统计显著性",
        caption = "误差棒表示95% Bootstrap置信区间"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "5b_Controlled_Comparison_Sensitivity_Gain.pdf"),
      plot = sens_gain_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 5c. Combined Gain Comparison (Scatter Plot) ---
    cat("  - 正在创建综合增益对比散点图...\n")
    
    # 筛选指定的模型
    combined_gain_data <- controlled_comp_data %>%
      filter(ModelName %in% models_for_controlled_comparison) %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = " - "),
        Both_Significant = if_else(
          Specificity_Gain_Q_Value < 0.05 & Sensitivity_Gain_Q_Value < 0.05,
          "均显著", 
          if_else(
            Specificity_Gain_Q_Value < 0.05 | Sensitivity_Gain_Q_Value < 0.05,
            "部分显著",
            "均不显著"
          )
        ),
        ModelName = factor(ModelName, levels = models_for_controlled_comparison)
      )
    
    combined_gain_plot <- ggplot(combined_gain_data, 
                                  aes(x = Specificity_Gain, y = Sensitivity_Gain,
                                      xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper,
                                      ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper,
                                      color = Both_Significant, shape = ModelName)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      # 添加横向error bar（Specificity增益的CI）
      geom_errorbarh(aes(xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper),
                     height = 0, alpha = 0.5, size = 0.5) +
      # 添加纵向error bar（Sensitivity增益的CI）
      geom_errorbar(aes(ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper),
                    width = 0, alpha = 0.5, size = 0.5) +
      # 点在error bar上层
      geom_point(size = 4, alpha = 0.8) +
      facet_wrap(~ Combined_Label, scales = "free") +
      scale_color_manual(
        values = c("均显著" = "#2ca02c", "部分显著" = "#ff7f0e", "均不显著" = "grey60")
      ) +
      labs(
        title = "控制变量比较：Sensitivity与Specificity增益的综合对比",
        subtitle = "相对于传统方法(logFC>0.5 & FDR<0.05)的性能增益（含95% CI error bars）",
        x = "Specificity增益（固定Sensitivity时）",
        y = "Sensitivity增益（固定Specificity时）",
        color = "显著性状态",
        shape = "模型类型",
        caption = "右上象限表示两个指标均优于传统方法；误差棒表示95% Bootstrap置信区间"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim, ylim = controlled_comparison_ylim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(breaks = seq(controlled_comparison_ylim[1], controlled_comparison_ylim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      filename = file.path(plots_dir, "5c_Controlled_Comparison_Combined_Gain.pdf"),
      plot = combined_gain_plot,
      width = 16, height = 12, bg = "white"
    )
    
    cat("控制变量比较的可视化图表已保存。\n")
  } else {
    cat("警告: Controlled_Variable_Comparison.xlsx 为空，跳过相关图表生成。\n")
  }
} else {
  cat("警告: 未找到 Controlled_Variable_Comparison.xlsx，跳过控制变量比较的图表生成。\n")
}
} else {
  cat("\n控制变量比较增益图生成已禁用，跳过该图表。\n")
}


# --- 7. Visualization for Controlled Variable Comparison vs 1D_logFC ---
controlled_comp_vs_logFC_path <- file.path(output_dir, "Controlled_Variable_Comparison_vs_1D_logFC.xlsx")
if (enable_plot_controlled_comparison_vs_logFC) {
  cat("\n正在为与1D_logFC模型的控制变量比较创建增益图...\n")

if (file.exists(controlled_comp_vs_logFC_path)) {
  controlled_comp_vs_logFC_data <- readxl::read_excel(controlled_comp_vs_logFC_path)
  
  if (nrow(controlled_comp_vs_logFC_data) > 0) {
    
    # --- 7a. Fixed Sensitivity - Specificity Gain Plot (vs 1D_logFC) ---
    cat("  - 正在创建固定Sensitivity下的Specificity增益图（vs 1D_logFC）...\n")
    
    spec_gain_logFC_data <- controlled_comp_vs_logFC_data %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Specificity_Gain_Q_Value < 0.05, "显著 (q<0.05)", "不显著"),
        ModelName = factor(ModelName, levels = rev(unique(ModelName)))
      )
    
    spec_gain_logFC_plot <- ggplot(spec_gain_logFC_data, aes(y = ModelName, x = Specificity_Gain, 
                                                   xmin = Specificity_Gain_CI_Lower, 
                                                   xmax = Specificity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("不显著" = "black", "显著 (q<0.05)" = "#d62728")) +
      labs(
        title = "控制变量比较（vs 1D_logFC）：固定Sensitivity下的Specificity增益",
        subtitle = "在1D_logFC模型的Sensitivity水平下，其他模型的Specificity增益（正值表示改进）",
        x = "Specificity增益（相对于1D_logFC模型）",
        y = "模型",
        color = "统计显著性",
        caption = "误差棒表示95% Bootstrap置信区间"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "7a_Controlled_Comparison_vs_1D_logFC_Specificity_Gain.pdf"),
      plot = spec_gain_logFC_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 7b. Fixed Specificity - Sensitivity Gain Plot (vs 1D_logFC) ---
    cat("  - 正在创建固定Specificity下的Sensitivity增益图（vs 1D_logFC）...\n")
    
    sens_gain_logFC_data <- controlled_comp_vs_logFC_data %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = "\n"),
        Significant = if_else(Sensitivity_Gain_Q_Value < 0.05, "显著 (q<0.05)", "不显著"),
        ModelName = factor(ModelName, levels = rev(unique(ModelName)))
      )
    
    sens_gain_logFC_plot <- ggplot(sens_gain_logFC_data, aes(y = ModelName, x = Sensitivity_Gain, 
                                                   xmin = Sensitivity_Gain_CI_Lower, 
                                                   xmax = Sensitivity_Gain_CI_Upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 1) +
      geom_point(aes(color = Significant), size = 3) +
      geom_errorbarh(aes(color = Significant), height = 0.2) +
      facet_wrap(~ Combined_Label, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("不显著" = "black", "显著 (q<0.05)" = "#d62728")) +
      labs(
        title = "控制变量比较（vs 1D_logFC）：固定Specificity下的Sensitivity增益",
        subtitle = "在1D_logFC模型的Specificity水平下，其他模型的Sensitivity增益（正值表示改进）",
        x = "Sensitivity增益（相对于1D_logFC模型）",
        y = "模型",
        color = "统计显著性",
        caption = "误差棒表示95% Bootstrap置信区间"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.05)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(plots_dir, "7b_Controlled_Comparison_vs_1D_logFC_Sensitivity_Gain.pdf"),
      plot = sens_gain_logFC_plot,
      width = 16, height = 12, bg = "white"
    )
    
    # --- 7c. Combined Gain Comparison (Scatter Plot vs 1D_logFC) ---
    cat("  - 正在创建综合增益对比散点图（vs 1D_logFC）...\n")
    
    combined_gain_logFC_data <- controlled_comp_vs_logFC_data %>%
      mutate(
        Combined_Label = paste(Df_Name, Analysis_Group, sep = " - "),
        Both_Significant = if_else(
          Specificity_Gain_Q_Value < 0.05 & Sensitivity_Gain_Q_Value < 0.05,
          "均显著", 
          if_else(
            Specificity_Gain_Q_Value < 0.05 | Sensitivity_Gain_Q_Value < 0.05,
            "部分显著",
            "均不显著"
          )
        ),
        ModelName = factor(ModelName, levels = unique(ModelName))
      )
    
    combined_gain_logFC_plot <- ggplot(combined_gain_logFC_data, 
                                  aes(x = Specificity_Gain, y = Sensitivity_Gain,
                                      xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper,
                                      ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper,
                                      color = Both_Significant, shape = ModelName)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      # 添加横向error bar（Specificity增益的CI）
      geom_errorbarh(aes(xmin = Specificity_Gain_CI_Lower, xmax = Specificity_Gain_CI_Upper),
                     height = 0, alpha = 0.5, size = 0.5) +
      # 添加纵向error bar（Sensitivity增益的CI）
      geom_errorbar(aes(ymin = Sensitivity_Gain_CI_Lower, ymax = Sensitivity_Gain_CI_Upper),
                    width = 0, alpha = 0.5, size = 0.5) +
      # 点在error bar上层
      geom_point(size = 4, alpha = 0.8) +
      facet_wrap(~ Combined_Label, scales = "free") +
      scale_color_manual(
        values = c("均显著" = "#2ca02c", "部分显著" = "#ff7f0e", "均不显著" = "grey60")
      ) +
      labs(
        title = "控制变量比较（vs 1D_logFC）：Sensitivity与Specificity增益的综合对比",
        subtitle = "相对于1D_logFC模型的性能增益（含95% CI error bars）",
        x = "Specificity增益（固定Sensitivity时）",
        y = "Sensitivity增益（固定Specificity时）",
        color = "显著性状态",
        shape = "模型类型",
        caption = "右上象限表示两个指标均优于1D_logFC模型；误差棒表示95% Bootstrap置信区间"
      ) +
      coord_cartesian(xlim = controlled_comparison_xlim, ylim = controlled_comparison_ylim) +
      scale_x_continuous(breaks = seq(controlled_comparison_xlim[1], controlled_comparison_xlim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(breaks = seq(controlled_comparison_ylim[1], controlled_comparison_ylim[2], by = 0.1),
                         labels = scales::number_format(accuracy = 0.1)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      filename = file.path(plots_dir, "7c_Controlled_Comparison_vs_1D_logFC_Combined_Gain.pdf"),
      plot = combined_gain_logFC_plot,
      width = 16, height = 12, bg = "white"
    )
    
    cat("与1D_logFC模型的控制变量比较可视化图表已保存。\n")
  } else {
    cat("警告: Controlled_Variable_Comparison_vs_1D_logFC.xlsx 为空，跳过相关图表生成。\n")
  }
} else {
  cat("警告: 未找到 Controlled_Variable_Comparison_vs_1D_logFC.xlsx，跳过与1D_logFC的控制变量比较图表生成。\n")
}
} else {
  cat("\n与1D_logFC的控制变量增益图生成已禁用，跳过该图表。\n")
}


# --- 6. Visualization for Candidate Protein Localization Distribution ---
candidates_path <- file.path(output_dir, "detailed_candidates_by_model.xlsx")
if (enable_plot_localization_distribution) {
  cat("\n正在为候选蛋白的定位分布创建堆叠柱状图...\n")

# 读取候选蛋白数据

if (file.exists(candidates_path)) {
  candidates_data <- readxl::read_excel(candidates_path)
  
  if (nrow(candidates_data) > 0) {
    
    # --- 6a. 准备数据：统计每个模型+组别的定位分布 ---
    cat("  - 正在统计候选蛋白的定位分布...\n")
    
    # 只保留通过阈值的候选蛋白（PassThreshold == "Pass"）
    candidates_passed <- candidates_data %>%
      filter(PassThreshold == "Pass")
    
    # 统计每个模型、组别、定位类别的蛋白数量
    localization_summary <- candidates_passed %>%
      group_by(DataAnalysisMethod, ExperimentGroup, Model, MultiBait_Localization) %>%
      summarise(
        Count = n(),
        .groups = "drop"
      ) %>%
      # 计算每个模型+组别组合内的百分比
      group_by(DataAnalysisMethod, ExperimentGroup, Model) %>%
      mutate(
        Total_Count = sum(Count),
        Percentage = Count / Total_Count * 100,
        # 创建组合标签用于绘图
        Model_Group = paste(Model, ExperimentGroup, sep = " - ")
      ) %>%
      ungroup()
    
    # 导出详细统计数据
    localization_summary_export <- localization_summary %>%
      select(
        DataAnalysisMethod,
        ExperimentGroup, 
        Model,
        MultiBait_Localization,
        Count,
        Total_Count,
        Percentage
      ) %>%
      arrange(DataAnalysisMethod, ExperimentGroup, Model, desc(Count))
    
    writexl::write_xlsx(
      list(Localization_Distribution = localization_summary_export),
      path = file.path(output_dir, "Candidate_Localization_Distribution.xlsx")
    )
    cat("  - 定位分布统计数据已导出到 Candidate_Localization_Distribution.xlsx\n")
    
    # --- 6b. 为每个DataAnalysisMethod创建堆叠柱状图 ---
    
    # 定义一致的颜色方案（根据示例图）
    localization_colors <- c(
      "Nuclear&SGs" = "#4472C4",           # 蓝色
      "Nuclear&Cytosol&SGs" = "#70AD47",  # 青色
      "SGs&Other" = "#FFC000",             # 黄色
      "Cytosol&SGs" = "#C00000"            # 红色
    )
    
    # 获取所有唯一的DataAnalysisMethod
    data_methods <- unique(localization_summary$DataAnalysisMethod)
    
    for (data_method in data_methods) {
      cat(paste("  - 正在为", data_method, "创建定位分布图...\n"))
      
      # 筛选当前data method的数据
      plot_data_loc <- localization_summary %>%
        filter(DataAnalysisMethod == data_method) %>%
        # 按照选定组别的顺序排列
        mutate(
          ExperimentGroup = factor(ExperimentGroup, levels = selected_groups),
          # 确保定位类别有序（从底部到顶部）
          # 顺序：红色Cytosol&SGs(底) → 黄色SGs&Other → 绿色Nuclear&Cytosol&SGs → 蓝色Nuclear&SGs(顶)
          MultiBait_Localization = factor(
            MultiBait_Localization,
            levels = c("Cytosol&SGs", "SGs&Other", "Nuclear&Cytosol&SGs", "Nuclear&SGs")
          )
        )
      
      # 为x轴创建有序的标签
      # 格式：Model (Group)
      plot_data_loc <- plot_data_loc %>%
        mutate(
          X_Label = paste0(Model, "\n(", ExperimentGroup, ")"),
          # 创建排序键：先按ExperimentGroup再按Model
          Sort_Key = paste(
            sprintf("%02d", as.numeric(ExperimentGroup)),
            Model
          )
        ) %>%
        arrange(Sort_Key) %>%
        mutate(
          X_Label = factor(X_Label, levels = unique(X_Label))
        )
      
      # 创建堆叠柱状图
      loc_barplot <- ggplot(plot_data_loc, aes(x = X_Label, y = Percentage, fill = MultiBait_Localization)) +
        geom_col(width = 0.8, color = "white", size = 0.3) +
        # 在柱子顶部添加蛋白数量标注
        geom_text(
          data = plot_data_loc %>% 
            group_by(X_Label) %>% 
            summarise(
              Total_Count = first(Total_Count),
              Y_Position = 101,
              .groups = "drop"
            ),
          aes(x = X_Label, y = Y_Position, label = Total_Count),
          inherit.aes = FALSE,
          vjust = 0,
          size = 3,
          fontface = "bold"
        ) +
        # 在每个堆叠段内显示百分比（只显示>5%的）
        geom_text(
          aes(label = ifelse(Percentage > 5, paste0(round(Percentage), "%"), "")),
          position = position_stack(vjust = 0.5),
          size = 2.5,
          color = "white",
          fontface = "bold"
        ) +
        scale_fill_manual(
          values = localization_colors,
          name = "Localization",
          breaks = c("Nuclear&SGs", "Nuclear&Cytosol&SGs", "SGs&Other", "Cytosol&SGs")
        ) +
        labs(
          title = paste("候选蛋白定位分布 -", data_method),
          subtitle = paste("按模型和实验组统计通过阈值的", myTP_vector[1], "相关蛋白的定位分布"),
          x = "模型 (实验组)",
          y = "相对丰度 (%)",
          caption = "柱子顶部数字表示该组合的候选蛋白总数"
        ) +
        scale_y_continuous(
          limits = c(0, 105),
          breaks = seq(0, 100, 25),
          expand = c(0, 0)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        )
      
      # 保存图表
      ggsave(
        filename = file.path(plots_dir, paste0("6_Localization_Distribution_", data_method, ".pdf")),
        plot = loc_barplot,
        width = max(14, length(unique(plot_data_loc$X_Label)) * 0.8),
        height = 10,
        bg = "white"
      )
    }
    
    # --- 6c. 创建对比图：比较不同模型的定位偏好 ---
    cat("  - 正在创建模型间定位分布对比图...\n")
    
    # 计算每个模型在所有组别上的平均定位分布
    model_avg_distribution <- localization_summary %>%
      group_by(DataAnalysisMethod, Model, MultiBait_Localization) %>%
      summarise(
        Avg_Percentage = mean(Percentage, na.rm = TRUE),
        Total_Proteins = sum(Count),
        .groups = "drop"
      ) %>%
      mutate(
        # 确保定位类别有序（从底部到顶部）
        # 顺序：红色Cytosol&SGs(底) → 黄色SGs&Other → 绿色Nuclear&Cytosol&SGs → 蓝色Nuclear&SGs(顶)
        MultiBait_Localization = factor(
          MultiBait_Localization,
          levels = c("Cytosol&SGs", "SGs&Other", "Nuclear&Cytosol&SGs", "Nuclear&SGs")
        )
      )
    
    for (data_method in data_methods) {
      plot_data_model_comp <- model_avg_distribution %>%
        filter(DataAnalysisMethod == data_method) %>%
        # 按模型维度排序
        mutate(
          Model = factor(Model, levels = c(
            "logFC>0.5 & FDR<0.05",
            setdiff(unique(Model), "logFC>0.5 & FDR<0.05")
          ))
        )
      
      model_comp_plot <- ggplot(plot_data_model_comp, 
                                 aes(x = Model, y = Avg_Percentage, fill = MultiBait_Localization)) +
        geom_col(width = 0.7, color = "white", size = 0.3) +
        geom_text(
          data = plot_data_model_comp %>% 
            group_by(Model) %>% 
            summarise(Total_Proteins = sum(Total_Proteins), .groups = "drop"),
          aes(x = Model, y = 102, label = Total_Proteins),
          inherit.aes = FALSE,
          vjust = 0,
          size = 3.5,
          fontface = "bold"
        ) +
        geom_text(
          aes(label = ifelse(Avg_Percentage > 5, paste0(round(Avg_Percentage), "%"), "")),
          position = position_stack(vjust = 0.5),
          size = 3,
          color = "white",
          fontface = "bold"
        ) +
        scale_fill_manual(
          values = localization_colors,
          name = "Localization",
          breaks = c("Nuclear&SGs", "Nuclear&Cytosol&SGs", "SGs&Other", "Cytosol&SGs")
        ) +
        labs(
          title = paste("模型间定位偏好对比 -", data_method),
          subtitle = "各模型在所有实验组上的平均定位分布",
          x = "模型",
          y = "平均相对丰度 (%)",
          caption = "数字表示该模型在所有组别中的候选蛋白总数"
        ) +
        scale_y_continuous(
          limits = c(0, 110),
          breaks = seq(0, 100, 25),
          expand = c(0, 0)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        )
      
      ggsave(
        filename = file.path(plots_dir, paste0("6_Model_Localization_Comparison_", data_method, ".pdf")),
        plot = model_comp_plot,
        width = 12,
        height = 8,
        bg = "white"
      )
    }
    
    cat("候选蛋白定位分布可视化已完成。\n")
    
  } else {
    cat("警告: detailed_candidates_by_model.xlsx 为空，跳过定位分布分析。\n")
  }
} else {
  cat("警告: 未找到 detailed_candidates_by_model.xlsx，跳过定位分布分析。\n")
}
} else {
  cat("\n候选蛋白定位分布图生成已禁用，跳过该图表。\n")
}


cat("\n--- 可视化阶段完成 ---\n")

} else {
  cat("\n可视化生成已禁用，跳过图表创建。\n")
}


# --- New Section: Cross-Group (Inter-Method) Comparison ---
if (enable_group_comparison) {
  cat("\n开始执行跨邻近方法比较分析（重复部分）...\n")

# Define the baseline model for this comparison
baseline_model_for_comparison <- "2D_Abundance_logFC"

# --- FIX: Create a dedicated deep copy for this section to avoid side effects ---
all_results_for_comparison <- rlang::duplicate(all_results, shallow = FALSE)

# Group results by the original dataframe name (e.g., "ForStep19_APEX", "ForStep19_BioID2")
results_by_df <- all_results_for_comparison %>% 
  purrr::keep(~!is.null(.x)) %>% # Filter out failed runs
  split(., purrr::map_chr(., "df_name"))

# List to store the comparison results for each dataframe
all_group_comparisons <- list()

for (df_name in names(results_by_df)) {
  
  # Extract the ROC objects for the baseline model from each group
  roc_objects_for_comparison <- results_by_df[[df_name]] %>% 
    purrr::map(~.x$roc_objects[[baseline_model_for_comparison]]) %>% 
    # Name the list elements by their group name for clarity
    setNames(purrr::map_chr(results_by_df[[df_name]], "group_name")) %>% 
    # Remove any NULL entries if a model failed for a specific group
    purrr::compact()

  if (length(roc_objects_for_comparison) < 2) {
    warning(paste("Skipping group comparison for", df_name, "as it has fewer than 2 valid groups for the model", baseline_model_for_comparison))
    next
  }
  
  # --- Perform pairwise DeLong's test ---
  group_names <- names(roc_objects_for_comparison)
  p_value_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                           dimnames = list(group_names, group_names))
  
  group_pairs <- combn(group_names, 2, simplify = FALSE)
  
  p_values <- purrr::map_dbl(group_pairs, ~ {
    tryCatch({
      test_result <- roc.test(roc_objects_for_comparison[[.x[1]]], roc_objects_for_comparison[[.x[2]]], method = "delong")
      test_result$p.value
    }, error = function(e) { NA })
  })
  
  for (i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    p_value_matrix[pair[1], pair[2]] <- p_values[i]
    p_value_matrix[pair[2], pair[1]] <- p_values[i]
  }
  diag(p_value_matrix) <- 1.0

  # --- Consolidate results into a summary table ---
  # Get AUC and CIs for each group
  auc_ci_summary <- purrr::map_df(roc_objects_for_comparison, ~{
    ci_vals <- pROC::ci.auc(.x, method = "delong")
    tibble::tibble(AUC = as.numeric(ci_vals[2]), AUC_Lower_CI = as.numeric(ci_vals[1]), AUC_Upper_CI = as.numeric(ci_vals[3]))
  }, .id = "GroupName")
  
  # Convert p-value matrix to a long format and apply BH correction
  p_values_long_temp <- as.data.frame(p_value_matrix) %>% 
    tibble::rownames_to_column("Group1") %>% 
    tidyr::pivot_longer(-Group1, names_to = "Group2", values_to = "p_value") %>% 
    filter(Group1 < Group2) # Avoid duplicates and self-comparisons
  
  p_values_long <- p_values_long_temp %>%
    mutate(
      q_value = if (nrow(p_values_long_temp) == 1) {
        # 单次组间比较时，q值等于p值
        p_value
      } else {
        # 多次组间比较时，应用BH校正
        p.adjust(p_value, method = "BH")
      },
      single_group_comparison = (nrow(p_values_long_temp) == 1)
    )

  # Store results
  all_group_comparisons[[df_name]] <- list(
    auc_ci_summary = auc_ci_summary,
    pairwise_tests = p_values_long
  )
}

# --- Export Cross-Group Comparison Results ---
if (length(all_group_comparisons) > 0) {
  # Create a single Excel file with sheets for each type of result
  export_list_group_comp <- list()
  
  # Add AUC/CI summaries for each dataset
  auc_sheets <- purrr::map(all_group_comparisons, "auc_ci_summary")
  names(auc_sheets) <- paste("AUCs", names(auc_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, auc_sheets)
  
  # Add pairwise test results for each dataset
  pairwise_sheets <- purrr::map(all_group_comparisons, "pairwise_tests")
  names(pairwise_sheets) <- paste("P-Values", names(pairwise_sheets), sep = "_")
  export_list_group_comp <- c(export_list_group_comp, pairwise_sheets)
  
  writexl::write_xlsx(export_list_group_comp, path = file.path(output_dir, "Group_Comparison_Summary.xlsx"))
  cat("\n跨邻近方法比较结果已导出到 Group_Comparison_Summary.xlsx\n")
} else {
  cat("\n未能生成任何跨邻近方法比较结果。\n")
}

} else {
  cat("\n组间比较分析已禁用，跳过重复的组间比较分析。\n")
}

# -------------------------------------------------------------------
# 7. Shutdown
# -------------------------------------------------------------------

# --- (REMOVED) Save Phase 1 Results to RData ---
# cat("\n正在将 Phase 1 的关键结果保存到 RData 文件中...\n")
# save(
#   all_results, 
#   all_group_comparisons, 
#   final_auc_summary,
#   file = file.path(output_dir, "analysis_phase1_output.RData")
# )
# cat("Phase 1 结果已成功保存到", file.path(output_dir, "analysis_phase1_output.RData"), "\n")

# Gracefully terminate the parallel backend
plan(sequential)
cat("\n并行会话已终止。脚本执行完毕。\n")


save.image(file = "Add_Phase1_2_alldata.RData")
print("工作环境已成功保存到 Add_Phase1_2_alldata.RData")
