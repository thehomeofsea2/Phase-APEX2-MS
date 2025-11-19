# ============================================================================
# 质谱分析流程 - 主控脚本
# ============================================================================
# 说明：
# 1. 逐个运行模块，每个模块完成后检查输出
# 2. 模块间通过 RData 文件传递数据
# 3. 首次运行前确保 Reference/ 和 Rawdata/ 目录存在
# ============================================================================

# ============================================================================
# 全局配置
# ============================================================================

# 工作目录（根据实际情况修改）
setwd("D:/Bioinfomatics/MS/CodeNorm_SG_Batch2")

# 加载必需包
required_packages <- c(
  "readr", "tidyverse", "readxl", "purrr", "ggplot2", 
  "ggsci", "ggthemes", "gridExtra", "pheatmap", "pROC", 
  "openxlsx", "preprocessCore", "limma", "scales", "ggrepel"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("✓ 包加载完成\n")

# 目录配置
dir_config <- list(
  root = getwd(),
  module = file.path(getwd(), "Module"),
  dev = file.path(getwd(), "Dev"),
  output = file.path(getwd(), "Output"),
  reference = file.path(getwd(), "Reference"),
  rawdata = file.path(getwd(), "Rawdata")
)

# 数据读取配置
data_config <- list(
  # TSV文件名匹配模式（正则表达式）
  file_pattern = "SG_batch2*\\.tsv$"
  # 可根据实际文件名修改，例如：
  # file_pattern = ".*\\.tsv$"  # 所有TSV文件
  # file_pattern = "sg_.*\\.tsv$"  # sg_开头的TSV文件
)

# ============================================================================
# Module 1: 环境初始化
# ============================================================================
cat("\n========================================\n")
cat("Module 1: 环境初始化\n")
cat("========================================\n")

source(file.path(dir_config$module, "module01_setup.R"))
config <- module01_setup(dir_config)

# 保存Module01环境（包含所有数据）
save.image(file = "Module01_workspace.RData")
cat("✓ 已保存: Module01_workspace.RData (工作目录)\n")
cat("  包含：dir_config, config\n")

cat("✓ Module 1 完成\n")

# ============================================================================
# Module 2: 数据读取与分组表
# ============================================================================
cat("\n========================================\n")
cat("Module 2: 数据读取与分组表\n")
cat("========================================\n")

source(file.path(dir_config$module, "module02_data_import.R"))

# 步骤1: 生成模板
result <- module02_read_and_generate_template(dir_config, file_pattern = data_config$file_pattern)

# 步骤2: 等待用户填写
cat("\n⚠ 请完成以下步骤后继续:\n")
cat("  1. 打开 Module02_sampleGroup_template.csv\n")
cat("  2. 填写所有必需列（bioGroup必填）\n")
cat("  3. 可选：修改Order列调整数据列顺序\n")
cat("  4. 保存 Module02_sampleGroup_template.csv\n")
cat("  5. 按回车继续（程序将自动生成Module02_sampleGroup.csv）...\n")
readline()

# 步骤3: 处理用户填写的数据
result <- module02_process_samplegroup(dir_config, result$data)

# 将返回的数据赋值给全局变量
data_raw <- result$data
sampleGroup <- result$sampleGroup

# 保存Module02环境（包含Module01所有数据 + Module02数据）
save.image(file = "Module02_workspace.RData")
cat("✓ 已保存: Module02_workspace.RData (工作目录)\n")
cat("  包含：Module01所有数据 + data_raw + sampleGroup\n")

cat("✓ Module 2 完成\n")

# ============================================================================
# Module 3: 注释系统
# ============================================================================
cat("\n========================================\n")
cat("Module 3: 注释系统\n")
cat("========================================\n")

source(file.path(dir_config$module, "module03_annotation.R"))

# 选择注释模式：
#   - "SG": 使用HaloMap/GO/MultiBait SGs预设
#   - "Nucleolus": 使用HPA Nucleolus + CLL核仁参考
annotation_mode <- "SG"

# 高级：如需自定义TP注释列，可在此提供列表（设置后将忽略annotation_mode）
custom_annotations_override <- NULL
# 示例：
# custom_annotations_override <- list(
#   list(
#     column_name = "Custom_Localization",
#     TP_source = "GO_SGs",
#     TP_column = "Gene",
#     TP_label = "SGs"
#   )
# )

# 可选：在选择SG/Nucleolus后追加额外注释列（不会覆盖原有配置）
additional_annotations <- list()
# 示例：
# additional_annotations <- list(
#   list(
#     column_name = "MyTurbo_Localization",
#     TP_source = "MyTurboSet",
#     TP_column = "Gene",
#     TP_label = "TurboTP"
#   )
# )

# 可选：声明自定义TP参考文件（可放在Reference目录或提供绝对路径）
custom_tp_sources <- list()
# 示例：
# custom_tp_sources <- list(
#   list(
#     source_name = "MyNucleolusSet",
#     file = "MyNucleolus.csv",   # 默认相对Reference目录
#     file_type = "csv",          # 可省略，自动根据扩展名判断
#     gene_column = "Gene"        # 文件中基因列名称
#   )
# )

result <- module03_annotation(
  dir_config,
  data_raw,
  sampleGroup,
  custom_annotations = custom_annotations_override,
  annotation_mode = annotation_mode,
  custom_tp_sources = custom_tp_sources,
  additional_annotations = additional_annotations
)

# 如果只需要基础注释（Cytosol/Nuclear/Mitochondrion），将custom_annotations_override保持为NULL，
# 并在module03_annotation内部选择"SG"/"Nucleolus"模式即可

# 将返回的数据赋值给全局变量
data_annotated <- result$data_annotated
annotation_references <- result$annotation_references

# 保存Module03环境（包含Module01-02所有数据 + Module03数据）
save.image(file = "Module03_workspace.RData")
cat("✓ 已保存: Module03_workspace.RData (工作目录)\n")
cat("  包含：Module01-02所有数据 + data_annotated + annotation_references\n")

cat("✓ Module 3 完成\n")

# ============================================================================
# Module 4: 标准化
# ============================================================================
cat("\n========================================\n")
cat("Module 4: 标准化\n")
cat("========================================\n")

source(file.path(dir_config$module, "module04_standardization.R"))

# 配置标准化类型
# 可选值："noNorm", "Global_QNorm", "Global_MNorm", "Local_QNorm", "Local_MNorm"
# 注意：必须至少包含一个全局标准化（Global_QNorm或Global_MNorm）
norm_types <- c("noNorm","Global_QNorm", "Local_QNorm")

result <- module04_standardization(dir_config, data_annotated, sampleGroup, norm_types)

# 将返回的数据赋值给全局变量
standardized_data_list <- result$standardized_data_list
norm_types_used <- result$norm_types_used

# 保存Module04环境（包含Module01-03所有数据 + Module04数据）
save.image(file = "Module04_workspace.RData")
cat("✓ 已保存: Module04_workspace.RData (工作目录)\n")
cat("  包含：Module01-03所有数据 + standardized_data_list + norm_types_used\n")

cat("✓ Module 4 完成\n")

# ============================================================================
# Module 5: 缺失值填补
# ============================================================================
cat("\n========================================\n")
cat("Module 5: 缺失值填补\n")
cat("========================================\n")

source(file.path(dir_config$module, "module05_imputation.R"))

# 配置填补参数
# impute_cat_mean: 是否对Cat组的n_valid=2情况用平均值填补（默认FALSE）
# random_seed: 随机数种子，用于保证Perseus填补的可重复性
impute_cat_mean <- FALSE  # 可改为TRUE以启用Cat组填补
random_seed <- 123        # 可改为其他数值以改变随机结果

result <- module05_imputation(dir_config, standardized_data_list, sampleGroup, 
                              impute_cat_mean, random_seed)

# 将返回的数据赋值给全局变量
imputed_data_list <- result$imputed_data_list
imputation_params <- result$imputation_params

# 保存Module05环境（包含Module01-04所有数据 + Module05数据）
save.image(file = "Module05_workspace.RData")
cat("✓ 已保存: Module05_workspace.RData (工作目录)\n")
cat("  包含：Module01-04所有数据 + imputed_data_list + imputation_params\n")

cat("✓ Module 5 完成\n")

# ============================================================================
# Module 6: 热图系统
# ============================================================================
cat("\n========================================\n")
cat("Module 6: 热图系统\n")
cat("========================================\n")

source(file.path(dir_config$module, "module06_heatmap.R"))

# 配置热图参数
# selected_versions: 选择用于绘制热图的标准化版本（支持向量）
# heatmap_types: 需要绘制的热图类型
# correlation_config: 相关性热图配置
# localization_columns: 用于分类的注释列（NULL表示自动检测）
# color_params: 热图颜色参数（用于AllLocalization和by_localization）

selected_versions <- c("noNorm_Imputed","Global_QNorm_Imputed", "Local_QNorm_Imputed")  # 支持多个版本
#selected_versions <- c("Global_QNorm_Imputed", "Local_QNorm_Imputed")  # 支持多个版本
## 可选值："noNorm_Imputed", "Global_QNorm_Imputed", "Global_MNorm_Imputed", "Local_QNorm_Imputed", "Local_MNorm _Imputed"
heatmap_types <- c("all", "correlation")
#heatmap_types <- c("all", "correlation", "by_localization")
# 相关性热图配置
correlation_config <- list(
  # 相关性范围和颜色
  corr_min = 0.8,
  corr_max = 1,
  corr_center = 0.9,
  col_low = "#4D97CD",
  col_center = "white",
  col_high = "#DB6968",
  
  # 方式1：使用排除规则（排除不需要的样本）
  exclude_context = NULL,       # 例如：c("Control") 排除所有Control组
  exclude_pltype = NULL,         # 例如：c("PL") 排除所有PL组
  exclude_catalytic = c("NoCat"),      # 例如：c("NoCat") 排除所有NoCat组
  
  # 方式2：直接指定包含的bioGroup（优先于排除规则）
  include_biogroups = NULL       # 例如：c("Light_Experiment", "H2O2_Experiment")
)

# 如果要使用排除规则示例：
# correlation_config$exclude_context <- c("Control")  # 只保留Experiment组
# correlation_config$exclude_pltype <- c("PL")        # 排除PL组

# 如果要直接指定bioGroup示例：
# correlation_config$include_biogroups <- c("Light_Experiment", "H2O2_Experiment")

localization_columns <- NULL  # NULL表示自动检测所有Localization列

# AllLocalization和by_localization热图的颜色参数
color_params <- list(
  custom_min = -4,
  custom_max = 4,
  custom_center = 0,
  col_low = "#4D97CD",
  col_center = "white",
  col_high = "#DB6968"
)

result <- module06_heatmap(dir_config, imputed_data_list, sampleGroup,
                          selected_versions, heatmap_types, 
                          correlation_config, localization_columns, color_params)

# 将返回的信息赋值给全局变量
heatmap_info <- result

# 保存Module06环境（包含Module01-05所有数据 + Module06数据）
save.image(file = "Module06_workspace.RData")
cat("✓ 已保存: Module06_workspace.RData (工作目录)\n")
cat("  包含：Module01-05所有数据 + heatmap_info\n")

cat("✓ Module 6 完成\n")

# ============================================================================
# Module 7: 第一次差异分析
# ============================================================================
cat("\n========================================\n")
cat("Module 7: 第一次差异分析\n")
cat("========================================\n")

source(file.path(dir_config$module, "module07_diff_analysis1.R"))

# 配置差异分析参数
# Module 7 会根据sampleGroup$FirstROCgroup自动构建对比组
# FirstROCgroup逻辑：
#   - 包含相同元素的组为同一组（如A、A/B、A&B都属于A组）
#   - 在同一组内，Context为Experiment的bioGroup vs Context为Control的bioGroup
#
# selected_versions: 需要分析的数据版本（NULL表示所有版本）
#   - 可选: c("noNorm_Imputed", "Local_QNorm_Imputed", "Global_QNorm_Imputed", etc.)

# 选择分析版本（NULL表示所有版本）
selected_versions_diff <- c("Local_QNorm_Imputed")

result <- module07_diff_analysis1(dir_config, imputed_data_list, sampleGroup,
                                   selected_versions_diff)

# 将返回的信息赋值给全局变量
diff_results1 <- result$diff_results1
comparisons_used <- result$comparisons_used
comparison_info <- result$comparison_info

# 保存Module07环境（包含Module01-06所有数据 + Module07数据）
save.image(file = "Module07_workspace.RData")
cat("✓ 已保存: Module07_workspace.RData (工作目录)\n")
cat("  包含：Module01-06所有数据 + diff_results1, comparisons_used, comparison_info\n")

cat("✓ Module 7 完成\n")

# ============================================================================
# Module 8: 第一次ROC分析
# ============================================================================
cat("\n========================================\n")
cat("Module 8: 第一次ROC分析\n")
cat("========================================\n")

source(file.path(dir_config$module, "module08_roc_analysis1.R"))

# 配置ROC分析参数
# SubMito转化：将注释中的Mitochondrion替换为MitoCarta3的亚定位（MIM, Matrix, MOM, IMS）
# ROC分析：使用pROC包，基于logFC和注释列计算ROC曲线和最佳阈值
#
# selected_versions: 需要分析的数据版本（NULL表示所有版本）
#   - 可选: c("noNorm_Imputed", "Local_QNorm_Imputed", "Global_QNorm_Imputed", etc.)
# roc_annotation_column: 用于ROC分析的注释列（默认"GO_Localization"）
# tp_label: True Positive标签（默认"SGs"）
# fp_label: False Positive标签（默认"Matrix"）
# enable_submito: 是否启用SubMito转化（默认TRUE）
# submito_annotation_columns: SubMito转化的注释列（NULL表示自动检测所有Localization列）
# min_tp: 最小TP阈值（默认0.3，用于筛选最佳阈值）

# 选择分析版本（NULL表示所有版本）
selected_versions_roc <- c("Local_QNorm_Imputed")

# ROC分析参数
roc_annotation_column <- "GO_Localization"  # 用于ROC分析的注释列
tp_label <- "SGs"                            # True Positive标签
fp_label <- "Matrix"                         # False Positive标签（SubMito转化后）
enable_submito <- TRUE                       # 启用SubMito转化
submito_annotation_columns <- NULL           # 自动检测所有Localization列
min_tp <- 0.3                                # 最小TP阈值

result <- module08_roc_analysis1(dir_config, diff_results1, annotation_references,
                                  selected_versions_roc,
                                  roc_annotation_column,
                                  tp_label,
                                  fp_label,
                                  enable_submito,
                                  submito_annotation_columns,
                                  min_tp)

# 将返回的信息赋值给全局变量
roc_results1 <- result$roc_results
roc_thresholds1 <- result$all_thresholds
data_with_submito <- result$data_with_submito
expr_fdr_df_list <- result$expr_fdr_df_list

# 保存Module08环境（包含Module01-07所有数据 + Module08数据）
save.image(file = "Module08_workspace.RData")
cat("✓ 已保存: Module08_workspace.RData (工作目录)\n")
cat("  包含：Module01-07所有数据 + roc_results1, roc_thresholds1, data_with_submito, expr_fdr_df_list\n")

cat("✓ Module 8 完成\n")

# ============================================================================
# Module 9: 背景扣除
# ============================================================================
cat("\n========================================\n")
cat("Module 9: 背景扣除\n")
cat("========================================\n")

source(file.path(dir_config$module, "module09_background_subtraction.R"))

# 配置背景扣除参数
# 自动生成Method A和Method B两种结果：
# - Method A：对指定comparison不使用FDR过滤
# - Method B：对所有comparison使用FDR过滤
# - Context为"Experiment"的bioGroup：使用阈值过滤
# - Context为"Spatial"的bioGroup：不使用阈值过滤
# - Context为"Control"的bioGroup：不操作
#
# selected_versions: 需要分析的数据版本（NULL表示所有版本）
# fdr_threshold: FDR阈值（默认0.05）
# no_fdr_comparisons: Method A中不使用FDR过滤的comparison向量（默认NULL，即按Method B处理）
# use_roc_threshold: 是否使用ROC threshold（TRUE）还是固定FC阈值（FALSE）
# fixed_fc_threshold: 如果不使用ROC threshold，使用的固定FC阈值
# min_valid_lfq: 最小有效LFQ值个数（默认2，用于催化组过滤）
# annotation_column: 用于分组和作图的注释列（默认"GO_Localization"）
# plot_all_annotations: 是否显示所有注释（TRUE）还是只显示TP（FALSE）
# tp_label: TP标签（默认"SGs"）
# tp_color: TP颜色（默认"#DB6968"）

# 选择分析版本（NULL表示所有版本）
selected_versions_bg <- c("Local_QNorm_Imputed")

# 背景扣除参数
fdr_threshold <- 0.05
# Method A：指定不使用FDR过滤的comparison向量（如果为NULL，Method A等同于Method B）
no_fdr_comparisons <- c("K69A1B3_Light_vs_A1B3_Light","K69C3_Light_vs_C3_Light")  # 可设为NULL使用默认Method B
use_roc_threshold <- TRUE
fixed_fc_threshold <- NULL  # 如果use_roc_threshold=FALSE，设置固定阈值如2.0
min_valid_lfq <- 2
annotation_column <- "GO_Localization"
plot_all_annotations <- TRUE  # FALSE只显示TP，TRUE显示所有注释
tp_label <- "SGs"
tp_color <- "#DB6968"

result <- module09_background_subtraction(
  dir_config, sampleGroup, diff_results1, roc_thresholds1, comparison_info,
  data_with_submito,
  expr_fdr_df_list,
  selected_versions_bg, fdr_threshold, no_fdr_comparisons,
  use_roc_threshold, fixed_fc_threshold, min_valid_lfq,
  annotation_column, plot_all_annotations, tp_label, tp_color
)

# 将返回的信息赋值给全局变量
filtered_data_A <- result$filtered_data_A
filtered_data_B <- result$filtered_data_B
merged_data_A <- result$merged_data_A
merged_data_B <- result$merged_data_B

# 保存Module09环境（包含Module01-08所有数据 + Module09数据）
save.image(file = "Module09_workspace.RData")
cat("✓ 已保存: Module09_workspace.RData (工作目录)\n")
cat("  包含：Module01-08所有数据 + filtered_data_A, filtered_data_B, merged_data_A, merged_data_B\n")

cat("✓ Module 9 完成\n")

# ============================================================================
# Module 10: 数据替换
# ============================================================================
cat("\n========================================\n")
cat("Module 10: 数据替换\n")
cat("========================================\n")

source(file.path(dir_config$module, "module10_data_replacement.R"))

# 使用 Module 5 的 imputed_data_list，优先 Global_QNorm_Imputed（无则用 Global_MNorm_Imputed）
# 对 Module 9 的 merged 数据（AfterROC_merged）进行定点替换
# 替换逻辑：保留所有NA，所有非NA值用Global标准化结果替换
result <- module10_data_replacement(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  merged_data_A = merged_data_A,        # 使用merged数据（不是filtered_data）
  merged_data_B = merged_data_B,        # 使用merged数据（不是filtered_data）
  imputed_data_list = imputed_data_list,
  selected_versions = selected_versions_bg,
  prefer_version = "Global_QNorm_Imputed",
  fallback_version = "Global_MNorm_Imputed",
  test_n = 20
)

# 将返回的数据赋值给全局变量
replaced_data_A <- result$replaced_data_A
replaced_data_B <- result$replaced_data_B

# 保存Module10环境（包含Module01-09所有数据 + Module10数据）
save.image(file = "Module10_workspace.RData")
cat("✓ 已保存: Module10_workspace.RData (工作目录)\n")
cat("  包含：Module01-09所有数据 + replaced_data_A, replaced_data_B\n")

cat("✓ Module 10 完成\n")


load("Module10_workspace.RData")
# ============================================================================
# Module 11: 第二次差异分析
# ============================================================================
cat("\n========================================\n")
cat("Module 11: 第二次差异分析\n")
cat("========================================\n")

source(file.path(dir_config$module, "module11_diff_analysis2.R"))

# 对 Module 10 的数据替换结果进行第二次差异分析
# 只选择 Experiment 和 Spatial 样本，构建 Exp vs Exp 和 Exp vs Spatial 的比较
result <- module11_diff_analysis2(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  replaced_data_A = replaced_data_A,
  replaced_data_B = replaced_data_B,
  selected_versions = selected_versions_bg
)

# 将返回的数据赋值给全局变量
FDR_combined_df_list_2nd <- result$FDR_combined_df_list_2nd
sample_info_2nd <- result$sample_info
comparisons_list_2nd <- result$comparisons_list
bioGroups_selected_2nd <- result$bioGroups_selected
bioGroup_info_2nd <- result$bioGroup_info

# 保存Module11环境（包含Module01-10所有数据 + Module11数据）
save.image(file = "Module11_workspace.RData")
cat("✓ 已保存: Module11_workspace.RData (工作目录)\n")
cat("  包含：Module01-10所有数据 + FDR_combined_df_list_2nd等\n")

cat("✓ Module 11 完成\n")

# ============================================================================
# Module 12: 第二次ROC分析
# ============================================================================
cat("\n========================================\n")
cat("Module 12: 第二次ROC分析\n")
cat("========================================\n")

source(file.path(dir_config$module, "module12_second_roc.R"))

# 准备表达矩阵列表（Method A/B 合并，并对名称进行清理）
expr_data_list_2nd <- list()
if (length(replaced_data_A) > 0) {
  expr_data_list_2nd <- c(
    expr_data_list_2nd,
    setNames(replaced_data_A, paste0(names(replaced_data_A), "_A"))
  )
}
if (length(replaced_data_B) > 0) {
  expr_data_list_2nd <- c(
    expr_data_list_2nd,
    setNames(replaced_data_B, paste0(names(replaced_data_B), "_B"))
  )
}

name_cleanup <- function(x) stringr::str_replace(x, "_Imputed", "_New")
names(expr_data_list_2nd) <- name_cleanup(names(expr_data_list_2nd))
names(FDR_combined_df_list_2nd) <- name_cleanup(names(FDR_combined_df_list_2nd))

names(FDR_combined_df_list_2nd) # 查看可选版本

# 配置第二次ROC分析接口（无SubMito）
# second_roc_versions: 指定要进入第二次ROC的版本（名称需与FDR列表匹配）
#   例如 c("noNorm_New_A","Local_QNorm_New_A","noNorm_New_B","Local_QNorm_New_B")
#   设为 NULL 表示自动包含所有版本
second_roc_versions <- c(
  "Local_QNorm_New_A"
)
second_roc_annotation_column <- "GO_Localization"
second_roc_tp_label <- "SGs"
second_roc_fp_label <- "Cytosol"
second_roc_min_tp <- 0.3

second_roc_order <- second_roc_versions
if (is.null(second_roc_order)) {
  second_roc_order <- names(FDR_combined_df_list_2nd)
}
second_roc_order <- second_roc_order[second_roc_order %in% names(FDR_combined_df_list_2nd)]

result <- module12_second_roc(
  dir_config = dir_config,
  FDR_combined_df_list_2nd = FDR_combined_df_list_2nd,
  expr_data_list = expr_data_list_2nd,
  comparisons_list = if (exists("comparisons_list_2nd")) comparisons_list_2nd else NULL,
  annotation_column = second_roc_annotation_column,
  tp_label = second_roc_tp_label,
  fp_label = second_roc_fp_label,
  desired_order = if (length(second_roc_order) > 0) second_roc_order else names(FDR_combined_df_list_2nd),
  min_tp = second_roc_min_tp,
  sample_columns = if (exists("sample_info_2nd")) sample_info_2nd$SampleName else NULL,
  youden_ylim = c(0, 0.8)
)

FDR_combined_df_list_2nd <- result$FDR_combined_df_list_2nd
Expr_FDR_df_list_2nd <- result$Expr_FDR_df_list_2nd
all_desired_thresholds_2nd <- result$all_desired_thresholds_2nd

save.image(file = "Module12_workspace.RData")
cat("✓ 已保存: Module12_workspace.RData (工作目录)\n")
cat("  包含：Module01-11所有数据 + second ROC结果\n")

cat("✓ Module 12 完成\n")

# ============================================================================
# Module 13: Sub注释（SG / Nucleolus）
# ============================================================================
cat("\n========================================\n")
cat("Module 13: Sub注释\n")
cat("========================================\n")

source(file.path(dir_config$module, "module13_subsg_annotation.R"))

# 配置 Sub注释模式与列（subSG 默认写入MultiBait，subNucleolus写入Sub_HPA）
subsg_mode <- "subSG"  # 可选："subSG" 或 "subNucleolus"
subsg_target_column <- if (subsg_mode == "subNucleolus") "Sub_HPA_Localization" else "MultiBait_Localization"
subsg_levels_order <- NULL  # 使用默认顺序，如需自定义可赋值为字符向量
names(Expr_FDR_df_list_2nd)
subsg_forstep16_versions <- c(
  "Local_QNorm_New_A"
)

result <- module13_subsg_annotation(
  dir_config = dir_config,
  expr_fdr_df_list_2nd = Expr_FDR_df_list_2nd,
  annotation_references = annotation_references,
  mode = subsg_mode,
  target_column = subsg_target_column,
  levels_order = subsg_levels_order,
  forstep16_versions = subsg_forstep16_versions
)

Expr_FDR_df_list_2nd_SubSGs <- result$Expr_FDR_df_list_2nd_SubSGs
ForStep16 <- result$ForStep16
subsg_output_files <- result$output_files
names(Expr_FDR_df_list_2nd_SubSGs)
save.image(file = "Module13_workspace.RData")
cat("✓ 已保存: Module13_workspace.RData (工作目录)\n")
cat("  包含：Module01-12所有数据 + Expr_FDR_df_list_2nd_SubSGs\n")

cat("✓ Module 13 完成\n")
load("Module13_workspace.RData")
# ============================================================================
# Module 14: 火山图
# ============================================================================
cat("\n========================================\n")
cat("Module 14: 火山图\n")
cat("========================================\n")

source(file.path(dir_config$module, "module14_volcano_plots.R"))

if (!exists("Expr_FDR_df_list_2nd_SubSGs")) {
  stop("✗ 错误：未检测到 Expr_FDR_df_list_2nd_SubSGs，请先完成 Module 13")
}
if (length(Expr_FDR_df_list_2nd_SubSGs) == 0) {
  stop("✗ 错误：Expr_FDR_df_list_2nd_SubSGs 为空，无法绘制火山图")
}
cat("✓ 检测到 Expr_FDR_df_list_2nd_SubSGs，可继续配置 Module 14 火山图参数\n")

module14_comparison_metadata <- NULL
if (exists("comparisons_list_2nd") && exists("bioGroup_info_2nd")) {
  module14_comparison_metadata <- module14_build_comparison_metadata(
    comparisons_list = comparisons_list_2nd,
    bioGroup_info = bioGroup_info_2nd
  )
}
if (is.null(module14_comparison_metadata)) {
  cat("⚠ 提示：未能获取比较分类信息（comparisons_list_2nd 或 bioGroup_info_2nd 缺失）\n")
} else {
  cat("✓ 已根据 Module11 信息构建比较分类（Exp_vs_Exp / Exp_vs_Spatial）\n")
}

# 配置火山图接口
# - volcano_versions: 需进入火山图的版本（NULL = 使用全部）
# - volcano_annotation_column: 选择用于着色的注释列
# - volcano_annotation_color_map: 注释 -> 颜色映射，可按需增删
# - volcano_use_threshold_colors: TRUE 时按阈值着色，忽略注释颜色
# - volcano_logfc_threshold / volcano_fdr_threshold: 阈值
# - volcano_label_mode: c("with","without") 表示同时输出带/不带标签 PDF
# - volcano_label_annotations: 哪些注释类别允许显示标签
# - volcano_comparison_sets: 定义 logFC & FDR 列组合
volcano_versions <- c("Local_QNorm_New_A") # 可选：names(Expr_FDR_df_list_2nd_SubSGs)，NULL=全部
volcano_annotation_column <- "HaloMap_Localization" # 可选：上方打印的 *_Localization 列
volcano_annotation_color_map <- c(
  "SGs" = "red",
  "Mitochondrion" = "black",
  "Nuclear" = "blue"
) # 未列出的注释会自动回退为灰色
volcano_use_threshold_colors <- FALSE # TRUE=按阈值着色并覆盖注释颜色
volcano_logfc_threshold <- 0.5       # 阈值逻辑 |logFC| ≥ 该值
volcano_fdr_threshold <- 0.05        # 阈值逻辑 FDR ≤ 该值
volcano_threshold_color_above <- "red"   # 阈值着色：满足条件
volcano_threshold_color_below <- "grey70" # 阈值着色：不满足条件
volcano_label_mode <- c("without") # 允许 "with" 或 "without"，可同时填写
volcano_label_annotations <- c("SGs", "Mitochondrion") # 仅这些注释允许打标签
volcano_label_size <- 3
volcano_label_max_overlaps <- 15
# 可选：c("Exp_vs_Exp","Exp_vs_Spatial")，或单独选择某一类，NULL=全部
volcano_comparison_categories <- c("Exp_vs_Spatial")
volcano_xlim_override <- c(-3, 6)
volcano_ylim_override <- c(0, 15)
# 自定义示例（如需替换默认值）：
# volcano_comparison_sets <- list(
#   list(
#     name = "Custom1",
#     fc_cols = c("K69A1B3_vs_K69C3_logFC", "K69A1B3_vs_C3_logFC"),
#     fdr_cols = c("K69A1B3_vs_K69C3_adj.P.Val", "K69A1B3_vs_C3_adj.P.Val"),
#     xlim = c(-2, 6),
#     ylim = c(0, 13)
#   ),
#   list(
#     name = "Custom2",
#     fc_cols = c("K69C3_vs_K20_logFC"),
#     fdr_cols = c("K69C3_vs_K20_adj.P.Val")
#   )
# )
volcano_comparison_sets <- module14_default_comparison_sets(module14_comparison_metadata)

module14_versions <- if (is.null(volcano_versions)) NULL else volcano_versions

volcano_result <- module14_volcano_plots(
  dir_config = dir_config,
  expr_fdr_df_list_2nd_subsgs = Expr_FDR_df_list_2nd_SubSGs,
  versions = module14_versions,
  annotation_column = volcano_annotation_column,
  annotation_color_map = volcano_annotation_color_map,
  use_threshold_colors = volcano_use_threshold_colors,
  logfc_threshold = volcano_logfc_threshold,
  fdr_threshold = volcano_fdr_threshold,
  threshold_color_above = volcano_threshold_color_above,
  threshold_color_below = volcano_threshold_color_below,
  label_mode = volcano_label_mode,
  label_annotations = volcano_label_annotations,
  label_size = volcano_label_size,
  label_max_overlaps = volcano_label_max_overlaps,
  comparison_sets = volcano_comparison_sets,
  xlim_override = volcano_xlim_override,
  ylim_override = volcano_ylim_override,
  comparison_metadata = module14_comparison_metadata,
  comparison_categories = volcano_comparison_categories
)

save.image(file = "Module14_workspace.RData")
cat("✓ 已保存: Module14_workspace.RData (工作目录)\n")
cat("  包含：Module01-13所有数据 + 火山图 PDF\n")

cat("✓ Module 14 完成\n")

# ============================================================================
cat("\n========================================\n")
cat("所有模块执行完成！\n")
cat("========================================\n")
load("Module14_workspace.RData")
forstep16_r_object <- file.path(dir_config$output, "Module14_ForStep16.rds")
saveRDS(ForStep16, forstep16_r_object)
cat(sprintf("✓ 已导出 ForStep16 R对象: %s\n", forstep16_r_object))

forstep16_excel <- file.path(dir_config$output, "Module14_ForStep16.xlsx")
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  warning("⚠ 提示：缺少 openxlsx 包，无法导出 ForStep16 Excel，请先安装 openxlsx")
} else {
  wb_forstep16 <- openxlsx::createWorkbook()
  existing_names <- character()
  for (sheet_name in names(ForStep16)) {
    df <- ForStep16[[sheet_name]]
    if (!is.data.frame(df)) {
      df <- as.data.frame(df)
    }
    sanitized <- gsub("[^A-Za-z0-9]+", "_", sheet_name)
    if (nchar(sanitized) == 0) sanitized <- "Sheet"
    base_name <- substr(sanitized, 1, 31)
    unique_name <- base_name
    suffix <- 1
    while (unique_name %in% existing_names) {
      suffix_label <- paste0("_", suffix)
      max_base_len <- max(31 - nchar(suffix_label), 1)
      unique_name <- paste0(substr(base_name, 1, max_base_len), suffix_label)
      suffix <- suffix + 1
    }
    existing_names <- c(existing_names, unique_name)
    openxlsx::addWorksheet(wb_forstep16, unique_name)
    openxlsx::writeData(wb_forstep16, unique_name, df)
  }
  openxlsx::saveWorkbook(wb_forstep16, forstep16_excel, overwrite = TRUE)
  cat(sprintf("✓ 已导出 ForStep16 多sheet Excel: %s\n", forstep16_excel))
}
# 接下来可以运行Base_model.R 脚本来构建1D/2D过滤模型以及输出Final蛋白列表

