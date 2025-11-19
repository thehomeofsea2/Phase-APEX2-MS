# Pipeline 模块接口详解文档

本文档详细解读 Pipeline 中 14 个模块的函数接口设计。

重点说明：函数的输入参数（Arguments）、期望的数据结构以及输出对象。这将帮助您快速理解如何调用、配置以及调试每个步骤。

## Phase I: 背景扣除与数据预处理 (Background Subtraction)

### Module 01: 环境初始化 (Setup)

* **文件** : `module01_setup.R`
* **主函数** : `module01_setup(dir_config)`
* **功能** : 检查必需目录，创建输出目录，保存运行环境元数据。

| **参数名** | **类型** | **默认值** | **说明**                                                                               |
| ---------------- | -------------- | ---------------- | -------------------------------------------------------------------------------------------- |
| `dir_config`   | List           | (必需)           | 包含 `root`,`rawdata`,`reference`,`output`等路径的列表。由 `main_pipeline.R`定义。 |

* **输出 (Return)** :
* `config` (List): 包含运行时间戳、R版本、系统平台信息。
* **副作用** : 自动创建 `Output/`, `Module/`, `Dev/` 文件夹。

### Module 02: 数据导入 (Data Import)

* **文件** : `module02_data_import.R`
* **主函数** : `module02_data_import(dir_config, file_pattern, auto_process)`
* **功能** : 读取原始 TSV，并采用交互式或自动模式生成实验分组表 (`sampleGroup`)。

| **参数名** | **类型** | **默认值**       | **说明**                                                                                               |
| ---------------- | -------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------ |
| `dir_config`   | List           | (必需)                 | 目录配置。                                                                                                   |
| `file_pattern` | String         | `"_matrix.*\\.tsv$"` | 用于匹配 Rawdata 目录中文件名的正则表达式。                                                                  |
| `auto_process` | Bool           | `FALSE`              | `FALSE`: 生成模板后暂停，等待用户填写 CSV。`TRUE`: 假设 CSV 已填写好，直接读取并处理（用于自动化脚本）。 |

* **输出 (Return)** :
* List:
  * `data`: 清洗列名后的原始表达矩阵 (DataFrame)。
  * `sampleGroup`: 整理并验证后的样本分组信息 (DataFrame)。
* **文件输出** : `Module02_sampleGroup_template.csv`, `Module02_sampleGroup.csv`, `Output/Module02_data_raw.csv`。

### Module 03: 注释系统 (Annotation)

* **文件** : `module03_annotation.R`
* **主函数** : `module03_annotation(...)`
* **功能** : 加载 HPA, MitoCarta 等数据库，对蛋白进行亚细胞定位注释。

| **参数名**           | **类型** | **默认值** | **说明**                                                                |
| -------------------------- | -------------- | ---------------- | ----------------------------------------------------------------------------- |
| `dir_config`             | List           | (必需)           | 目录配置。                                                                    |
| `data_raw`               | DataFrame      | (必需)           | 原始数据。                                                                    |
| `sampleGroup`            | DataFrame      | `NULL`         | 样本信息（可选）。                                                            |
| `custom_annotations`     | List           | `NULL`         | **高级接口** 。完全覆盖默认注释配置。                                   |
| `annotation_mode`        | String         | `"SG"`         | `"SG"`(应激颗粒模式) 或 `"Nucleolus"`(核仁模式)。决定默认加载哪些参考集。 |
| `custom_tp_sources`      | List           | `NULL`         | 用于加载用户自定义的外部参考文件路径和读取方式。                              |
| `additional_annotations` | List           | `NULL`         | 在默认配置基础上**追加**新的注释列配置。                                |

* **输出 (Return)** :
* List:
  * `data_annotated`: 带注释列的数据框。
  * `annotation_references`: 包含所有加载的参考数据库（HPA, MitoCarta 等）的大列表。
* **文件输出** : `Output/Module03_data_annotated.csv`。

### Module 04: 数据标准化 (Standardization)

* **文件** : `module04_standardization.R`
* **主函数** : `module04_standardization(...)`
* **功能** : 执行 Log2 转化及多种归一化策略。

| **参数名**   | **类型** | **默认值**                   | **说明**                                                |
| ------------------ | -------------- | ---------------------------------- | ------------------------------------------------------------- |
| `dir_config`     | List           | (必需)                             | -                                                             |
| `data_annotated` | DataFrame      | (必需)                             | 来自 Module 03 的数据。                                       |
| `sampleGroup`    | DataFrame      | (必需)                             | 需包含 `bioGroup`列用于局部标准化。                         |
| `norm_types`     | Vector         | `c("noNorm", "Global_QNorm"...)` | 指定要执行的标准化方法。**必须包含至少一种全局方法** 。 |

* **输出 (Return)** :
* List:
  * `standardized_data_list`: 列表，键为标准化方法名，值为 DataFrame。
  * `norm_types_used`: 实际执行的方法名列表。
* **文件输出** : 各版本的 CSV 文件及由 Boxplot 组成的 PDF 质控图。

### Module 05: 缺失值填补 (Imputation)

* **文件** : `module05_imputation.R`
* **主函数** : `module05_imputation(...)`
* **功能** : 基于 Perseus 策略（左偏正态分布）填补缺失值。

| **参数名**           | **类型** | **默认值** | **说明**                                                        |
| -------------------------- | -------------- | ---------------- | --------------------------------------------------------------------- |
| `dir_config`             | List           | (必需)           | -                                                                     |
| `standardized_data_list` | List           | (必需)           | 来自 Module 04 的列表。                                               |
| `sampleGroup`            | DataFrame      | (必需)           | 需包含 `CatalyticGroup`列。                                         |
| `impute_cat_mean`        | Bool           | `FALSE`        | 针对 `Cat`(阳性) 组，如果 n_valid=2，是否用均值填补（默认不填补）。 |
| `random_seed`            | Integer        | `123`          | 随机种子，确保填补结果可复现。                                        |

* **关键逻辑** : `NoCat` 组（Context=Control）会强制填补，模拟背景噪音。
* **输出 (Return)** :
* List:
  * `imputed_data_list`: 列表，键名添加 `_Imputed` 后缀。
  * `imputation_params`: 记录填补使用的参数。
* **文件输出** : 填补前后分布对比的 Boxplot PDF。

### Module 06: 热图质控 (Heatmap)

* **文件** : `module06_heatmap.R`
* **主函数** : `module06_heatmap(...)`
* **功能** : 绘制样本相关性热图和聚类热图。

| **参数名**         | **类型** | **默认值**                | **说明**                                                               |
| ------------------------ | -------------- | ------------------------------- | ---------------------------------------------------------------------------- |
| `imputed_data_list`    | List           | (必需)                          | -                                                                            |
| `sampleGroup`          | DataFrame      | (必需)                          | -                                                                            |
| `selected_versions`    | Vector         | `NULL`                        | 指定用哪个版本的数据画图 (默认自动选择第一个)。                              |
| `heatmap_types`        | Vector         | `c("all", "by_localization")` | 绘图类型：`correlation`,`all`,`by_localization`。                      |
| `correlation_config`   | List           | `NULL`                        | **重要** 。配置相关性热图的细节，如 `exclude_context`,`corr_min`。 |
| `localization_columns` | Vector         | `NULL`                        | 指定用于分类绘图的注释列名（NULL则自动检测）。                               |
| `color_params`         | List           | `NULL`                        | 自定义热图颜色范围 (`custom_min`,`custom_max`)。                         |

* **输出 (Return)** :
* `heatmap_info` (List): 包含绘图参数信息的列表。
* **文件输出** : 相关性热图 PDF, 全局/分定位聚类热图 PDF。

### Module 07: 差异分析 I (Diff Analysis - Background)

* **文件** : `module07_diff_analysis1.R`
* **主函数** : `module07_diff_analysis1(...)`
* **功能** : 基于 `FirstROCgroup` 进行去背景差异分析 (Exp vs Control)。

| **参数名**      | **类型** | **默认值** | **说明**                           |
| --------------------- | -------------- | ---------------- | ---------------------------------------- |
| `imputed_data_list` | List           | (必需)           | -                                        |
| `sampleGroup`       | DataFrame      | (必需)           | 需包含 `FirstROCgroup`和 `Context`。 |
| `selected_versions` | Vector         | `NULL`         | 指定分析的数据版本（NULL则分析全部）。   |

* **核心逻辑** : 读取 `sampleGroup$FirstROCgroup`。同组内的 `Context=Experiment` 样本 vs `Context=Control` 样本。
* **输出 (Return)** :
* List:
  * `diff_results1`: 核心结果列表。包含 `combined` (大宽表), `expr_matrix` (表达矩阵), `comparisons` (对比组信息)。
  * `comparisons_used`: 实际执行的对比组列表。
  * `comparison_info`: 详细的对比组元数据。
* **文件输出** : `_DiffAnalysis.csv` (合并结果), `_DiffAnalysis_Full.xlsx` (含 topTable 原始结果).

### Module 08: ROC 分析 I (ROC Analysis - Background)

* **文件** : `module08_roc_analysis1.R`
* **主函数** : `module08_roc_analysis1(...)`
* **功能** : 计算去背景的最佳 LogFC 阈值。

| **参数名**          | **类型** | **默认值**      | **说明**                                                   |
| ------------------------- | -------------- | --------------------- | ---------------------------------------------------------------- |
| `diff_results1`         | List           | (必需)                | 来自 Module 07。                                                 |
| `annotation_references` | List           | (必需)                | 来自 Module 03，需包含 MitoCarta。                               |
| `selected_versions`     | Vector         | `NULL`              | -                                                                |
| `roc_annotation_column` | String         | `"GO_Localization"` | ROC 的“真值”参考列。                                           |
| `tp_label`              | String         | `"SGs"`             | 真阳性标签 (Target)。                                            |
| `fp_label`              | String         | `"Matrix"`          | 假阳性标签 (Background)，通常由 SubMito 转化而来。               |
| `enable_submito`        | Bool           | `TRUE`              | 是否启用 SubMito 转化（将 Mitochondrion 拆分为 Matrix/IMS 等）。 |
| `min_tp`                | Numeric        | `0.3`               | 计算阈值时要求的最小 TP 比例 (Sensitivity)。                     |

* **输出 (Return)** :
* List:
  * `roc_results`: 包含 ROC 对象和数据框。
  * `all_thresholds`: 计算出的最佳 LogFC 阈值列表。
  * `data_with_submito`: 经过 SubMito 转化后的数据。
  * `expr_fdr_df_list`: 集成表达量、FDR 和原始注释的中间数据（供 Module 09 使用）。
* **文件输出** : ROC 曲线 PDF, Youden Index PDF, 阈值 Excel。

### Module 09: 背景扣除 (Background Subtraction)

* **文件** : `module09_background_subtraction.R`
* **主函数** : `module09_background_subtraction(...)`
* **功能** : 执行最终过滤，生成“真实”互作蛋白列表。

| **参数名**                    | **类型** | **默认值**        | **说明**                                                                        |
| ----------------------------------- | -------------- | ----------------------- | ------------------------------------------------------------------------------------- |
| `diff_results`,`roc_thresholds` | List           | (必需)                  | -                                                                                     |
| `fdr_threshold`                   | Numeric        | `0.05`                | 通用 FDR 阈值。                                                                       |
| `no_fdr_comparisons`              | Vector         | `NULL`                | **Method A 核心参数** 。指定哪些对比组**不**使用 FDR 过滤（仅用 LogFC）。 |
| `use_roc_threshold`               | Bool           | `TRUE`                | 是否使用 Module 08 计算的动态阈值。                                                   |
| `fixed_fc_threshold`              | Numeric        | `NULL`                | 若 `use_roc_threshold=FALSE`，则使用此固定值。                                      |
| `min_valid_lfq`                   | Integer        | `2`                   | 催化组最小有效值数量要求。                                                            |
| `tp_label`,`tp_color`           | String         | `"SGs"`,`"#DB6968"` | 用于绘图的 TP 标签和颜色。                                                            |

* **输出 (Return)** :
* List:
  * `filtered_data_A/B`: 过滤后的列表 (按 BioGroup 分表)。
  * `merged_data_A/B`: 过滤后的列表 (合并去重，用于 Module 10)。
* **文件输出** : `FilteredData.xlsx`, 堆积柱状图 PDF (展示定位比例)。

## Phase II: 生物学探索 (Biological Discovery)

### Module 10: 数据替换 (Data Replacement)

* **文件** : `module10_data_replacement.R`
* **主函数** : `module10_data_replacement(...)`
* **功能** : 将筛选出的蛋白回填定量数值，为组间比较做准备。

| **参数名**      | **类型** | **默认值**           | **说明**                  |
| --------------------- | -------------- | -------------------------- | ------------------------------- |
| `merged_data_A/B`   | List           | (必需)                     | 来自 Module 09 的“掩码”列表。 |
| `imputed_data_list` | List           | (必需)                     | 来自 Module 05 的“数值源”。   |
| `selected_versions` | Vector         | `NULL`                   | -                               |
| `prefer_version`    | String         | `"Global_QNorm_Imputed"` | 优先使用的数值版本。            |
| `fallback_version`  | String         | `"Global_MNorm_Imputed"` | 备选数值版本。                  |
| `test_n`            | Integer        | `20`                     | 抽样检查的基因数量。            |

* **输出 (Return)** :
* List:
  * `replaced_data_A`, `replaced_data_B`: 完成数值回填的矩阵列表。
  * `base_version_used`: 实际使用的基准版本。
* **文件输出** : `Output/Check/` 下的校验报告 CSV/TXT, `_AfterFirstROC_Intersect.xlsx`.

### Module 11: 差异分析 II (Diff Analysis - Biological)

* **文件** : `module11_diff_analysis2.R`
* **主函数** : `module11_diff_analysis2(...)`
* **功能** : 比较 Experiment 组之间或 Experiment vs Spatial 的差异。

| **参数名**      | **类型** | **默认值** | **说明**                             |
| --------------------- | -------------- | ---------------- | ------------------------------------------ |
| `replaced_data_A/B` | List           | (必需)           | 来自 Module 10。                           |
| `sampleGroup`       | DataFrame      | (必需)           | 用于解析 `SecondROCgroup`和 `PLtype`。 |
| `selected_versions` | Vector         | `NULL`         | -                                          |

* **核心逻辑** :
* 仅保留 `Context %in% c("Experiment", "Spatial")`。
* 自动构建 `Exp vs Exp` (组间) 和 `Exp vs Spatial` (组与参考) 的比较对。
* **输出 (Return)** :
* List:
  * `FDR_combined_df_list_2nd`: 包含第二轮 LogFC 和 FDR 的汇总表。
  * `comparisons_list`: 生成的对比组列表。
  * `bioGroup_info`: 参与分析的样本组信息。
* **文件输出** : `Module11_FC_FDR_2nd_combined.xlsx` (核心结果), `_Raw_FDR_test_list_*.xlsx` (原始数据).

### Module 12: ROC 分析 II (ROC Analysis - Biological)

* **文件** : `module12_second_roc.R`
* **主函数** : `module12_second_roc(...)`
* **功能** : 针对特定的生物学比较计算富集阈值。

| **参数名**             | **类型** | **默认值**      | **说明**                         |
| ---------------------------- | -------------- | --------------------- | -------------------------------------- |
| `FDR_combined_df_list_2nd` | List           | (必需)                | 来自 Module 11。                       |
| `expr_data_list`           | List           | (必需)                | 表达矩阵列表（通常即 replaced_data）。 |
| `annotation_column`        | String         | `"GO_Localization"` | ROC 真值列。                           |
| `tp_label`                 | String         | `"SGs"`             | 真阳性 (如应激颗粒)。                  |
| `fp_label`                 | String         | `"Cytosol"`         | 假阳性 (如胞浆)。                      |
| `min_tp`                   | Numeric        | `0.3`               | 最小 TP 比例。                         |
| `youden_ylim`              | Vector         | `c(0, 0.7)`         | Youden Index 图的 Y 轴范围。           |

* **输出 (Return)** :
* List:
  * `all_desired_thresholds_2nd`: 第二轮最佳 LogFC 阈值。
  * `Expr_FDR_df_list_2nd`: 集成了表达量、FDR 和注释的完整数据对象。
* **文件输出** : ROC 曲线 PDF, Youden Index PDF, 阈值 Excel, 数据 Excel。

### Module 13: 精细注释 (Sub-Annotation)

* **文件** : `module13_subsg_annotation.R`
* **主函数** : `module13_subsg_annotation(...)`
* **功能** : 利用逻辑运算生成复杂的亚细胞定位标签 (如 "Nuclear&SGs")。

| **参数名**          | **类型** | **默认值** | **说明**                                                          |
| ------------------------- | -------------- | ---------------- | ----------------------------------------------------------------------- |
| `expr_fdr_df_list_2nd`  | List           | (必需)           | 来自 Module 12（或 11）。                                               |
| `annotation_references` | List           | (必需)           | 需包含 HPA/HaloMap 参考。                                               |
| `mode`                  | String         | `"subSG"`      | `"subSG"`(SGs分析) 或 `"subNucleolus"`(核仁分析)。                  |
| `target_column`         | String         | `NULL`         | 输出列名 (如 `"MultiBait_Localization"`)。NULL 则自动根据 mode 决定。 |
| `levels_order`          | Vector         | `NULL`         | 指定输出 Factor 的 Level 顺序 (影响绘图颜色顺序)。                      |

* **输出 (Return)** :
* List:
  * `Expr_FDR_df_list_2nd_SubSGs`: 最终用于绘图的完整数据集 (含复杂注释)。
  * `ForStep16`: 供后续步骤使用的数据子集。
* **文件输出** : `Module13_*.csv`。

### Module 14: 火山图 (Volcano Plots)

* **文件** : `module14_volcano_plots.R`
* **主函数** : `module14_volcano_plots(...)`
* **功能** : 生成最终的可视化结果。

| **参数名**                | **类型** | **默认值**              | **说明**                                                                              |
| ------------------------------- | -------------- | ----------------------------- | ------------------------------------------------------------------------------------------- |
| `expr_fdr_df_list_2nd_subsgs` | List           | (必需)                        | 来自 Module 13。                                                                            |
| `comparison_sets`             | List           | `NULL`                      | **核心配置** 。定义哪些对比画在一起，以及对应的 xlim/ylim。若为 NULL 则使用默认配置。 |
| `annotation_column`           | String         | `"GO_Localization"`         | 用于着色的列 (通常是 Module 13 生成的列)。                                                  |
| `annotation_color_map`        | Vector         | (默认红/蓝/黑)                | 颜色映射表。                                                                                |
| `use_threshold_colors`        | Bool           | `FALSE`                     | 是否忽略注释，仅按是否过阈值着色 (红/灰)。                                                  |
| `logfc_threshold`             | Numeric        | `0.5`                       | 划线阈值。                                                                                  |
| `fdr_threshold`               | Numeric        | `0.05`                      | 划线阈值。                                                                                  |
| `label_mode`                  | Vector         | `c("with",`"without")       | 同时生成带标签和不带标签的 PDF。                                                            |
| `label_annotations`           | Vector         | `c("SGs", "Mitochondrion")` | 指定哪些注释类别的蛋白需要显示基因名标签。                                                  |

* **输出 (Return)** :
* List: 包含生成的 PDF 文件路径列表 (`pdf_files`) 和源数据文件路径 (`source_data_file`)。
* **文件输出** : 多个 PDF 火山图, `_VolcanoSource.xlsx` (绘图源数据)。



# 亚细胞组分多维模型筛选与验证系统 (Step29) 项目指南

## 1. 项目概述与设计框架

### 1.1 项目背景

本项目模块（Step29）紧接 `main_pipeline.R` 的 Module14，输入数据为 `ForStep16`。其核心目标是突破传统的单维（仅 logFC 或 仅丰度）筛选限制，通过整合 **质谱定量数据（丰度/差异）** 与 **蛋白互作网络（PPI）** 信息，构建多维广义线性模型（GLM），以识别具有特定亚细胞定位属性（如核仁、应激颗粒）的候选蛋白。

### 1.2 核心设计框架

脚本采用 **三阶段 (Three-Phase)** 设计模式，实现了从计算到可视化再到评估的完整闭环：

* **Phase 1: 特征工程与模型构建 (Data Calculation & Modeling)**
  * 数据预处理：整合本地 STRING 数据库，计算全局和局部 PPI 分数。
  * 模型训练：基于已知定位（Ground Truth）训练 Logistic Regression 模型。
  * 打分与预测：利用优化的 Youden Index 阈值对全蛋白组进行打分和分类。
  * *输出：RData 对象与 Cytoscape 网络文件。*
* **Phase 2: 模型评估与可视化 (Visualization)**
  * 性能评估：ROC 曲线、Youden 指数图。
  * 特征空间可视化：决策边界图、3D 散点图。
  * 网络密度分析：PPI 热图、联合概率（JP）矩阵分析（核心创新点，用于评估蛋白在网络中的聚集紧密度）。
* **Phase 3: 候选比较与最终筛选 (Comparison & Selection)**
  * 方法对标：将机器学习模型结果与传统筛选（logFC/FDR）及空间对照组进行对比。
  * 成分分析：分析候选列表的亚细胞定位组成（MultiBait Analysis）。
  * 富集度分析：TP（真阳性）百分比曲线，评估排序靠前的候选蛋白的准确性。

---

## 2. 统计思想与数学原理

### 2.1 贝叶斯先验思想的融合

系统不仅仅依赖实验数据的后验分布（logFC），还引入了“先验知识”（PPI 网络）。假设： **定位在同一无膜细胞器中的蛋白倾向于在 PPI 网络中彼此互作** 。

### 2.2 广义线性模型 (GLM)

使用 Logistic Regression（逻辑回归）来预测蛋白属于目标定位（TP）的概率 $P(Y=1|X)$：

$$
\ln(\frac{P}{1-P}) = \beta_0 + \beta_1 X_{abundance} + \beta_2 X_{logFC} + \beta_3 X_{PPI}
$$

其中，$X$ 为输入特征，$\beta$ 为模型学习到的权重。

### 2.3 概率质量与空间密度 (Spatial Density)

在 Phase 2 中，引入了基于 Kernel Density Estimation (KDE) 的互作密度分析。

* **思想** ：如果一组候选蛋白是真实的复合物组分，它们在 PPI 矩阵中应呈现“聚集”状态，而非随机分布。
* **实现** ：通过 Joint Probability (JP) 和 Probability Mass (PM) 矩阵，计算蛋白对在排序空间中的互作概率密度。

---

## 3. 接口配置说明 (Configuration)

在运行脚本前，需检查以下关键接口：

### 3.1 输入与分组 (`group_info`)

定义实验组结构，包含重复样本列名、logFC 列和 FDR 列。

**R**

```
group_info <- list(
  E7A2B4 = list(samples = c("..."), logFC = "...", logFC_FDR="..."),
  ...
)
```

### 3.2 真阳性定义 (`myTP_vector` & `myTP_vector2`)

定义哪些亚细胞定位标签被视为“金标准”（Ground Truth），用于训练模型。

* **当前配置** ：侧重于核仁 (Nucleolus) 相关组分。
* **注意** ：`myTP_vector2` 专门用于 Phase 3 的 TPR 曲线计算，逻辑更加宽泛（包含所有含 "Nucleolus" 的字符）。

### 3.3 筛选参数 (`PARAM_*`)

* `PARAM_MIN_LFQ_NONNA`: 建模前要求定量值的最小非NA样本数（默认2）。
* `PARAM_MODELS_TO_PLOT`: Phase 3 中需要展示和对比的模型列表。
* `PARAM_LOCALIZATION_COLORS`: 绘图颜色映射，必须覆盖数据中出现的所有定位类别。

---

## 4. 模块详细设计

### Phase 1: 数据计算与导出

1. **PPI 预处理** ：

* 加载本地 STRING 数据库（`9606.protein.links`）。
* **Local PPI Score** ：计算每个蛋白与其在丰度排序上邻近蛋白（Window Size=10）的互作强度。这是为了捕捉在梯度上共迁移的互作蛋白。
* **Global PPI Score** ：计算蛋白在整个网络中的总互作强度。

1. **模型训练循环** ：

* 针对每个实验组，构建多种模型（1D, 2D, 3D）。
* 使用 `pROC` 包计算 AUC，并利用 Youden Index (`sensitivity + specificity - 1`) 确定最佳分类阈值。
* 输出：`step29_Analysis_Phase1_DataExport` 目录下的 RData 和 Cytoscape CSV 文件。

### Phase 2: 可视化系统

此阶段生成一系列 PDF 用于评估模型：

1. **ROC/Youden (文件夹 1 & 2)** ：评估模型的分类效能。
2. **决策边界 (文件夹 3)** ：

* 在二维坐标系（如 Abundance vs logFC）上画出模型的分割线。
* **意义** ：直观展示模型是如何“切割”数据的（例如，是否同时要求高丰度和高变化倍数）。

1. **3D 散点图 (文件夹 4)** ：

* 使用颜色代表第三维（PPI Score），展示数据分布。

1. **PPI 热图与密度 (文件夹 5, 7, 13, 15)** ：

* **Candidate PPI Heatmaps** ：仅展示通过筛选的候选蛋白之间的互作。
* **JP/PM Heatmaps (高级)** ：将蛋白按丰度排序映射到 **$N \times N$** 网格。颜色越红/亮，表示该区域（通常是左下角，即高丰度区域）的互作密度越高。这是验证无膜细胞器相分离特性的强力证据。

### Phase 3: 比较与决策

此阶段用于回答“哪个模型最好？”和“筛选出了什么？”：

1. **详细列表导出** ：生成包含所有模型候选蛋白的 Excel 表。
2. **成分堆积图 (Plot_Count/Abundance)** ：

* 展示不同方法筛选出的候选蛋白中，TP（目标定位）与其他定位的比例。
* **意义** ：如果某个模型筛选出的 TP 比例显著高于 Control 或传统方法，说明模型特异性好。

1. **MultiBait 分析** ：

* 针对多诱饵实验设计的细化分析，展示更精细的定位分类（如 Nuclear&Nucleolus vs Cytosol&Nucleolus）。

1. **TP 百分比曲线 (图9)** ：

* X轴为排名百分比，Y轴为 TP 占比。
* **意义** ：曲线越靠左上，说明模型能将真正的目标蛋白排在列表的最前面。这是评估排序质量的核心指标。

---

## 5. 图表意义速查表

| **图表类型**   | **文件路径关键词**  | **核心意义**       | **理想形态**                                          |
| -------------------- | ------------------------- | ------------------------ | ----------------------------------------------------------- |
| **ROC 曲线**   | `1_ROC_Curves`          | 模型区分 TP/Other 的能力 | AUC 接近 1，曲线拱向左上角                                  |
| **决策边界**   | `3_Decision_Boundaries` | 模型筛选的截止线         | 绿线（2D模型）能比红/蓝虚线（1D）更灵活地包围目标群         |
| **PPI 热图**   | `7_Candidate_PPI`       | 候选蛋白内部的互作网络   | 出现明显的红色色块，表示候选蛋白间有强互作                  |
| **JP/PM 热图** | `13_JP`,`15_JPw`      | 蛋白互作的空间分布密度   | **左下角** （高排名区域）出现高亮色块，且显著高于背景 |
| **成分堆积图** | `Plot_AbundancePercent` | 候选集的纯度             | 目标定位（如红色 TP）的色块占比最高                         |
| **TP 曲线**    | `TP_Percent_Curve`      | 排序的准确性             | 曲线在 X 轴早期（Top 10-20%）保持极高的 Y 值（>80%）        |

---

## 6. 下一步建议 (Next Step)

作为科研工作者，您在运行完此脚本后，建议执行以下操作：

1. **检查 Phase 3 的 `TP_Percent_Curve`** ：选择曲线下面积最大或起始点最高的模型作为最终的候选筛选标准。
2. **查看 `Step29_Analysis_Phase3_Comparison` 下的 Excel 表** ：提取最佳模型对应的 `Detailed_Candidates` Sheet。
3. **Cytoscape 可视化** ：利用 Phase 1 导出的 `_edges_ModelScores.csv` 和 `_nodes.csv`，在 Cytoscape 中绘制带权重的互作网络，权重可直接使用模型预测概率（Prob）。



# 亚细胞组分模型比较与稳健性评估系统 (Step29_ModelComparison) 项目指南

## 1. 项目概述与设计框架

### 1.1 项目背景

本项目是亚细胞定位分析流程的高级评估模块。在构建了多种预测模型（1D/2D/3D）后，必须回答以下关键科学问题：

1. **模型优越性** ：引入 PPI 或丰度数据后，模型性能是否显著优于传统的 logFC 筛选？
2. **方法学比较** ：不同的邻近标记方法（如 BioID, APEX 等，即不同的 Group）在捕捉特定细胞器组分时的表现有何差异？
3. **稳健性 (Robustness)** ：模型的高性能是来自于特定样本的偶然性，还是具有广泛的稳定性？

### 1.2 核心设计框架

脚本采用 **并行计算 + 引导抽样 (Bootstrap) + 非参数检验** 的高强度统计框架：

* **配置与初始化 (Configuration)** ：高度参数化的接口，允许用户灵活选择对比的组别、模型及可视化开关。
* **Phase 1: 核心计算 (Core Calculation)**
  * **PPI 整合** ：加载本地 STRING 数据库。
  * **并行模型训练** ：对每个组别训练全套 GLM 模型。
  * **统计检验** ：计算 AUC，执行 DeLong 检验（对比传统方法），执行 Bootstrap 检验（计算特异性/敏感性增益）。
* **Phase 1.1: 汇总与稳健性分析 (Consolidation & Robustness)**
  * **Type 1 稳健性** ：样本交集公平性分析（Union vs Intersection）。
  * **Type 2 稳健性** ：模型稳定性分析（AUC 的标准差）。
  * **跨组比较** ：使用 DeLong 检验对比不同实验组（邻近标记方法）的 ROC 曲线。
* **Phase 1.2: 可视化系统 (Visualization)**
  * 生成 7 类关键图表（森林图、热图、哑铃图、散点图、增益图等），全方位展示评估结果。

---

## 2. 统计思想与数学原理

### 2.1 ROC 曲线比较 (DeLong's Test)

脚本不只看 AUC 数值大小，还使用 **DeLong's Test** 来计算两个 ROC 曲线之间差异的统计显著性（P-value）。

* **应用场景** ：
* 模型 vs 传统方法（Model vs `logFC>0.5 & FDR<0.05`）。
* 模型 vs 1D_logFC 模型。
* 实验组 vs 实验组（如 Group A vs Group B）。

### 2.2 控制变量比较 (Controlled Variable Comparison via Bootstrap)

这是一个高级统计分析，用于证明“引入新变量（如 PPI）带来的性能提升不是偶然的”。

* **原理** ：

1. **固定 Sensitivity** ：在 ROC 曲线上找到与对照模型（如 logFC）相同的敏感性点，比较此时的特异性（Specificity）差异。
2. **固定 Specificity** ：同理，比较敏感性（Sensitivity）差异。
3. **Bootstrap** ：重采样 500-1000 次，构建差异分布，计算 95% 置信区间 (CI) 和 P 值。

* **意义** ：证明了模型在“不牺牲查全率的情况下提高了查准率”（或反之）。

### 2.3 多重检验校正 (BH Correction)

由于脚本涉及大量两两比较，所有输出的 P 值均经过 **Benjamini-Hochberg (BH)** 校正，生成 Q 值 (q-value)，以控制假阳性率。

---

## 3. 接口配置说明 (User Interface)

位于脚本开头的 `配置接口 - 用户可选设置` 区域是控制分析流向的关键。

### 3.1 核心选择

* `selected_groups`: 指定参与分析的实验组（如 `c("E7A2B4", "C2")`）。设为 `NULL` 则分析所有组。
* `selected_models`: 指定参与对比的模型列表。建议包含 `1D_Abundance`, `1D_logFC`, `2D_Abundance_logFC` 等关键模型。

### 3.2 可视化开关 (`enable_*`)

一套布尔值 (`TRUE`/`FALSE`) 开关，用于精细控制输出内容：

* `enable_cytoscape_export`: 是否输出网络文件。
* `enable_robustness_analysis`: 是否跑耗时的交集分析。
* `enable_plot_auc_summary`: 是否画森林图（图1）。
* `enable_plot_controlled_comparison`: 是否画 Bootstrap 增益图（图5）。
* ... (涵盖图1至图7)

### 3.3 统计参数

* `n_bootstrap`: 控制变量比较的重采样次数。**建议：测试时设 200，正式出图设 1000，严格发表设 2000+。**
* `parallel_cores`: 并行核心数（默认 20），根据服务器性能调整。

---

## 4. 模块详细设计与图表意义

### 4.1 Phase 1: 数据计算 (`process_group` 函数)

* **功能** ：这是并行计算的最小单元。它在一个隔离的环境中加载数据、计算局部 PPI、训练 GLM 模型，并立即运行 DeLong 和 Bootstrap 检验。
* **输出** ：不仅返回 AUC，还返回完整的 `roc_objects` 和 `controlled_comparison` 数据框。

### 4.2 Phase 1.1: 汇总与稳健性

* **Group Comparison (跨组比较)** ：
* 选取一个基准模型（默认为 `2D_Abundance_logFC`）。
* 在不同实验组间进行两两 DeLong 检验。
* **输出** ：`Group_Comparison_Summary.xlsx`。
* **Intersection Analysis (交集分析)** ：
* **思想** ：如果只看所有组都鉴定到的那些“共有蛋白”，各组的表现差异如何？这排除了因不同方法鉴定深度不同带来的偏差。
* **输出** ：`Group_Comparison_Intersection_Summary.xlsx`。

### 4.3 Phase 1.2: 可视化系统 (Plots)

| **图表编号** | **图表名称**                | **文件名关键词**                | **统计意义与解读**                                                                              |
| ------------------ | --------------------------------- | ------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| **图 1**     | **分组森林图**              | `1_AUC_Faceted_Forest`              | **模型优劣概览** 。展示各模型 AUC 及 95% CI。红色点表示显著优于传统方法。                       |
| **图 1b**    | **vs logFC 森林图**         | `1b_AUC_..._vs_1D_logFC`            | **增量价值证明** 。专门展示各模型是否显著优于仅使用 logFC 的模型。                              |
| **图 2a/b**  | **组间比较图**              | `2a_Bar`,`2b_Heatmap`             | **方法学评估** 。哪种邻近标记方法（Group）的 AUC 最高？热图展示两两差异的显著性。               |
| **图 3**     | **哑铃图**                  | `3_Intersection_Dumbbell`           | **数据公平性** 。对比“全数据 AUC”与“交集数据 AUC”。若差异巨大，说明方法依赖于特有蛋白。     |
| **图 4**     | **性能-稳定性图**           | `4_Performance_vs_Stability`        | **模型选择指南** 。X轴为 AUC 标准差（稳定性），Y轴为 AUC 均值。 **优选左上角模型** 。     |
| **图 5a/b**  | **增益图 (vs Traditional)** | `5a_Specificity`,`5b_Sensitivity` | **临床/应用价值** 。展示在固定传统方法的灵敏度下，特异性提升了多少。误差棒不跨过 0 线即为显著。 |
| **图 5c**    | **综合增益散点图**          | `5c_Combined_Gain`                  | **综合评估** 。右上象限的点表示双重指标均优于传统方法。                                         |
| **图 6**     | **定位分布堆叠图**          | `6_Localization_Distribution`       | **生物学验证** 。展示筛选出的候选蛋白在已知定位（MultiBait）中的分布比例。                      |
| **图 7**     | **增益图 (vs logFC)**       | `7a/b/c_..._vs_1D_logFC`            | 同图 5，但参照系换为 1D_logFC 模型，更严格地证明多维模型的价值。                                      |

---

## 5. 关键注意事项

1. **`myTP_vector` 的作用域** ：

* 在脚本中，`myTP_vector` 在 `process_group` 函数内部被重新定义了一次。这是为了防止并行计算环境（`future` 包）中全局变量无法传递的问题。如果您修改了 TP 的定义， **必须同时修改函数内外的定义** 。

1. **Bootstrap 耗时** ：

* `n_bootstrap` 设置为 1000 时，计算量巨大（模型数 × 组数 × 1000 次重采样）。在调试代码时，建议将其临时改为 50 或 100。

1. **Cytoscape 导出** ：

* 脚本会为每个组别导出 `_node_attributes.csv` 和 `_edge_attributes.csv`。
* **Node Attributes** 包含了所有模型的预测概率 (`prob_...`)，可在 Cytoscape 中通过映射颜色来直观展示不同模型的预测差异。

## 6. 下一步行动建议 (Next Step)

1. **运行脚本** ：首先使用默认参数运行，检查生成的 `Plot 1`（森林图）和 `Plot 4`（性能-稳定性图）。
2. **选择最佳模型** ：根据 Plot 4，选择位于 **左上角** （高 AUC 均值，低标准差）的模型。通常 `2D_Abundance_logFC` 是一个很好的平衡点。
3. **查阅增益图 (Plot 7c)** ：确认您选择的最佳模型是否位于右上象限（即在灵敏度和特异性上均优于 1D_logFC），并将此图作为论文中的核心证据。
4. **导出候选列表** ：去 `Step29_ModelComparision` 文件夹下查看 `detailed_candidates_by_model.xlsx`，获取最终筛选的蛋白列表。
