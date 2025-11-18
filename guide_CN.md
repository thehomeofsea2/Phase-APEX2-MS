# 质谱分析流程模块化项目指南

## 项目概述

将 CleanCode.R (3322行) 重构为14个独立可调试模块，实现手动分步执行的质谱数据分析流程。

**核心特点：**
- 模块间通过 RData 传递数据
- 每个模块可独立调试
- 保留原代码所有分析细节

## 目录结构

```
CodeNorm/
├── main_pipeline.R          # 主控脚本
├── guide.md                 # 本文档
├── CleanCode.R              # 原始代码（参考）
├── Module/                  # 模块代码
│   ├── module01_setup.R
│   ├── module02_data_import.R
│   └── ...
├── Dev/                     # 调试代码
├── Output/                  # 输出结果
├── Reference/               # 注释文件（必需）
├── Rawdata/                 # 原始数据（必需）
└── *.RData                  # 模块间数据传递
```


## 模块列表

### Module 1: 环境初始化 ✓

**文件：** `Module/module01_setup.R`

**功能：**
1. 检查必需目录（Reference, Rawdata）- 不存在则报错
2. 自动创建可选目录（Output, Module, Dev）
3. 保存配置信息到 RData

**输入：**
- `dir_config` - 目录配置列表

**输出：**
- `Module01_workspace.RData` - 所有环境变量（工作目录）
  - 包含：dir_config, config

**函数：**
```r
module01_setup(dir_config)
```

**使用示例：**
```r
# 在 main_pipeline.R 中
source(file.path(dir_config$module, "module01_setup.R"))
module01_setup(dir_config)
load("Module01_workspace.RData")
```

**测试方法：**
1. 确保 Reference/ 和 Rawdata/ 目录存在
2. 在 main_pipeline.R 中运行 Module 1 部分
3. 检查是否生成 Module01_workspace.RData
4. 检查 Output/, Module/, Dev/ 目录是否创建

**注意事项：**
- Reference 和 Rawdata 为必需目录，缺失会报错终止
- 其他目录会自动创建

---

**以下模块待完成，将在用户确认后逐步更新**

### Module 2: 数据读取与分组表 ✓

**文件：** `Module/module02_data_import.R`

**功能：**
1. 读取TSV数据文件（支持pattern匹配）
2. 自动清理列名（去除共同前缀）
3. 生成sampleGroup模板（含所有可选值提示）
4. 读取用户填写的sampleGroup
5. 生成FinalName并重命名数据列

**输入：**
- `Module01_workspace.RData` - Module01环境变量（通过load加载）
- TSV数据文件（在Rawdata/目录）
- 用户填写的`Module02_sampleGroup_template.csv`

**输出：**
- `Module02_sampleGroup_template.csv` - 供用户填写的模板（工作目录，用户编辑）
- `Module02_sampleGroup.csv` - 验证后的最终分组信息（工作目录，程序生成）
- `Module02_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01所有 + data_raw, sampleGroup
- `Output/Module02_data_raw.csv` - CSV格式数据（Output目录）

**关键函数：**
```r
# 步骤1: 读取数据并生成模板
result <- module02_read_and_generate_template(dir_config, file_pattern = "_matrix.*\\.tsv$")

# 步骤2: 等待用户填写 Module02_sampleGroup_template.csv

# 步骤3: 读取template，验证，生成Module02_sampleGroup.csv
result <- module02_process_samplegroup(dir_config, result$data)
```

**sampleGroup表结构：**
- OriginalName: 清理后的原始列名
- bioGroup: 生物学分组名（用户填写，必填）
- CatalyticGroup: Cata/NoCat（模板已填充示例值，如Cata/Cata/NoCat/NoCat...）
- PLtype: Light/H2O2/PL（模板已填充示例值，如Light/Light/H2O2/H2O2/PL...）
- Context: Experiment/Control/Spatial（模板已填充示例值）
- replicate: 1/2/3（模板已填充，循环123）
- FirstROCgroup: A/B/C/A&B/NA（模板已填充示例值）
- SecondROCgroup: A/B/C/D/E/F（模板已填充示例值）
- **Order: 整数（1/2/3...），控制数据列展示顺序，1最先出现**（模板默认1,2,3...）
- FinalName: 自动生成（bioGroup_LFQ_replicate）

**模板示例值逻辑：**
- 所有可选元素至少出现一次
- 智能均匀分配（如10个样本分配3个选项：A/A/A/A/B/B/B/C/C/C）
- replicate自动循环1/2/3
- Order默认为1,2,3...（按原始顺序）

**Order列说明：**
- 用户可修改Order值来调整数据列的展示顺序
- 数据处理流程：**Order排序 → 列重排 → FinalName重命名**
- Order=1的列最先出现，Order=2次之，依此类推
- Gene列始终保持在第一列

**使用流程：**
1. 运行main_pipeline.R中Module 2部分
2. 打开生成的`Module02_sampleGroup_template.csv`
3. 填写所有必需列（参考提示值，bioGroup必填）
4. **可选：修改Order列调整数据列顺序**（例如将某些样本移到前面）
5. **保存`Module02_sampleGroup_template.csv`**（不需要另存为）
6. 按回车继续，程序将自动验证并生成`Module02_sampleGroup.csv`

**测试方法：**
1. 检查Rawdata/目录中是否有TSV文件
2. 运行Module 2步骤1，检查是否生成`Module02_sampleGroup_template.csv`
3. 填写template（可测试修改Order值，如将第5列改为Order=1）
4. 保存`Module02_sampleGroup_template.csv`并按回车继续
5. 程序自动读取template，验证数据，生成`Module02_sampleGroup.csv`
6. 验证数据列是否按Order重新排列
7. 检查`Module02_workspace.RData`是否生成
8. 验证列名是否正确重命名为FinalName格式

**配置说明：**
- 在`main_pipeline.R`中配置`data_config$file_pattern`修改文件匹配模式
- 默认：`"_matrix.*\\.tsv$"`匹配包含"_matrix"的TSV文件

**注意事项：**
- 模板已智能填充示例值，可直接参考或修改
- FinalName会自动生成，模板中留空即可
- bioGroup和Order必须填写，否则报错
- Order必须为整数，用于控制数据列展示顺序
- 数据处理顺序：**先按Order排序重排列，再按FinalName重命名**
- **工作流程**：用户填写template → 保存 → 程序读取验证 → 自动生成Module02_sampleGroup.csv
- **文件命名**：所有输出文件带Module02_前缀
- RData文件保存在工作目录，CSV输出到Output/目录

---

### Module 3: 注释系统 ✓

**文件：** `Module/module03_annotation.R`

**功能：**
1. 读取HPA数据库，生成Cytosol/Nuclear/Nuclear_Cytosol注释
2. 读取MitoCarta数据，生成Mitochondrion注释
3. 读取SGs参考数据（GO_SGs, HaloMap_SGs, HaloMap_DifMethods）
4. 对数据进行多种注释（默认3种：HaloMap_Localization, GO_Localization, MultiBait_Localization）
5. 支持用户自定义额外注释

**输入：**
- `Module02_workspace.RData` - Module01-02所有环境变量（包含data_raw, sampleGroup等）
- Reference/ 目录中的注释文件：
  - `proteinatlas.tsv` - HPA数据库
  - `MitoCarta3.0.csv` - 线粒体蛋白
  - `CytoSGs GO.tsv` - GO注释的SGs
  - `HaloMap SGs reference.csv` - HaloMap SGs参考
  - `HaloMap differentMethods SGs.csv` - 不同方法的SGs

**输出：**
- `Module03_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-02所有 + data_annotated, annotation_references
- `Output/Module03_data_annotated.csv` - CSV格式输出（Output目录）

**关键函数：**
```r
# 加载注释参考
HPA_anno <- load_HPA_annotations(dir_config)
MitoCarta_anno <- load_MitoCarta_annotations(dir_config)
SGs_refs <- load_SGs_references(dir_config)

# 默认注释
module03_annotation(dir_config, data_raw, sampleGroup)

# 自定义注释
custom_annotations <- list(
  list(
    column_name = "Custom_Localization",
    TP_source = "GO_SGs",
    TP_column = "Gene",
    TP_label = "SGs"
  )
)
module03_annotation(dir_config, data_raw, sampleGroup, custom_annotations)
```

**注释优先级：**
1. TP（用户指定的True Positive，如SGs）
2. Nuclear / Cytosol / Nuclear_Cytosol（HPA注释）
3. Mitochondrion（MitoCarta注释）
4. Other（其他未注释）

**HPA注释规则：**
- **Cytosol**: 包含"yto/crotubule/ctin"，排除"ucleol/ucleoplasm/uclear/itochond"
- **Nuclear**: 包含"ucleol/ucleoplasm/uclear"，排除"yto/crotubule/ctin/itochond/Plasma/Vesicles/Golgi/Endoplasmic"
- **Nuclear_Cytosol**: 同时包含细胞质和核内特征，排除"itochond"

**默认注释列：**
- **HaloMap_Localization**: 使用HaloMap SGs reference作为TP
- **GO_Localization**: 使用GO注释的SGs作为TP
- **MultiBait_Localization**: 使用HaloMap不同方法的SGs作为TP

**测试方法：**
1. 确保Reference/目录包含所有必需注释文件
2. 运行Module 3
3. 检查Module03_workspace.RData是否生成
4. 验证注释列是否正确添加（HaloMap_Localization, GO_Localization, MultiBait_Localization）
5. 检查各注释类别的蛋白数量
6. 确认Output/目录中生成Module03_data_annotated.csv

**注意事项：**
- HPA数据会过滤Reliability为"Uncertain"的条目
- 各注释集合间可能有交集，按优先级处理
- 自定义注释会追加到默认注释之后
- **文件命名**：所有输出文件带Module03_前缀

---

### Module 4: 数据标准化 ✓

**文件：** `Module/module04_standardization.R`

**功能：**
1. Log2转化：对数据进行log2转换
2. 全局标准化（Global Normalization）：
   - Global_QNorm：Quantile标准化
   - Global_MNorm：Median标准化
3. 局部标准化（Local Normalization）：按bioGroup分组标准化
   - Local_QNorm：组内Quantile标准化
   - Local_MNorm：组内Median标准化
4. Boxplot可视化：展示不同标准化方法的效果，按bioGroup着色

**输入：**
- `Module03_workspace.RData` - Module01-03所有环境变量（包含data_annotated, sampleGroup等）
- `norm_types` - 指定需要的标准化类型（可配置）

**输出：**
- `Module04_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-03所有 + standardized_data_list, norm_types_used
- CSV文件（Output目录）：
  - `Module04_log2.csv` - Log2转化后的数据
  - `Module04_Global_QNorm.csv` - 全局Quantile标准化
  - `Module04_Global_MNorm.csv` - 全局Median标准化
  - `Module04_Local_QNorm.csv` - 局部Quantile标准化
  - `Module04_Local_MNorm.csv` - 局部Median标准化
- PDF文件（Output目录）：
  - `Module04_log2_boxplot.pdf` - Log2转化后的boxplot
  - `Module04_Global_Norm_boxplot.pdf` - 全局标准化对比
  - `Module04_Local_Norm_boxplot.pdf` - 局部标准化对比
  - `Module04_All_Norm_comparison_boxplot.pdf` - 所有方法综合对比

**关键函数：**
```r
# 使用所有标准化方法
norm_types <- c("noNorm", "Global_QNorm", "Global_MNorm", "Local_QNorm", "Local_MNorm")
result <- module04_standardization(dir_config, data_annotated, sampleGroup, norm_types)

# 只使用全局标准化
norm_types <- c("noNorm", "Global_QNorm", "Global_MNorm")
result <- module04_standardization(dir_config, data_annotated, sampleGroup, norm_types)

# 只使用Median标准化方法
norm_types <- c("noNorm", "Global_MNorm", "Local_MNorm")
result <- module04_standardization(dir_config, data_annotated, sampleGroup, norm_types)
```

**标准化方法说明：**

**1. noNorm（仅Log2转化）**
- 对原始数据进行log2转换
- 不进行任何标准化处理
- 用于对比标准化效果

**2. Global_QNorm（全局Quantile标准化）**
- 使用`preprocessCore::normalize.quantiles()`
- 对所有样本统一进行分位数标准化
- 使各样本的数据分布一致

**3. Global_MNorm（全局Median标准化）**
- 计算每个样本的中位数
- 计算所有样本的全局中位数
- 用中位数比值对每个样本进行缩放

**4. Local_QNorm（局部Quantile标准化）**
- 按`sampleGroup$bioGroup`分组
- 在每个bioGroup内部分别进行Quantile标准化
- 然后合并所有组的结果

**5. Local_MNorm（局部Median标准化）**
- 按`sampleGroup$bioGroup`分组
- 在每个bioGroup内部分别进行Median标准化
- 然后合并所有组的结果

**Boxplot可视化特点：**
- 不同bioGroup使用不同颜色
- 自动生成颜色调色板（≤12个组用hue_pal，>12个组用rainbow）
- 每个图添加图例标识bioGroup
- 轴标签字体缩小（cex.axis=0.4）便于查看大量样本

**配置说明：**
在`main_pipeline.R`中配置`norm_types`：
```r
# 必须至少包含一个全局标准化（Global_QNorm或Global_MNorm）
norm_types <- c("noNorm", "Global_QNorm", "Global_MNorm", "Local_QNorm", "Local_MNorm")
```

**测试方法：**
1. 确保Module 3已完成，生成`Module03_workspace.RData`
2. 配置需要的`norm_types`
3. 运行Module 4
4. 检查Output/目录中生成的CSV文件（数量应与norm_types匹配）
5. 检查PDF文件，验证boxplot按bioGroup着色
6. 验证不同标准化方法的效果差异
7. 检查`Module04_workspace.RData`是否生成
8. 验证`standardized_data_list`包含所有标准化版本

**注意事项：**
- **必须包含至少一个全局标准化**（Global_QNorm或Global_MNorm）
- 需要安装`preprocessCore`包：`BiocManager::install("preprocessCore")`
- Local标准化需要sampleGroup中的bioGroup信息
- 数据列和注释列自动识别（注释列需包含"Localization"后缀）
- Boxplot使用log="y"刻度，便于观察数据分布
- 颜色映射：同一bioGroup的样本使用相同颜色
- **文件命名**：所有输出文件带Module04_前缀
- RData文件保存在工作目录，CSV/PDF输出到Output/目录

---

### Module 5: 缺失值填补 ✓

**文件：** `Module/module05_imputation.R`

**功能：**
1. Perseus方法缺失值填补：基于正态分布随机抽样
2. 分组智能填补：区分NoCat组和Cat组
3. n_valid分类填补：根据有效值数量采用不同策略
4. Boxplot可视化：填补前后对比，按bioGroup着色

**输入：**
- `Module04_workspace.RData` - Module01-04所有环境变量（包含standardized_data_list, sampleGroup等）
- `impute_cat_mean` - 是否对Cat组填补（默认FALSE）
- `random_seed` - 随机数种子（默认123）

**输出：**
- `Module05_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-04所有 + imputed_data_list, imputation_params
- CSV文件（Output目录）：
  - `Module05_noNorm_Imputed.csv` - noNorm填补后数据
  - `Module05_Global_QNorm_Imputed.csv` - Global_QNorm填补后数据
  - `Module05_Global_MNorm_Imputed.csv` - Global_MNorm填补后数据
  - `Module05_Local_QNorm_Imputed.csv` - Local_QNorm填补后数据
  - `Module05_Local_MNorm_Imputed.csv` - Local_MNorm填补后数据
- PDF文件（Output目录）：
  - `Module05_Imputation_Before_After_boxplot.pdf` - 所有数据集的填补前后对比

**关键函数：**
```r
# 默认配置：只填补NoCat组，不填补Cat组
impute_cat_mean <- FALSE
random_seed <- 123
result <- module05_imputation(dir_config, standardized_data_list, sampleGroup, 
                              impute_cat_mean, random_seed)

# 启用Cat组填补（n_valid=2时用平均值填补）
impute_cat_mean <- TRUE
result <- module05_imputation(dir_config, standardized_data_list, sampleGroup, 
                              impute_cat_mean, random_seed)

# 更改随机种子
random_seed <- 456
result <- module05_imputation(dir_config, standardized_data_list, sampleGroup, 
                              impute_cat_mean, random_seed)
```

**填补策略详解：**

**1. Perseus参数计算**
对每个数据列（样本）计算：
- `imputation_mean = col_mean - 1.8 * col_sd`
- `imputation_sd = 0.3 * col_sd`

**2. NoCat组填补规则**（必填）
根据`sampleGroup$CatalyticGroup == "NoCat"`识别：
- **n_valid == 2**：用bioGroup内两个有效值的平均值填补
- **n_valid < 2**：用Perseus参数随机抽样填补
  - 从正态分布`N(imputation_mean, imputation_sd)`中抽样
  - 保证结果可重复（使用random_seed）

**3. Cat组填补规则**（可选）
根据`sampleGroup$CatalyticGroup == "Cat"`识别：
- **默认（impute_cat_mean=FALSE）**：不填补，保留NA
- **启用（impute_cat_mean=TRUE）**：
  - n_valid == 2：用bioGroup内平均值填补
  - n_valid < 2：不填补，保留NA

**4. 填补逻辑说明**
```r
for (每个bioGroup) {
  if (组为NoCat) {
    if (n_valid == 2) {
      填补值 = mean(有效值)
    } else if (n_valid < 2) {
      填补值 = rnorm(n, imputation_mean, imputation_sd)
    }
  } else if (组为Cat && impute_cat_mean == TRUE) {
    if (n_valid == 2) {
      填补值 = mean(有效值)
    }
  }
}
```

**Boxplot可视化特点：**
- 每个标准化版本生成一对boxplot（Before vs After）
- 不同bioGroup使用不同颜色
- 自动生成颜色调色板
- 添加图例标识bioGroup
- 主标题显示数据集名称

**配置说明：**
在`main_pipeline.R`中配置填补参数：
```r
# impute_cat_mean: 是否对Cat组的n_valid=2情况用平均值填补
impute_cat_mean <- FALSE  # 默认不填补Cat组

# random_seed: 随机数种子，保证Perseus填补可重复
random_seed <- 123        # 可改为其他数值
```

**测试方法：**
1. 确保Module 4已完成，生成`Module04_workspace.RData`
2. 配置`impute_cat_mean`和`random_seed`
3. 运行Module 5
4. 检查Output/目录中生成的CSV文件（数量应与Module 4的标准化版本匹配）
5. 检查PDF文件，验证填补前后的差异
6. 验证注释列完整性（应与原始数据一致）
7. 检查`Module05_workspace.RData`是否生成
8. 验证`imputed_data_list`包含所有填补版本
9. 使用相同random_seed重复运行，验证结果一致性

**注意事项：**
- **必须先完成Module 2**（需要sampleGroup的CatalyticGroup和bioGroup信息）
- **必须先完成Module 4**（需要standardized_data_list）
- Perseus填补使用正态分布随机抽样，需要random_seed保证可重复性
- 注释列不会被填补，填补前后注释信息保持一致
- NoCat组的缺失值会被强制填补（这是必需的预处理步骤）
- Cat组的缺失值填补是可选的（根据impute_cat_mean参数）
- **n_valid定义**：每个蛋白质在同一bioGroup内的有效（非NA）样本数量
- **文件命名**：所有输出文件带Module05_前缀
- RData文件保存在工作目录，CSV/PDF输出到Output/目录
- 填补后的数据列数和行数与原始数据完全一致
- 使用`set.seed(random_seed)`确保Perseus填补的可重复性

---

### Module 6: 热图系统 ✓

**文件：** `Module/module06_heatmap.R`

**功能：**
1. 相关性热图（Correlation Heatmap）：样本间相关性分析
   - 高度可定制的样本筛选
   - PLtype×Context组合注释
   - 智能颜色映射规则
2. 全局表达热图（AllLocalization Heatmap）：所有蛋白质表达模式
3. 分类表达热图（by_localization Heatmap）：按Localization分类展示
4. 支持多版本数据同时绘制

**输入：**
- `Module05_workspace.RData` - Module01-05所有环境变量（包含imputed_data_list, sampleGroup等）
- `selected_versions` - 选择用于绘制的标准化版本（向量，如`c("noNorm_Imputed", "Local_QNorm_Imputed")`）
- `heatmap_types` - 需要绘制的热图类型（all/correlation/by_localization）
- `correlation_config` - 相关性热图详细配置（见下文）
- `localization_columns` - 用于分类的注释列（NULL表示自动检测）
- `color_params` - 表达热图颜色参数

**输出：**
- `Module06_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-05所有 + heatmap_info
- PDF文件（Output目录）：
  - `Module06_{version}_correlation_heatmap.pdf` - 相关性热图
  - `Module06_{version}_AllLocalization_heatmap.pdf` - 全局表达热图
  - `Module06_{version}_{localization_name}_heatmap.pdf` - 分类表达热图

**关键函数：**
```r
# 基本用法
result <- module06_heatmap(dir_config, imputed_data_list, sampleGroup,
                          selected_versions = c("noNorm_Imputed", "Local_QNorm_Imputed"),
                          heatmap_types = c("all", "correlation", "by_localization"))

# 高级配置
correlation_config <- list(
  # 相关性范围和颜色
  corr_min = 0.8, corr_max = 1, corr_center = 0.9,
  col_low = "#4D97CD", col_center = "white", col_high = "#DB6968",
  
  # 样本筛选（方式1：排除规则）
  exclude_context = NULL,      # 例如：c("Control") 排除所有Control组
  exclude_pltype = NULL,       # 例如：c("PL") 排除所有PL组
  exclude_catalytic = NULL,    # 例如：c("NoCat") 排除所有NoCat组
  
  # 样本筛选（方式2：直接指定，优先级高）
  include_biogroups = NULL     # 例如：c("Light_Experiment", "H2O2_Experiment")
)
```

**热图类型详解：**

**1. 相关性热图（correlation）**
- **功能**：展示样本间的Pearson相关系数
- **特点**：
  - 高度可定制的样本筛选
  - 三层列注释（PLtype、Context、Sample_Color）
  - 智能颜色映射（见下文"相关性热图颜色规则"）
  - 支持排除/包含特定样本类型

**2. 全局表达热图（all）**
- **功能**：展示所有蛋白质在所有样本中的表达
- **特点**：
  - 行注释：所有Localization列（如HaloMap_Localization）
  - 列注释：bioGroup（自动着色）
  - 自动过滤NA行（所有样本都是NA的蛋白质）

**3. 分类表达热图（by_localization）**
- **功能**：按Localization分类绘制独立热图
- **特点**：
  - 为每个Localization类别（SGs/Nuclear/Cytosol等）生成单独热图
  - 只包含该类别的蛋白质
  - 列注释：bioGroup
  - 自动命名：`Module06_{version}_{localization_name}_heatmap.pdf`

**相关性热图配置详解：**

**样本筛选：**
```r
# 方式1：排除规则（适合简单筛选）
correlation_config <- list(
  exclude_context = c("Control"),   # 只保留Experiment和Spatial组
  exclude_pltype = c("PL"),         # 排除所有PL实验
  exclude_catalytic = c("NoCat")    # 只保留催化实验组
)

# 方式2：直接指定bioGroup（适合精确控制，优先级高）
correlation_config <- list(
  include_biogroups = c("K69A1B3_Light", "K69C3_Light", "C3_H2O2")
)
```

**相关性热图颜色规则：**

**列注释结构：**
1. **PLtype行**：显示实验类型（Light/H2O2/PL）
2. **Context行**：显示实验环境（Experiment/Control/Spatial）
3. **Sample_Color行**：PLtype×Context组合颜色（主要可视化指标）

**Sample_Color颜色映射规则（优先级从高到低）：**

1. **Spatial优先级最高**：
   - 任何bioGroup包含"Spatial"关键词 → **黑色**
   
2. **Light实验（蓝色系）**：
   - Light + Experiment → **深蓝色 (#0066CC)**
   - Light + Control → **浅蓝色 (#87CEEB)**
   
3. **H2O2实验（橘黄色系）**：
   - H2O2 + Experiment → **深橙色 (#FF8C00)**
   - H2O2 + Control → **浅橙色/金色 (#FFD700)**
   
4. **PL实验（绿色系）**：
   - PL + Experiment → **深绿色 (#228B22)**
   - PL + Control → **浅绿色 (#90EE90)**
   
5. **其他实验（灰色系）**：
   - 其他 + Experiment → **深灰色 (#696969)**
   - 其他 + Control → **浅灰色 (#D3D3D3)**

**表达热图颜色配置：**
```r
# AllLocalization和by_localization热图的颜色参数
color_params <- list(
  custom_min = -4,           # 低值阈值
  custom_max = 4,            # 高值阈值
  custom_center = 0,         # 中心值
  col_low = "#4D97CD",       # 低值颜色（蓝色）
  col_center = "white",      # 中心颜色（白色）
  col_high = "#DB6968"       # 高值颜色（红色）
)
```

**使用示例：**

**示例1：绘制所有类型热图**
```r
selected_versions <- c("noNorm_Imputed", "Local_QNorm_Imputed")
heatmap_types <- c("all", "correlation", "by_localization")

correlation_config <- list(
  corr_min = 0.8, corr_max = 1, corr_center = 0.9,
  col_low = "#4D97CD", col_center = "white", col_high = "#DB6968",
  exclude_context = NULL, exclude_pltype = NULL, 
  exclude_catalytic = c("NoCat"),  # 只分析催化实验组
  include_biogroups = NULL
)

result <- module06_heatmap(dir_config, imputed_data_list, sampleGroup,
                          selected_versions, heatmap_types, 
                          correlation_config, NULL, color_params)
```

**示例2：只绘制相关性热图，排除对照组**
```r
selected_versions <- c("Local_QNorm_Imputed")
heatmap_types <- c("correlation")

correlation_config <- list(
  corr_min = 0.8, corr_max = 1, corr_center = 0.9,
  col_low = "#4D97CD", col_center = "white", col_high = "#DB6968",
  exclude_context = c("Control"),  # 排除所有对照组
  exclude_pltype = NULL, exclude_catalytic = NULL,
  include_biogroups = NULL
)

result <- module06_heatmap(dir_config, imputed_data_list, sampleGroup,
                          selected_versions, heatmap_types, 
                          correlation_config, NULL, NULL)
```

**示例3：直接指定特定bioGroup进行相关性分析**
```r
correlation_config <- list(
  corr_min = 0.8, corr_max = 1, corr_center = 0.9,
  col_low = "#4D97CD", col_center = "white", col_high = "#DB6968",
  exclude_context = NULL, exclude_pltype = NULL, exclude_catalytic = NULL,
  include_biogroups = c("K69A1B3_Light", "K69C3_Light", "K20_Light")  # 直接指定
)
```

**示例4：自定义表达热图颜色（蓝-白-红）**
```r
color_params <- list(
  custom_min = -3, custom_max = 3, custom_center = 0,
  col_low = "#0000FF",       # 纯蓝色
  col_center = "#FFFFFF",    # 纯白色
  col_high = "#FF0000"       # 纯红色
)
```

**测试方法：**
1. 确保Module 5已完成，生成`Module05_workspace.RData`
2. 配置`selected_versions`（至少一个标准化版本）
3. 配置`heatmap_types`（至少一种类型）
4. 配置`correlation_config`（如果绘制相关性热图）
5. 运行Module 6
6. 检查Output/目录中生成的PDF文件：
   - 相关性热图：检查列注释颜色是否符合规则
   - 全局表达热图：检查行列注释是否完整
   - 分类表达热图：检查是否为每个Localization生成独立PDF
7. 检查`Module06_workspace.RData`是否生成
8. 验证多个`selected_versions`时是否为每个版本生成独立热图

**技术细节：**

**聚类方法：**
- 聚类方法：`ward.D2`（Ward层次聚类）
- 距离度量：`correlation`（相关性距离）
- 适用范围：所有热图类型

**数据过滤：**
- 自动过滤所有样本都是NA的行
- 不显示行名（gene names），避免图片拥挤
- 列名自动旋转45度（便于阅读）

**注释自动检测：**
- 如果`localization_columns = NULL`，自动检测所有包含"Localization"的列
- 自动为每个Localization类别生成独立热图

**颜色生成：**
- bioGroup颜色：使用`scales::hue_pal()`或`rainbow()`自动生成
- Sample_Color：基于PLtype×Context组合的固定规则
- 表达热图：支持自定义三色渐变（low-center-high）

**配置说明：**
在`main_pipeline.R`中配置热图参数：
```r
# 选择标准化版本（支持向量，会为每个版本生成独立热图）
selected_versions <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# 选择热图类型
heatmap_types <- c("all", "correlation", "by_localization")

# 相关性热图配置（详细）
correlation_config <- list(
  corr_min = 0.8, corr_max = 1, corr_center = 0.9,
  col_low = "#4D97CD", col_center = "white", col_high = "#DB6968",
  exclude_context = NULL, exclude_pltype = NULL, 
  exclude_catalytic = c("NoCat"),
  include_biogroups = NULL
)

# Localization列（NULL=自动检测）
localization_columns <- NULL

# 表达热图颜色参数
color_params <- list(
  custom_min = -4, custom_max = 4, custom_center = 0,
  col_low = "#4D97CD", col_center = "white", col_high = "#DB6968"
)
```

**注意事项：**
- **必须先完成Module 5**（需要imputed_data_list）
- **相关性热图的样本筛选**：`include_biogroups`优先级高于排除规则
- **Spatial组识别**：任何bioGroup包含"Spatial"（不区分大小写）都会被标记为黑色
- **TP注释识别**：在Localization列中的"SGs"类别会在行注释中显示为红色
- **颜色规则的逻辑**：Spatial > PLtype×Context > 其他
- **Context颜色规则**：Experiment组颜色深，Control组颜色浅
- **需要安装pheatmap包**：`install.packages("pheatmap")`
- **需要安装scales包**：`install.packages("scales")`
- **PDF尺寸**：
  - 相关性热图：20×18英寸
  - AllLocalization热图：28×24英寸
  - by_localization热图：动态调整（基于蛋白数量）
- **文件命名**：所有输出文件带Module06_前缀和版本名称
- RData文件保存在工作目录，PDF输出到Output/目录
- 热图自动聚类，同类样本会自动聚集在一起
- 相关性热图只计算数值列，自动排除注释列

---

### Module 7: 第一次差异分析 ✓

**文件：** `Module/module07_diff_analysis1.R`

**功能：**
1. **自动构建对比组**：基于sampleGroup$FirstROCgroup自动解析并构建对比
2. 基于limma包进行差异表达分析
3. 使用经验贝叶斯方法（eBayes）调整统计量
4. 计算logFC（log2 Fold Change）和adj.P.Val（FDR校正后的p值）
5. 输出合并的差异分析结果、表达矩阵和LogFC宽表

**FirstROCgroup逻辑：**
- **同组识别**：包含相同元素的为同一组
  - 例如："A"、"A/B"、"A&B"都属于A组
- **对比构建**：在同一组内，Context为"Experiment"的bioGroup vs Context为"Control"的bioGroup
- **示例**：
  - FirstROCgroup=A, Context=Experiment: K69A1B3_Light
  - FirstROCgroup=A, Context=Control: K69A1B3_noLight, A1B3_Light
  - FirstROCgroup=A&B, Context=Control: K69_Light（同时属于A组和B组）
  - **自动生成的对比**：
    * K69A1B3_Light vs K69A1B3_noLight
    * K69A1B3_Light vs K69_Light
    * K69A1B3_Light vs A1B3_Light

**输入：**
- `Module06_workspace.RData` - Module01-06所有环境变量（包含imputed_data_list, sampleGroup等）
- `sampleGroup`必需列：FinalName, bioGroup, FirstROCgroup, Context
- `selected_versions` - 需要分析的数据版本（默认NULL表示所有版本）

**输出：**
- `Module07_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-06所有 + diff_results1, comparisons_used, comparison_info
- CSV文件（Output目录）：
  - `Module07_{version}_DiffAnalysis.csv` - 合并的差异分析结果（所有对比组的logFC和adj.P.Val）
  - `Module07_{version}_ExprMatrix.csv` - 表达矩阵（Gene × Samples）
  - `Module07_{version}_LogFC_Wide.csv` - LogFC宽表（Gene × Comparisons）
- XLSX文件（Output目录）：
  - `Module07_{version}_DiffAnalysis_Full.xlsx` - 完整结果（多sheet）

**关键函数：**
```r
# 基本用法（对比组自动从sampleGroup$FirstROCgroup构建）
result <- module07_diff_analysis1(dir_config, imputed_data_list, sampleGroup,
                                   selected_versions)

# 选择分析版本
selected_versions <- c("noNorm_Imputed", "Local_QNorm_Imputed")
```

**差异分析流程：**
```
1. 解析FirstROCgroup
   └─> 识别所有ROC组元素（A, B, C等）
   └─> 分割复合标记（如A&B、A/B）
   └─> 为每个ROC组找出Experiment和Control的bioGroup

2. 自动构建对比组
   └─> 在每个ROC组内：Experiment bioGroup vs 每个Control bioGroup
   └─> 示例：A组内 K69A1B3_Light vs {K69A1B3_noLight, K69_Light, A1B3_Light}

3. 数据准备
   └─> 提取表达矩阵（只保留数值列）
   └─> 行名设为Gene symbol
   └─> 确保sampleGroup顺序与数据列一致

4. 构建设计矩阵（Design Matrix）
   └─> 使用 model.matrix(~ 0 + Group)
   └─> Group来自sampleGroup$bioGroup

5. 拟合线性模型
   └─> 使用 limma::lmFit(expr_matrix, design)

6. 定义对比矩阵（Contrast Matrix）
   └─> 使用 limma::makeContrasts()
   └─> 每个对比：GroupA - GroupB

7. 计算对比并应用经验贝叶斯调整
   └─> 使用 limma::contrasts.fit()
   └─> 使用 limma::eBayes()

8. 提取结果
   └─> 使用 limma::topTable(n = Inf, adjust.method = "BH")
   └─> 提取logFC和adj.P.Val
   └─> 合并所有对比结果
   └─> 生成LogFC宽表和表达矩阵
```

**输出文件结构：**

**CSV文件列结构（合并结果）：**
```
Gene  |  Comparison1_logFC  |  Comparison1_adj.P.Val  |  Comparison2_logFC  |  Comparison2_adj.P.Val  |  ...  |  [注释列]
```

**CSV文件列结构（ExprMatrix）：**
```
Gene  |  Sample1  |  Sample2  |  Sample3  |  ...
```

**CSV文件列结构（LogFC_Wide）：**
```
Gene  |  Comparison1  |  Comparison2  |  Comparison3  |  ...
```

**XLSX文件结构（多sheet）：**
- **Combined** sheet：合并结果（与DiffAnalysis.csv相同）
- **LogFC_Wide** sheet：LogFC宽表
- **ExprMatrix** sheet：表达矩阵
- **Comparison1** sheet：对比组1的完整topTable结果
  - 列：Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B
- **Comparison2** sheet：对比组2的完整topTable结果
- ...

**使用示例：**

**示例1：分析所有标准化版本（对比组自动构建）**
```r
# Module 7会自动从sampleGroup$FirstROCgroup构建对比组
result <- module07_diff_analysis1(dir_config, imputed_data_list, sampleGroup,
                                   selected_versions = NULL)
```

**示例2：只分析特定版本**
```r
# 只分析两个版本
selected_versions <- c("noNorm_Imputed", "Local_QNorm_Imputed")

result <- module07_diff_analysis1(dir_config, imputed_data_list, sampleGroup,
                                   selected_versions)
```

**示例3：FirstROCgroup配置示例**

| bioGroup | FirstROCgroup | Context | 说明 |
|----------|---------------|---------|------|
| K69A1B3_Light | A | Experiment | A组实验组 |
| K69A1B3_noLight | A | Control | A组对照1 |
| A1B3_Light | A | Control | A组对照2 |
| K69_Light | A&B | Control | 同时属于A组和B组的对照 |
| K69C3_Light | B | Experiment | B组实验组 |
| K69C3_noLight | B | Control | B组对照1 |
| C3_Light | B | Control | B组对照2 |
| C3_H2O2 | C | Experiment | C组实验组 |
| C3_noH2O2 | C | Control | C组对照 |
| K20_Light | - | Spatial | 不参与第一次差异分析 |

**自动生成的对比：**
- A组：K69A1B3_Light vs K69A1B3_noLight
- A组：K69A1B3_Light vs A1B3_Light
- A组：K69A1B3_Light vs K69_Light
- B组：K69C3_Light vs K69C3_noLight
- B组：K69C3_Light vs C3_Light
- B组：K69C3_Light vs K69_Light
- C组：C3_H2O2 vs C3_noH2O2

**测试方法：**
1. 确保Module 5已完成，生成`Module05_workspace.RData`（包含imputed_data_list）
2. 确保sampleGroup包含必需列：FinalName, bioGroup, FirstROCgroup, Context
3. 检查FirstROCgroup配置：
   - 确保有Experiment和Control组
   - 使用&或/分隔符标记多组成员（如"A&B"）
4. 配置`selected_versions`（NULL表示所有版本）
5. 运行Module 7
6. 检查控制台输出：
   - 查看检测到的ROC组
   - 查看自动生成的对比组
7. 检查Output/目录中生成的文件：
   - `*_DiffAnalysis.csv`：合并的差异分析结果
   - `*_ExprMatrix.csv`：表达矩阵
   - `*_LogFC_Wide.csv`：LogFC宽表
   - `*_DiffAnalysis_Full.xlsx`：完整结果（多sheet）
8. 检查`Module07_workspace.RData`是否生成
9. 验证`diff_results1`列表结构：
   ```r
   names(diff_results1)  # 应包含所有分析的版本
   names(diff_results1[[1]])  # 应包含: combined, logfc_wide, expr_matrix, raw_results, comparisons, comparison_info, contrast_names
   ```
10. 检查自动生成的对比：
   ```r
   # 查看自动生成的对比组
   comparison_info  # 查看每个对比的ROC组、Experiment组、Control组信息
   ```
11. 检查统计显著性：
   ```r
   # 查看显著差异蛋白数量（adj.P.Val < 0.05, |logFC| > 1）
   sig_proteins <- diff_results1$noNorm_Imputed$combined %>%
     filter(abs(K69A1B3_Light_vs_K69A1B3_noLight_logFC) > 1 &
            K69A1B3_Light_vs_K69A1B3_noLight_adj.P.Val < 0.05)
   nrow(sig_proteins)
   ```

**技术细节：**

**limma差异分析统计学：**
- **Linear Model**: 对每个蛋白质拟合线性模型
- **Empirical Bayes**: 借用所有蛋白质的信息，改进单个蛋白质的方差估计
- **Multiple Testing Correction**: 使用Benjamini-Hochberg（BH）方法控制FDR

**topTable输出列说明：**
- **logFC**: log2 fold change（log2倍数变化）
- **AveExpr**: average log2-expression（平均log2表达量）
- **t**: moderated t-statistic（调整后的t统计量）
- **P.Value**: raw p-value（原始p值）
- **adj.P.Val**: adjusted p-value（FDR校正后的p值）
- **B**: log-odds that the gene is differentially expressed（差异表达的对数优势比）

**对比组命名规则：**
- 格式：`GroupA_vs_GroupB`
- 特殊字符会被`make.names()`转换（如空格转为`.`）
- logFC > 0：GroupA中表达量高于GroupB
- logFC < 0：GroupA中表达量低于GroupB

**配置说明：**
在`main_pipeline.R`中配置差异分析参数：
```r
# 选择分析版本（NULL表示所有版本）
# 对比组会自动从sampleGroup$FirstROCgroup构建
selected_versions_diff <- c("noNorm_Imputed", "Local_QNorm_Imputed")
```

**注意事项：**
- **必须先完成Module 2**（需要sampleGroup，包含FirstROCgroup和Context列）
- **必须先完成Module 5**（需要imputed_data_list）
- **sampleGroup必需列**：FinalName, bioGroup, FirstROCgroup, Context
- **FirstROCgroup格式**：
  - 单组：直接使用字母（如"A"、"B"、"C"）
  - 多组：使用&或/分隔（如"A&B"、"A/B"）
  - 不参与：留空或NA（如Spatial组）
- **Context取值**：必须为"Experiment"或"Control"
- **对比方向**：Experiment - Control，正值表示Experiment组高表达
- **需要安装tidyr包**：`install.packages("tidyr")`
- **需要安装limma包**：`BiocManager::install("limma")`
- **需要安装openxlsx包**：`install.packages("openxlsx")`
- **Excel sheet名称限制**：最长31字符，超过会被截断
- **文件命名**：所有输出文件带Module07_前缀和版本名称
- RData文件保存在工作目录，CSV/XLSX输出到Output/目录
- 差异分析会自动排除注释列，只分析数值列
- 结果中的Gene列来自原始数据，如果数据中没有Gene列会使用行号
- **adj.P.Val阈值建议**：< 0.05（控制FDR在5%以下）
- **logFC阈值建议**：|logFC| > 1（表示至少2倍变化）
- **输出文件数量**：每个版本生成3个CSV文件和1个XLSX文件
- **自动化程度高**：无需手动指定对比组，完全基于FirstROCgroup自动构建

**技术要点：对比名称匹配机制**
- **limma生成的对比系数名称**：格式为 `GroupA - GroupB`（带空格和减号）
- **模块内部友好名称**：格式为 `GroupA_vs_GroupB`（下划线格式，用于文件输出）
- **关键实现**：
  - 从 `colnames(contrast_matrix)` 获取limma生成的实际系数名称
  - 使用实际系数名称调用 `limma::topTable(..., coef = actual_coef_name)`
  - 使用友好名称保存结果和生成输出文件（便于文件系统和R对象命名）
- **错误处理**：三层检查机制确保稳定性
  1. 检查bioGroup是否在设计矩阵中（排除数据中不存在的组）
  2. 检查系数名称是否在contrast_matrix中（排除无效对比）
  3. 使用tryCatch包裹topTable调用（捕获个别对比失败）
- **"Partial NA coefficients"警告**：正常现象，limma对某些蛋白质的系数估计不完整，不影响整体结果

---

### Module 8: 第一次ROC分析

**功能：**
1. **SubMito定位转化**：将注释中的"Mitochondrion"替换为MitoCarta3的亚定位（MIM, Matrix, MOM, IMS）
2. **ROC分析**：基于logFC和注释列进行ROC曲线分析，识别最佳阈值
3. **输出可视化**：绘制ROC曲线和Youden Index图
4. **保存数据**：输出ROC数据、最佳阈值、转化后的注释数据

**输入：**
- `Module07_workspace.RData` - Module01-07所有环境变量（包含diff_results1, annotation_references等）
- MitoCarta3.0注释数据（用于SubMito转化）

**输出：**
- `Module08_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-07所有 + roc_results1, roc_thresholds1, data_with_submito, expr_fdr_df_list
- `Output/Module08_{version}_ROC_curves.pdf` - ROC曲线图
- `Output/Module08_{version}_Youden_Index.pdf` - Youden Index图
- `Output/Module08_{version}_ROC_data.xlsx` - ROC数据表（多sheet）
- `Output/Module08_{version}_data_with_SubMito.csv` - SubMito转化后的数据
- `Output/Module08_all_thresholds.xlsx` - 所有版本的最佳阈值
- `Output/Module08_Expr_FDR_df_list.xlsx` - 每个版本一张sheet，用于Step13 A/B过滤的输入底表（不含SubMito注释）

**Expr_FDR_df_list 构建逻辑（对齐 CleanCode.R 1728-1739）：**
1. `mydata1` = 标准化且补缺失值后的表达矩阵（即进入热图的数据，等同于本项目的 `expr_matrix`/`mylist2`）
2. `mydata2` = 差异分析宽表（FC/FDR 等统计），只选择“Gene 到注释列之间”的列，排除 SubMito 注释列
3. 结果 = `mydata1` 的 `Gene + 数据列 + 原始注释列` 再 left_join `mydata2` 的 `FC/FDR` 列
4. 用途 = 作为 Step13 A/B 背景扣除的输入数据库；该表“不带 SubMito 注释”

**关键函数：**
```r
# SubMito转化
data_transformed <- apply_submito_transformation(
  data = diff_data,
  mitocarta_anno = mitocarta_anno,
  annotation_columns = NULL,  # NULL表示自动检测所有Localization列
  enable_submito = TRUE
)

# ROC分析
roc_analysis <- perform_roc_analysis(
  data = data_transformed,
  logfc_columns = logfc_columns,
  annotation_column = "GO_Localization",
  tp_label = "SGs",
  fp_label = "Matrix",
  direction = "<"
)

# 计算最佳阈值
thresholds <- calculate_optimal_thresholds(roc_dfs, min_tp = 0.3)

# 完整模块调用
result <- module08_roc_analysis1(
  dir_config, 
  diff_results1, 
  annotation_references,
  selected_versions = c("noNorm_Imputed", "Local_QNorm_Imputed"),
  roc_annotation_column = "GO_Localization",
  tp_label = "SGs",
  fp_label = "Matrix",
  enable_submito = TRUE,
  submito_annotation_columns = NULL,
  min_tp = 0.3
)
```

**SubMito转化逻辑：**
1. 从MitoCarta3.0提取Symbol和MitoCarta3.0_SubMitoLocalization
2. 过滤出四种亚定位类别：MIM（线粒体内膜）, Matrix（基质）, MOM（线粒体外膜）, IMS（膜间隙）
3. 通过left_join将亚定位信息添加到数据中
4. 对所有Localization列，如果Gene在MitoCarta3中有亚定位信息，则替换原有的"Mitochondrion"值
5. 移除辅助列，保持数据整洁

**ROC分析逻辑：**
1. **数据过滤**：只保留TP（如SGs）和FP（如Matrix）两类样本
2. **样本量检查**：TP和FP各至少5个样本，否则跳过
3. **pROC调用**：
   - `levels = c(FP_label, TP_label)` - 设置阴性和阳性标签
   - `direction = "<"` - 值越小越可能是FP
4. **提取ROC数据**：
   - Threshold（阈值）
   - TP（True Positive Rate，即Sensitivity）
   - FP（False Positive Rate，即1-Specificity）
   - TP_FP（Youden's Index = TP - FP）
5. **最佳阈值计算**：
   - 找出TP_FP最大的行（Youden's Index最大）
   - 如果这些行中所有TP < 0.3，则设阈值为0（不满足要求）
   - 否则，在TP >= 0.3的候选行中，找出最小的Threshold

**ROC曲线图说明：**
- **X轴**：FPR（False Positive Rate，假阳性率）
- **Y轴**：TPR（True Positive Rate，真阳性率，即Sensitivity）
- **AUC**：Area Under Curve，曲线下面积（0.5-1.0，越大越好）
- **最佳阈值**：在ROC曲线上标注，通常是Youden's Index最大的点

**Youden Index图说明：**
- **X轴**：Log2 FC（对比组的logFC值）
- **Y轴**：TPR - FPR（Youden's Index）
- **最优点**：黑色菱形标注，对应最大的Youden's Index
- **垂直线**：灰色虚线，标注最佳阈值位置

**输出文件结构：**

**ROC_data.xlsx（多sheet）：**
- **roc1** sheet：第一个对比组的ROC数据
  - 列：Comparison1_logFC_Threshold, Comparison1_logFC_TP, Comparison1_logFC_FP, Comparison1_logFC_TP_FP
- **roc2** sheet：第二个对比组的ROC数据
- ...

**all_thresholds.xlsx（多sheet）：**
- **version1** sheet：第一个数据版本的最佳阈值
  - 列：Comparison（对比组）, Threshold（最佳阈值）
- **version2** sheet：第二个数据版本的最佳阈值
- ...

**测试方法：**
1. 完成Module 1-7
2. 运行Module 8
3. 检查控制台输出，确认：
   - SubMito转化完成（显示转化的注释列数量）
   - MitoCarta3亚定位类别（MIM, Matrix, MOM, IMS）
   - ROC分析结果（TP和FP数量、AUC值）
   - 最佳阈值
4. 检查Output/目录输出文件：
   - PDF文件（ROC曲线和Youden Index图）
   - XLSX文件（ROC数据和阈值表）
   - CSV文件（SubMito转化后的数据）
5. 检查`Module08_workspace.RData`是否生成
6. 验证`roc_results1`列表结构：
   ```r
   names(roc_results1)  # 应包含所有分析的版本
   names(roc_results1[[1]])  # 应包含: roc_objects, roc_dfs, thresholds, data_transformed
   ```
7. 验证SubMito转化：
   ```r
   # 查看转化后的注释列
   unique(data_with_submito$noNorm_Imputed$GO_Localization)
   # 应包含：SGs, Matrix, MIM, MOM, IMS, Nuclear, Cytosol, Other等
   ```
8. 查看ROC曲线AUC值：
   ```r
   # 查看所有ROC对象的AUC
   lapply(roc_results1$noNorm_Imputed$roc_objects, pROC::auc)
   ```
9. 查看最佳阈值：
   ```r
   # 查看各版本的最佳阈值
   roc_thresholds1
   ```

**配置说明：**
在`main_pipeline.R`中配置ROC分析参数：
```r
# 选择分析版本（NULL表示所有版本）
selected_versions_roc <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# ROC分析参数
roc_annotation_column <- "GO_Localization"  # 用于ROC分析的注释列
tp_label <- "SGs"                            # True Positive标签
fp_label <- "Matrix"                         # False Positive标签（SubMito转化后）
enable_submito <- TRUE                       # 启用SubMito转化
submito_annotation_columns <- NULL           # 自动检测所有Localization列
min_tp <- 0.3                                # 最小TP阈值
```

**注意事项：**
- **必须先完成Module 3**（需要annotation_references，包含MitoCarta注释）
- **必须先完成Module 7**（需要diff_results1，包含logFC数据）
- **SubMito转化默认启用**，如果不需要可设置`enable_submito = FALSE`
- **注释列自动检测**：如果`submito_annotation_columns = NULL`，会自动检测所有包含"Localization"的列
- **ROC分析要求**：TP和FP各至少5个样本，否则跳过该版本的分析
- **最佳阈值筛选**：要求TP >= 0.3（即True Positive Rate >= 30%），避免漏检太多真阳性
- **需要安装pROC包**：`install.packages("pROC")`
- **需要安装openxlsx包**：`install.packages("openxlsx")`
- **文件命名**：所有输出文件带Module08_前缀和版本名称
- RData文件保存在工作目录，PDF/XLSX/CSV输出到Output/目录
- **FP标签选择**：第一次ROC分析通常使用"Matrix"作为FP（线粒体基质），因为它是非SG的参考定位
- **direction参数**：`direction = "<"`表示logFC值越小越可能是FP，这适用于SG富集分析（SG蛋白logFC通常较高）

**技术细节：**

**Youden's Index（约登指数）：**
- **定义**：J = Sensitivity + Specificity - 1 = TPR - FPR
- **范围**：0-1，越大越好
- **最优阈值**：使Youden's Index最大的阈值
- **生物学意义**：在保证一定灵敏度（不漏掉太多真阳性）的前提下，最大化特异性（减少假阳性）

**ROC曲线统计学：**
- **AUC解释**：
  - AUC = 0.5：随机猜测，无区分能力
  - AUC = 0.6-0.7：较差的区分能力
  - AUC = 0.7-0.8：可接受的区分能力
  - AUC = 0.8-0.9：优秀的区分能力
  - AUC = 0.9-1.0：极好的区分能力
- **置信区间**：pROC可以计算AUC的95%置信区间（需要设置`ci = TRUE`）

**SubMito转化的必要性：**
- **原因**：线粒体蛋白（Mitochondrion）包含多种亚定位，如基质（Matrix）、内膜（MIM）、外膜（MOM）、膜间隙（IMS）
- **目的**：区分不同亚定位的线粒体蛋白，提高ROC分析的特异性
- **第一次ROC分析**：使用Matrix作为FP，因为Matrix蛋白不应富集在SG中
- **第二次ROC分析**：通常使用Cytosol作为FP，进一步区分SG和细胞质蛋白

---

### Module 9: 背景扣除

**功能：**
1. **自动生成两种方法**：同时生成Method A和Method B两种过滤结果
   - **Method A**：对指定comparison不使用FDR过滤（FDR设为1）
   - **Method B**：对所有comparison使用FDR过滤
2. **阈值过滤**：根据Module 8的ROC最佳阈值，对每个bioGroup单独进行数据过滤
3. **分组处理**（基于sampleGroup的Context列）：
   - Context为"Experiment"的bioGroup：应用logFC和FDR阈值过滤
   - Context为"Spatial"的bioGroup：不使用阈值过滤（直接选列和去除NA）
   - Context为"Control"的bioGroup：不操作（跳过）
4. **催化组过滤**：确保每个蛋白至少有指定数量的有效LFQ值
5. **统计汇总**：按注释列（如GO_Localization）分组，计算count, mean, sum等统计指标
6. **数据合并**：生成merged数据（full_join所有bioGroup的过滤数据）
7. **可视化**：绘制堆积条形图（数量比例和丰度比例），展示不同定位的分布

**输入：**
- `Module08_workspace.RData` - Module01-08所有环境变量（包含diff_results1, roc_thresholds1, comparison_info等）
- `sampleGroup` - 必须包含Context列（标识Experiment/Spatial/Control）
- `expr_fdr_df_list`（可选，推荐）- 来自Module 8，用作过滤底表（Gene+数据列+原始注释+FC/FDR，不含SubMito注释）
- `data_with_submito`（可选）- 来自Module 8，仅在未提供 `expr_fdr_df_list` 时作为备选输入

**输出：**
- `Module09_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-08所有 + filtered_data_A, filtered_data_B, merged_data_A, merged_data_B
- `Output/Module09_{method}_{version}_FilteredData.xlsx` - 过滤后的数据（多sheet）
  - 先Exp组（数据sheet + _Summarise汇总sheet）
  - 再Spatial组（数据sheet + _Summarise汇总sheet）
- `Output/Module09_{method}_{version}_AfterROC_data_summarise.xlsx` - 专用于检查的工作簿：
  - sheet 顺序：先所有 Exp“数据” → 所有 Spatial“数据” → 所有 Exp“_Summarise” → 所有 Spatial“_Summarise”
- `Output/Module09_{method}_{version}_AfterROC_data_merged.xlsx` - 合并后的 AfterROC_merged（单sheet）
- `Output/Module09_{method}_{version}_CountPercent_StackBar.pdf` - 数量比例堆积条形图
- `Output/Module09_{method}_{version}_AbundancePercent_StackBar.pdf` - 丰度比例堆积条形图

**关键函数：**
```r
# 完整模块调用（自动生成Method A和Method B两种结果）
result <- module09_background_subtraction(
  dir_config,
  sampleGroup,  # 必须包含Context列
  diff_results1,
  roc_thresholds1,
  comparison_info,
  data_with_submito,   # 可选
  expr_fdr_df_list,    # 可选，优先使用
  selected_versions = c("noNorm_Imputed", "Local_QNorm_Imputed"),
  fdr_threshold = 0.05,
  no_fdr_comparisons = c("K69A1B3_Light_vs_A1B3_Light"),  # Method A专用，NULL表示Method A=Method B
  use_roc_threshold = TRUE,
  fixed_fc_threshold = NULL,  # 如果use_roc_threshold=FALSE，设置如2.0
  min_valid_lfq = 2,
  annotation_column = "GO_Localization",
  plot_all_annotations = FALSE,  # FALSE只显示TP，TRUE显示所有注释
  tp_label = "SGs",
  tp_color = "#DB6968"
)

# 提取结果
filtered_data_A <- result$filtered_data_A  # Method A的过滤数据
filtered_data_B <- result$filtered_data_B  # Method B的过滤数据
merged_data_A <- result$merged_data_A      # Method A的合并数据
merged_data_B <- result$merged_data_B      # Method B的合并数据
```

**使用示例：**
```r
# 示例1：标准用法（同时生成Method A和Method B）
result <- module09_background_subtraction(
  dir_config, sampleGroup, diff_results1, roc_thresholds1, comparison_info,
  selected_versions = c("noNorm_Imputed"),
  fdr_threshold = 0.05,
  no_fdr_comparisons = c("K69A1B3_Light_vs_A1B3_Light", "K69C3_Light_vs_C3_Light"),  # 多个comparison
  use_roc_threshold = TRUE,
  min_valid_lfq = 2,
  annotation_column = "GO_Localization",
  plot_all_annotations = FALSE,
  tp_label = "SGs"
)

# 示例2：不指定no_fdr_comparisons（Method A和Method B完全相同）
result <- module09_background_subtraction(
  dir_config, sampleGroup, diff_results1, roc_thresholds1, comparison_info,
  selected_versions = c("noNorm_Imputed"),
  fdr_threshold = 0.05,
  no_fdr_comparisons = NULL,  # Method A = Method B
  use_roc_threshold = TRUE,
  min_valid_lfq = 2,
  annotation_column = "GO_Localization",
  plot_all_annotations = FALSE,
  tp_label = "SGs"
)

# 示例3：使用固定FC阈值而非ROC阈值
result <- module09_background_subtraction(
  dir_config, sampleGroup, diff_results1, roc_thresholds1, comparison_info,
  selected_versions = c("noNorm_Imputed"),
  fdr_threshold = 0.05,
  no_fdr_comparisons = NULL,
  use_roc_threshold = FALSE,
  fixed_fc_threshold = 2.0,  # 使用固定阈值2.0
  min_valid_lfq = 2,
  annotation_column = "GO_Localization",
  plot_all_annotations = FALSE,
  tp_label = "SGs"
)

# 示例4：显示所有注释（而非只显示TP）
result <- module09_background_subtraction(
  dir_config, sampleGroup, diff_results1, roc_thresholds1, comparison_info,
  selected_versions = c("noNorm_Imputed"),
  fdr_threshold = 0.05,
  no_fdr_comparisons = c("K69A1B3_Light_vs_A1B3_Light"),
  use_roc_threshold = TRUE,
  min_valid_lfq = 2,
  annotation_column = "GO_Localization",
  plot_all_annotations = TRUE,  # 显示所有注释
  tp_label = "SGs"
)
```

**背景扣除逻辑：**

**过滤流程：**
1. **bioGroup识别**：从sampleGroup中提取所有bioGroup
2. **Context判断**：从sampleGroup的Context列确定每个bioGroup的类型（Experiment/Spatial/Control）
3. **bioGroup排序**：先Exp组，再Spatial组，Control组跳过
4. **过滤执行**：
   - **Experiment组**：
     - 找出该bioGroup对应的所有comparisons（从comparison_info中）
     - 对每个comparison：
       - 应用logFC阈值（ROC阈值或固定阈值）
       - 应用FDR阈值：
         - Method A：如果comparison在no_fdr_comparisons中，FDR=1（不过滤）
         - Method B：所有comparison都用fdr_threshold
     - 选择列：Gene, LFQ列, logFC列, FDR列, 所有注释列
     - 催化组过滤：rowSums(!is.na(LFQ_cols)) >= min_valid_lfq
   - **Spatial组**：
     - 不应用阈值过滤
     - 选择列：Gene, LFQ列, 所有注释列
     - 去除NA：na.omit()
   - **Control组**：
     - 不操作（跳过）

**统计汇总逻辑：**
```r
# 按注释列分组统计
summarise_data <- filtered_data %>%
  group_by(annotation_column) %>%
  summarise(
    count = n(),
    # 每个LFQ列的mean和sum
    across(all_of(lfq_cols), 
           list(mean = ~mean(., na.rm = TRUE), 
                sum = ~sum(., na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  mutate(
    # 数量百分比
    Percent = count / sum(count),
    # 丰度百分比（每个LFQ列的sum百分比）
    across(ends_with("_sum"), 
           ~. / sum(., na.rm = TRUE), 
           .names = "{.col}_percent"),
    # MeanSum：所有sum_percent的平均值
    MeanSum = rowMeans(select(., ends_with("_sum_percent")), na.rm = TRUE)
  )
```

**堆积条形图特性：**
1. **X轴**：bioGroup（每个bioGroup一列）
2. **Y轴**：Percent（数量比例）或MeanSum（丰度比例）
3. **填充色**：按注释列分组
   - TP（如SGs）：红色（#DB6968），永远在底部
   - 其他注释：
     - `plot_all_annotations = TRUE`：使用调色板
     - `plot_all_annotations = FALSE`：灰色
4. **标注**：TP在顶部标注count和百分比
5. **标题和副标题**：
   - 标题：版本名称 + Method + FDR阈值
   - 副标题：每个bioGroup使用的阈值（每行一个comparison）
6. **图例**：显示注释类型和对应颜色

**结果列表结构（AfterROC_Summarise_list顺序）：**
1. 先放Exp组：
   - `bioGroup1` - 过滤后的数据
   - `bioGroup1_Summarise` - 统计汇总
   - `bioGroup2` - 过滤后的数据
   - `bioGroup2_Summarise` - 统计汇总
   - ...
2. 再放Spatial组：
   - `bioGroup_spatial1` - 过滤后的数据
   - `bioGroup_spatial1_Summarise` - 统计汇总
   - ...

**merged数据生成逻辑（AfterROC_merged）：**
1. 提取所有过滤后的数据（不包含_Summarise）
2. 使用`reduce(data_list, full_join, by="Gene")`合并所有bioGroup
3. 只保留Gene和所有LFQ列
4. 从原始combined_data中提取注释列
5. 使用`left_join`添加注释列到merged数据

**Method A vs Method B：**

| 特性 | Method A | Method B |
|------|----------|----------|
| 主要用途 | 对某些低置信度comparison放宽FDR要求 | 严格控制所有comparison的FDR |
| FDR过滤 | 部分comparison不用FDR（FDR设为1） | 所有comparison都用FDR |
| 配置参数 | 指定`no_fdr_comparisons`向量 | `no_fdr_comparisons`设为NULL |
| 适用场景 | 探索性分析，避免漏掉潜在候选 | 保守分析，确保高置信度结果 |
| 自动生成 | 每次运行都自动生成 | 每次运行都自动生成 |

**阈值显示示例（副标题）：**
```
K69A1B3: K69A1B3_Light_vs_K69A1B3_noLight: logFC>2.35, FDR<0.05; K69A1B3_Light_vs_K69_Light: logFC>1.87, FDR<0.05
K69C3: K69C3_Light_vs_K69C3_noLight: logFC>2.12, FDR<0.05; K69C3_Light_vs_K69_Light: logFC>1.65, no FDR filter
K73: No threshold (Spatial)
K20: No threshold (Control)
```

**测试方法：**
1. 完成Module 1-8
2. 确保sampleGroup包含Context列
3. 运行Module 9
4. 检查控制台输出，确认：
   - Experiment/Spatial/Control组的数量
   - Method A和Method B分别处理的结果
   - 每个bioGroup的过滤结果（保留蛋白数量）
   - 每个bioGroup的最小logFC值（验证过滤成功）
   - 统计汇总的分组数量
   - merged数据的蛋白数量和列数
5. 检查Output/目录输出文件：
   - `Module09_A_{version}_FilteredData.xlsx` - Method A结果
   - `Module09_B_{version}_FilteredData.xlsx` - Method B结果
   - `Module09_A|B_{version}_AfterROC_data_summarise.xlsx` - 检查版工作簿（固定sheet顺序：Exp数据→Spatial数据→Exp_Summarise→Spatial_Summarise）
   - `Module09_A|B_{version}_AfterROC_data_merged.xlsx` - 合并数据
   - PDF文件（两种方法各有数量和丰度两个图）
6. 检查`Module09_workspace.RData`是否生成
7. 验证结果列表结构：
   ```r
   # Method A结果
   names(filtered_data_A)  # 应包含所有分析的版本
   names(filtered_data_A$noNorm_Imputed)  # 应包含所有bioGroup和_Summarise
   head(filtered_data_A$noNorm_Imputed$K69A1B3)  # 查看过滤后的数据
   filtered_data_A$noNorm_Imputed$K69A1B3_Summarise  # 查看统计汇总
   
   # Method B结果
   names(filtered_data_B)  # 应包含所有分析的版本
   names(filtered_data_B$noNorm_Imputed)  # 应包含所有bioGroup和_Summarise
   
   # Merged数据
   head(merged_data_A$noNorm_Imputed)  # 查看合并后的数据
   colnames(merged_data_A$noNorm_Imputed)  # 应包含Gene, 所有LFQ列, 注释列
   ```
8. 检查过滤是否成功：
   ```r
   # 查看最小logFC值（应接近或大于阈值）
   min(filtered_data_A$noNorm_Imputed$K69A1B3$K69A1B3_Light_vs_K69A1B3_noLight_logFC, na.rm = TRUE)
   
   # 查看FDR值分布
   summary(filtered_data_A$noNorm_Imputed$K69A1B3$K69A1B3_Light_vs_K69A1B3_noLight_adj.P.Val)
   ```
9. 比较Method A和Method B的差异：
   ```r
   # Method A保留的蛋白数量
   nrow(filtered_data_A$noNorm_Imputed$K69A1B3)
   
   # Method B保留的蛋白数量
   nrow(filtered_data_B$noNorm_Imputed$K69A1B3)
   
   # 差异（Method A通常>=Method B，因为部分comparison不用FDR）
   nrow(filtered_data_A$noNorm_Imputed$K69A1B3) - nrow(filtered_data_B$noNorm_Imputed$K69A1B3)
   ```
10. 检查结果列表顺序：
    ```r
    # 应该是：Exp组（数据+汇总）在前，Spatial组（数据+汇总）在后
    names(filtered_data_A$noNorm_Imputed)
    ```

**配置说明：**
在`main_pipeline.R`中配置背景扣除参数：
```r
# 选择分析版本（NULL表示所有版本）
selected_versions_bg <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# 背景扣除参数（自动生成Method A和Method B）
fdr_threshold <- 0.05
# Method A：指定不使用FDR过滤的comparison向量（NULL表示Method A = Method B）
no_fdr_comparisons <- c("K69A1B3_Light_vs_A1B3_Light", "K69C3_Light_vs_C3_Light")
use_roc_threshold <- TRUE
fixed_fc_threshold <- NULL  # 如果use_roc_threshold=FALSE，设置固定阈值如2.0
min_valid_lfq <- 2
annotation_column <- "GO_Localization"
plot_all_annotations <- FALSE  # FALSE只显示TP，TRUE显示所有注释
tp_label <- "SGs"
tp_color <- "#DB6968"
```

**注意事项：**
- **必须先完成Module 8**（需要roc_thresholds1和comparison_info）
- **必须先完成Module 7**（需要diff_results1和comparison_info）
- **sampleGroup必须包含Context列**：标识每个bioGroup的类型（Experiment/Spatial/Control）
- **自动生成两种方法**：每次运行都会同时生成Method A和Method B结果
- **no_fdr_comparisons**：支持向量输入，可同时指定多个comparison不使用FDR过滤
- **Method A vs Method B**：Method A更宽松（探索性），Method B更严格（保守）
- **固定阈值vs ROC阈值**：推荐使用ROC阈值（数据驱动），固定阈值仅用于特殊情况
- **min_valid_lfq**：默认2，确保催化组有足够的有效重复
- **annotation_column**：必须在数据中存在，用于分组统计和作图
- **plot_all_annotations**：FALSE时图更清晰，只突出显示TP；TRUE时显示完整的定位分布
- **Spatial组**：从Context列识别，不使用阈值过滤，只选择列和去除NA
- **Control组**：从Context列识别，不操作（跳过）
- **过滤验证**：通过检查最小logFC值确认过滤成功（应接近阈值）
- **结果列表顺序**：先Exp组（数据+汇总），再Spatial组（数据+汇总）
- **需要安装ggplot2包**：`install.packages("ggplot2")`
- **需要安装dplyr包**：`install.packages("dplyr")`
- **需要安装openxlsx包**：`install.packages("openxlsx")`
- **文件命名**：所有输出文件带Module09_前缀、方法（A/B）和版本名称
- RData文件保存在工作目录，XLSX/PDF输出到Output/目录

**技术细节：**

**催化组有效值过滤的必要性：**
- **原因**：LFQ定量可能存在缺失值（missing values），如果一个蛋白在大部分重复中缺失，其定量不可靠
- **min_valid_lfq = 2**：要求至少2个重复有有效值（假设3个重复的实验）
- **实现**：`rowSums(!is.na(select(., all_of(lfq_cols)))) >= min_valid_lfq`

**Spatial组的特殊处理：**
- **原因**：Spatial组通常是特定亚细胞定位的参考样本（如K20, K73），不用于差异分析
- **处理方式**：不使用logFC和FDR阈值，直接选择Gene, LFQ列和注释列
- **数据完整性**：使用`na.omit()`确保数据完整

**丰度百分比计算逻辑：**
1. 对每个LFQ列，计算每个注释类别的sum
2. 对每个LFQ列，计算sum的百分比：`sum / total_sum`
3. 计算所有LFQ列的百分比平均值：`MeanSum = mean(all_sum_percents)`
4. **生物学意义**：反映不同定位蛋白的丰度贡献（而非数量贡献）

**阈值信息的提取和显示：**
- 从`comparison_info`中获取每个bioGroup对应的comparisons
- 从`roc_thresholds1`中获取每个comparison的最佳阈值
- 根据`filter_method`和`no_fdr_comparisons`确定FDR文本
- 在堆积条形图的副标题中显示所有阈值信息

---

### Module 10: 数据替换 ✓

**文件：** `Module/module10_data_replacement.R`

**功能：**
1. 对Module 9生成的**merged数据**（AfterROC_merged）进行定点替换
2. 替换逻辑：保留所有NA值，所有非NA值用Module 5的Global_QNorm_Imputed对应值替换
3. 如果没有Global_QNorm_Imputed，使用Global_MNorm_Imputed代替
4. 生成抽样校验报告，验证替换是否成功
5. Check文件单独存放在**Output/Check/**子目录

**输入：**
- `Module09_workspace.RData` - Module01-09所有环境变量（包含merged_data_A, merged_data_B等）
- `imputed_data_list` - Module 5的全局标准化+填补结果
- `sampleGroup` - 用于识别LFQ列
- `selected_versions` - 需要处理的版本

**输出：**
- `Module10_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-09所有 + replaced_data_A, replaced_data_B
- `Output/Module10_{method}_{version}_AfterFirstROC_Intersect.xlsx` - 替换后的数据
- `Output/Check/Module10_{method}_{version}_replacement_check.csv` - 抽样校验表
- `Output/Check/Module10_{method}_{version}_replacement_summary.txt` - 替换汇总

**关键函数：**
```r
# 完整模块调用
result <- module10_data_replacement(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  merged_data_A = merged_data_A,        # Module 9的merged数据
  merged_data_B = merged_data_B,        # Module 9的merged数据
  imputed_data_list = imputed_data_list,
  selected_versions = selected_versions_bg,
  prefer_version = "Global_QNorm_Imputed",
  fallback_version = "Global_MNorm_Imputed",
  test_n = 20
)

# 提取结果
replaced_data_A <- result$replaced_data_A  # Method A的替换数据
replaced_data_B <- result$replaced_data_B  # Method B的替换数据
```

**替换逻辑（对齐CleanCode.R 2303-2361）：**
```r
# 对每个数据单元格：
if (is.na(原始值)) {
  保留 NA
} else {
  替换为 Global_QNorm_Imputed 的对应值
}
```

**详细流程：**
1. **选择基准版本**：
   - 优先使用`Global_QNorm_Imputed`
   - 如果不存在，使用`Global_MNorm_Imputed`
   - 都不存在则报错
   
2. **创建Check子目录**：
   - 在`Output/`下创建`Check/`子目录
   - 所有校验文件存放在该目录

3. **对每个版本的merged数据执行替换**：
   - Method A和Method B分别处理
   - 识别所有LFQ列（从sampleGroup）
   - 按Gene列对齐merged数据和基准数据
   - 执行单元格级替换（非NA替换，NA保留）

4. **抽样校验**：
   - 随机选择最多20个基因
   - 选择最多3个LFQ列
   - 对比原始值、替换后值、参考值
   - 判断CellMatch（非NA的替换后值应等于参考值）
   - 导出CSV到Check/子目录

5. **统计汇总**：
   - 记录替换的单元格数量
   - 记录保留的NA数量
   - 导出TXT到Check/子目录

6. **导出结果**：
   - 每个version的每个method生成一个XLSX
   - 文件名：`Module10_{method}_{version}_AfterFirstROC_Intersect.xlsx`

**使用示例：**

**示例1：标准用法（使用Global_QNorm_Imputed）**
```r
result <- module10_data_replacement(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  merged_data_A = merged_data_A,
  merged_data_B = merged_data_B,
  imputed_data_list = imputed_data_list,
  selected_versions = c("noNorm_Imputed", "Local_QNorm_Imputed"),
  prefer_version = "Global_QNorm_Imputed",
  fallback_version = "Global_MNorm_Imputed",
  test_n = 20
)
```

**示例2：使用Global_MNorm_Imputed作为首选**
```r
result <- module10_data_replacement(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  merged_data_A = merged_data_A,
  merged_data_B = merged_data_B,
  imputed_data_list = imputed_data_list,
  selected_versions = c("noNorm_Imputed"),
  prefer_version = "Global_MNorm_Imputed",  # 改为首选Median标准化
  fallback_version = "Global_QNorm_Imputed",
  test_n = 20
)
```

**示例3：增加校验样本数**
```r
result <- module10_data_replacement(
  dir_config = dir_config,
  sampleGroup = sampleGroup,
  merged_data_A = merged_data_A,
  merged_data_B = merged_data_B,
  imputed_data_list = imputed_data_list,
  selected_versions = c("noNorm_Imputed", "Local_QNorm_Imputed"),
  prefer_version = "Global_QNorm_Imputed",
  fallback_version = "Global_MNorm_Imputed",
  test_n = 50  # 增加到50个基因
)
```

**测试方法：**
1. 确保Module 9已完成，生成`Module09_workspace.RData`（包含merged_data_A和merged_data_B）
2. 确保Module 5已完成，生成`imputed_data_list`（包含Global标准化版本）
3. 运行Module 10
4. 检查控制台输出，确认：
   - 选择的基准版本（Global_QNorm_Imputed或Global_MNorm_Imputed）
   - 每个版本的原始数据规模
   - 替换的单元格数量
   - 保留的NA数量
5. 检查Output/目录输出文件：
   - `Module10_A_{version}_AfterFirstROC_Intersect.xlsx` - Method A结果
   - `Module10_B_{version}_AfterFirstROC_Intersect.xlsx` - Method B结果
6. 检查Output/Check/目录校验文件：
   - `Module10_A_{version}_replacement_check.csv` - 抽样校验
   - `Module10_A_{version}_replacement_summary.txt` - 替换汇总
   - Method B的对应文件
7. 验证替换成功：
   ```r
   # 检查校验文件
   check_A <- read.csv("Output/Check/Module10_A_noNorm_Imputed_replacement_check.csv")
   
   # 查看CellMatch列（TRUE表示替换成功）
   table(check_A$CellMatch)
   # 应该全部为TRUE（非NA的单元格）
   
   # 查看替换汇总
   readLines("Output/Check/Module10_A_noNorm_Imputed_replacement_summary.txt")
   ```
8. 验证NA保留：
   ```r
   # 读取替换后的数据
   library(openxlsx)
   replaced_A <- read.xlsx("Output/Module10_A_noNorm_Imputed_AfterFirstROC_Intersect.xlsx")
   
   # 检查NA数量（应与原始merged数据一致）
   sum(is.na(replaced_A))
   
   # 对比原始merged数据的NA数量
   sum(is.na(merged_data_A$noNorm_Imputed))
   ```
9. 验证非NA值替换：
   ```r
   # 读取基准数据
   base_data <- imputed_data_list$Global_QNorm_Imputed
   
   # 随机选择一个基因和一个LFQ列
   test_gene <- replaced_A$Gene[1]
   test_col <- "K69A1B3_Light_LFQ_1"
   
   # 如果merged数据中该单元格非NA，替换后的值应等于base_data中的值
   merged_val <- merged_data_A$noNorm_Imputed[merged_data_A$noNorm_Imputed$Gene == test_gene, test_col]
   replaced_val <- replaced_A[replaced_A$Gene == test_gene, test_col]
   base_val <- base_data[base_data$Gene == test_gene, test_col]
   
   if (!is.na(merged_val)) {
     # 应该相等
     print(all.equal(replaced_val, base_val))
   }
   ```
10. 检查`Module10_workspace.RData`是否生成

**配置说明：**
在`main_pipeline.R`中配置替换参数：
```r
# 优先版本和备选版本
prefer_version = "Global_QNorm_Imputed"
fallback_version = "Global_MNorm_Imputed"

# 抽样校验的基因数量
test_n = 20  # 可增加到50或更多
```

**注意事项：**
- **必须先完成Module 9**（需要merged_data_A和merged_data_B）
- **必须先完成Module 5**（需要imputed_data_list，包含Global标准化版本）
- **使用merged数据**：Module 10处理的是merged矩阵，不是filtered_data的individual bioGroup表
- **Check文件位置**：所有校验文件存放在`Output/Check/`子目录
- **替换逻辑**：只替换非NA值，保留所有NA值不变
- **NA保留验证**：替换前后NA的位置和数量应完全一致
- **基准版本选择**：优先使用Global_QNorm_Imputed，确保全局一致性
- **抽样校验**：默认检查20个基因×3个LFQ列，可通过test_n参数调整
- **CellMatch判断**：对非NA的原始值，替换后应等于参考值（容差1e-9）
- **需要安装dplyr包**：`install.packages("dplyr")`
- **需要安装openxlsx包**：`install.packages("openxlsx")`
- **文件命名**：所有输出文件带Module10_前缀、方法（A/B）和版本名称
- RData文件保存在工作目录，XLSX输出到Output/目录，Check文件到Output/Check/目录

**技术细节：**

**替换算法实现：**
```r
# 伪代码
for (每个LFQ列) {
  v_old <- merged数据[, 列]
  v_ref <- Global_QNorm_Imputed[, 列]
  
  for (每个单元格) {
    if (!is.na(v_old[i])) {
      v_new[i] <- v_ref[i]  # 替换为参考值
    } else {
      v_new[i] <- NA  # 保留NA
    }
  }
}
```

**对齐CleanCode.R实现：**
- CleanCode.R (2303-2361)：
  ```r
  replace_cols <- names(AfterROC_merged_list_AB[[1]])[2:16]
  ref_df <- noMBR_QNorm_Imputed %>% select(Gene, all_of(replace_cols))
  
  replace_nonNA_cellwise <- function(df, ref_df, cols) {
    merged <- left_join(df, ref_df, by = "Gene", suffix = c("", ".ref"))
    for (col in cols) {
      ref_col <- paste0(col, ".ref")
      merged[[col]] <- mapply(function(x, y) {
        if (!is.na(x)) y else x
      }, merged[[col]], merged[[ref_col]])
    }
    merged %>% select(-ends_with(".ref"))
  }
  ```
- 本模块使用相同逻辑，确保结果一致

**校验CSV列说明：**
- **Gene**：基因名
- **Sample**：样本列名
- **Original**：替换前的值（来自merged数据）
- **Replaced**：替换后的值（执行替换后）
- **Imputed**：参考值（来自Global_QNorm_Imputed）
- **CellMatch**：是否匹配（对非NA的Original，Replaced应等于Imputed）

**汇总TXT内容：**
```
版本：noNorm_Imputed - Method A
原始数据：1234 蛋白, 20 列
替换列数：15 LFQ列
替换结果：12345 个非NA单元格被替换
NA保留：6789/19035 单元格
基准版本：Global_QNorm_Imputed
```

---

### Module 11: 第二次差异分析 ✓

**文件：** `Module/module11_diff_analysis2.R`

**功能：**
1. 对Module 10的数据替换结果进行第二次差异分析
2. 只选择**Context为Experiment或Spatial**的样本（去除Control样本）
3. 根据Context和PLtype构建比较组：
   - **Exp vs Exp**：所有Experiment组之间的比较（不限PLtype）
   - **Exp vs Spatial**：Experiment组与Spatial组的比较（**必须PLtype相同**）
4. 使用limma进行差异分析，得到logFC和adj.P.Val
5. 生成FDR_combined_df_list_2nd（包含logFC、adj.P.Val和注释列）

**输入：**
- `Module10_workspace.RData` - Module01-10所有环境变量
- `replaced_data_A/B` - Module 10的数据替换结果
- `sampleGroup` - 用于识别Context和SecondROCgroup
- `selected_versions` - 需要处理的版本

**输出：**
- `Module11_workspace.RData` - 所有环境变量（工作目录）
  - 包含：Module01-10所有 + FDR_combined_df_list_2nd, sample_info_2nd等
- `Output/Module11_Raw_FDR_test_list_{version}_{method}.xlsx` - 原始topTable结果（每个比较一个sheet）
- `Output/Module11_FC_FDR_2nd_combined.xlsx` - 合并的logFC、adj.P.Val和注释列（所有版本和方法）

**关键对象：**
- `FDR_combined_df_list_2nd` - 包含所有版本的差异分析结果（命名格式：`{version}_{method}`）
- `sample_info_2nd` - 第二次分析的样本分组信息
- `comparisons_list_2nd` - 比较组列表
- `bioGroups_selected_2nd` - 选择的bioGroup

**差异分析流程：**
1. **样本筛选**：从sampleGroup中提取Context为Experiment或Spatial的样本
2. **比较矩阵构建**：
   - 提取每个bioGroup的Context和SecondROCgroup
   - 构建Exp vs Exp和Exp vs Spatial的比较组
3. **limma分析**：
   - 创建设计矩阵（`design <- model.matrix(~ 0 + Group)`）
   - 拟合线性模型（`lmFit`）
   - 构建对比矩阵（`makeContrasts`）
   - 经验贝叶斯调整（`eBayes`）
   - 提取topTable结果（`topTable(..., adjust.method = "BH")`）
4. **结果合并**：
   - 提取每个比较的Gene, logFC, adj.P.Val
   - 使用full_join按Gene合并所有比较
   - 添加注释列（最后3列）

**配置参数（在main_pipeline.R中）：**
```r
result <- module11_diff_analysis2(
  dir_config = dir_config,              # 目录配置
  sampleGroup = sampleGroup,            # 样本分组表
  replaced_data_A = replaced_data_A,    # Method A数据
  replaced_data_B = replaced_data_B,    # Method B数据
  selected_versions = selected_versions_bg  # 要处理的版本（NULL表示全部）
)
```

**比较组构建逻辑（示例）：**

假设有以下样本：
- K69A1B3_Light (Experiment, PLtype=Light, SecondROCgroup=A)
- K69C3_Light (Experiment, PLtype=Light, SecondROCgroup=B)
- C3_H2O2 (Experiment, PLtype=H2O2, SecondROCgroup=C)
- K20_Light (Spatial, PLtype=Light, SecondROCgroup=A&B)
- K73_H2O2 (Spatial, PLtype=H2O2, SecondROCgroup=空)

构建的比较组（6个）：
1. **K69A1B3_vs_K69C3** (Exp vs Exp, 都是Light) ✓
2. **K69A1B3_vs_C3** (Exp vs Exp, Light vs H2O2) ✓
3. **K69C3_vs_C3** (Exp vs Exp, Light vs H2O2) ✓
4. **K69A1B3_vs_K20** (Exp vs Spatial, 都是Light) ✓
5. **K69C3_vs_K20** (Exp vs Spatial, 都是Light) ✓
6. **C3_vs_K73** (Exp vs Spatial, 都是H2O2) ✓

**说明：**
- Exp vs Exp：所有3个Experiment组两两比较，共3对（不限PLtype）
- Exp vs Spatial：只比较PLtype相同的组
  - K69A1B3(Light) 和 K69C3(Light) 可以和 K20(Light) 比较
  - C3(H2O2) 可以和 K73(H2O2) 比较

**输出示例：**

FDR_combined_df_list_2nd结构（每个版本一个数据框）：
```
$noNorm_Imputed_A
  Gene  K69A1B3_vs_K69C3_logFC  K69A1B3_vs_K69C3_adj.P.Val  K69A1B3_vs_C3_logFC  K69A1B3_vs_C3_adj.P.Val  ...  MainLocalization  MitoGene  MainFunction
  P12345  1.23  0.001  0.87  0.012  ...  Cytoplasm  No  Transport
  Q67890  -0.45  0.023  -0.32  0.045  ...  Nucleus  No  Signaling
  ...
```

**列结构：**
- 第1列：Gene（基因名）
- 第2-13列：6个比较的logFC和adj.P.Val（共12列）
  - K69A1B3_vs_K69C3_logFC / adj.P.Val
  - K69A1B3_vs_C3_logFC / adj.P.Val
  - K69C3_vs_C3_logFC / adj.P.Val
  - K69A1B3_vs_K20_logFC / adj.P.Val
  - K69C3_vs_K20_logFC / adj.P.Val
  - C3_vs_K73_logFC / adj.P.Val
- 最后3列：注释列（MainLocalization, MitoGene, MainFunction）

**技术细节：**
1. **设计矩阵**：使用`~ 0 + Group`，列名直接对应bioGroup名称
2. **对比矩阵**：使用makeContrasts动态构建，处理特殊字符（`make.names`）
3. **topTable**：使用Benjamini-Hochberg方法进行FDR校正
4. **Excel导出**：Sheet名称最多31个字符，自动截断

**使用示例：**
```r
# 运行Module 11
source("main_pipeline.R")

# 查看比较组
names(comparisons_list_2nd)

# 查看某个版本的差异分析结果
head(FDR_combined_df_list_2nd$noNorm_Imputed_A)

# 筛选显著差异蛋白（以第一个比较为例）
comp_name <- names(comparisons_list_2nd)[1]
logFC_col <- paste0(comp_name, "_logFC")
fdr_col <- paste0(comp_name, "_adj.P.Val")

significant <- FDR_combined_df_list_2nd$noNorm_Imputed_A %>%
  filter(abs(.data[[logFC_col]]) > 1 & .data[[fdr_col]] < 0.05)
```

**调试指南：**
```r
# 加载Module11环境
load("Module11_workspace.RData")

# 检查选择的样本
print(sample_info_2nd)

# 检查比较组
print(comparisons_list_2nd)

# 检查bioGroup选择
print(bioGroups_selected_2nd)

# 查看某个版本的结果
View(FDR_combined_df_list_2nd$noNorm_Imputed_A)
```

**故障排除：**
1. **错误：缺少样本列** → 检查replaced_data中的列名是否与sampleGroup$FinalName匹配
2. **错误：没有Experiment或Spatial样本** → 检查sampleGroup的Context列是否正确填写
3. **警告：跳过版本** → 某些版本缺少必需的样本列，检查数据完整性

**测试结果示例：**
```
=== Module 11: 第二次差异分析 ===

✓ 准备分析 4 个数据版本
  版本：noNorm_Imputed_A, noNorm_Imputed_B, Local_QNorm_Imputed_A, Local_QNorm_Imputed_B

--- 步骤 1: 提取样本分组信息 ---
✓ 找到 15 个 Experiment/Spatial 样本
✓ 选择的 bioGroup (5):
  1. K69A1B3_Light (Context: Experiment, SecondROCgroup: A)
  2. K69C3_Light (Context: Experiment, SecondROCgroup: B)
  3. K20_Light (Context: Spatial, SecondROCgroup: A&B)
  4. C3_H2O2 (Context: Experiment, SecondROCgroup: C)
  5. K73_H2O2 (Context: Spatial, SecondROCgroup: 空)
✓ 构建 sample_info，共 15 个样本

--- 步骤 2: 构建比较矩阵 ---
✓ Experiment bioGroups: K69A1B3_Light, K69C3_Light, C3_H2O2
✓ Spatial bioGroups: K20_Light, K73_H2O2

✓ 构建 6 个比较组：
  1. K69A1B3_Light_vs_K69C3_Light
  2. K69A1B3_Light_vs_C3_H2O2
  3. K69C3_Light_vs_C3_H2O2
  4. K69A1B3_Light_vs_K20_Light
  5. K69C3_Light_vs_K20_Light
  6. C3_H2O2_vs_K73_H2O2

--- 步骤 3: 执行差异分析 ---

[1/4] 处理：noNorm_Imputed_A
  ✓ 提取数据矩阵：1234 genes × 15 samples
  ✓ 创建设计矩阵
  ✓ 拟合线性模型
  ✓ 创建对比矩阵
  ✓ 经验贝叶斯调整
  ✓ 提取 6 个比较的差异分析结果
  ✓ 导出原始 topTable：Module11_noNorm_Imputed_A_Raw_topTable.xlsx
  ✓ 合并结果：1234 genes × 15 columns

...

--- 步骤 4: 导出合并结果 ---
✓ 导出合并的差异分析结果：
  Module11_FC_FDR_2nd_combined.xlsx
  包含 4 个数据版本

=== Module 11 完成 ===
✓ 分析了 4 个数据版本
✓ 构建了 6 个比较组
✓ 输出文件数量：4 个原始 topTable + 1 个合并结果
```

---

### Module 12: 第二次ROC分析(见后面)
### Module 13: SubSG注释(见后面)
### Module 14: 火山图(见后面)

## 文件命名规范总结

### 工作目录文件
- **环境数据**：`ModuleXX_workspace.RData` - 累积保存所有环境变量
- **分组表**：
  - `Module02_sampleGroup_template.csv` - 用户填写
  - `Module02_sampleGroup.csv` - 程序生成

### Output目录文件
- **CSV**：`ModuleXX_功能描述.csv`
- **PDF**：`ModuleXX_功能描述.pdf`

### Module目录文件
- **模块代码**：`moduleXX_功能描述.R`


## 配置接口总结

### 数据读取配置（Module 2）
在 `main_pipeline.R` 中修改：
```r
data_config <- list(
  file_pattern = "_matrix.*\\.tsv$"  # TSV文件名匹配模式
)
```

### 注释配置（Module 3）
可添加自定义注释：
```r
custom_annotations <- list(
  list(
    column_name = "Custom_Localization",
    TP_source = "GO_SGs",
    TP_column = "Gene",
    TP_label = "SGs"
  )
)
```

### 标准化配置（Module 4）
选择需要的标准化方法：
```r
# 必须至少包含一个全局标准化（Global_QNorm或Global_MNorm）
norm_types <- c("noNorm", "Global_QNorm", "Global_MNorm", "Local_QNorm", "Local_MNorm")

# 或只使用部分方法
norm_types <- c("noNorm", "Global_QNorm")
```

### 缺失值填补配置（Module 5）
配置填补策略和随机种子：
```r
# impute_cat_mean: 是否对Cat组的n_valid=2情况用平均值填补
impute_cat_mean <- FALSE  # 默认不填补Cat组

# random_seed: 随机数种子，保证Perseus填补可重复
random_seed <- 123        # 可改为其他数值以获得不同的随机填补结果

# 启用Cat组填补
impute_cat_mean <- TRUE
```

### 热图系统配置（Module 6）
配置热图版本、类型和样式：
```r
# 选择标准化版本（支持向量，会为每个版本生成独立热图）
selected_versions <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# 选择热图类型
heatmap_types <- c("all", "correlation", "by_localization")

# 相关性热图配置（高度可定制）
correlation_config <- list(
  # 相关性范围和颜色
  corr_min = 0.8,                # 最低相关性值
  corr_max = 1,                  # 最高相关性值
  corr_center = 0.9,             # 中心值
  col_low = "#4D97CD",           # 低值颜色（蓝色）
  col_center = "white",          # 中心颜色（白色）
  col_high = "#DB6968",          # 高值颜色（红色）
  
  # 方式1：使用排除规则（排除不需要的样本）
  exclude_context = NULL,        # 例如：c("Control") 排除所有Control组
  exclude_pltype = NULL,         # 例如：c("PL") 排除所有PL组
  exclude_catalytic = c("NoCat"), # 例如：c("NoCat") 排除所有NoCat组
  
  # 方式2：直接指定包含的bioGroup（优先于排除规则）
  include_biogroups = NULL       # 例如：c("Light_Experiment", "H2O2_Experiment")
)

# Localization列（NULL=自动检测）
localization_columns <- NULL

# 表达热图颜色参数
color_params <- list(
  custom_min = -4,               # 低值阈值
  custom_max = 4,                # 高值阈值
  custom_center = 0,             # 中心值
  col_low = "#4D97CD",           # 低值颜色（蓝色）
  col_center = "white",          # 中心颜色（白色）
  col_high = "#DB6968"           # 高值颜色（红色）
)
```

### 差异分析配置（Module 7）
配置分析版本（对比组自动从FirstROCgroup构建）：
```r
# 选择分析版本（NULL表示所有版本）
selected_versions_diff <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# 对比组由sampleGroup$FirstROCgroup自动构建
# 确保sampleGroup包含以下必需列：
# - FirstROCgroup: A, B, C, A&B, A/B等（不参与的留空或NA）
# - Context: Experiment, Control, Spatial
# - bioGroup: 每个样本组的唯一标识
```

### ROC分析配置（Module 8）
配置ROC分析参数：
```r
# 选择分析版本（NULL表示所有版本）
selected_versions_roc <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# ROC分析参数
roc_annotation_column <- "GO_Localization"  # 用于ROC分析的注释列
tp_label <- "SGs"                            # True Positive标签
fp_label <- "Matrix"                         # False Positive标签（SubMito转化后）
enable_submito <- TRUE                       # 启用SubMito转化
submito_annotation_columns <- NULL           # 自动检测所有Localization列
min_tp <- 0.3                                # 最小TP阈值（保证至少30%灵敏度）

# 修改TP标签（例如使用HaloMap注释）
roc_annotation_column <- "HaloMap_Localization"
tp_label <- "SGs"
fp_label <- "Matrix"

# 禁用SubMito转化（保留原始Mitochondrion注释）
enable_submito <- FALSE

# 指定特定注释列进行SubMito转化
submito_annotation_columns <- c("GO_Localization", "HaloMap_Localization")
```

### 背景扣除配置（Module 9）
配置过滤方法和参数（自动生成Method A和Method B）：
```r
# 选择分析版本（NULL表示所有版本）
selected_versions_bg <- c("noNorm_Imputed", "Local_QNorm_Imputed")

# 背景扣除参数
fdr_threshold <- 0.05

# Method A：指定不使用FDR过滤的comparison向量（NULL表示Method A = Method B）
no_fdr_comparisons <- c("K69A1B3_Light_vs_A1B3_Light", "K69C3_Light_vs_C3_Light")

# 阈值类型
use_roc_threshold <- TRUE
fixed_fc_threshold <- NULL  # 如果use_roc_threshold=FALSE，设置如2.0

# 催化组过滤
min_valid_lfq <- 2

# 注释和作图
annotation_column <- "GO_Localization"
plot_all_annotations <- FALSE  # FALSE只显示TP，TRUE显示所有注释
tp_label <- "SGs"
tp_color <- "#DB6968"
```

**常用配置示例：**

```r
# 示例1：只分析催化实验组的相关性
correlation_config$exclude_catalytic <- c("NoCat")

# 示例2：只分析实验组，排除对照组
correlation_config$exclude_context <- c("Control")

# 示例3：直接指定特定bioGroup
correlation_config$include_biogroups <- c("K69A1B3_Light", "K69C3_Light", "C3_H2O2")

# 示例4：使用蓝-白-红配色方案（表达热图）
color_params <- list(
  custom_min = -3, custom_max = 3, custom_center = 0,
  col_low = "#0000FF", col_center = "#FFFFFF", col_high = "#FF0000"
)

# 示例5：不指定no_fdr_comparisons（Method A = Method B）
no_fdr_comparisons <- NULL

# 示例6：使用固定FC阈值而非ROC阈值
use_roc_threshold <- FALSE
fixed_fc_threshold <- 2.0

# 示例7：指定多个comparison不使用FDR过滤
no_fdr_comparisons <- c("K69A1B3_Light_vs_A1B3_Light", 
                        "K69C3_Light_vs_C3_Light",
                        "C3_H2O2_vs_noH2O2")
```

### Module 12: 第二次ROC分析 ✓

**文件：** `Module/module12_second_roc.R`

**功能：**
1. 基于 Module 11 输出的 `FDR_combined_df_list_2nd` 逐版本生成第二次 ROC 曲线
2. 默认使用 `GO_Localization` 中的 `SGs`（TP）与 `Cytosol`（FP），并输出 ROC 数据 + 图形
3. 计算 6 组比较的最优 log2FC 阈值（满足 TPR ≥ 0.3，否则阈值记为 0）
4. 构建 `Expr_FDR_df_list_2nd`：表达矩阵（Method A/B）左连第二次差异分析结果

**输入：**
- `FDR_combined_df_list_2nd`
- `expr_data_list`（由 `replaced_data_A/B` 合并并追加 `_A/_B` 后的列表）
- `sample_columns`（可选）用于限制表达矩阵只保留 Experiment / Spatial 样本列（main_pipeline 默认传入 `sample_info_2nd$SampleName`）
- `dir_config$output`
- 可选参数：`desired_order`、`annotation_column`、`comparison_cols`

**输出：**
- `Output/Module12_Step16_{version}_{FP}_{TP}_ROCdata.xlsx`
- `Output/Module12_Step16_1_{version}_GO_CytosolROC.pdf`
- `Output/Module12_Step16_1_{version}_Youden_Index_Plots_BaseR.pdf`
- `Output/Module12_Step16_all_desired_thresholds_2nd.xlsx`（每个版本的 sheet 中包含 Comparison、对应 logFC 列名、注释列、TP/FP 标签与阈值，便于追溯）
- `Output/Module12_Step16_Expr_FDR_df_list_2nd.xlsx`
- `Module12_workspace.RData`

**调用示例：**
```r
result <- module12_second_roc(
  dir_config = dir_config,
  FDR_combined_df_list_2nd = FDR_combined_df_list_2nd,
  expr_data_list = expr_data_list_2nd,
  desired_order = if (length(second_roc_order) > 0) second_roc_order else names(FDR_combined_df_list_2nd)
)
```

**提示：**
- 主流程默认用 `stringr::str_replace("_Imputed", "_New")` 统一命名，便于与 CleanCode 对齐
- `expr_data_list_2nd` 由 `replaced_data_A/B` 组合，名称以 `_A/_B` 标识 Method
- 所有输出均添加 `Module12_` 前缀，保持 `main_pipeline` 的命名习惯，同时保留 Step16 信息

---

### Module 13: SubSG注释 ✓

**文件：** `Module/module13_subsg_annotation.R`

**功能：**
1. 按 CleanCode.R (2708-2742) 逻辑，在 `Expr_FDR_df_list_2nd` 中追加 SubSG 注释
2. 默认对 `MultiBait_Localization` 列进行重写；若该列缺失则在最后一个 `_Localization` 列之后自动新建
3. 使用 HaloMap + HPA + MitoCarta 注释组合生成 9 种分类（Nuclear&Cytosol&SGs、Nuclear&SGs、Cytosol&SGs 等）
4. 结果导出为 `Module13_{version}_SubSGs.csv`，并返回 `Expr_FDR_df_list_2nd_SubSGs`
5. 可配置 `forstep16_versions`，用于生成与 CleanCode `ForStep16` 相同的子集（默认等于全部 SubSG 结果）

**输入：**
- `expr_fdr_df_list_2nd`（Module 12 输出）
- `annotation_references`（Module 3 输出，提供 HaloMap/HPA/MitoCarta 数据）
- 可选参数：`target_column`、`levels_order`、`forstep16_versions`、`output_prefix`

**输出：**
- `Expr_FDR_df_list_2nd_SubSGs`
- `ForStep16`（可用于与 Spatial 组比较的子集）
- `Output/Module13_{version}_SubSGs.csv`
- `Module13_workspace.RData`

**调用示例：**
```r
result <- module13_subsg_annotation(
  dir_config = dir_config,
  expr_fdr_df_list_2nd = Expr_FDR_df_list_2nd,
  annotation_references = annotation_references,
  target_column = "MultiBait_Localization",
  forstep16_versions = c("noMBR_New_A","noMBR_QNorm_New_A","noMBR_Local_QNorm_New_A")
)
Expr_FDR_df_list_2nd_SubSGs <- result$Expr_FDR_df_list_2nd_SubSGs
ForStep16 <- result$ForStep16
```

---

### Module 14: 火山图 ✓

**文件：** `Module/module14_volcano_plots.R`

**功能：**
1. 读取 `Expr_FDR_df_list_2nd_SubSGs`，自动打印可用的版本、注释列（*_Localization）及其独立元素，帮助确认可选项
2. 检索所有 `_logFC` 与 `_adj.P.Val` 成对列，列出可用 comparison 名称，并根据 Module11 的 `comparisons_list` + `bioGroup_info` 自动判断每个 comparison 属于 Exp_vs_Exp 还是 Exp_vs_Spatial
3. 可通过 `comparison_categories` 接口选择只绘制某一类（如只看 Exp_vs_Exp）
4. 根据配置的版本、注释列、颜色映射或阈值着色方案批量生成火山图 PDF；`label_mode` 支持 with/without 同时输出
5. 不会修改 `Expr_FDR_df_list_2nd_SubSGs` 原始数据

**输入：**
- `expr_fdr_df_list_2nd_subsgs`（Module 13 输出）
- `dir_config$output`
- 关键参数：
  - `versions`：进入绘图的版本（默认全部）
  - `annotation_column`：用来映射颜色的注释列（如 `GO_Localization` / `HaloMap_Localization`）
  - `annotation_color_map`：注释值到颜色的命名向量（如 `c("SGs"="red","Nuclear"="blue")`）
  - `use_threshold_colors`：TRUE 时按照 `logFC`/`FDR` 阈值着色，并覆盖注释颜色
  - `logfc_threshold`、`fdr_threshold`：显著性判定阈值（绝对值比较）
  - `threshold_color_above` / `threshold_color_below`：阈值模式下的颜色
  - `label_mode`：`"with"`、`"without"` 或向量（如 `c("with","without")` 表示同时导出）
  - `label_annotations`：允许打标签的注释类别
  - `comparison_sets`：自定义 logFC/FDR 列组合、坐标轴范围与输出名称
  - `comparison_categories`：`c("Exp_vs_Exp","Exp_vs_Spatial")` 或其子集，决定绘图时保留哪些比较；`NULL` 表示全部

**输出（每个 comparison_set × label_mode 各生成一个 PDF）：**
- `Output/Module14_Step20_{annotation}_{comparison}_{labelMode}.pdf`
- `Output/Module14_Step20_VolcanoSource.xlsx`（包含每个版本×comparison 的源数据，含 logFC/FDR/注释/阈值标记）
- `Module14_workspace.RData`
- 返回列表：`pdf_files`（字典形式，便于后续汇报）

**调用示例（main_pipeline 中默认配置）：**
```r
volcano_result <- module14_volcano_plots(
  dir_config = dir_config,
  expr_fdr_df_list_2nd_subsgs = Expr_FDR_df_list_2nd_SubSGs,
  versions = volcano_versions,
  annotation_column = volcano_annotation_column,
  annotation_color_map = volcano_annotation_color_map,
  use_threshold_colors = volcano_use_threshold_colors,
  logfc_threshold = volcano_logfc_threshold,
  fdr_threshold = volcano_fdr_threshold,
  threshold_color_above = volcano_threshold_color_above,
  threshold_color_below = volcano_threshold_color_below,
  label_mode = volcano_label_mode,
  label_annotations = volcano_label_annotations,
  comparison_sets = volcano_comparison_sets,
  comparison_categories = volcano_comparison_categories
)
```

---
