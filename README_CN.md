# 1. 解压

文件包解压即用，不要改变目录结构，将需要分析的数据整理为以下格式:

| Genes(列名固定)      | Sample1(列名可以不固定) | Sample2(列名可以不固定) | …… |
| -------------------- | ----------------------- | ----------------------- | ---- |
| Official_ProteinName |                         |                         |      |
| Official_ProteinName |                         |                         |      |

放进Rawdata文件夹中

# 2. 主流程
按照guide.md说明，在R(R 4.3.3)中运行main_pipeline.R

分析思路为：
1. 配置环境
2. 读取数据
3. 使用Human Protein Atlas (HPA)数据库和MitoCarta3.0数据库完成基本蛋白质定位注释。TP集合使用对应的无膜细胞器的数据集注释。
4. 数据标准化(不标准化，全局标准化，局部标准化三种情况)
    全局标准化：催化组与对照组等全部样本一起完成中位数标准化或Quantile标准化。必须至少完成全局quantile标准化。
    局部标准化：每个分组内(biogroup)，所有技术/生物学重复组内进行标准化，实现样品的中位数或分位数对齐，以减少后面第一次ROC分析扣除非催化组背景时出现的假阴性问题。使用此参数需要严格控制，Steptavidin-pull down实验时所有样品（催化组和对照组）的上样量和实验操作流程要尽量一致。
5. 缺失值填补（只对非催化组进行缺失值补充，按照Perseus方法根据单个样品内的数据分布模拟一个随机较低值补充，用于支持Limma包在第一次组间背景扣除时计算组间差异FC/FDR值）

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
  
6. 相关性热图+注释蛋白定位的聚类热图，观察评估数据质量
7. 第一次差异分析（催化vs非催化）。使用Limma包完成组间比较，给出组间比较的FC和FDR值并汇总到表达矩阵中。
8. 扣除背景的第一次ROC分析。使用pROC包完成ROC分析，根据指定的TP/FP注释，无监督寻找最适用于分隔背景蛋白的FC(催化vs非催化)阈值。输出每种比较的最适阈值和ROC分析图表。
9. 正式的第一次背景扣除，使用上述ROC分析中给出的最佳阈值，在每个biogroup中去除低于阈值的蛋白，且默认去除组内（生物学重复）有效值小于2的蛋白。此处有A/B两种方案，A可以指定对特定的组间比较不使用FDR过滤只使用最适FC，用于解决一些数据集中高背景的问题。B方案默认对所有组间比较实用FDR和最适FC过滤。
10. 从步骤9中获得每个biogroup的新蛋白列表，再合并为新的数据矩阵。只保留数据矩阵中的蛋白列表和NA结构，所有非0值替换为全局quantile标准化的对应值 (只有全局标准化的数据才能支持催化组之间的横向比较，如不同邻近标记反应数据集；选择替换而不是重新标准化是为了避免蛋白数量差异大时产生的系统偏差)
11. 第二次差异分析。使用步骤10产生的新数据矩阵(已经过滤了背景且替换为全局标准化后的数据)，用Limma包计算组间差异(催化组vsSpatial control)，给出LogFC和FDR值。
12. 第二次ROC分析(指定TP和FP完成ROC分析，看与spatial control相比组分差异和最适分隔阈值；此步骤只用于观察数据质量不进入后续计算)
13. 新创建或替换MultiBait_Localization注释列，在SG注释中加入SG_subset注释，以判断SG_ref中的潜在系统偏差；或加入Nucleolus_subset注释列。合并数据集输出ForStep16对象由于后续模块计算。
14. 火山图绘制（不用于后续计算）


# 3. 1D/2D模型构建与PPI分析
使用main_pipeline.R输出的模块14.RData作为251117_BaseModels_Concise.R的输入，运行251117_BaseModels_Concise.R脚本。脚本中为模版，需要根据使用情况更改各种接口。脚本依赖从Module14.RData中继承的ForStep19列表，不要改动此名称，进入Step29 多模型复杂需求分析部分时需要ForStep19对象.

# 4. 模型性能比较
使用main_pipeline.R输出的模块14.RData作为251117 ModelComparision.R的输入，同样需要改接口，但是保留最后的ForStep19对象。
