# Pipeline Module Interface Documentation

This document details the function interfaces for the 14 modules in the Pipeline.

**Focus**: Function input arguments (Arguments), expected data structures, and output objects. This will help you quickly understand how to call, configure, and debug each step.

## Phase I: Background Subtraction & Data Preprocessing

### Module 01: Setup

* **File**: `module01_setup.R`
* **Main Function**: `module01_setup(dir_config)`
* **Function**: Checks necessary directories, creates output directories, and saves runtime environment metadata.

| **Argument** | **Type** | **Default** | **Description**                                                                                        |
| :----------------- | :------------- | :---------------- | :----------------------------------------------------------------------------------------------------------- |
| `dir_config`     | List           | (Required)        | List containing paths like `root`, `rawdata`, `reference`, `output`. Defined by `main_pipeline.R`. |

* **Return**:
  * `config` (List): Contains runtime timestamp, R version, and system platform information.
* **Side Effects**: Automatically creates `Output/`, `Module/`, `Dev/` folders.

### Module 02: Data Import

* **File**: `module02_data_import.R`
* **Main Function**: `module02_data_import(dir_config, file_pattern, auto_process)`
* **Function**: Reads original TSV and generates the experiment grouping table (`sampleGroup`) using interactive or automatic modes.

| **Argument** | **Type** | **Default**      | **Description**                                                                                                                                          |
| :----------------- | :------------- | :--------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dir_config`     | List           | (Required)             | Directory configuration.                                                                                                                                       |
| `file_pattern`   | String         | `"_matrix.*\\.tsv$"` | Regex to match filenames in the Rawdata directory.                                                                                                             |
| `auto_process`   | Bool           | `FALSE`              | `FALSE`: Pauses after generating the template, waiting for user to fill CSV. `TRUE`: Assumes CSV is filled, reads and processes directly (for automation). |

* **Return**:
  * List:
    * `data`: Original expression matrix with cleaned column names (DataFrame).
    * `sampleGroup`: Organized and verified sample grouping information (DataFrame).
* **File Output**: `Module02_sampleGroup_template.csv`, `Module02_sampleGroup.csv`, `Output/Module02_data_raw.csv`.

### Module 03: Annotation

* **File**: `module03_annotation.R`
* **Main Function**: `module03_annotation(...)`
* **Function**: Loads databases like HPA and MitoCarta to annotate subcellular protein localization.

| **Argument**         | **Type** | **Default** | **Description**                                                                                                      |
| :------------------------- | :------------- | :---------------- | :------------------------------------------------------------------------------------------------------------------------- |
| `dir_config`             | List           | (Required)        | Directory configuration.                                                                                                   |
| `data_raw`               | DataFrame      | (Required)        | Raw data.                                                                                                                  |
| `sampleGroup`            | DataFrame      | `NULL`          | Sample information (optional).                                                                                             |
| `custom_annotations`     | List           | `NULL`          | **Advanced Interface**. Completely overrides default annotation configuration.                                       |
| `annotation_mode`        | String         | `"SG"`          | `"SG"` (Stress Granule mode) or `"Nucleolus"` (Nucleolus mode). Determines which reference sets are loaded by default. |
| `custom_tp_sources`      | List           | `NULL`          | Used to load user-defined external reference file paths and reading methods.                                               |
| `additional_annotations` | List           | `NULL`          | **Appends** new annotation column configurations on top of the default configuration.                                |

* **Return**:
  * List:
    * `data_annotated`: Data frame with annotation columns.
    * `annotation_references`: Large list containing all loaded reference databases (HPA, MitoCarta, etc.).
* **File Output**: `Output/Module03_data_annotated.csv`.

### Module 04: Standardization

* **File**: `module04_standardization.R`
* **Main Function**: `module04_standardization(...)`
* **Function**: Performs Log2 transformation and various normalization strategies.

| **Argument** | **Type** | **Default**                  | **Description**                                                                         |
| :----------------- | :------------- | :--------------------------------- | :-------------------------------------------------------------------------------------------- |
| `dir_config`     | List           | (Required)                         | -                                                                                             |
| `data_annotated` | DataFrame      | (Required)                         | Data from Module 03.                                                                          |
| `sampleGroup`    | DataFrame      | (Required)                         | Must contain `bioGroup` column for local normalization.                                     |
| `norm_types`     | Vector         | `c("noNorm", "Global_QNorm"...)` | Specifies normalization methods to execute.**Must include at least one global method**. |

* **Return**:
  * List:
    * `standardized_data_list`: List where keys are method names and values are DataFrames.
    * `norm_types_used`: List of actual method names executed.
* **File Output**: CSV files for each version and PDF QC charts (Boxplots).

### Module 05: Imputation

* **File**: `module05_imputation.R`
* **Main Function**: `module05_imputation(...)`
* **Function**: Fills missing values based on the Perseus strategy (left-skewed normal distribution).

| **Argument**         | **Type** | **Default** | **Description**                                                                               |
| :------------------------- | :------------- | :---------------- | :-------------------------------------------------------------------------------------------------- |
| `dir_config`             | List           | (Required)        | -                                                                                                   |
| `standardized_data_list` | List           | (Required)        | List from Module 04.                                                                                |
| `sampleGroup`            | DataFrame      | (Required)        | Must contain `CatalyticGroup` column.                                                             |
| `impute_cat_mean`        | Bool           | `FALSE`         | For `Cat` (Positive) group, if n_valid=2, whether to impute with mean (default is no imputation). |
| `random_seed`            | Integer        | `123`           | Random seed to ensure reproducibility.                                                              |

* **Key Logic**: `NoCat` group (Context=Control) is forcibly imputed to simulate background noise.
* **Return**:
  * List:
    * `imputed_data_list`: List, keys have `_Imputed` suffix added.
    * `imputation_params`: Records parameters used for imputation.
* **File Output**: Boxplot PDF comparing distributions before and after imputation.

### Module 06: Heatmap

* **File**: `module06_heatmap.R`
* **Main Function**: `module06_heatmap(...)`
* **Function**: Plots sample correlation heatmaps and clustering heatmaps.

| **Argument**       | **Type** | **Default**               | **Description**                                                                                |
| :----------------------- | :------------- | :------------------------------ | :--------------------------------------------------------------------------------------------------- |
| `imputed_data_list`    | List           | (Required)                      | -                                                                                                    |
| `sampleGroup`          | DataFrame      | (Required)                      | -                                                                                                    |
| `selected_versions`    | Vector         | `NULL`                        | Specifies which data version to plot (defaults to selecting the first one).                          |
| `heatmap_types`        | Vector         | `c("all", "by_localization")` | Plot types:`correlation`, `all`, `by_localization`.                                            |
| `correlation_config`   | List           | `NULL`                        | **Important**. Configure correlation heatmap details, e.g., `exclude_context`, `corr_min`. |
| `localization_columns` | Vector         | `NULL`                        | Specifies annotation column names for categorical plotting (NULL detects automatically).             |
| `color_params`         | List           | `NULL`                        | Custom heatmap color range (`custom_min`, `custom_max`).                                         |

* **Return**:
  * `heatmap_info` (List): List containing plotting parameter information.
* **File Output**: Correlation heatmap PDF, Global/By-Localization clustering heatmap PDFs.

### Module 07: Diff Analysis I (Background)

* **File**: `module07_diff_analysis1.R`
* **Main Function**: `module07_diff_analysis1(...)`
* **Function**: Background subtraction differential analysis (Exp vs Control) based on `FirstROCgroup`.

| **Argument**    | **Type** | **Default** | **Description**                                   |
| :-------------------- | :------------- | :---------------- | :------------------------------------------------------ |
| `imputed_data_list` | List           | (Required)        | -                                                       |
| `sampleGroup`       | DataFrame      | (Required)        | Must contain `FirstROCgroup` and `Context`.         |
| `selected_versions` | Vector         | `NULL`          | Specifies data versions to analyze (NULL analyzes all). |

* **Core Logic**: Reads `sampleGroup$FirstROCgroup`. Within the same group, compares `Context=Experiment` samples vs `Context=Control` samples.
* **Return**:
  * List:
    * `diff_results1`: Core result list. Contains `combined` (large wide table), `expr_matrix` (expression matrix), `comparisons` (comparison group info).
    * `comparisons_used`: List of comparison groups actually executed.
    * `comparison_info`: Detailed comparison group metadata.
* **File Output**: `_DiffAnalysis.csv` (merged results), `_DiffAnalysis_Full.xlsx` (containing raw topTable results).

### Module 08: ROC Analysis I (Background)

* **File**: `module08_roc_analysis1.R`
* **Main Function**: `module08_roc_analysis1(...)`
* **Function**: Calculates the optimal LogFC threshold for background removal.

| **Argument**        | **Type** | **Default**     | **Description**                                                                |
| :------------------------ | :------------- | :-------------------- | :----------------------------------------------------------------------------------- |
| `diff_results1`         | List           | (Required)            | From Module 07.                                                                      |
| `annotation_references` | List           | (Required)            | From Module 03, must contain MitoCarta.                                              |
| `selected_versions`     | Vector         | `NULL`              | -                                                                                    |
| `roc_annotation_column` | String         | `"GO_Localization"` | "Truth" reference column for ROC.                                                    |
| `tp_label`              | String         | `"SGs"`             | True Positive Label (Target).                                                        |
| `fp_label`              | String         | `"Matrix"`          | False Positive Label (Background), usually converted from SubMito.                   |
| `enable_submito`        | Bool           | `TRUE`              | Whether to enable SubMito conversion (splitting Mitochondrion into Matrix/IMS etc.). |
| `min_tp`                | Numeric        | `0.3`               | Minimum TP proportion (Sensitivity) required when calculating thresholds.            |

* **Return**:
  * List:
    * `roc_results`: Contains ROC objects and data frames.
    * `all_thresholds`: List of calculated optimal LogFC thresholds.
    * `data_with_submito`: Data after SubMito conversion.
    * `expr_fdr_df_list`: Intermediate data integrating expression, FDR, and original annotations (for Module 09).
* **File Output**: ROC Curve PDF, Youden Index PDF, Threshold Excel.

### Module 09: Background Subtraction

* **File**: `module09_background_subtraction.R`
* **Main Function**: `module09_background_subtraction(...)`
* **Function**: Performs final filtering to generate the "True" interactor list.

| **Argument**                   | **Type** | **Default**        | **Description**                                                                                           |
| :----------------------------------- | :------------- | :----------------------- | :-------------------------------------------------------------------------------------------------------------- |
| `diff_results`, `roc_thresholds` | List           | (Required)               | -                                                                                                               |
| `fdr_threshold`                    | Numeric        | `0.05`                 | General FDR threshold.                                                                                          |
| `no_fdr_comparisons`               | Vector         | `NULL`                 | **Method A Core Parameter**. Specifies which comparisons do **not** use FDR filtering (LogFC only). |
| `use_roc_threshold`                | Bool           | `TRUE`                 | Whether to use dynamic thresholds calculated in Module 08.                                                      |
| `fixed_fc_threshold`               | Numeric        | `NULL`                 | If `use_roc_threshold=FALSE`, use this fixed value.                                                           |
| `min_valid_lfq`                    | Integer        | `2`                    | Minimum valid value count required for the catalytic group.                                                     |
| `tp_label`, `tp_color`           | String         | `"SGs"`, `"#DB6968"` | TP label and color for plotting.                                                                                |

* **Return**:
  * List:
    * `filtered_data_A/B`: Filtered list (separated by BioGroup).
    * `merged_data_A/B`: Filtered list (merged and deduplicated, for Module 10).
* **File Output**: `FilteredData.xlsx`, Stacked Bar Chart PDF (showing localization proportions).

## Phase II: Biological Discovery

### Module 10: Data Replacement

* **File**: `module10_data_replacement.R`
* **Main Function**: `module10_data_replacement(...)`
* **Function**: Backfills quantitative values for selected proteins to prepare for inter-group comparison.

| **Argument**    | **Type** | **Default**          | **Description**                   |
| :-------------------- | :------------- | :------------------------- | :-------------------------------------- |
| `merged_data_A/B`   | List           | (Required)                 | "Mask" list from Module 09.             |
| `imputed_data_list` | List           | (Required)                 | "Value Source" from Module 05.          |
| `selected_versions` | Vector         | `NULL`                   | -                                       |
| `prefer_version`    | String         | `"Global_QNorm_Imputed"` | Preferred value version.                |
| `fallback_version`  | String         | `"Global_MNorm_Imputed"` | Fallback value version.                 |
| `test_n`            | Integer        | `20`                     | Number of genes to sample for checking. |

* **Return**:
  * List:
    * `replaced_data_A`, `replaced_data_B`: Matrix lists with values backfilled.
    * `base_version_used`: The base version actually used.
* **File Output**: Validation reports CSV/TXT in `Output/Check/`, `_AfterFirstROC_Intersect.xlsx`.

### Module 11: Diff Analysis II (Biological)

* **File**: `module11_diff_analysis2.R`
* **Main Function**: `module11_diff_analysis2(...)`
* **Function**: Compares differences between Experiment groups or Experiment vs Spatial.

| **Argument**    | **Type** | **Default** | **Description**                            |
| :-------------------- | :------------- | :---------------- | :----------------------------------------------- |
| `replaced_data_A/B` | List           | (Required)        | From Module 10.                                  |
| `sampleGroup`       | DataFrame      | (Required)        | Used to parse `SecondROCgroup` and `PLtype`. |
| `selected_versions` | Vector         | `NULL`          | -                                                |

* **Core Logic**:
  * Retains only `Context %in% c("Experiment", "Spatial")`.
  * Automatically constructs `Exp vs Exp` (inter-group) and `Exp vs Spatial` (group vs reference) comparison pairs.
* **Return**:
  * List:
    * `FDR_combined_df_list_2nd`: Summary table containing second-round LogFC and FDR.
    * `comparisons_list`: List of generated comparisons.
    * `bioGroup_info`: Info of sample groups involved in analysis.
* **File Output**: `Module11_FC_FDR_2nd_combined.xlsx` (core results), `_Raw_FDR_test_list_*.xlsx` (raw data).

### Module 12: ROC Analysis II (Biological)

* **File**: `module12_second_roc.R`
* **Main Function**: `module12_second_roc(...)`
* **Function**: Calculates enrichment thresholds for specific biological comparisons.

| **Argument**           | **Type** | **Default**     | **Description**                               |
| :--------------------------- | :------------- | :-------------------- | :-------------------------------------------------- |
| `FDR_combined_df_list_2nd` | List           | (Required)            | From Module 11.                                     |
| `expr_data_list`           | List           | (Required)            | Expression matrix list (usually `replaced_data`). |
| `annotation_column`        | String         | `"GO_Localization"` | ROC "Truth" column.                                 |
| `tp_label`                 | String         | `"SGs"`             | True Positive (e.g., Stress Granules).              |
| `fp_label`                 | String         | `"Cytosol"`         | False Positive (e.g., Cytosol).                     |
| `min_tp`                   | Numeric        | `0.3`               | Minimum TP proportion.                              |
| `youden_ylim`              | Vector         | `c(0, 0.7)`         | Y-axis range for Youden Index plot.                 |

* **Return**:
  * List:
    * `all_desired_thresholds_2nd`: Second-round optimal LogFC thresholds.
    * `Expr_FDR_df_list_2nd`: Complete data object integrating expression, FDR, and annotation.
* **File Output**: ROC Curve PDF, Youden Index PDF, Threshold Excel, Data Excel.

### Module 13: Sub-Annotation

* **File**: `module13_subsg_annotation.R`
* **Main Function**: `module13_subsg_annotation(...)`
* **Function**: Generates complex subcellular localization labels (e.g., "Nuclear&SGs") using logical operations.

| **Argument**        | **Type** | **Default** | **Description**                                                                             |
| :------------------------ | :------------- | :---------------- | :------------------------------------------------------------------------------------------------ |
| `expr_fdr_df_list_2nd`  | List           | (Required)        | From Module 12 (or 11).                                                                           |
| `annotation_references` | List           | (Required)        | Must contain HPA/HaloMap references.                                                              |
| `mode`                  | String         | `"subSG"`       | `"subSG"` (SGs analysis) or `"subNucleolus"` (Nucleolus analysis).                            |
| `target_column`         | String         | `NULL`          | Output column name (e.g.,`"MultiBait_Localization"`). NULL decides automatically based on mode. |
| `levels_order`          | Vector         | `NULL`          | Specifies Factor Level order (affects plot color order).                                          |

* **Return**:
  * List:
    * `Expr_FDR_df_list_2nd_SubSGs`: Final complete dataset for plotting (contains complex annotations).
    * `ForStep16`: Data subset for subsequent steps.
* **File Output**: `Module13_*.csv`.

### Module 14: Volcano Plots

* **File**: `module14_volcano_plots.R`
* **Main Function**: `module14_volcano_plots(...)`
* **Function**: Generates final visualization results.

| **Argument**              | **Type** | **Default**             | **Description**                                                                                                                    |
| :------------------------------ | :------------- | :---------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------- |
| `expr_fdr_df_list_2nd_subsgs` | List           | (Required)                    | From Module 13.                                                                                                                          |
| `comparison_sets`             | List           | `NULL`                      | **Core Config**. Defines which comparisons are plotted together, and corresponding xlim/ylim. If NULL, uses default configuration. |
| `annotation_column`           | String         | `"GO_Localization"`         | Column used for coloring (usually generated in Module 13).                                                                               |
| `annotation_color_map`        | Vector         | (Default Red/Blue/Black)      | Color map.                                                                                                                               |
| `use_threshold_colors`        | Bool           | `FALSE`                     | Whether to ignore annotation and color only by threshold crossing (Red/Gray).                                                            |
| `logfc_threshold`             | Numeric        | `0.5`                       | Line drawing threshold.                                                                                                                  |
| `fdr_threshold`               | Numeric        | `0.05`                      | Line drawing threshold.                                                                                                                  |
| `label_mode`                  | Vector         | `c("with", "without")`      | Generates both labeled and unlabeled PDFs.                                                                                               |
| `label_annotations`           | Vector         | `c("SGs", "Mitochondrion")` | Specifies which annotation categories require gene name labels.                                                                          |

* **Return**:
  * List: Contains list of generated PDF file paths (`pdf_files`) and source data file path (`source_data_file`).
* **File Output**: Multiple PDF Volcano plots, `_VolcanoSource.xlsx` (Plotting source data).

---

# Subcellular Component Multidimensional Model Selection & Validation System (Step29) Project Guide

## 1. Project Overview & Design Framework

### 1.1 Project Background

This project module (Step29) follows Module 14 of `main_pipeline.R`, taking `ForStep16` as input. Its core objective is to break through the limitations of traditional single-dimension (LogFC only or Abundance only) screening. By integrating **Mass Spectrometry Quantitative Data (Abundance/Difference)** with **Protein-Protein Interaction (PPI)** information, it constructs Multidimensional Generalized Linear Models (GLM) to identify candidate proteins with specific subcellular localization attributes (e.g., Nucleolus, Stress Granules).

### 1.2 Core Design Framework

The script adopts a **Three-Phase** design pattern, implementing a closed loop from calculation to visualization and evaluation:

* **Phase 1: Feature Engineering & Model Construction (Data Calculation & Modeling)**
  * Data Preprocessing: Integrate local STRING database, calculate Global and Local PPI scores.
  * Model Training: Train Logistic Regression models based on known localization (Ground Truth).
  * Scoring & Prediction: Score and classify the whole proteome using optimized Youden Index thresholds.
  * *Output: RData objects and Cytoscape network files.*
* **Phase 2: Model Evaluation & Visualization (Visualization)**
  * Performance Evaluation: ROC Curves, Youden Index plots.
  * Feature Space Visualization: Decision Boundary plots, 3D Scatter plots.
  * Network Density Analysis: PPI Heatmaps, Joint Probability (JP) Matrix analysis (Core innovation, used to evaluate protein clustering tightness in the network).
* **Phase 3: Candidate Comparison & Final Selection (Comparison & Selection)**
  * Method Benchmarking: Compare machine learning model results with traditional screening (logFC/FDR) and spatial controls.
  * Composition Analysis: Analyze subcellular localization composition of candidate lists (MultiBait Analysis).
  * Enrichment Analysis: TP (True Positive) Percentage Curve, evaluating the accuracy of top-ranked candidates.

---

## 2. Statistical Concepts & Mathematical Principles

### 2.1 Integration of Bayesian Prior Concepts

The system relies not only on the posterior distribution of experimental data (logFC) but also introduces "Prior Knowledge" (PPI Network). Assumption: **Proteins localized in the same membraneless organelle tend to interact with each other in the PPI network.**

### 2.2 Generalized Linear Model (GLM)

Use Logistic Regression to predict the probability $P(Y=1|X)$ of a protein belonging to the target localization (TP):

$$
\ln(\frac{P}{1-P}) = \beta_0 + \beta_1 X_{abundance} + \beta_2 X_{logFC} + \beta_3 X_{PPI}
$$

Where $X$ are input features, and $\beta$ are weights learned by the model.

### 2.3 Probability Mass & Spatial Density

In Phase 2, interaction density analysis based on Kernel Density Estimation (KDE) is introduced.

* **Concept**: If a group of candidate proteins are true complex components, they should appear "clustered" in the PPI matrix rather than randomly distributed.
* **Implementation**: Calculate the interaction probability density of protein pairs in the ranking space using Joint Probability (JP) and Probability Mass (PM) matrices.

---

## 3. Interface Configuration (Configuration)

Check the following key interfaces before running the script:

### 3.1 Input & Grouping (`group_info`)

Define experiment group structure, including replicate column names, logFC column, and FDR column.

**R**

```r
group_info <- list(
  E7A2B4 = list(samples = c("..."), logFC = "...", logFC_FDR="..."),
  ...
)
```


### 3.2 True Positive Definition (`myTP_vector` & `myTP_vector2`)

Define which subcellular localization labels are considered "Gold Standards" (Ground Truth) for model training.

* **Current Config** : Focuses on Nucleolus related components.
* **Note** : `myTP_vector2` is specifically for Phase 3 TPR curve calculation, with broader logic (includes all strings containing "Nucleolus").

### 3.3 Screening Parameters (`PARAM_*`)

* `PARAM_MIN_LFQ_NONNA`: Minimum number of non-NA samples required for quantitation before modeling (Default 2).
* `PARAM_MODELS_TO_PLOT`: List of models to display and compare in Phase 3.
* `PARAM_LOCALIZATION_COLORS`: Plotting color map, must cover all localization categories appearing in the data.

---

## 4. Detailed Module Design

### Phase 1: Data Calculation & Export

1. **PPI Preprocessing** :

* Load local STRING database (`9606.protein.links`).
* **Local PPI Score** : Calculate interaction strength of each protein with its neighbors in abundance ranking (Window Size=10). This captures interacting proteins co-migrating in gradients.
* **Global PPI Score** : Calculate total interaction strength of the protein in the entire network.

1. **Model Training Loop** :

* Build various models (1D, 2D, 3D) for each experimental group.
* Use `pROC` package to calculate AUC and determine optimal classification threshold using Youden Index (`sensitivity + specificity - 1`).
* Output: RData and Cytoscape CSV files in `step29_Analysis_Phase1_DataExport` directory.

### Phase 2: Visualization System

This phase generates a series of PDFs for model evaluation:

1. **ROC/Youden (Folders 1 & 2)** : Evaluate model classification efficiency.
2. **Decision Boundaries (Folder 3)** :

* Draw model separation lines in 2D coordinate systems (e.g., Abundance vs logFC).
* **Significance** : Visually demonstrate how the model "cuts" the data (e.g., whether it requires both high abundance and high fold change).

1. **3D Scatter Plots (Folder 4)** :

* Use color to represent the third dimension (PPI Score), showing data distribution.

1. **PPI Heatmap & Density (Folders 5, 7, 13, 15)** :

* **Candidate PPI Heatmaps** : Display interactions only among candidates passing the screen.
* **JP/PM Heatmaps (Advanced)** : Map proteins to an **$N \times N$** grid by abundance ranking. Redder/brighter colors indicate higher interaction density in that region (usually bottom-left, i.e., high abundance region). This is strong evidence for verifying phase separation properties of membraneless organelles.

### Phase 3: Comparison & Decision

This phase answers "Which model is best?" and "What was selected?":

1. **Detailed List Export** : Generate Excel tables containing candidate proteins from all models.
2. **Composition Stacked Plot (Plot_Count/Abundance)** :

* Show the proportion of TP (target localization) vs other localizations in candidates selected by different methods.
* **Significance** : If the TP proportion of a model is significantly higher than Control or traditional methods, it indicates good model specificity.

1. **MultiBait Analysis** :

* Refined analysis for multi-bait experimental designs, showing finer localization classifications (e.g., Nuclear&Nucleolus vs Cytosol&Nucleolus).

1. **TP Percentage Curve (Figure 9)** :

* X-axis is ranking percentage, Y-axis is TP proportion.
* **Significance** : Curves closer to the top-left indicate the model ranks true target proteins at the very top of the list. This is a core metric for evaluating ranking quality.

---

## 5. Chart Significance Cheat Sheet

| **Chart Type**        | **File Path Keyword** | **Core Significance**                  | **Ideal Form**                                                                              |
| --------------------------- | --------------------------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| **ROC Curve**         | `1_ROC_Curves`            | Ability to distinguish TP/Other              | AUC close to 1, curve arches to top-left                                                          |
| **Decision Boundary** | `3_Decision_Boundaries`   | Model screening cutoff line                  | Green line (2D model) encloses target group more flexibly than red/blue dashed lines (1D)         |
| **PPI Heatmap**       | `7_Candidate_PPI`         | Interaction network within candidates        | Distinct red blocks appear, indicating strong interaction among candidates                        |
| **JP/PM Heatmap**     | `13_JP`,`15_JPw`        | Spatial distribution density of interactions | **Bottom-left**(High rank region) shows bright blocks, significantly higher than background |
| **Composition Plot**  | `Plot_AbundancePercent`   | Purity of candidate set                      | Target localization (e.g., Red TP) block has highest proportion                                   |
| **TP Curve**          | `TP_Percent_Curve`        | Ranking accuracy                             | Curve maintains extremely high Y value (>80%) in early X axis (Top 10-20%)                        |

---

## 6. Next Step Recommendations

As a researcher, after running this script, it is recommended to:

1. **Check `TP_Percent_Curve` in Phase 3** : Select the model with the largest area under the curve or highest starting point as the final candidate screening standard.
2. **Review Excel tables in `Step29_Analysis_Phase3_Comparison`** : Extract the `Detailed_Candidates` Sheet corresponding to the best model.
3. **Cytoscape Visualization** : Use `_edges_ModelScores.csv` and `_nodes.csv` exported in Phase 1 to draw weighted interaction networks in Cytoscape, where weights can directly use model prediction probabilities (Prob).

---

# Subcellular Component Model Comparison & Robustness Evaluation System (Step29_ModelComparison) Project Guide

## 1. Project Overview & Design Framework

### 1.1 Project Background

This project is the advanced evaluation module of the subcellular localization analysis workflow. After building multiple prediction models (1D/2D/3D), specific key scientific questions must be answered:

1. **Model Superiority** : Does model performance significantly improve over traditional logFC screening after introducing PPI or abundance data?
2. **Methodological Comparison** : How do different proximity labeling methods (e.g., BioID, APEX, etc., i.e., different Groups) differ in capturing specific organelle components?
3. **Robustness** : Is high model performance due to chance in specific samples, or is it broadly stable?

### 1.2 Core Design Framework

The script adopts a high-intensity statistical framework of  **Parallel Computing + Bootstrap + Non-parametric Tests** :

* **Configuration & Initialization** : Highly parameterized interface allowing flexible selection of groups, models, and visualization switches for comparison.
* **Phase 1: Core Calculation**
  * **PPI Integration** : Load local STRING database.
  * **Parallel Model Training** : Train full set of GLM models for each group.
  * **Statistical Testing** : Calculate AUC, perform DeLong Test (compare vs traditional), perform Bootstrap Test (calculate specificity/sensitivity gain).
* **Phase 1.1: Consolidation & Robustness**
  * **Type 1 Robustness** : Sample intersection fairness analysis (Union vs Intersection).
  * **Type 2 Robustness** : Model stability analysis (Standard Deviation of AUC).
  * **Cross-Group Comparison** : Compare ROC curves of different experimental groups (proximity labeling methods) using DeLong Test.
* **Phase 1.2: Visualization System**
  * Generates 7 types of key charts (Forest plots, Heatmaps, Dumbbell plots, Scatter plots, Gain plots, etc.) to display evaluation results comprehensively.

---

## 2. Statistical Concepts & Mathematical Principles

### 2.1 ROC Curve Comparison (DeLong's Test)

The script looks not only at AUC values but also uses **DeLong's Test** to calculate the statistical significance (P-value) of the difference between two ROC curves.

* **Application Scenarios** :
* Model vs Traditional Method (Model vs `logFC>0.5 & FDR<0.05`).
* Model vs 1D_logFC Model.
* Exp Group vs Exp Group (e.g., Group A vs Group B).

### 2.2 Controlled Variable Comparison via Bootstrap

This is an advanced statistical analysis used to prove "Performance improvement from introducing new variables (like PPI) is not accidental."

* **Principle** :

1. **Fixed Sensitivity** : Find the point on the ROC curve with the same sensitivity as the control model (e.g., logFC) and compare the Specificity difference at this point.
2. **Fixed Specificity** : Similarly, compare Sensitivity difference.
3. **Bootstrap** : Resample 500-1000 times to build difference distribution, calculate 95% Confidence Interval (CI) and P-value.

* **Significance** : Proves the model "improves precision without sacrificing recall" (or vice versa).

### 2.3 Multiple Testing Correction (BH Correction)

Since the script involves many pairwise comparisons, all output P-values are corrected using **Benjamini-Hochberg (BH)** to generate Q-values (q-value) to control the False Discovery Rate.

---

## 3. Interface Configuration (User Interface)

The `Configuration Interface - User Optional Settings` area at the top of the script is key to controlling analysis flow.

### 3.1 Core Selections

* `selected_groups`: Specify experimental groups for analysis (e.g., `c("E7A2B4", "C2")`). If `NULL`, analyzes all groups.
* `selected_models`: Specify list of models for comparison. Suggest including key models like `1D_Abundance`, `1D_logFC`, `2D_Abundance_logFC`.

### 3.2 Visualization Switches (`enable_*`)

A set of Boolean (`TRUE`/`FALSE`) switches to finely control output content:

* `enable_cytoscape_export`: Whether to output network files.
* `enable_robustness_analysis`: Whether to run time-consuming intersection analysis.
* `enable_plot_auc_summary`: Whether to draw Forest plot (Figure 1).
* `enable_plot_controlled_comparison`: Whether to draw Bootstrap Gain plot (Figure 5).
* ... (Covers Figure 1 to Figure 7)

### 3.3 Statistical Parameters

* `n_bootstrap`: Resampling count for controlled variable comparison. **Recommendation: Set to 200 for testing, 1000 for formal plots, 2000+ for rigorous publication.**
* `parallel_cores`: Number of parallel cores (default 20), adjust according to server performance.

---

## 4. Detailed Module Design & Chart Significance

### 4.1 Phase 1: Data Calculation (`process_group` function)

* **Function** : This is the smallest unit of parallel computing. It loads data, calculates local PPI, trains GLM models, and immediately runs DeLong and Bootstrap tests in an isolated environment.
* **Output** : Returns not only AUC but also complete `roc_objects` and `controlled_comparison` data frames.

### 4.2 Phase 1.1: Consolidation & Robustness

* **Group Comparison** :
* Select a baseline model (default `2D_Abundance_logFC`).
* Perform pairwise DeLong tests between different experimental groups.
* **Output** : `Group_Comparison_Summary.xlsx`.
* **Intersection Analysis** :
* **Concept** : If we only look at "common proteins" identified by all groups, how do performance differences look? This excludes bias caused by different identification depths of different methods.
* **Output** : `Group_Comparison_Intersection_Summary.xlsx`.

### 4.3 Phase 1.2: Visualization System (Plots)

| **Chart No.** | **Chart Name**                 | **Filename Keyword**            | **Statistical Significance & Interpretation**                                                                                                     |
| ------------------- | ------------------------------------ | ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Fig 1**     | **Group Forest Plot**          | `1_AUC_Faceted_Forest`              | **Model Overview** . Shows AUC & 95% CI for each model. Red dots indicate significantly better than traditional methods.                          |
| **Fig 1b**    | **vs logFC Forest Plot**       | `1b_AUC_..._vs_1D_logFC`            | **Incremental Value** . Specifically shows if models are significantly better than the logFC-only model.                                          |
| **Fig 2a/b**  | **Group Comparison**           | `2a_Bar`,`2b_Heatmap`             | **Methodology Eval** . Which proximity labeling method (Group) has the highest AUC? Heatmap shows pairwise significance.                          |
| **Fig 3**     | **Dumbbell Plot**              | `3_Intersection_Dumbbell`           | **Data Fairness** . Compares "Full Data AUC" vs "Intersection Data AUC". Huge difference implies method relies on unique proteins.                |
| **Fig 4**     | **Performance-Stability**      | `4_Performance_vs_Stability`        | **Model Selection Guide** . X-axis is AUC SD (Stability), Y-axis is AUC Mean. **Prefer top-left models** .                                  |
| **Fig 5a/b**  | **Gain Plot (vs Traditional)** | `5a_Specificity`,`5b_Sensitivity` | **Clinical/App Value** . Shows Specificity gain at fixed sensitivity of traditional method. Error bars not crossing 0 line indicate significance. |
| **Fig 5c**    | **Combined Gain Scatter**      | `5c_Combined_Gain`                  | **Comprehensive Eval** . Dots in top-right quadrant indicate superiority over traditional methods in both metrics.                                |
| **Fig 6**     | **Localization Distribution**  | `6_Localization_Distribution`       | **Biological Validation** . Shows distribution proportion of screened candidates in known localizations (MultiBait).                              |
| **Fig 7**     | **Gain Plot (vs logFC)**       | `7a/b/c_..._vs_1D_logFC`            | Same as Fig 5, but reference is 1D_logFC model, proving value of multidimensional models more strictly.                                                 |

---

## 5. Key Precautions

1. **`myTP_vector` Scope** :

* In the script, `myTP_vector` is redefined once inside the `process_group` function. This prevents issues with global variables not being passed in the parallel computing environment (`future` package). If you modify the TP definition,  **you must modify it both inside and outside the function** .

1. **Bootstrap Time Consumption** :

* When `n_bootstrap` is set to 1000, calculation volume is huge (Models × Groups × 1000 resamples). For debugging, temporarily change it to 50 or 100.

1. **Cytoscape Export** :

* The script exports `_node_attributes.csv` and `_edge_attributes.csv` for each group.
* **Node Attributes** contain prediction probabilities (`prob_...`) for all models, which can be visualized in Cytoscape by mapping colors to show prediction differences between models.

## 6. Next Step Recommendations

1. **Run Script** : First run with default parameters, check generated `Plot 1` (Forest Plot) and `Plot 4` (Performance-Stability Plot).
2. **Select Best Model** : Based on Plot 4, select the model in the **Top-Left** (High Mean AUC, Low SD). `2D_Abundance_logFC` is usually a good balance point.
3. **Review Gain Plot (Plot 7c)** : Confirm if your chosen best model is in the top-right quadrant (better than 1D_logFC in both Sensitivity and Specificity), and use this plot as core evidence in your paper.
4. **Export Candidate List** : Go to `Step29_ModelComparision` folder and view `detailed_candidates_by_model.xlsx` to get the final screened protein list.
