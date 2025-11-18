# 1. Unpacking

Extract the archive and keep the directory structure unchanged. Organize the data you want to analyze using the following format:

| Genes (column name fixed) | Sample1 (column name may vary) | Sample2 (column name may vary) | … |
| ------------------------ | ------------------------------ | ------------------------------ | --- |
| Official_ProteinName     |                                |                                |     |
| Official_ProteinName     |                                |                                |     |

Place the formatted file(s) inside the `Rawdata/` folder.

# 2. Main workflow
Follow the instructions in `guide.md` and run `main_pipeline.R` in R (tested with R 4.3.3).

The analysis flow is:
1. Configure the environment
2. Load the data
3. Annotate protein localization using the Human Protein Atlas (HPA) and MitoCarta 3.0 databases, assigning TP sets with the corresponding acellular organelle datasets.
4. Normalize the data (three modes: no normalization, global normalization, local normalization)
   - **Global normalization**: all samples from catalytic and control groups undergo median or quantile normalization together. At least global quantile normalization must be completed.
   - **Local normalization**: within each `bioGroup`, normalize every technical/biological replicate set to align medians or quantiles, which helps avoid false negatives when subtracting the background from the first ROC analysis. Use this option cautiously—during Streptavidin pull-down experiments the loading amount and workflow for all samples (both catalytic and control) must be as consistent as possible.
5. Impute missing values (only for the non-catalytic group). Missing values are filled via Perseus-style simulation: for each sample a random low value is drawn from the sample-specific distribution to support Limma in calculating FC/FDR during the first background subtraction.

**1. Perseus parameter calculation**
For each column (sample) compute:
- `imputation_mean = col_mean - 1.8 * col_sd`
- `imputation_sd = 0.3 * col_sd`

**2. NoCat group imputation rules** (mandatory)
Based on `sampleGroup$CatalyticGroup == "NoCat"`:
- **n_valid == 2**: impute with the mean of the two valid values within the same `bioGroup`
- **n_valid < 2**: draw a sample from `N(imputation_mean, imputation_sd)`
  - sampling is deterministic using the provided `random_seed`

**3. Cat group imputation rules** (optional)
Based on `sampleGroup$CatalyticGroup == "Cat"`:
- **Default (`impute_cat_mean = FALSE`)**: leave `NA`
- **If `impute_cat_mean = TRUE`**:
  - `n_valid == 2`: impute with the group mean
  - `n_valid < 2`: keep as `NA`

6. Generate correlation and clustering heatmaps (annotated with localization data) to assess data quality
7. First differential analysis (catalytic vs. non-catalytic) with the Limma package, outputting FC and FDR added into the expression matrix
8. First ROC analysis for background subtraction using the `pROC` package. Given the specified TP/FP annotation, find the FC threshold that best separates background proteins (catalytic vs. non-catalytic) in an unsupervised manner and export the optimal thresholds plus ROC plots per comparison.
9. Apply the first explicit background subtraction: remove proteins with values below the optimal thresholds (from step 8) within each `bioGroup`, and by default drop proteins with fewer than 2 valid values per biological replicate. There are two strategies (A/B):
   - Strategy A allows overriding FDR usage for specific group comparisons (only the optimal FC is used)
   - Strategy B applies both FDR and optimal FC filtering to all comparisons
10. Build a new protein list from step 9 (per `bioGroup`) and merge it into a new data matrix. Keep only the protein list and the `NA` structure; replace every non-zero value with the globally quantile-normalized counterpart (only globally normalized data supports cross-catalytic comparison; replacing instead of renormalizing avoids biases caused by large protein count differences).
11. Second differential analysis using the cleaned, globally normalized matrix (from step 10). Run Limma to compare catalytic groups vs. Spatial control and report LogFC/FDR.
12. Second ROC analysis (observe how the catalytic groups compare to spatial control and log the optimal separating threshold; this step is for diagnostics only and does not feed into later calculations)
13. Create or refresh the `MultiBait_Localization` annotation column, add `SG_subset` annotations to assess potential biases in `SG_ref`, or add `Nucleolus_subset`. Merge datasets and export the `ForStep16` object needed by later modules.
14. Draw volcano plots (informational only).

# 3. 1D/2D modeling and PPI analysis
Use the `Module14.RData` output from `main_pipeline.R` as input to `251117_BaseModels_Concise.R`. The script is a template—adapt the interfaces based on your data. It depends on the `ForStep19` list saved inside `Module14.RData`; do not rename it, because Step 29 (multi-model analysis) expects `ForStep19`.

# 4. Model performance comparison
Also use the `Module14.RData` output as input to `251117 ModelComparision.R`. Update any script interfaces as needed, but keep the final `ForStep19` object unchanged.
