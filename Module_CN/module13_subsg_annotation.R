## ============================================================================
## Module 13: Sub注释模块（SubSG / SubNucleolus）
## 对齐 CleanCode.R (2708-2742) 与 20250725.R (3530-3593) 逻辑
## ============================================================================

module13_subsg_annotation <- function(
    dir_config,
    expr_fdr_df_list_2nd,
    annotation_references,
    mode = c("subSG", "subNucleolus"),
    target_column = NULL,
    output_prefix = "Module13",
    levels_order = NULL,
    forstep16_versions = NULL
) {

  mode <- match.arg(mode)
  mode_label <- ifelse(mode == "subNucleolus", "SubNucleolus", "SubSG")
  cat(sprintf("\n=== Module 13: %s 注释 ===\n\n", mode_label))

  required_pkgs <- c("dplyr", "stringr", "purrr", "readr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("✗ 错误：需要安装 %s 包", pkg))
    }
  }

  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)

  if (!exists("load_HPA_annotations")) {
    module03_path <- file.path(dir_config$module, "module03_annotation.R")
    if (!file.exists(module03_path)) {
      stop("✗ 错误：无法找到 module03_annotation.R 以加载HPA逻辑")
    }
    source(module03_path)
  }

  if (length(expr_fdr_df_list_2nd) == 0) {
    stop("✗ 错误：expr_fdr_df_list_2nd 为空，无法进行 SubSG 注释")
  }
  if (is.null(annotation_references)) {
    stop("✗ 错误：annotation_references 为空，请先运行 Module 3")
  }

  if (is.null(target_column)) {
    target_column <- if (mode == "subNucleolus") {
      "Sub_HPA_Localization"
    } else {
      "MultiBait_Localization"
    }
  }

  if (is.null(levels_order)) {
    levels_order <- if (mode == "subNucleolus") {
      c(
        "Other",
        "Mitochondrion",
        "Cytosol",
        "Nuclear_Cytosol",
        "Nuclear",
        "Cytosol&Nucleolus",
        "Nuclear&Cytosol&Nucleolus",
        "Other&Nucleolus",
        "Nuclear&Nucleolus",
        "Nucleolus"
      )
    } else {
      c(
        "Other",
        "Mitochondrion",
        "Nuclear",
        "Nuclear_Cytosol",
        "Cytosol",
        "Nuclear&SGs",
        "Nuclear&Cytosol&SGs",
        "SGs&Other",
        "Cytosol&SGs"
      )
    }
  }

  if (mode == "subSG") {
    halo_refs <- annotation_references$SGs_refs
    if (is.null(halo_refs$HaloMap_SGs)) {
      stop("✗ 错误：未找到 HaloMap_SGs 参考数据，请检查 Reference 目录或 Module 3 配置")
    }
  }

  halomap_sg_genes <- if (mode == "subSG") halo_refs$HaloMap_SGs$Known.SG.Reference else NULL

  hpa_mode <- if (mode == "subNucleolus") "Nucleolus" else "SG"
  hpa_sets <- load_HPA_annotations(dir_config, mode = hpa_mode)

  nuclear_genes <- hpa_sets$Nuclear$Gene
  cytosol_genes <- hpa_sets$Cytosol$Gene
  nuclear_cytosol_genes <- hpa_sets$Nuclear_Cytosol$Gene
  mito_genes <- annotation_references$MitoCarta_anno$Symbol

  subn_sets <- NULL
  if (mode == "subNucleolus") {
    subn_sets <- build_subnucleolus_sets(dir_config, hpa_sets)
  }

  annotated_list <- list()
  output_files <- list()

  for (dataset_name in names(expr_fdr_df_list_2nd)) {
    cat(sprintf("[Module13] 处理数据集：%s\n", dataset_name))
    df <- expr_fdr_df_list_2nd[[dataset_name]]

    if (!("Gene" %in% colnames(df))) {
      warning(sprintf("⚠ 警告：%s 缺少 Gene 列，跳过", dataset_name))
      next
    }

    target_exists <- target_column %in% colnames(df)

    mutated_df <- if (mode == "subNucleolus") {
      annotate_subnucleolus_df(
        df,
        target_column,
        subn_sets,
        mito_genes
      )
    } else {
      annotate_subsg_df(
        df,
        target_column,
        halomap_sg_genes,
        nuclear_genes,
        cytosol_genes,
        nuclear_cytosol_genes,
        mito_genes
      )
    }

    mutated_df[[target_column]] <- factor(
      mutated_df[[target_column]],
      levels = levels_order
    )

    if (!target_exists) {
      loc_cols <- grep("_Localization$", colnames(mutated_df), value = TRUE)
      loc_cols <- loc_cols[loc_cols != target_column]

      if (length(loc_cols) > 0) {
        last_loc <- loc_cols[length(loc_cols)]
        mutated_df <- mutated_df %>%
          relocate(all_of(target_column), .after = last_loc)
      } else {
        mutated_df <- mutated_df %>%
          relocate(all_of(target_column), .after = "Gene")
      }
    }

    annotated_list[[dataset_name]] <- mutated_df

    output_file <- file.path(
      dir_config$output,
      sprintf(
        "%s_%s_%s.csv",
        output_prefix,
        dataset_name,
        ifelse(mode == "subNucleolus", "SubNucleolus", "SubSGs")
      )
    )
    write.csv(mutated_df, output_file, row.names = FALSE)
    output_files[[dataset_name]] <- output_file
    cat(sprintf("  ✓ 导出 %s 注释结果：%s\n", mode_label, basename(output_file)))
  }

  if (length(annotated_list) == 0) {
    warning("⚠ 提示：未生成任何注释结果")
  }

  ForStep16 <- annotated_list
  if (!is.null(forstep16_versions)) {
    existing <- names(annotated_list)
    keep <- forstep16_versions[forstep16_versions %in% existing]
    missing <- setdiff(forstep16_versions, existing)
    if (length(missing) > 0) {
      warning(sprintf("⚠ 警告：以下数据集不存在，无法加入 ForStep16：%s",
                      paste(missing, collapse = ", ")))
    }
    if (length(keep) > 0) {
      ForStep16 <- annotated_list[keep]
    } else {
      ForStep16 <- list()
    }
  }
  if (length(ForStep16) == 0) {
    ForStep16 <- annotated_list
  }

  cat("\n=== Module 13 完成 ===\n")
  cat(sprintf("✓ 处理数据集数量：%d\n", length(annotated_list)))

  return(list(
    Expr_FDR_df_list_2nd_SubSGs = annotated_list,
    ForStep16 = ForStep16,
    output_files = output_files,
    mode = mode,
    target_column = target_column
  ))
}


annotate_subsg_df <- function(df,
                              target_column,
                              halomap_sg_genes,
                              nuclear_genes,
                              cytosol_genes,
                              nuclear_cytosol_genes,
                              mito_genes) {
  df %>%
    mutate(
      !!target_column := case_when(
        Gene %in% halomap_sg_genes & Gene %in% nuclear_cytosol_genes ~ "Nuclear&Cytosol&SGs",
        Gene %in% halomap_sg_genes & Gene %in% nuclear_genes ~ "Nuclear&SGs",
        Gene %in% halomap_sg_genes & Gene %in% cytosol_genes ~ "Cytosol&SGs",
        Gene %in% halomap_sg_genes ~ "SGs&Other",
        Gene %in% nuclear_cytosol_genes ~ "Nuclear_Cytosol",
        Gene %in% nuclear_genes ~ "Nuclear",
        Gene %in% cytosol_genes ~ "Cytosol",
        Gene %in% mito_genes ~ "Mitochondrion",
        TRUE ~ "Other"
      )
    )
}


build_subnucleolus_sets <- function(dir_config, hpa_sets) {
  hpa_file <- file.path(dir_config$reference, "proteinatlas.tsv")
  if (!file.exists(hpa_file)) {
    stop(sprintf("✗ 错误：未找到HPA文件 - %s", hpa_file))
  }

  cat("  - 读取HPA数据构建SubNucleolus参考\n")
  hpa_data <- readr::read_tsv(hpa_file, show_col_types = FALSE)
  if (!("Reliability (IF)" %in% names(hpa_data))) {
    stop("✗ 错误：HPA数据缺少 'Reliability (IF)' 列")
  }

  hpa_data <- hpa_data %>%
    filter(!`Reliability (IF)` == "Uncertain")

  Independent_Nuclear_Cytosol <- hpa_data %>%
    filter(str_detect(`Subcellular location`, "itochond|yto|crotubule|ctin|Intermediate|Plasma|Vesicles|Golgi|Endoplasmic") &
             str_detect(`Subcellular location`, "ucleoplasm|uclear")) %>%
    filter(!str_detect(`Subcellular main location`, "membrane"))

  Independent_Nuclear <- hpa_data %>%
    filter(str_detect(`Subcellular main location`, "ucleoplasm|uclear")) %>%
    filter(!str_detect(`Subcellular location`, "itochond|yto|crotubule|ctin|Intermediate")) %>%
    filter(!str_detect(`Subcellular location`, "Plasma|Vesicles|Golgi|Endoplasmic")) %>%
    filter(!str_detect(`Subcellular main location`, "membrane"))

  Independent_Cytosol <- hpa_data %>%
    filter(str_detect(`Subcellular main location`, "itochond|yto|crotubule|ctin|Intermediate|Plasma|Vesicles|Golgi|Endoplasmic")) %>%
    filter(!str_detect(`Subcellular location`, "ucleoplasm")) %>%
    filter(!str_detect(`Subcellular location`, "uclear")) %>%
    filter(!str_detect(`Subcellular main location`, "membrane"))

  nucleolus_df <- hpa_sets$Nucleolus
  nucleolus_main_df <- hpa_sets$Nucleolus_main
  if (is.null(nucleolus_df) || is.null(nucleolus_main_df)) {
    stop("✗ 错误：未在 annotation_references 中找到 Nucleolus 注释，请在 Module 3 中启用相应模式")
  }

  list(
    nucleolus_genes = unique(nucleolus_df$Gene),
    nucleolus_main_genes = unique(nucleolus_main_df$Gene),
    independent_nuclear_cytosol = unique(Independent_Nuclear_Cytosol$Gene),
    independent_nuclear = unique(Independent_Nuclear$Gene),
    independent_cytosol = unique(Independent_Cytosol$Gene)
  )
}


annotate_subnucleolus_df <- function(df,
                                     target_column,
                                     subn_sets,
                                     mito_genes) {
  nucleolus_genes <- subn_sets$nucleolus_genes
  nucleolus_main_genes <- subn_sets$nucleolus_main_genes
  indep_nc_genes <- subn_sets$independent_nuclear_cytosol
  indep_nuclear <- subn_sets$independent_nuclear
  indep_cytosol <- subn_sets$independent_cytosol

  df %>%
    mutate(
      !!target_column := case_when(
        Gene %in% nucleolus_genes & Gene %in% indep_nc_genes ~ "Nuclear&Cytosol&Nucleolus",
        Gene %in% nucleolus_genes & Gene %in% indep_nuclear ~ "Nuclear&Nucleolus",
        Gene %in% nucleolus_genes & Gene %in% indep_cytosol ~ "Cytosol&Nucleolus",
        Gene %in% nucleolus_main_genes ~ "Nucleolus",
        Gene %in% nucleolus_genes ~ "Other&Nucleolus",
        Gene %in% indep_nc_genes ~ "Nuclear_Cytosol",
        Gene %in% indep_nuclear ~ "Nuclear",
        Gene %in% indep_cytosol ~ "Cytosol",
        Gene %in% mito_genes ~ "Mitochondrion",
        TRUE ~ "Other"
      )
    )
}

