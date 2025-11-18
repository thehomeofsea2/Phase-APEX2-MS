## ============================================================================
## Module 13: SubSG 注释模块
## 对齐 CleanCode.R (2708-2742) 逻辑
## ============================================================================

module13_subsg_annotation <- function(
    dir_config,
    expr_fdr_df_list_2nd,
    annotation_references,
    target_column = "MultiBait_Localization",
    output_prefix = "Module13",
    levels_order = c(
      "Other",
      "Mitochondrion",
      "Nuclear",
      "Nuclear_Cytosol",
      "Cytosol",
      "Nuclear&SGs",
      "Nuclear&Cytosol&SGs",
      "SGs&Other",
      "Cytosol&SGs"
    ),
    forstep16_versions = NULL
) {

  cat("\n=== Module 13: SubSG 注释 ===\n\n")

  required_pkgs <- c("dplyr", "stringr", "purrr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("✗ 错误：需要安装 %s 包", pkg))
    }
  }

  library(dplyr)
  library(stringr)
  library(purrr)

  if (length(expr_fdr_df_list_2nd) == 0) {
    stop("✗ 错误：expr_fdr_df_list_2nd 为空，无法进行 SubSG 注释")
  }
  if (is.null(annotation_references)) {
    stop("✗ 错误：annotation_references 为空，请先运行 Module 3")
  }

  halo_refs <- annotation_references$SGs_refs
  if (is.null(halo_refs$HaloMap_SGs)) {
    stop("✗ 错误：未找到 HaloMap_SGs 参考数据，请检查 Reference 目录或 Module 3 配置")
  }

  halomap_sg_genes <- halo_refs$HaloMap_SGs$Known.SG.Reference
  nuclear_genes <- annotation_references$HPA_anno$Nuclear$Gene
  cytosol_genes <- annotation_references$HPA_anno$Cytosol$Gene
  nuclear_cytosol_genes <- annotation_references$HPA_anno$Nuclear_Cytosol$Gene
  mito_genes <- annotation_references$MitoCarta_anno$Symbol

  Expr_FDR_df_list_2nd_SubSGs <- list()
  output_files <- list()

  for (dataset_name in names(expr_fdr_df_list_2nd)) {
    cat(sprintf("[Module13] 处理数据集：%s\n", dataset_name))
    df <- expr_fdr_df_list_2nd[[dataset_name]]

    if (!("Gene" %in% colnames(df))) {
      warning(sprintf("⚠ 警告：%s 缺少 Gene 列，跳过", dataset_name))
      next
    }

    target_exists <- target_column %in% colnames(df)

    mutated_df <- df %>%
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

    Expr_FDR_df_list_2nd_SubSGs[[dataset_name]] <- mutated_df

    output_file <- file.path(
      dir_config$output,
      sprintf("%s_%s_SubSGs.csv", output_prefix, dataset_name)
    )
    write.csv(mutated_df, output_file, row.names = FALSE)
    output_files[[dataset_name]] <- output_file
    cat(sprintf("  ✓ 导出 SubSG 注释结果：%s\n", basename(output_file)))
  }

  if (length(Expr_FDR_df_list_2nd_SubSGs) == 0) {
    warning("⚠ 提示：未生成任何 SubSG 注释结果")
  }

  ForStep16 <- Expr_FDR_df_list_2nd_SubSGs
  if (!is.null(forstep16_versions)) {
    existing <- names(Expr_FDR_df_list_2nd_SubSGs)
    keep <- forstep16_versions[forstep16_versions %in% existing]
    missing <- setdiff(forstep16_versions, existing)
    if (length(missing) > 0) {
      warning(sprintf("⚠ 警告：以下数据集不存在，无法加入 ForStep16：%s",
                      paste(missing, collapse = ", ")))
    }
    if (length(keep) > 0) {
      ForStep16 <- Expr_FDR_df_list_2nd_SubSGs[keep]
    } else {
      ForStep16 <- list()
    }
  }
  if (length(ForStep16) == 0) {
    ForStep16 <- Expr_FDR_df_list_2nd_SubSGs
  }

  cat("\n=== Module 13 完成 ===\n")
  cat(sprintf("✓ 处理数据集数量：%d\n", length(Expr_FDR_df_list_2nd_SubSGs)))

  return(list(
    Expr_FDR_df_list_2nd_SubSGs = Expr_FDR_df_list_2nd_SubSGs,
    ForStep16 = ForStep16,
    output_files = output_files
  ))
}

