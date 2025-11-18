# ============================================================================
# Module 3: 注释系统
# ============================================================================
# 功能：
#   1. 读取HPA数据库，生成Cytosol/Nuclear/Nuclear_Cytosol注释
#   2. 读取MitoCarta数据，生成Mitochondrion注释
#   3. 读取用户指定的TP参考数据（如GO_SGs, HaloMap等）
#   4. 根据用户配置添加TP注释列（所有TP注释列需手动配置）
#   5. 如未配置TP注释，则只保留基础注释（Cytosol/Nuclear/Mitochondrion）
# 
# 输入：dir_config, data_raw, sampleGroup, custom_annotations
# 输出：Module03_workspace.RData, Module03_data_annotated.csv
# ============================================================================

#' 读取并处理HPA数据
#' 
#' @param dir_config 目录配置
#' @return 包含Cytosol/Nuclear/Nuclear_Cytosol注释的列表
load_HPA_annotations <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("步骤1: 读取HPA数据库\n")
  cat("----------------------------------------\n")
  
  hpa_file <- file.path(dir_config$reference, "proteinatlas.tsv")
  
  if (!file.exists(hpa_file)) {
    stop(sprintf("✗ 错误：未找到HPA文件 - %s", hpa_file))
  }
  
  HPA_DATA <- read_tsv(hpa_file, show_col_types = FALSE)
  cat(sprintf("✓ 已读取HPA数据: %d 行\n", nrow(HPA_DATA)))
  
  # 选择需要的列并过滤Uncertain
  HPA_DATA <- HPA_DATA %>% 
    select(1:11,
           `Reliability (IF)`,
           `Subcellular location`,
           `Subcellular main location`,
           `Subcellular additional location`) %>% 
    filter(!`Reliability (IF)` == "Uncertain")
  
  cat(sprintf("✓ 过滤后: %d 行\n", nrow(HPA_DATA)))
  
  # --------------------------------------------------------------------------
  # 生成Cytosol注释
  # --------------------------------------------------------------------------
  # 规则：包含yto/crotubule/ctin，排除ucleol/ucleoplasm/uclear/itochond
  Cytosol_HPA_anno <- HPA_DATA %>% 
    filter(str_detect(`Subcellular main location`, "yto|crotubule|ctin")) %>% 
    filter(!str_detect(`Subcellular location`, "ucleol")) %>% 
    filter(!str_detect(`Subcellular location`, "ucleoplasm")) %>% 
    filter(!str_detect(`Subcellular location`, "uclear")) %>% 
    filter(!str_detect(`Subcellular location`, "itochond"))
  
  cat(sprintf("✓ Cytosol注释: %d 个蛋白\n", nrow(Cytosol_HPA_anno)))
  
  # --------------------------------------------------------------------------
  # 生成Nuclear注释
  # --------------------------------------------------------------------------
  # 规则：包含ucleol/ucleoplasm/uclear，排除yto/crotubule/ctin/itochond/Plasma/Vesicles/Golgi/Endoplasmic
  Nuclear_HPA_anno <- HPA_DATA %>% 
    filter(str_detect(`Subcellular main location`, "ucleol|ucleoplasm|uclear")) %>% 
    filter(!str_detect(`Subcellular location`, "yto|crotubule|ctin")) %>% 
    filter(!str_detect(`Subcellular location`, "itochond")) %>% 
    filter(!str_detect(`Subcellular location`, "Plasma|Vesicles|Golgi|Endoplasmic"))
  
  cat(sprintf("✓ Nuclear注释: %d 个蛋白\n", nrow(Nuclear_HPA_anno)))
  
  # --------------------------------------------------------------------------
  # 生成Nuclear_Cytosol注释
  # --------------------------------------------------------------------------
  # 规则：同时包含yto/crotubule/ctin 和 ucleol/ucleoplasm/uclear，排除itochond
  Nuclear_Cytosol_HPA_anno <- HPA_DATA %>% 
    filter(str_detect(`Subcellular location`, "yto|crotubule|ctin") &
             str_detect(`Subcellular location`, "ucleol|ucleoplasm|uclear")) %>% 
    filter(!str_detect(`Subcellular location`, "itochond"))
  
  cat(sprintf("✓ Nuclear_Cytosol注释: %d 个蛋白\n", nrow(Nuclear_Cytosol_HPA_anno)))
  
  # 检查交集
  cyto_nuclear_intersect <- intersect(Cytosol_HPA_anno$Gene, Nuclear_HPA_anno$Gene)
  cyto_nc_intersect <- intersect(Cytosol_HPA_anno$Gene, Nuclear_Cytosol_HPA_anno$Gene)
  nuclear_nc_intersect <- intersect(Nuclear_HPA_anno$Gene, Nuclear_Cytosol_HPA_anno$Gene)
  
  cat(sprintf("\n交集检查:\n"))
  cat(sprintf("  Cytosol ∩ Nuclear: %d\n", length(cyto_nuclear_intersect)))
  cat(sprintf("  Cytosol ∩ Nuclear_Cytosol: %d\n", length(cyto_nc_intersect)))
  cat(sprintf("  Nuclear ∩ Nuclear_Cytosol: %d\n", length(nuclear_nc_intersect)))
  
  return(list(
    Cytosol = Cytosol_HPA_anno,
    Nuclear = Nuclear_HPA_anno,
    Nuclear_Cytosol = Nuclear_Cytosol_HPA_anno
  ))
}


#' 读取MitoCarta数据
#' 
#' @param dir_config 目录配置
#' @return MitoCarta注释数据框
load_MitoCarta_annotations <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("步骤2: 读取MitoCarta数据\n")
  cat("----------------------------------------\n")
  
  mito_file <- file.path(dir_config$reference, "MitoCarta3.0.csv")
  
  if (!file.exists(mito_file)) {
    stop(sprintf("✗ 错误：未找到MitoCarta文件 - %s", mito_file))
  }
  
  MitoCarta_anno <- read.csv(mito_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ 已读取MitoCarta数据: %d 个蛋白\n", nrow(MitoCarta_anno)))
  
  return(MitoCarta_anno)
}


#' 读取SGs参考数据
#' 
#' @param dir_config 目录配置
#' @return 包含各种SGs参考的列表
load_SGs_references <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("步骤3: 读取SGs参考数据\n")
  cat("----------------------------------------\n")
  
  sgs_refs <- list()
  
  # GO SGs
  go_file <- file.path(dir_config$reference, "CytoSGs GO.tsv")
  if (file.exists(go_file)) {
    sgs_refs$GO_SGs <- read_tsv(go_file, show_col_types = FALSE)
    cat(sprintf("✓ GO_SGs: %d 个蛋白\n", nrow(sgs_refs$GO_SGs)))
  } else {
    warning(sprintf("⚠ 未找到GO SGs文件: %s", go_file))
  }
  
  # HaloMap SGs reference
  halomap_file <- file.path(dir_config$reference, "HaloMap SGs reference.csv")
  if (file.exists(halomap_file)) {
    sgs_refs$HaloMap_SGs <- read.csv(halomap_file, stringsAsFactors = FALSE)
    cat(sprintf("✓ HaloMap_SGs: %d 个蛋白\n", nrow(sgs_refs$HaloMap_SGs)))
  } else {
    warning(sprintf("⚠ 未找到HaloMap SGs文件: %s", halomap_file))
  }
  
  # HaloMap different methods
  halomap_methods_file <- file.path(dir_config$reference, "HaloMap differentMethods SGs.csv")
  if (file.exists(halomap_methods_file)) {
    sgs_refs$HaloMap_DifMethods <- read.csv(halomap_methods_file, stringsAsFactors = FALSE)
    cat(sprintf("✓ HaloMap_DifMethods: %d 个蛋白\n", nrow(sgs_refs$HaloMap_DifMethods)))
  } else {
    warning(sprintf("⚠ 未找到HaloMap不同方法文件: %s", halomap_methods_file))
  }
  
  return(sgs_refs)
}


#' 对数据进行注释
#' 
#' @param data 待注释的数据
#' @param HPA_anno HPA注释列表
#' @param MitoCarta_anno MitoCarta注释
#' @param SGs_refs SGs参考列表
#' @param annotation_configs 注释配置列表，每个元素格式为：
#'   list(column_name = "注释列名", TP_source = "TP数据源", TP_column = "TP列名", TP_label = "TP标签")
#' @return 注释后的数据
annotate_data <- function(data, HPA_anno, MitoCarta_anno, SGs_refs, annotation_configs) {
  
  cat("\n----------------------------------------\n")
  cat("步骤4: 数据注释\n")
  cat("----------------------------------------\n")
  
  data_annotated <- data
  
  for (config in annotation_configs) {
    col_name <- config$column_name
    tp_source <- config$TP_source
    tp_column <- config$TP_column
    tp_label <- config$TP_label
    
    cat(sprintf("\n添加注释列: %s\n", col_name))
    cat(sprintf("  TP来源: %s$%s\n", tp_source, tp_column))
    
    # 获取TP基因列表
    if (tp_source == "GO_SGs") {
      tp_genes <- SGs_refs$GO_SGs[[tp_column]]
    } else if (tp_source == "HaloMap_SGs") {
      tp_genes <- SGs_refs$HaloMap_SGs[[tp_column]]
    } else if (tp_source == "HaloMap_DifMethods") {
      tp_genes <- SGs_refs$HaloMap_DifMethods[[tp_column]]
    } else {
      stop(sprintf("✗ 错误：未知的TP来源 - %s", tp_source))
    }
    
    # 进行注释（优先级：TP > Nuclear/Cytosol/Nuclear_Cytosol > Mitochondrion > Other）
    data_annotated <- data_annotated %>%
      mutate(!!col_name := case_when(
        Gene %in% tp_genes ~ tp_label,
        Gene %in% HPA_anno$Nuclear$Gene ~ "Nuclear",
        Gene %in% HPA_anno$Cytosol$Gene ~ "Cytosol",
        Gene %in% HPA_anno$Nuclear_Cytosol$Gene ~ "Nuclear_Cytosol",
        Gene %in% MitoCarta_anno$Symbol ~ "Mitochondrion",
        TRUE ~ "Other"
      ))
    
    # 统计
    anno_counts <- table(data_annotated[[col_name]])
    cat(sprintf("  注释统计:\n"))
    for (label in names(anno_counts)) {
      cat(sprintf("    %s: %d\n", label, anno_counts[label]))
    }
  }
  
  cat("\n✓ 数据注释完成\n")
  return(data_annotated)
}


#' Module 3 主函数
#' 
#' @param dir_config 目录配置
#' @param data_raw 原始数据
#' @param sampleGroup 分组信息（可选，用于后续分析）
#' @param custom_annotations 用户自定义注释配置（可选）
module03_annotation <- function(dir_config, data_raw, sampleGroup = NULL, custom_annotations = NULL) {
  
  cat("\n========================================\n")
  cat("Module 3: 注释系统\n")
  cat("========================================\n")
  
  # --------------------------------------------------------------------------
  # 加载参考数据
  # --------------------------------------------------------------------------
  HPA_anno <- load_HPA_annotations(dir_config)
  MitoCarta_anno <- load_MitoCarta_annotations(dir_config)
  SGs_refs <- load_SGs_references(dir_config)
  
  # --------------------------------------------------------------------------
  # 用户自定义注释配置（所有TP注释列需手动配置）
  # --------------------------------------------------------------------------
  if (is.null(custom_annotations) || length(custom_annotations) == 0) {
    cat("⚠ 未配置TP注释，将只使用基础注释（Cytosol/Nuclear/Mitochondrion）\n")
    annotation_configs <- list()
  } else {
    annotation_configs <- custom_annotations
    cat(sprintf("✓ 已配置 %d 个TP注释列\n", length(annotation_configs)))
  }
  
  # --------------------------------------------------------------------------
  # 执行注释
  # --------------------------------------------------------------------------
  data_annotated <- annotate_data(data_raw, HPA_anno, MitoCarta_anno, SGs_refs, annotation_configs)
  
  # --------------------------------------------------------------------------
  # 保存数据
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("步骤5: 保存数据\n")
  cat("----------------------------------------\n")
  
  setwd(dir_config$root)
  
  # 保存参考数据（供后续模块使用）
  annotation_references <- list(
    HPA_anno = HPA_anno,
    MitoCarta_anno = MitoCarta_anno,
    SGs_refs = SGs_refs
  )
  
  # 保存CSV到Output目录
  output_file <- file.path(dir_config$output, "Module03_data_annotated.csv")
  write.csv(data_annotated, output_file, row.names = FALSE)
  cat(sprintf("✓ 已保存: %s (Output目录)\n", output_file))
  
  cat("\n✓ Module 3 数据处理完成，返回数据到主程序\n")
  
  return(list(
    data_annotated = data_annotated,
    annotation_references = annotation_references
  ))
}

