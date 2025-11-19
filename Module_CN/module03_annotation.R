# ============================================================================
# Module 3: 注释系统
# ============================================================================
# 功能：
#   1. 读取HPA数据库，生成Cytosol/Nuclear/Nuclear_Cytosol注释
#   2. 读取MitoCarta数据，生成Mitochondrion注释
#   3. 读取用户指定的TP参考数据（如GO_SGs, HaloMap等）
#   4. 支持额外的核仁与自定义TP参考文件
#   5. 根据预设模式或用户配置添加TP注释列
#   6. 如未配置TP注释，则只保留基础注释（Cytosol/Nuclear/Mitochondrion）
#
# 输入：dir_config, data_raw, sampleGroup, custom_annotations, custom_tp_sources
# 输出：Module03_workspace.RData, Module03_data_annotated.csv
# ============================================================================

#' 读取并处理HPA数据
#' 
#' @param dir_config 目录配置
#' @param mode 注释模式（"SG" 或 "Nucleolus"），影响HPA筛选逻辑
#' @return 包含Cytosol/Nuclear/Nuclear_Cytosol注释的列表
load_HPA_annotations <- function(dir_config, mode = c("SG", "Nucleolus")) {
  mode <- match.arg(mode)
  
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
  
  if (mode == "Nucleolus") {
    # ------------------------------------------------------------------------
    # 核仁模式：严格遵循20250725.R的定义
    # ------------------------------------------------------------------------
    Cytosol_HPA_anno <- HPA_DATA %>% 
      filter(str_detect(`Subcellular main location`, "yto|crotubule|ctin|Intermediate")) %>% 
      filter(!str_detect(`Subcellular location`, "ucleol")) %>% 
      filter(!str_detect(`Subcellular location`, "ucleoplasm")) %>% 
      filter(!str_detect(`Subcellular location`, "uclear")) %>% 
      filter(!str_detect(`Subcellular location`, "itochond"))
    
    Nuclear_HPA_anno <- HPA_DATA %>% 
      filter(str_detect(`Subcellular main location`, "ucleoplasm|uclear")) %>% 
      filter(!str_detect(`Subcellular location`, "yto|crotubule|ctin|Intermediate")) %>% 
      filter(!str_detect(`Subcellular location`, "itochond")) %>% 
      filter(!str_detect(`Subcellular location`, "Plasma|Vesicles|Golgi|Endoplasmic")) %>% 
      filter(!str_detect(`Subcellular location`, "ucleol"))
    
    Nuclear_Cytosol_HPA_anno <- HPA_DATA %>% 
      filter(str_detect(`Subcellular location`, "yto|crotubule|ctin|Intermediate") &
               str_detect(`Subcellular location`, "ucleoplasm|uclear")) %>% 
      filter(!str_detect(`Subcellular location`, "itochond"))
    
  } else {
    # ------------------------------------------------------------------------
    # SG模式：保留原module03的默认定义
    # ------------------------------------------------------------------------
    Cytosol_HPA_anno <- HPA_DATA %>% 
      filter(str_detect(`Subcellular main location`, "yto|crotubule|ctin")) %>% 
      filter(!str_detect(`Subcellular location`, "ucleol")) %>% 
      filter(!str_detect(`Subcellular location`, "ucleoplasm")) %>% 
      filter(!str_detect(`Subcellular location`, "uclear")) %>% 
      filter(!str_detect(`Subcellular location`, "itochond"))
    
    Nuclear_HPA_anno <- HPA_DATA %>% 
      filter(str_detect(`Subcellular main location`, "ucleol|ucleoplasm|uclear")) %>% 
      filter(!str_detect(`Subcellular location`, "yto|crotubule|ctin")) %>% 
      filter(!str_detect(`Subcellular location`, "itochond")) %>% 
      filter(!str_detect(`Subcellular location`, "Plasma|Vesicles|Golgi|Endoplasmic"))
    
    Nuclear_Cytosol_HPA_anno <- HPA_DATA %>% 
      filter(str_detect(`Subcellular location`, "yto|crotubule|ctin") &
               str_detect(`Subcellular location`, "ucleol|ucleoplasm|uclear")) %>% 
      filter(!str_detect(`Subcellular location`, "itochond"))
  }
  
  cat(sprintf("✓ Cytosol注释: %d 个蛋白\n", nrow(Cytosol_HPA_anno)))
  cat(sprintf("✓ Nuclear注释: %d 个蛋白\n", nrow(Nuclear_HPA_anno)))
  cat(sprintf("✓ Nuclear_Cytosol注释: %d 个蛋白\n", nrow(Nuclear_Cytosol_HPA_anno)))
  
  # --------------------------------------------------------------------------
  # 生成Nucleolus注释（main location + 所有location）
  # --------------------------------------------------------------------------
  Nucleolus_main_Anno <- HPA_DATA %>%
    filter(str_detect(`Subcellular main location`, "ucleol"))
  cat(sprintf("✓ Nucleolus(main location)注释: %d 个蛋白\n", nrow(Nucleolus_main_Anno)))
  
  Nucleolus_Anno <- HPA_DATA %>%
    filter(str_detect(`Subcellular location`, "ucleol"))
  cat(sprintf("✓ Nucleolus(all location)注释: %d 个蛋白\n", nrow(Nucleolus_Anno)))
  
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
    Nuclear_Cytosol = Nuclear_Cytosol_HPA_anno,
    Nucleolus = Nucleolus_Anno,
    Nucleolus_main = Nucleolus_main_Anno
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


#' 读取核仁参考数据（CLL等）
#'
#' @param dir_config 目录配置
#' @return 包含核仁相关参考的列表
load_nucleolus_references <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("附加参考: 读取核仁注释（CLL等）\n")
  cat("----------------------------------------\n")
  
  candidate_files <- c("CLL_Nucleolus_subNucleolus.csv", "CLL_Nucleolus.csv")
  nucleolus_refs <- list()
  cll_file <- NULL
  
  for (fname in candidate_files) {
    file_path <- file.path(dir_config$reference, fname)
    if (file.exists(file_path)) {
      cll_file <- file_path
      break
    }
  }
  
  if (is.null(cll_file)) {
    warning("⚠ 未找到CLL核仁参考文件（CLL_Nucleolus*.csv），将跳过CLL_Localization注释")
    return(nucleolus_refs)
  }
  
  cll_anno <- read.csv(cll_file, stringsAsFactors = FALSE, check.names = FALSE)
  gene_col <- intersect(c("Gene", "CLL_Nucleolus"), names(cll_anno))
  if (length(gene_col) == 0) {
    warning(sprintf("⚠ CLL核仁参考缺少Gene列: %s", basename(cll_file)))
  } else {
    if (!"Gene" %in% names(cll_anno)) {
      names(cll_anno)[names(cll_anno) == gene_col[1]] <- "Gene"
    }
    nucleolus_refs$CLL_Nucleolus <- cll_anno
    cat(sprintf("✓ CLL_Nucleolus: %d 个蛋白\n", nrow(cll_anno)))
  }
  
  return(nucleolus_refs)
}


#' 读取自定义TP参考文件
#'
#' @param custom_tp_sources 列表，每个元素包含source_name、file、file_type、gene_column
#' @param dir_config 目录配置
#' @return 命名列表（source_name -> 数据框）
load_custom_tp_sources <- function(custom_tp_sources, dir_config) {
  
  if (is.null(custom_tp_sources) || length(custom_tp_sources) == 0) {
    return(list())
  }
  
  cat("\n----------------------------------------\n")
  cat("附加参考: 读取自定义TP文件\n")
  cat("----------------------------------------\n")
  
  custom_refs <- list()
  
  for (i in seq_along(custom_tp_sources)) {
    cfg <- custom_tp_sources[[i]]
    
    if (is.null(cfg$source_name) || cfg$source_name == "") {
      stop(sprintf("✗ 错误：第%d个自定义TP缺少source_name", i))
    }
    if (is.null(cfg$file) || cfg$file == "") {
      stop(sprintf("✗ 错误：自定义TP %s 缺少file路径", cfg$source_name))
    }
    
    source_name <- cfg$source_name
    file_path <- cfg$file
    
    if (!file.exists(file_path)) {
      candidate <- file.path(dir_config$reference, file_path)
      if (file.exists(candidate)) {
        file_path <- candidate
      }
    }
    
    if (!file.exists(file_path)) {
      stop(sprintf("✗ 错误：自定义TP %s 的文件不存在 - %s", source_name, cfg$file))
    }
    
    file_type <- cfg$file_type
    if (is.null(file_type) || file_type == "") {
      ext <- tolower(tools::file_ext(file_path))
      file_type <- ifelse(ext == "", "csv", ext)
    }
    file_type <- tolower(file_type)
    
    gene_column <- if (!is.null(cfg$gene_column)) cfg$gene_column else "Gene"
    
    data <- switch(
      file_type,
      csv = read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE),
      tsv = read_tsv(file_path, show_col_types = FALSE),
      txt = read_tsv(file_path, show_col_types = FALSE),
      xlsx = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("✗ 错误：需要readxl包以读取xlsx文件")
        }
        readxl::read_excel(file_path)
      },
      xls = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("✗ 错误：需要readxl包以读取xls文件")
        }
        readxl::read_excel(file_path)
      },
      stop(sprintf("✗ 错误：不支持的文件类型 %s（source: %s）", file_type, source_name))
    )
    
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    
    if (!(gene_column %in% colnames(data))) {
      stop(sprintf("✗ 错误：自定义TP %s 缺少基因列 %s", source_name, gene_column))
    }
    
    cat(sprintf("✓ 自定义TP: %s (%d 行, 文件: %s)\n", source_name, nrow(data), basename(file_path)))
    custom_refs[[source_name]] <- data
  }
  
  return(custom_refs)
}


#' 内置注释配置
#'
#' @param mode 模式：SG 或 Nucleolus
#' @return 注释配置列表
get_preset_annotation_configs <- function(mode = c("SG", "Nucleolus")) {
  mode <- match.arg(mode)
  
  if (mode == "SG") {
    return(list(
      list(
        column_name = "HaloMap_Localization",
        TP_source = "HaloMap_SGs",
        TP_column = "Known.SG.Reference",
        TP_label = "SGs"
      ),
      list(
        column_name = "GO_Localization",
        TP_source = "GO_SGs",
        TP_column = "Gene",
        TP_label = "SGs"
      ),
      list(
        column_name = "MultiBait_Localization",
        TP_source = "HaloMap_DifMethods",
        TP_column = "Multi.Bait.APEX.Marmor.Kollet.et..al.2020",
        TP_label = "SGs"
      )
    ))
  }
  
  if (mode == "Nucleolus") {
    return(list(
      list(
        column_name = "Nucleolus_Localization",
        TP_source = "HPA_Nucleolus",
        TP_column = "Gene",
        TP_label = "Nucleolus"
      ),
      list(
        column_name = "Nucleolus_main_Localization",
        TP_source = "HPA_Nucleolus_Main",
        TP_column = "Gene",
        TP_label = "Nucleolus"
      ),
      list(
        column_name = "CLL_Localization",
        TP_source = "CLL_Nucleolus",
        TP_column = "Gene",
        TP_label = "Nucleolus"
      )
    ))
  }
  
  return(list())
}


#' 构建TP查找表
#' 
#' @param SGs_refs SGs参考列表
#' @param HPA_anno HPA注释
#' @param nucleolus_refs 核仁参考
#' @return 名称到数据框的映射
build_tp_lookup <- function(SGs_refs, HPA_anno, nucleolus_refs, custom_tp_refs = list()) {
  tp_lookup <- list()
  
  add_dataset <- function(name, dataset) {
    if (!is.null(dataset) && NROW(dataset) > 0) {
      tp_lookup[[name]] <<- dataset
    }
  }
  
  if (!is.null(SGs_refs) && length(SGs_refs) > 0) {
    for (nm in names(SGs_refs)) {
      add_dataset(nm, SGs_refs[[nm]])
    }
  }
  
  add_dataset("HPA_Nucleolus", HPA_anno$Nucleolus)
  add_dataset("HPA_Nucleolus_Main", HPA_anno$Nucleolus_main)
  
  if (!is.null(nucleolus_refs) && length(nucleolus_refs) > 0) {
    for (nm in names(nucleolus_refs)) {
      add_dataset(nm, nucleolus_refs[[nm]])
    }
  }
  
  if (!is.null(custom_tp_refs) && length(custom_tp_refs) > 0) {
    for (nm in names(custom_tp_refs)) {
      add_dataset(nm, custom_tp_refs[[nm]])
    }
  }
  
  return(tp_lookup)
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
annotate_data <- function(data, HPA_anno, MitoCarta_anno, annotation_configs, tp_lookup = list()) {
  
  cat("\n----------------------------------------\n")
  cat("步骤4: 数据注释\n")
  cat("----------------------------------------\n")
  
  data_annotated <- data
  
  get_tp_genes <- function(tp_source, tp_column) {
    tp_data <- tp_lookup[[tp_source]]
    if (is.null(tp_data)) {
      warning(sprintf("⚠ 未找到TP来源：%s，列将仅包含基础定位", tp_source))
      return(character(0))
    }
    if (!(tp_column %in% colnames(tp_data))) {
      warning(sprintf("⚠ TP来源 %s 缺少列 %s，列将仅包含基础定位", tp_source, tp_column))
      return(character(0))
    }
    genes <- tp_data[[tp_column]]
    genes <- genes[!is.na(genes)]
    return(unique(as.character(genes)))
  }
  
  for (config in annotation_configs) {
    col_name <- config$column_name
    tp_source <- config$TP_source
    tp_column <- config$TP_column
    tp_label <- config$TP_label
    
    cat(sprintf("\n添加注释列: %s\n", col_name))
    cat(sprintf("  TP来源: %s$%s\n", tp_source, tp_column))
    
    tp_genes <- get_tp_genes(tp_source, tp_column)
    
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
#' @param custom_annotations 用户自定义注释配置（可选，非空时覆盖预设）
#' @param annotation_mode 预设注释模式（"SG" 或 "Nucleolus"）
#' @param custom_tp_sources 自定义TP参考文件配置（可选）
#' @param additional_annotations 附加的注释配置（在预设/自定义基础上追加）
module03_annotation <- function(dir_config,
                                data_raw,
                                sampleGroup = NULL,
                                custom_annotations = NULL,
                                annotation_mode = c("SG", "Nucleolus"),
                                custom_tp_sources = NULL,
                                additional_annotations = NULL) {
  
  cat("\n========================================\n")
  cat("Module 3: 注释系统\n")
  cat("========================================\n")
  annotation_mode <- match.arg(annotation_mode)
  cat(sprintf("当前注释模式: %s\n", annotation_mode))
  
  # --------------------------------------------------------------------------
  # 加载参考数据
  # --------------------------------------------------------------------------
  HPA_anno <- load_HPA_annotations(dir_config, annotation_mode)
  MitoCarta_anno <- load_MitoCarta_annotations(dir_config)
  SGs_refs <- load_SGs_references(dir_config)
  custom_tp_refs <- load_custom_tp_sources(custom_tp_sources, dir_config)
  
  # --------------------------------------------------------------------------
  # 注释模式与配置
  # --------------------------------------------------------------------------
  use_custom <- !is.null(custom_annotations) && length(custom_annotations) > 0
  if (use_custom) {
    annotation_configs <- custom_annotations
    cat(sprintf("✓ 使用自定义注释配置（%d 列）\n", length(annotation_configs)))
  } else {
    annotation_configs <- get_preset_annotation_configs(annotation_mode)
    if (length(annotation_configs) == 0) {
      cat("⚠ 未配置TP注释，将只使用基础注释（Cytosol/Nuclear/Mitochondrion）\n")
    } else {
      cat(sprintf("✓ 已加载预设注释配置（模式: %s, 列数: %d）\n",
                  annotation_mode, length(annotation_configs)))
    }
  }
  
  has_additional <- !is.null(additional_annotations) && length(additional_annotations) > 0
  if (has_additional) {
    cat(sprintf("✓ 追加自定义注释列（%d 列）\n", length(additional_annotations)))
    annotation_configs <- c(annotation_configs, additional_annotations)
  }
  
  needs_cll <- any(vapply(annotation_configs, function(cfg) {
    identical(cfg$TP_source, "CLL_Nucleolus")
  }, logical(1)))
  if (!needs_cll) {
    nucleolus_refs <- list()
  } else {
    nucleolus_refs <- load_nucleolus_references(dir_config)
  }
  
  tp_lookup <- build_tp_lookup(SGs_refs, HPA_anno, nucleolus_refs, custom_tp_refs)
  
  # --------------------------------------------------------------------------
  # 用户自定义注释配置（所有TP注释列需手动配置）
  # --------------------------------------------------------------------------
  if (length(annotation_configs) == 0) {
    cat("  - 未提供注释配置，结果将保留原始列\n")
  }

  # --------------------------------------------------------------------------
  # 执行注释
  # --------------------------------------------------------------------------
  data_annotated <- annotate_data(
    data_raw,
    HPA_anno,
    MitoCarta_anno,
    annotation_configs,
    tp_lookup
  )
  
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
    MitoCarta = MitoCarta_anno,
    SGs_refs = SGs_refs,
    nucleolus_refs = nucleolus_refs,
    custom_tp_refs = custom_tp_refs,
    annotation_mode = if (use_custom) "custom" else annotation_mode,
    annotation_configs = annotation_configs
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

