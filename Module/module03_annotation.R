# ============================================================================
# Module 3: Annotation System
# ============================================================================
# Features:
#   1. Read HPA database to generate Cytosol/Nuclear/Nuclear_Cytosol annotations
#   2. Read MitoCarta data to generate Mitochondrion annotations
#   3. Read user-specified TP reference data (e.g., GO_SGs, HaloMap)
#   4. Support extra nucleolus and custom TP reference files
#   5. Add TP annotation columns based on preset mode or user config
#   6. If no TP annotation configured, keep only basic annotations (Cytosol/Nuclear/Mitochondrion)
#
# Input: dir_config, data_raw, sampleGroup, custom_annotations, custom_tp_sources
# Output: Module03_workspace.RData, Module03_data_annotated.csv
# ============================================================================

#' Read and process HPA data
#' 
#' @param dir_config Directory configuration
#' @param mode Annotation mode ("SG" or "Nucleolus"), affects HPA filters
#' @return List of Cytosol/Nuclear/Nuclear_Cytosol annotations
load_HPA_annotations <- function(dir_config, mode = c("SG", "Nucleolus")) {
  mode <- match.arg(mode)
  
  cat("\n----------------------------------------\n")
  cat("Step 1: Read HPA database\n")
  cat("----------------------------------------\n")
  
  hpa_file <- file.path(dir_config$reference, "proteinatlas.tsv")
  
  if (!file.exists(hpa_file)) {
    stop(sprintf("✗ Error: HPA file not found - %s", hpa_file))
  }
  
  HPA_DATA <- read_tsv(hpa_file, show_col_types = FALSE)
  cat(sprintf("✓ Loaded HPA data: %d rows\n", nrow(HPA_DATA)))
  
  # Select needed columns and filter out Uncertain
  HPA_DATA <- HPA_DATA %>% 
    select(1:11,
           `Reliability (IF)`,
           `Subcellular location`,
           `Subcellular main location`,
           `Subcellular additional location`) %>% 
    filter(!`Reliability (IF)` == "Uncertain")
  
  cat(sprintf("✓ After filtering: %d rows\n", nrow(HPA_DATA)))
  
  if (mode == "Nucleolus") {
    # ------------------------------------------------------------------------
    # Nucleolus mode: follow definitions from 20250725.R
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
    # SG mode: keep default definitions from original module03
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
  
  cat(sprintf("✓ Cytosol annotations: %d proteins\n", nrow(Cytosol_HPA_anno)))
  cat(sprintf("✓ Nuclear annotations: %d proteins\n", nrow(Nuclear_HPA_anno)))
  cat(sprintf("✓ Nuclear_Cytosol annotations: %d proteins\n", nrow(Nuclear_Cytosol_HPA_anno)))
  
  # --------------------------------------------------------------------------
  # Build Nucleolus annotations (main location + all locations)
  # --------------------------------------------------------------------------
  Nucleolus_main_Anno <- HPA_DATA %>%
    filter(str_detect(`Subcellular main location`, "ucleol"))
  cat(sprintf("✓ Nucleolus (main location) annotations: %d proteins\n", nrow(Nucleolus_main_Anno)))
  
  Nucleolus_Anno <- HPA_DATA %>%
    filter(str_detect(`Subcellular location`, "ucleol"))
  cat(sprintf("✓ Nucleolus (all locations) annotations: %d proteins\n", nrow(Nucleolus_Anno)))
  
  # Check overlaps
  cyto_nuclear_intersect <- intersect(Cytosol_HPA_anno$Gene, Nuclear_HPA_anno$Gene)
  cyto_nc_intersect <- intersect(Cytosol_HPA_anno$Gene, Nuclear_Cytosol_HPA_anno$Gene)
  nuclear_nc_intersect <- intersect(Nuclear_HPA_anno$Gene, Nuclear_Cytosol_HPA_anno$Gene)
  
  cat(sprintf("\nIntersection checks:\n"))
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


#' Read MitoCarta data
#' 
#' @param dir_config Directory configuration
#' @return MitoCarta annotation data frame
load_MitoCarta_annotations <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("Step 2: Read MitoCarta data\n")
  cat("----------------------------------------\n")
  
  mito_file <- file.path(dir_config$reference, "MitoCarta3.0.csv")
  
  if (!file.exists(mito_file)) {
    stop(sprintf("✗ Error: MitoCarta file not found - %s", mito_file))
  }
  
  MitoCarta_anno <- read.csv(mito_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ Loaded MitoCarta data: %d proteins\n", nrow(MitoCarta_anno)))
  
  return(MitoCarta_anno)
}


#' Read SGs reference data
#' 
#' @param dir_config Directory configuration
#' @return List containing SGs references
load_SGs_references <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("Step 3: Read SGs reference data\n")
  cat("----------------------------------------\n")
  
  sgs_refs <- list()
  
  # GO SGs
  go_file <- file.path(dir_config$reference, "CytoSGs GO.tsv")
  if (file.exists(go_file)) {
    sgs_refs$GO_SGs <- read_tsv(go_file, show_col_types = FALSE)
    cat(sprintf("✓ GO_SGs: %d proteins\n", nrow(sgs_refs$GO_SGs)))
  } else {
    warning(sprintf("⚠ GO SGs file not found: %s", go_file))
  }
  
  # HaloMap SGs reference
  halomap_file <- file.path(dir_config$reference, "HaloMap SGs reference.csv")
  if (file.exists(halomap_file)) {
    sgs_refs$HaloMap_SGs <- read.csv(halomap_file, stringsAsFactors = FALSE)
    cat(sprintf("✓ HaloMap_SGs: %d proteins\n", nrow(sgs_refs$HaloMap_SGs)))
  } else {
    warning(sprintf("⚠ HaloMap SGs file not found: %s", halomap_file))
  }
  
  # HaloMap different methods
  halomap_methods_file <- file.path(dir_config$reference, "HaloMap differentMethods SGs.csv")
  if (file.exists(halomap_methods_file)) {
    sgs_refs$HaloMap_DifMethods <- read.csv(halomap_methods_file, stringsAsFactors = FALSE)
    cat(sprintf("✓ HaloMap_DifMethods: %d proteins\n", nrow(sgs_refs$HaloMap_DifMethods)))
  } else {
    warning(sprintf("⚠ HaloMap different methods file not found: %s", halomap_methods_file))
  }
  
  return(sgs_refs)
}


#' Read nucleolus reference data (CLL, etc.)
#'
#' @param dir_config Directory configuration
#' @return List of nucleolus-related references
load_nucleolus_references <- function(dir_config) {
  
  cat("\n----------------------------------------\n")
  cat("Additional references: read nucleolus annotations (CLL, etc.)\n")
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
    warning("⚠ CLL nucleolus reference not found (CLL_Nucleolus*.csv); skipping CLL_Localization annotation")
    return(nucleolus_refs)
  }
  
  cll_anno <- read.csv(cll_file, stringsAsFactors = FALSE, check.names = FALSE)
  gene_col <- intersect(c("Gene", "CLL_Nucleolus"), names(cll_anno))
  if (length(gene_col) == 0) {
    warning(sprintf("⚠ CLL nucleolus reference missing Gene column: %s", basename(cll_file)))
  } else {
    if (!"Gene" %in% names(cll_anno)) {
      names(cll_anno)[names(cll_anno) == gene_col[1]] <- "Gene"
    }
    nucleolus_refs$CLL_Nucleolus <- cll_anno
    cat(sprintf("✓ CLL_Nucleolus: %d proteins\n", nrow(cll_anno)))
  }
  
  return(nucleolus_refs)
}


#' Read custom TP reference files
#'
#' @param custom_tp_sources List with elements source_name, file, file_type, gene_column
#' @param dir_config Directory configuration
#' @return Named list (source_name -> data frame)
load_custom_tp_sources <- function(custom_tp_sources, dir_config) {
  
  if (is.null(custom_tp_sources) || length(custom_tp_sources) == 0) {
    return(list())
  }
  
  cat("\n----------------------------------------\n")
  cat("Additional references: read custom TP files\n")
  cat("----------------------------------------\n")
  
  custom_refs <- list()
  
  for (i in seq_along(custom_tp_sources)) {
    cfg <- custom_tp_sources[[i]]
    
    if (is.null(cfg$source_name) || cfg$source_name == "") {
      stop(sprintf("✗ Error: custom TP #%d missing source_name", i))
    }
    if (is.null(cfg$file) || cfg$file == "") {
      stop(sprintf("✗ Error: custom TP %s missing file path", cfg$source_name))
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
      stop(sprintf("✗ Error: file for custom TP %s does not exist - %s", source_name, cfg$file))
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
          stop("✗ Error: readxl package required to read xlsx files")
        }
        readxl::read_excel(file_path)
      },
      xls = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("✗ Error: readxl package required to read xls files")
        }
        readxl::read_excel(file_path)
      },
      stop(sprintf("✗ Error: unsupported file type %s (source: %s)", file_type, source_name))
    )
    
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    
    if (!(gene_column %in% colnames(data))) {
      stop(sprintf("✗ Error: custom TP %s missing gene column %s", source_name, gene_column))
    }
    
    cat(sprintf("✓ Custom TP: %s (%d rows, file: %s)\n", source_name, nrow(data), basename(file_path)))
    custom_refs[[source_name]] <- data
  }
  
  return(custom_refs)
}


#' Built-in annotation configurations
#'
#' @param mode Mode: SG or Nucleolus
#' @return Annotation config list
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


#' Build TP lookup table
#' 
#' @param SGs_refs SGs reference list
#' @param HPA_anno HPA annotations
#' @param nucleolus_refs Nucleolus references
#' @return Mapping from name to data frame
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


#' Annotate data
#' 
#' @param data Data to annotate
#' @param HPA_anno HPA annotation list
#' @param MitoCarta_anno MitoCarta annotations
#' @param SGs_refs SGs references
#' @param annotation_configs Annotation configs: list(column_name, TP_source, TP_column, TP_label)
#' @return Annotated data
annotate_data <- function(data, HPA_anno, MitoCarta_anno, annotation_configs, tp_lookup = list()) {
  
  cat("\n----------------------------------------\n")
  cat("Step 4: Data annotation\n")
  cat("----------------------------------------\n")
  
  data_annotated <- data
  
  get_tp_genes <- function(tp_source, tp_column) {
    tp_data <- tp_lookup[[tp_source]]
    if (is.null(tp_data)) {
      warning(sprintf("⚠ TP source not found: %s; column will keep only basic localization", tp_source))
      return(character(0))
    }
    if (!(tp_column %in% colnames(tp_data))) {
      warning(sprintf("⚠ TP source %s missing column %s; column will keep only basic localization", tp_source, tp_column))
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
    
    cat(sprintf("\nAdd annotation column: %s\n", col_name))
    cat(sprintf("  TP source: %s$%s\n", tp_source, tp_column))
    
    tp_genes <- get_tp_genes(tp_source, tp_column)
    
    # Annotate (priority: TP > Nuclear/Cytosol/Nuclear_Cytosol > Mitochondrion > Other)
    data_annotated <- data_annotated %>%
      mutate(!!col_name := case_when(
        Gene %in% tp_genes ~ tp_label,
        Gene %in% HPA_anno$Nuclear$Gene ~ "Nuclear",
        Gene %in% HPA_anno$Cytosol$Gene ~ "Cytosol",
        Gene %in% HPA_anno$Nuclear_Cytosol$Gene ~ "Nuclear_Cytosol",
        Gene %in% MitoCarta_anno$Symbol ~ "Mitochondrion",
        TRUE ~ "Other"
      ))
    
    # Counts
    anno_counts <- table(data_annotated[[col_name]])
    cat(sprintf("  Annotation counts:\n"))
    for (label in names(anno_counts)) {
      cat(sprintf("    %s: %d\n", label, anno_counts[label]))
    }
  }
  
  cat("\n✓ Data annotation complete\n")
  return(data_annotated)
}


#' Module 3 main function
#' 
#' @param dir_config Directory configuration
#' @param data_raw Raw data
#' @param sampleGroup Group info (optional, for downstream analysis)
#' @param custom_annotations Custom annotation configs (optional, overrides presets)
#' @param annotation_mode Preset annotation mode ("SG" or "Nucleolus")
#' @param custom_tp_sources Custom TP reference configs (optional)
#' @param additional_annotations Extra annotation configs (appended to preset/custom)
module03_annotation <- function(dir_config,
                                data_raw,
                                sampleGroup = NULL,
                                custom_annotations = NULL,
                                annotation_mode = c("SG", "Nucleolus"),
                                custom_tp_sources = NULL,
                                additional_annotations = NULL) {
  
  cat("\n========================================\n")
  cat("Module 3: Annotation System\n")
  cat("========================================\n")
  annotation_mode <- match.arg(annotation_mode)
  cat(sprintf("Current annotation mode: %s\n", annotation_mode))
  
  # --------------------------------------------------------------------------
  # Load reference data
  # --------------------------------------------------------------------------
  HPA_anno <- load_HPA_annotations(dir_config, annotation_mode)
  MitoCarta_anno <- load_MitoCarta_annotations(dir_config)
  SGs_refs <- load_SGs_references(dir_config)
  custom_tp_refs <- load_custom_tp_sources(custom_tp_sources, dir_config)
  
  # --------------------------------------------------------------------------
  # Annotation mode and configs
  # --------------------------------------------------------------------------
  use_custom <- !is.null(custom_annotations) && length(custom_annotations) > 0
  if (use_custom) {
    annotation_configs <- custom_annotations
    cat(sprintf("✓ Using custom annotation configs (%d columns)\n", length(annotation_configs)))
  } else {
    annotation_configs <- get_preset_annotation_configs(annotation_mode)
    if (length(annotation_configs) == 0) {
      cat("⚠ No TP annotations configured; only basic annotations (Cytosol/Nuclear/Mitochondrion) will be used\n")
    } else {
      cat(sprintf("✓ Loaded preset annotation configs (mode: %s, columns: %d)\n",
                  annotation_mode, length(annotation_configs)))
    }
  }
  
  has_additional <- !is.null(additional_annotations) && length(additional_annotations) > 0
  if (has_additional) {
    cat(sprintf("✓ Appended custom annotation columns (%d columns)\n", length(additional_annotations)))
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
  # User custom annotation config (all TP columns must be configured manually)
  # --------------------------------------------------------------------------
  if (length(annotation_configs) == 0) {
    cat("  - No annotation configs provided; output will keep original columns\n")
  }

  # --------------------------------------------------------------------------
  # Run annotation
  # --------------------------------------------------------------------------
  data_annotated <- annotate_data(
    data_raw,
    HPA_anno,
    MitoCarta_anno,
    annotation_configs,
    tp_lookup
  )
  
  # --------------------------------------------------------------------------
  # Save data
  # --------------------------------------------------------------------------
  cat("\n----------------------------------------\n")
  cat("Step 5: Save data\n")
  cat("----------------------------------------\n")
  
  setwd(dir_config$root)
  
  # Save reference data for later modules
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
  
  # Save CSV to Output directory
  output_file <- file.path(dir_config$output, "Module03_data_annotated.csv")
  write.csv(data_annotated, output_file, row.names = FALSE)
  cat(sprintf("✓ Saved: %s (Output directory)\n", output_file))
  
  cat("\n✓ Module 3 processing complete; returning data to main pipeline\n")
  
  return(list(
    data_annotated = data_annotated,
    annotation_references = annotation_references
  ))
}

