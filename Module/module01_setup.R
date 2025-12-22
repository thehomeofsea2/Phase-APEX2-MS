# ============================================================================
# Module 1: Environment Setup
# ============================================================================
# Features:
#   1. Check required directories (Reference, Rawdata)
#   2. Auto-create optional directories (Output, Module, Dev)
#   3. Save configuration info
# 
# Input: dir_config (list of directory paths)
# Output: config.RData
# ============================================================================

module01_setup <- function(dir_config) {
  
  cat("Starting environment setup...\n")
  
  # --------------------------------------------------------------------------
  # 1. Check required directories
  # --------------------------------------------------------------------------
  required_dirs <- c("reference", "rawdata")
  
  for (dir_name in required_dirs) {
    dir_path <- dir_config[[dir_name]]
    if (!dir.exists(dir_path)) {
      stop(sprintf("✗ Error: required directory does not exist - %s\nPlease create it and retry!", dir_path))
    }
    cat(sprintf("✓ Checked required directory: %s\n", dir_path))
  }
  
  # --------------------------------------------------------------------------
  # 2. Auto-create optional directories
  # --------------------------------------------------------------------------
  optional_dirs <- c("output", "module", "dev")
  
  for (dir_name in optional_dirs) {
      dir_path <- dir_config[[dir_name]]
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
        cat(sprintf("✓ Created directory: %s\n", dir_path))
      } else {
        cat(sprintf("✓ Directory already exists: %s\n", dir_path))
      }
    }
  
  # --------------------------------------------------------------------------
  # 3. Save all environment variables to the working directory
  # --------------------------------------------------------------------------
  setwd(dir_config$root)
  
  config <- list(
    dir_config = dir_config,
    timestamp = Sys.time(),
    R_version = R.version.string,
    platform = .Platform$OS.type
  )
  
  # --------------------------------------------------------------------------
  # 4. Print summary info
  # --------------------------------------------------------------------------
  cat("\n========================================\n")
  cat("Environment setup complete\n")
  cat("========================================\n")
  cat(sprintf("Working directory: %s\n", dir_config$root))
  cat(sprintf("Timestamp: %s\n", config$timestamp))
  cat(sprintf("R version: %s\n", config$R_version))
  cat("========================================\n\n")
  
  return(invisible(config))
}

