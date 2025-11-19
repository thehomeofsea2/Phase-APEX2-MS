# ============================================================================
# Module 1: 环境初始化
# ============================================================================
# 功能：
#   1. 检查必需目录（Reference, Rawdata）
#   2. 自动创建可选目录（Output, Module, Dev）
#   3. 保存配置信息
# 
# 输入：dir_config（目录配置列表）
# 输出：config.RData
# ============================================================================

module01_setup <- function(dir_config) {
  
  cat("开始环境初始化...\n")
  
  # --------------------------------------------------------------------------
  # 1. 检查必需目录
  # --------------------------------------------------------------------------
  required_dirs <- c("reference", "rawdata")
  
  for (dir_name in required_dirs) {
    dir_path <- dir_config[[dir_name]]
    if (!dir.exists(dir_path)) {
      stop(sprintf("✗ 错误：必需目录不存在 - %s\n请创建该目录后重试！", dir_path))
    }
    cat(sprintf("✓ 检查必需目录: %s\n", dir_path))
  }
  
  # --------------------------------------------------------------------------
  # 2. 自动创建可选目录
  # --------------------------------------------------------------------------
  optional_dirs <- c("output", "module", "dev")
  
  for (dir_name in optional_dirs) {
    dir_path <- dir_config[[dir_name]]
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat(sprintf("✓ 创建目录: %s\n", dir_path))
    } else {
      cat(sprintf("✓ 目录已存在: %s\n", dir_path))
    }
  }
  
  # --------------------------------------------------------------------------
  # 3. 保存所有环境变量到工作目录
  # --------------------------------------------------------------------------
  setwd(dir_config$root)
  
  config <- list(
    dir_config = dir_config,
    timestamp = Sys.time(),
    R_version = R.version.string,
    platform = .Platform$OS.type
  )
  
  # --------------------------------------------------------------------------
  # 4. 打印摘要信息
  # --------------------------------------------------------------------------
  cat("\n========================================\n")
  cat("环境初始化完成\n")
  cat("========================================\n")
  cat(sprintf("工作目录: %s\n", dir_config$root))
  cat(sprintf("时间戳: %s\n", config$timestamp))
  cat(sprintf("R版本: %s\n", config$R_version))
  cat("========================================\n\n")
  
  return(invisible(config))
}

