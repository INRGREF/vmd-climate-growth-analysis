# =========================================================
# 02 - FULL VMD DECOMPOSITION PIPELINE
# Climate (Temperature + Precipitation)
# Growth (Tree-ring chronologies)
# Pinus halepensis study
# =========================================================

library(VMDecomp)
library(dplyr)
library(dplR)

# =========================================================
# 1. PROJECT STRUCTURE (GLOBAL STANDARD)
# =========================================================

results_dir <- "results"

climate_temp_dir <- file.path(results_dir, "VMD_climate_temperature")
climate_prec_dir <- file.path(results_dir, "VMD_climate_precipitation")
growth_dir       <- file.path(results_dir, "VMD_growth")
cycle_temp_dir <- file.path(climate_temp_dir, "Cycle_analysis")
cycle_prec_dir <- file.path(climate_prec_dir, "Cycle_analysis")

if(!dir.exists(cycle_temp_dir)) dir.create(cycle_temp_dir, recursive = TRUE)
if(!dir.exists(cycle_prec_dir)) dir.create(cycle_prec_dir, recursive = TRUE)
if(!dir.exists(climate_temp_dir)) dir.create(climate_temp_dir, recursive = TRUE)
if(!dir.exists(climate_prec_dir)) dir.create(climate_prec_dir, recursive = TRUE)
if(!dir.exists(growth_dir)) dir.create(growth_dir, recursive = TRUE)

# =========================================================
# 2. OPTIMAL K FUNCTION
# =========================================================

find_k_optimal <- function(series, max_k = 8, threshold = 0.06) {
  for (k in 2:max_k) {
    res <- vmd(series, alpha = 2000, tau = 0, K = k,
               DC = FALSE, init = 1, tol = 1e-7)
    omegas <- sort(res$omega[nrow(res$omega), ])
    if (min(diff(omegas)) < threshold) return(k - 1)
  }
  return(max_k)
}

# =========================================================
# 3. CLIMATE VMD FUNCTION
# =========================================================

run_vmd_climate <- function(file_path, output_dir, label) {
  
  data <- read.table(file_path, header = TRUE)
  monthly <- data[, 2:13]
  
  colnames(monthly) <- c(
    "Jan","Feb","Mar","Apr","May","Jun",
    "Jul","Aug","Sep","Oct","Nov","Dec"
  )
  
  vmd_results <- list()
  
  for (m in colnames(monthly)) {
    
    series <- monthly[[m]]
    k_opt <- find_k_optimal(series)
    
    res <- vmd(series, alpha = 2000, tau = 0,
               K = k_opt, DC = FALSE, init = 1, tol = 1e-7)
    
    vmd_results[[m]] <- list(
      modes = res$u,
      frequencies = res$omega[nrow(res$omega), ],
      K = k_opt
    )
  }
  
  # ---- periods ----
  maxK <- max(sapply(vmd_results, function(x) x$K))
  
  period_mat <- matrix(NA, 12, maxK)
  rownames(period_mat) <- names(vmd_results)
  colnames(period_mat) <- paste0("IMF_", 1:maxK)
  
  for (m in names(vmd_results)) {
    p <- sort(1 / vmd_results[[m]]$frequencies, decreasing = TRUE)
    period_mat[m, 1:length(p)] <- p
  }
  
  write.table(
    round(period_mat, 2),
    file = file.path(output_dir, paste0(label, "_vmd_periods_rounded.txt")),
    sep = "\t", quote = FALSE
  )
  
  # ---- export modes ----
  modes_dir <- file.path(output_dir, "modes")
  if(!dir.exists(modes_dir)) dir.create(modes_dir, recursive = TRUE)
  
  for (m in names(vmd_results)) {
    
    modes <- vmd_results[[m]]$modes
    if (nrow(modes) != nrow(data)) modes <- t(modes)
    
    df <- cbind(Year = data[,1], as.data.frame(modes))
    colnames(df) <- c("Year", paste0("IMF_", 1:ncol(modes)))
    
    write.csv(df,
              file.path(modes_dir, paste0("VMD_", m, ".csv")),
              row.names = FALSE)
  }
  
  return(vmd_results)
}

# =========================================================
# 4. CLIMATE DATA PROCESSING
# =========================================================

temp_file <- "raw_data/temperature.txt"
prec_file <- "raw_data/precipitation.txt"

temp_results <- run_vmd_climate(temp_file, climate_temp_dir, "temperature")
prec_results <- run_vmd_climate(prec_file, climate_prec_dir, "precipitation")

cat("\n✔ Climate VMD completed\n")

# =========================================================
# 5. GROWTH VMD (TREE-RINGS)
# =========================================================

tree_data <- rw_index

find_k_growth <- function(series, max_k = 7, threshold = 0.04) {
  for (k in 2:max_k) {
    res <- vmd(series, alpha = 2000, tau = 0,
               K = k, DC = FALSE, init = 1, tol = 1e-7)
    if (min(diff(sort(res$omega[nrow(res$omega), ]))) < threshold)
      return(k - 1)
  }
  return(max_k)
}

vmd_tree_results <- list()

for (i in 2:ncol(tree_data)) {
  
  id <- colnames(tree_data)[i]
  series <- tree_data[, i]
  years <- tree_data[, 1]
  
  idx <- which(!is.na(series))
  clean <- series[idx]
  
  if (length(clean) > 30) {
    
    k_opt <- find_k_growth(clean)
    
    res <- vmd(clean, alpha = 2000, tau = 0,
               K = k_opt, DC = FALSE, init = 1, tol = 1e-7)
    
    aligned <- matrix(NA, length(years), k_opt)
    modes <- res$u
    if (nrow(modes) != length(clean)) modes <- t(modes)
    
    aligned[idx, ] <- modes
    
    df <- cbind(Year = years, as.data.frame(aligned))
    colnames(df) <- c("Year", paste0("IMF_", 1:k_opt))
    
    write.table(df,
                file.path(growth_dir, paste0("VMD_", id, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    vmd_tree_results[[id]] <- list(
      frequencies = res$omega[nrow(res$omega), ],
      K = k_opt
    )
  }
}

# ---- growth periods ----
maxK <- max(sapply(vmd_tree_results, function(x) x$K))

growth_periods <- matrix(NA, length(vmd_tree_results), maxK)
rownames(growth_periods) <- names(vmd_tree_results)
colnames(growth_periods) <- paste0("IMF_", 1:maxK)

for (t in names(vmd_tree_results)) {
  p <- sort(1 / vmd_tree_results[[t]]$frequencies, decreasing = TRUE)
  growth_periods[t, 1:length(p)] <- p
}

write.table(
  round(growth_periods, 2),
  file = file.path(growth_dir, "tree_vmd_periods_rounded.txt"),
  sep = "\t", quote = FALSE, col.names = NA
)

cat("\n✔ Growth VMD completed\n")
