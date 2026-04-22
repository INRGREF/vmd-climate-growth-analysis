# =========================================================
# 02 - FULL VMD DECOMPOSITION PIPELINE (FINAL)
# Climate + Growth (Pinus halepensis, Tunisia)
# =========================================================

library(VMDecomp)
library(dplyr)
library(dplR)

# =========================================================
# 1. PROJECT STRUCTURE
# =========================================================

results_dir <- "results"

climate_temp_dir <- file.path(results_dir, "VMD_climate_temperature")
climate_prec_dir <- file.path(results_dir, "VMD_climate_precipitation")
growth_dir       <- file.path(results_dir, "VMD_growth")

cycle_temp_dir <- file.path(climate_temp_dir, "Cycle_analysis")
cycle_prec_dir <- file.path(climate_prec_dir, "Cycle_analysis")
cycle_growth_dir <- file.path(growth_dir, "Cycle_analysis")

dir.create(cycle_temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cycle_prec_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cycle_growth_dir, recursive = TRUE, showWarnings = FALSE)

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
# 3. VMD FUNCTION (CLIMATE + GROWTH)
# =========================================================

run_vmd <- function(file_path, output_dir, label, cycle_dir) {

  data <- read.table(file_path, header = TRUE)

  # =========================================================
  # CLIMATE FORMAT (12 months)
  # =========================================================
  monthly <- data[, 2:13]

  colnames(monthly) <- c(
    "Jan","Feb","Mar","Apr","May","Jun",
    "Jul","Aug","Sep","Oct","Nov","Dec"
  )

  vmd_results <- list()

  for (m in colnames(monthly)) {

    series <- monthly[[m]]
    k_opt <- find_k_optimal(series)

    res <- vmd(series,
               alpha = 2000,
               tau = 0,
               K = k_opt,
               DC = FALSE,
               init = 1,
               tol = 1e-7)

    vmd_results[[m]] <- list(
      modes = res$u,
      frequencies = res$omega[nrow(res$omega), ],
      K = k_opt
    )
  }

  # =========================================================
  # EXPORT PERIODS (FOR PCA + VARIANCE)
  # =========================================================

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
    file = file.path(cycle_dir,
                     paste0(label, "_vmd_periods.txt")),
    sep = "\t",
    quote = FALSE
  )

  # =========================================================
  # EXPORT MODES
  # =========================================================

  modes_dir <- file.path(output_dir, "modes")
  dir.create(modes_dir, recursive = TRUE, showWarnings = FALSE)

  for (m in names(vmd_results)) {

    modes <- vmd_results[[m]]$modes
    if (nrow(modes) != nrow(data)) modes <- t(modes)

    df <- cbind(Year = data[,1], as.data.frame(modes))
    colnames(df) <- c("Year", paste0("IMF_", 1:ncol(modes)))

    write.csv(df,
              file.path(modes_dir,
                        paste0("VMD_", m, ".csv")),
              row.names = FALSE)
  }

  # =========================================================
  # EXPORT CYCLES (CRITICAL FOR SCRIPT 3)
  # =========================================================

  dir.create(cycle_dir, recursive = TRUE, showWarnings = FALSE)

  for (m in names(vmd_results)) {

    modes <- vmd_results[[m]]$modes
    K <- ncol(modes)

    for (k in 1:K) {

      df_cycle <- data.frame(
        Year = data[,1],
        Cycle = modes[,k]
      )

      write.table(
        df_cycle,
        file = file.path(cycle_dir,
                         paste0("Cycle_", k, "_", label, "_", m, ".txt")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }

  return(vmd_results)
}

# =========================================================
# 4. GROWTH VMD FUNCTION (TREE-RINGS SPECIAL CASE)
# =========================================================

run_vmd_growth <- function(tree_data, output_dir, label, cycle_dir) {

  vmd_tree_results <- list()

  for (i in 2:ncol(tree_data)) {

    id <- colnames(tree_data)[i]
    series <- tree_data[, i]
    years <- tree_data[, 1]

    idx <- which(!is.na(series))
    clean <- series[idx]

    if (length(clean) > 30) {

      k_opt <- find_k_optimal(clean)

      res <- vmd(clean,
                 alpha = 2000,
                 tau = 0,
                 K = k_opt,
                 DC = FALSE,
                 init = 1,
                 tol = 1e-7)

      aligned <- matrix(NA, length(years), k_opt)
      modes <- res$u

      if (nrow(modes) != length(clean)) modes <- t(modes)

      aligned[idx, ] <- modes

      df <- cbind(Year = years, as.data.frame(aligned))
      colnames(df) <- c("Year", paste0("IMF_", 1:k_opt))

      write.table(df,
                  file.path(output_dir,
                            paste0("VMD_", id, ".txt")),
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)

      # store cycles
      for (k in 1:k_opt) {

        df_cycle <- data.frame(
          Year = years,
          Cycle = aligned[,k]
        )

        write.table(
          df_cycle,
          file = file.path(cycle_dir,
                           paste0("Cycle_", k, "_", label, "_", id, ".txt")),
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
      }

      vmd_tree_results[[id]] <- list(
        frequencies = res$omega[nrow(res$omega), ],
        K = k_opt
      )
    }
  }

  return(vmd_tree_results)
}

# =========================================================
# 5. RUN CLIMATE + GROWTH
# =========================================================

temp_file <- "raw_data/temperature.txt"
prec_file <- "raw_data/precipitation.txt"

run_vmd(temp_file, climate_temp_dir, "temperature", cycle_temp_dir)
run_vmd(prec_file, climate_prec_dir, "precipitation", cycle_prec_dir)

tree_data <- rw_index

run_vmd_growth(tree_data, growth_dir, "growth", cycle_growth_dir)

cat("\n✔ FULL VMD PIPELINE COMPLETED (CLIMATE + GROWTH)\n")
