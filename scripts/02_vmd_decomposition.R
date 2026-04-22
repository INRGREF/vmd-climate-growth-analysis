# =========================================================
# 02 - VMD FULL PIPELINE (PUBLICATION VERSION)
# Climate + Growth
# =========================================================

library(VMDecomp)
library(dplyr)
library(dplR)

# =========================================================
# 1. PATHS
# =========================================================

results_dir <- "results"

clim_t_dir <- file.path(results_dir, "VMD_climate_temperature")
clim_p_dir <- file.path(results_dir, "VMD_climate_precipitation")
growth_dir <- file.path(results_dir, "VMD_growth")

# STANDARD CYCLE OUTPUT (IMPORTANT)
cycle_t_dir <- file.path(clim_t_dir, "Cycle_analysis")
cycle_p_dir <- file.path(clim_p_dir, "Cycle_analysis")
cycle_g_dir <- file.path(growth_dir, "Cycle_analysis")

dir.create(cycle_t_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cycle_p_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cycle_g_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. OPTIMAL K
# =========================================================

find_k <- function(x, max_k = 8, thr = 0.06) {
  for (k in 2:max_k) {
    res <- vmd(x, alpha = 2000, tau = 0, K = k,
               DC = FALSE, init = 1, tol = 1e-7)
    om <- sort(res$omega[nrow(res$omega), ])
    if (min(diff(om)) < thr) return(k - 1)
  }
  max_k
}

# =========================================================
# 3. GENERIC VMD FUNCTION
# =========================================================

run_vmd <- function(file, out_dir, label) {

  df <- read.table(file, header = TRUE)
  clim <- df[, 2:13]

  colnames(clim) <- c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec")

  K_store <- list()

  for (m in names(clim)) {

    x <- clim[[m]]
    k <- find_k(x)

    res <- vmd(x, alpha = 2000, tau = 0, K = k,
               DC = FALSE, init = 1, tol = 1e-7)

    # save cycles per month
    df_out <- data.frame(
      Year = df[,1],
      res$u
    )

    colnames(df_out) <- c("Year", paste0("IMF_", 1:ncol(res$u)))

    write.table(df_out,
                file.path(out_dir, paste0(label, "_", m, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # cycle summary
    K_store[[m]] <- list(K = k,
                         freq = res$omega[nrow(res$omega), ])
  }

  # cycles table
  maxK <- max(sapply(K_store, `[[`, "K"))

  cyc_mat <- matrix(NA, 12, maxK)
  rownames(cyc_mat) <- names(K_store)
  colnames(cyc_mat) <- paste0("IMF_", 1:maxK)

  for (m in names(K_store)) {
    cyc_mat[m, 1:length(K_store[[m]]$freq)] <-
      sort(1 / K_store[[m]]$freq, decreasing = TRUE)
  }

  write.table(round(cyc_mat, 2),
              file.path(out_dir, paste0(label, "_Cycle_analysis.txt")),
              sep = "\t", quote = FALSE)

  return(K_store)
}

# =========================================================
# 4. RUN VMD
# =========================================================

temp <- run_vmd("raw_data/temperature.txt", cycle_t_dir, "Temperature")
prec <- run_vmd("raw_data/precipitation.txt", cycle_p_dir, "Precipitation")

# growth
tree <- rw_index

for (i in 2:ncol(tree)) {

  id <- colnames(tree)[i]
  x <- tree[, i]
  yrs <- tree[,1]

  ok <- which(!is.na(x))
  x_clean <- x[ok]

  if (length(x_clean) < 30) next

  k <- find_k(x_clean)

  res <- vmd(x_clean, alpha = 2000, tau = 0, K = k,
             DC = FALSE, init = 1, tol = 1e-7)

  mat <- matrix(NA, nrow(tree), k)
  mat[ok, ] <- t(res$u)

  df <- cbind(Year = yrs, mat)
  colnames(df) <- c("Year", paste0("IMF_", 1:k))

  write.table(df,
              file.path(cycle_g_dir, paste0("Growth_", id, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\n✔ VMD completed\n")
