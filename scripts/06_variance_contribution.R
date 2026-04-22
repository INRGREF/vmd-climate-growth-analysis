# =========================================================
# 06 - VARIANCE CONTRIBUTION (FINAL VERSION)
# =========================================================

library(dplyr)

# =========================================================
# PATHS
# =========================================================

growth_dir <- "results/VMD_growth"
cycle_dir  <- file.path(growth_dir, "Cycle_analysis")

out_dir <- file.path(growth_dir, "variance_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 1. MODE VARIANCE
# =========================================================

files <- list.files(cycle_dir, pattern = "Growth_.*\\.txt$", full.names = TRUE)

var_list <- list()

for (f in files) {

  df <- read.table(f, header = TRUE)
  df <- df[, -1]

  var_list[[f]] <- apply(df, 2, var, na.rm = TRUE)
}

mat <- do.call(rbind, var_list)

mean_var <- colMeans(mat, na.rm = TRUE)
contrib <- mean_var / sum(mean_var) * 100

res_modes <- data.frame(
  Mode = names(mean_var),
  Variance = mean_var,
  Contribution = contrib
)

write.csv(res_modes,
          file.path(out_dir, "variance_modes.csv"),
          row.names = FALSE)

# =========================================================
# 2. CYCLE VARIANCE
# =========================================================

files_c <- list.files(cycle_dir, pattern = "Cycle_.*\\.txt$", full.names = TRUE)

cycle_var <- sapply(files_c, function(f) {
  df <- read.table(f, header = TRUE)[,-1]
  mean(apply(df, 2, var, na.rm = TRUE))
})

res_cycles <- data.frame(
  Cycle = basename(files_c),
  Variance = cycle_var
)

res_cycles$Contribution <- res_cycles$Variance / sum(res_cycles$Variance) * 100

write.csv(res_cycles,
          file.path(out_dir, "variance_cycles.csv"),
          row.names = FALSE)

cat("\n✔ Variance analysis done\n")
