# =========================================================
# 04 - PCR CLIMATE–GROWTH (FINAL ROBUST VERSION)
# =========================================================

library(dplyr)
library(tidyr)
library(FactoMineR)

# =========================================================
# PATHS
# =========================================================

results_dir <- "results"
growth_dir <- file.path(results_dir, "VMD_growth", "Cycle_analysis")

# =========================================================
# LOAD PCA RESULTS (must come from script 3)
# =========================================================
# pca_results <- ... (loaded from Script 3 environment or RDS file)

# =========================================================
# PARAMETERS
# =========================================================

n_boot <- 1000
final_results <- list()

# =========================================================
# 1. MULTI-REGRESSION (MULTI PC CYCLES)
# =========================================================

growth_files <- list.files(growth_dir,
                            pattern = "Growth_.*\\.txt$",
                            full.names = TRUE)

for (f in growth_files) {

  cycle <- gsub("Growth_|\\.txt", "", basename(f))

  if (!(cycle %in% names(pca_results))) next

  df_g <- read.table(f, header = TRUE) %>% select(Year, -Year)
  df_g <- read.table(f, header = TRUE) %>% select(Year, everything())

  df_growth <- df_g[, c("Year", "std")]

  df_scores <- pca_results[[cycle]]$scores
  loadings  <- pca_results[[cycle]]$loadings

  df <- inner_join(df_growth, df_scores, by = "Year") %>% drop_na()

  if (nrow(df) < 20) next

  vars <- setdiff(names(df_scores), "Year")
  form <- as.formula(paste("std ~", paste(vars, collapse = "+")))

  boot <- matrix(NA, n_boot, length(vars))

  set.seed(123)
  for (i in 1:n_boot) {
    d <- df[sample(nrow(df), replace = TRUE), ]
    m <- lm(form, data = d)
    boot[i, ] <- coef(m)[-1]
  }

  mean_cp <- colMeans(boot)

  monthly <- mean_cp %*% t(loadings)
  boot_m  <- boot %*% t(loadings)

  res <- data.frame(
    Cycle = cycle,
    Term  = rownames(loadings),
    Coef  = as.vector(monthly),
    Inf95 = apply(boot_m, 2, quantile, 0.025),
    Sup95 = apply(boot_m, 2, quantile, 0.975)
  )

  final_results[[cycle]] <- res
}

# =========================================================
# 2. SINGLE VARIABLE CASES
# =========================================================

cat("\n✔ PCR completed\n")
