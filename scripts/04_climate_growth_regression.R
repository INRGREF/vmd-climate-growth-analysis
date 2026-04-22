# =========================================================
# 04 - CLIMATE–GROWTH REGRESSION (PCR Bootstrap)
# Multi-scale VMD-based climate–growth relationships
# Pinus halepensis (Tunisia)
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(FactoMineR)

# =========================================================
# 1. PATHS (HARMONIZED WITH PIPELINE)
# =========================================================

growth_dir <- "results/VMD_growth"

temp_dir   <- "results/VMD_temperature_modes/Cycle_analysis"
precip_dir <- "results/VMD_precipitation_modes/Cycle_analysis"

# =========================================================
# IMPORTANT: pca_results must come from Script 3
# =========================================================
# source("03_PCA_analysis.R")

# =========================================================
# 2. PCR BOOTSTRAP (MULTI-VARIABLE CYCLES)
# =========================================================

growth_files <- list.files(
  growth_dir,
  pattern = "VMD_.*\\.txt$",
  full.names = TRUE
)

final_results_list <- list()
n_boot <- 1000

for (f in growth_files) {

  cycle_raw <- gsub("VMD_|\\.txt", "", basename(f))
  cycle_name <- paste0("Cycle_", cycle_raw)

  if (!(cycle_name %in% names(pca_results))) next

  cat("\nPCR regression for:", cycle_name, "\n")

  df_growth <- read.table(f, header = TRUE) %>%
    select(Year)

  # growth proxy = mean of IMFs (robust option)
  df_growth$std <- rowMeans(df_growth[,-1], na.rm = TRUE)

  df_scores <- pca_results[[cycle_name]]$scores
  loadings  <- pca_results[[cycle_name]]$loadings

  df_pcr <- inner_join(df_growth, df_scores, by = "Year") %>% drop_na()
  if (nrow(df_pcr) < 20) next

  vars_cp <- setdiff(colnames(df_scores), "Year")
  formula <- as.formula(paste("std ~", paste(vars_cp, collapse = " + ")))

  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))

  set.seed(123)
  for (i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(formula, data = df_b)
    boot_coefs_cp[i, ] <- coef(mod)[-1]
  }

  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)

  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  boot_monthly  <- boot_coefs_cp %*% t(loadings)

  inf95 <- apply(boot_monthly, 2, quantile, 0.025, na.rm = TRUE)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975, na.rm = TRUE)

  res_df <- data.frame(
    Cycle = cycle_name,
    Term  = rownames(loadings),
    Coef  = as.vector(coefs_monthly),
    Inf95 = inf95,
    Sup95 = sup95
  )

  final_results_list[[cycle_name]] <- res_df
}

# =========================================================
# 3. SINGLE VARIABLE CASES
# =========================================================

results_solo <- list()

cycles_solo <- list(
  "Cycle_10"  = list(dir = temp_dir, var = "T"),
  "Cycle_12"  = list(dir = precip_dir, var = "P"),
  "Cycle_13"  = list(dir = precip_dir, var = "P"),
  "Cycle_15"  = list(dir = precip_dir, var = "P")
)

for (cyc in names(cycles_solo)) {

  var_code <- cycles_solo[[cyc]]$var

  suffix <- ifelse(
    var_code == "T",
    "_Temperature.transf.txt",
    "_Precipitation.transf.txt"
  )

  f_clim <- file.path(
    cycles_solo[[cyc]]$dir,
    paste0(cyc, suffix)
  )

  if (!file.exists(f_clim)) next

  data_clim <- read.table(f_clim, header = TRUE)
  data_input <- data_clim %>% select(-Year)
  data_input[is.na(data_input)] <- 0

  res.pca <- PCA(data_input, scale.unit = TRUE, ncp = 12, graph = FALSE)

  eigenvalues <- res.pca$eig[,1]
  nb_cp_keep <- max(1, sum(eigenvalues > mean(eigenvalues)))

  scores <- as.data.frame(res.pca$ind$coord[,1:nb_cp_keep])
  colnames(scores) <- paste0("Dim", 1:nb_cp_keep)
  scores$Year <- data_clim$Year

  loadings <- res.pca$var$coord[,1:nb_cp_keep]

  f_growth <- file.path(growth_dir, paste0("VMD_", cyc, ".txt"))
  if (!file.exists(f_growth)) next

  df_growth <- read.table(f_growth, header = TRUE)
  df_growth$std <- rowMeans(df_growth[,-1], na.rm = TRUE)

  df_pcr <- inner_join(df_growth, scores, by = "Year") %>% drop_na()
  if (nrow(df_pcr) < 20) next

  vars_cp <- setdiff(colnames(scores), "Year")

  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))

  set.seed(123)
  for (i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(std ~ ., data = df_b)
    boot_coefs_cp[i, ] <- coef(mod)[-1]
  }

  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)

  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  boot_monthly  <- boot_coefs_cp %*% t(loadings)

  inf95 <- apply(boot_monthly, 2, quantile, 0.025)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975)

  month_names <- paste0(var_code, "_", rownames(loadings))

  res_df <- data.frame(
    Cycle = cyc,
    Term  = month_names,
    Coef  = as.vector(coefs_monthly),
    Inf95 = inf95,
    Sup95 = sup95
  )

  results_solo[[cyc]] <- res_df
}

# =========================================================
# 4. FINAL MERGE
# =========================================================

all_results <- bind_rows(c(final_results_list, results_solo))

final_df <- all_results %>%
  mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", "")) %>%
  tidyr::separate(Term, into = c("Variable", "Month"),
                   sep = "_", fill = "right")

cat("\n✔ PCR analysis completed\n")
