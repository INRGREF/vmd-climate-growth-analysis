# =========================================================
# 04 - CLIMATE–GROWTH REGRESSION (PCR Bootstrap)
# Multi-scale VMD-based climate–growth relationships
# Pinus halepensis (Tunisia)
# =========================================================

# ---------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(FactoMineR)

# ---------------------------------------------------------
# 2. Global paths (HARMONIZED WITH VMD WORKFLOW)
# ---------------------------------------------------------

# Growth data (from VMD growth decomposition)
growth_dir <- "results/VMD_growth/CHRONOLOGIES_GROWTH_CYCLES"

# Climate PCA results (created in Script 3)
temp_dir   <- "results/VMD_temperature_modes/Cycle_analysis"
precip_dir <- "results/VMD_precipitation_modes/Cycle_analysis"

# ---------------------------------------------------------
# 3. PCR Bootstrap (multi-variable climate cycles)
# ---------------------------------------------------------

growth_files <- list.files(
  growth_dir,
  pattern = "Growth_Chronology_.*\\.csv$",
  full.names = TRUE
)

final_results_list <- list()
n_boot <- 1000

for (f in growth_files) {
  
  # -----------------------------
  # Cycle identification
  # -----------------------------
  cycle_raw <- gsub("Growth_Chronology_|\\.csv", "", basename(f))
  cycle_name <- ifelse(cycle_raw == "Trend", "Cycle_Trend", cycle_raw)
  
  # Check PCA availability
  if (!(cycle_name %in% names(pca_results))) next
  
  cat("\nPCR regression for:", cycle_name, "\n")
  
  # -----------------------------
  # Load growth + PCA scores
  # -----------------------------
  df_growth <- read.csv(f) %>% select(Year, std)
  df_scores <- pca_results[[cycle_name]]$scores
  loadings  <- pca_results[[cycle_name]]$loadings
  
  df_pcr <- inner_join(df_growth, df_scores, by = "Year") %>% drop_na()
  if (nrow(df_pcr) < 20) next
  
  # -----------------------------
  # Regression setup
  # -----------------------------
  vars_cp <- setdiff(colnames(df_scores), "Year")
  formula <- as.formula(paste("std ~", paste(vars_cp, collapse = " + ")))
  
  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))
  
  set.seed(123)
  for (i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(formula, data = df_b)
    boot_coefs_cp[i, ] <- coef(mod)[-1]
  }
  
  # -----------------------------
  # Back-transformation (CP → monthly climate)
  # -----------------------------
  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)
  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  
  boot_monthly <- boot_coefs_cp %*% t(loadings)
  
  inf95 <- apply(boot_monthly, 2, quantile, 0.025, na.rm = TRUE)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975, na.rm = TRUE)
  
  # -----------------------------
  # Output formatting
  # -----------------------------
  res_df <- data.frame(
    Cycle = cycle_name,
    Term  = rownames(loadings),
    Coef  = as.vector(coefs_monthly),
    Inf95 = inf95,
    Sup95 = sup95
  )
  
  final_results_list[[cycle_name]] <- res_df
}

# ---------------------------------------------------------
# 4. Single-variable cycles (T or P only)
# ---------------------------------------------------------

results_solo <- list()

cycles_solo <- list(
  "Cycle_10" = list(dir = temp_dir, var = "T"),
  "Cycle_12" = list(dir = precip_dir, var = "P"),
  "Cycle_13" = list(dir = precip_dir, var = "P"),
  "Cycle_15" = list(dir = precip_dir, var = "P")
)

for (cyc in names(cycles_solo)) {
  
  var_code <- cycles_solo[[cyc]]$var
  
  suffix <- ifelse(
    var_code == "T",
    "_Temperature.transf.txt",
    "_Precipitation.transf.txt"
  )
  
  f_clim <- file.path(cycles_solo[[cyc]]$dir, paste0(cyc, suffix))
  
  if (!file.exists(f_clim)) {
    cat("Missing file:", f_clim, "\n")
    next
  }
  
  cat("\nSolo PCR for:", cyc, "\n")
  
  data_clim <- read.table(f_clim, header = TRUE, check.names = FALSE)
  data_input <- data_clim %>% select(-Year)
  data_input[is.na(data_input)] <- 0
  
  res.pca <- PCA(data_input, scale.unit = TRUE, ncp = 12, graph = FALSE)
  
  eigenvalues <- res.pca$eig[, 1]
  nb_cp_keep <- sum(eigenvalues > mean(eigenvalues))
  if (nb_cp_keep < 1) nb_cp_keep <- 1
  
  scores <- as.data.frame(res.pca$ind$coord[, 1:nb_cp_keep, drop = FALSE])
  colnames(scores) <- paste0("Dim", 1:nb_cp_keep)
  scores$Year <- data_clim$Year
  
  loadings <- res.pca$var$coord[, 1:nb_cp_keep, drop = FALSE]
  
  # growth file
  f_growth <- file.path(growth_dir, paste0("Growth_Chronology_", cyc, ".csv"))
  if (!file.exists(f_growth)) next
  
  df_growth <- read.csv(f_growth) %>% select(Year, std)
  df_pcr <- inner_join(df_growth, scores, by = "Year") %>% drop_na()
  if (nrow(df_pcr) < 20) next
  
  vars_cp <- setdiff(colnames(scores), "Year")
  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))
  
  set.seed(123)
  for (i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(std ~ ., data = df_b %>% select(all_of(vars_cp), std))
    boot_coefs_cp[i, ] <- coef(mod)[-1]
  }
  
  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)
  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  
  boot_monthly <- boot_coefs_cp %*% t(loadings)
  
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

# ---------------------------------------------------------
# 5. Final merge
# ---------------------------------------------------------

all_results <- bind_rows(c(final_results_list, results_solo))

final_df <- all_results %>%
  mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", "")) %>%
  separate(Term, into = c("Variable", "Month"), sep = "_", fill = "right")

cat("\n✔ PCR analysis completed successfully.\n")
