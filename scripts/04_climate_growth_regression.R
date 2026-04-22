# =========================================================
# 04 - CLIMATE–GROWTH REGRESSION (PCR Bootstrap - FINAL)
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
# 2. PATHS
# ---------------------------------------------------------

growth_dir <- "results/VMD_growth/CHRONOLOGIES_GROWTH_CYCLES"

# IMPORTANT: must match Script 3
temp_dir   <- "results/VMD_climate_temperature/Cycle_analysis"
precip_dir <- "results/VMD_climate_precipitation/Cycle_analysis"

# =========================================================
# IMPORTANT: pca_results must exist from Script 3
# =========================================================

if (!exists("pca_results")) {
  stop("pca_results not found. Run Script 3 first.")
}

# ---------------------------------------------------------
# 3. MULTI-VARIABLE PCR
# ---------------------------------------------------------

growth_files <- list.files(
  growth_dir,
  pattern = "Chrono_.*\\.csv$",
  full.names = TRUE
)

final_results_list <- list()
n_boot <- 1000

for (f in growth_files) {
  
  # -----------------------------
  # Cycle identification (ROBUST)
  # -----------------------------
  cycle_raw <- gsub("Chrono_|\\.csv", "", basename(f))
  cycle_name <- ifelse(cycle_raw == "Trend",
                        "Cycle_Trend",
                        cycle_raw)
  
  # Skip if PCA missing
  if (!(cycle_name %in% names(pca_results))) {
    cat("Skipping missing PCA cycle:", cycle_name, "\n")
    next
  }
  
  cat("\nPCR regression for:", cycle_name, "\n")
  
  # -----------------------------
  # Load data
  # -----------------------------
  df_growth <- read.csv(f) %>% select(Year, std)
  
  df_scores <- pca_results[[cycle_name]]$scores
  loadings  <- pca_results[[cycle_name]]$loadings
  
  # safety check
  if (ncol(df_scores) < 2 || nrow(df_scores) < 10) next
  
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
  # Back-transform CP → monthly climate
  # -----------------------------
  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)
  
  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  boot_monthly  <- boot_coefs_cp %*% t(loadings)
  
  inf95 <- apply(boot_monthly, 2, quantile, 0.025, na.rm = TRUE)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975, na.rm = TRUE)
  
  # -----------------------------
  # Output
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
# 4. SINGLE VARIABLE CYCLES
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
  
  f_clim <- file.path(cycles_solo[[cyc]]$dir,
                      paste0(cyc, suffix))
  
  if (!file.exists(f_clim)) {
    cat("Missing file:", f_clim, "\n")
    next
  }
  
  cat("\nSolo PCR:", cyc, "\n")
  
  data_clim <- read.table(f_clim, header = TRUE, check.names = FALSE)
  data_input <- data_clim %>% select(-Year)
  data_input[is.na(data_input)] <- 0
  
  res.pca <- PCA(data_input,
                 scale.unit = TRUE,
                 ncp = min(12, ncol(data_input)),
                 graph = FALSE)
  
  eigenvalues <- res.pca$eig[,1]
  nb_cp <- max(1, sum(eigenvalues > mean(eigenvalues)))
  
  scores <- as.data.frame(res.pca$ind$coord[, 1:nb_cp, drop = FALSE])
  colnames(scores) <- paste0("Dim", 1:nb_cp)
  scores$Year <- data_clim$Year
  
  loadings <- res.pca$var$coord[, 1:nb_cp, drop = FALSE]
  
  f_growth <- file.path(growth_dir,
                        paste0("Growth_Chronology_", cyc, ".csv"))
  
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
  boot_monthly  <- boot_coefs_cp %*% t(loadings)
  
  inf95 <- apply(boot_monthly, 2, quantile, 0.025, na.rm = TRUE)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975, na.rm = TRUE)
  
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
# 5. FINAL MERGE
# ---------------------------------------------------------

all_results <- bind_rows(c(final_results_list, results_solo))

final_df <- all_results %>%
  mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", "")) %>%
  separate(Term, into = c("Variable", "Month"),
           sep = "_", fill = "right")

cat("\n✔ PCR analysis completed successfully\n")
