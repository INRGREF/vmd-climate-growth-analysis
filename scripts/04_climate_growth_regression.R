# =========================================================
# 04_climate_growth_regression.R
# Multi-scale Climate–Growth Relationships using PCR Bootstrap
# Pinus halepensis (Tunisia)
# =========================================================

# ---------------------------------------------------------
# 1. Load libraries
# ---------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------------------------
# 2. Define paths (aligned with project structure)
# ---------------------------------------------------------

# Growth chronologies (by cycle)
growth_dir <- "results/VMD_growth/CHRONOLOGIES_GROWTH_CYCLES"

# Climate PCA results (already in memory: liste_scores_acp)

# Climate cycle folders
temp_dir   <- "results/VMD_temperature_modes/Cycle_analysis"
precip_dir <- "results/VMD_precipitation_modes/Cycle_analysis"

# ---------------------------------------------------------
# 3. PCR Bootstrap on PCA scores (multi-variable cycles)
# ---------------------------------------------------------

growth_files <- list.files(growth_dir, pattern = "Chrono_.*\\.csv$", full.names = TRUE)

final_results_list <- list()
n_boot <- 1000

for (f in growth_files) {
  
  # Identify cycle
  cycle_raw <- gsub("Chrono_|.csv", "", basename(f))
  cycle_name <- ifelse(cycle_raw == "Tendance", "Cycle_Trend", cycle_raw)
  
  if (!(cycle_name %in% names(liste_scores_acp))) next
  
  cat("PCR regression for:", cycle_name, "\n")
  
  # Load data
  df_growth <- read.csv(f) %>% select(Year, std)
  df_scores <- liste_scores_acp[[cycle_name]]$scores
  loadings  <- liste_scores_acp[[cycle_name]]$loadings
  
  # Merge and clean
  df_pcr <- inner_join(df_growth, df_scores, by = "Year") %>% drop_na()
  if(nrow(df_pcr) < 20) next
  
  # Regression formula
  vars_cp <- setdiff(colnames(df_scores), "Year")
  formula <- as.formula(paste("std ~", paste(vars_cp, collapse = " + ")))
  
  # Bootstrap
  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))
  
  set.seed(123)
  for(i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(formula, data = df_b)
    boot_coefs_cp[i, ] <- coef(mod)[-1]
  }
  
  # Back-transformation to monthly scale
  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)
  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  
  boot_monthly <- boot_coefs_cp %*% t(loadings)
  inf95 <- apply(boot_monthly, 2, quantile, 0.025, na.rm = TRUE)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975, na.rm = TRUE)
  
  # Format output
  res_df <- data.frame(
    Cycle = cycle_name,
    Term = rownames(loadings),
    Coef = as.vector(coefs_monthly),
    Inf95 = inf95,
    Sup95 = sup95
  )
  
  final_results_list[[cycle_name]] <- res_df
}

# ---------------------------------------------------------
# 4. PCR for single-variable cycles 
# ---------------------------------------------------------

results_solo <- list()

cycles_solo <- list(
  "Cycle_10_ans" = list(dir = temp_dir, var = "T"),
  "Cycle_12_ans" = list(dir = precip_dir, var = "P"),
  "Cycle_13_ans" = list(dir = precip_dir, var = "P"),
  "Cycle_15_ans" = list(dir = precip_dir, var = "P")
)

for (cyc in names(cycles_solo)) {
  
  var_code <- cycles_solo[[cyc]]$var
  suffix <- ifelse(var_code == "T", "_Temperature.transf.txt", "_Precipitation.transf.txt")
  f_clim <- file.path(cycles_solo[[cyc]]$dir, paste0(cyc, suffix))
  
  if(!file.exists(f_clim)) {
    cat("Missing file:", f_clim, "\n")
    next
  }
  
  cat("Solo PCR for:", cyc, "\n")
  
  data_clim <- read.table(f_clim, header = TRUE, check.names = FALSE)
  data_input <- data_clim %>% select(-Year)
  data_input[is.na(data_input)] <- 0
  
  # PCA
  res.pca <- FactoMineR::PCA(data_input, scale.unit = TRUE, ncp = 12, graph = FALSE)
  
  eigenvalues <- res.pca$eig[,1]
  nb_cp_keep <- sum(eigenvalues > mean(eigenvalues))
  if(nb_cp_keep < 1) nb_cp_keep <- 1
  
  scores <- as.data.frame(res.pca$ind$coord[, 1:nb_cp_keep, drop=FALSE])
  colnames(scores) <- paste0("Dim", 1:nb_cp_keep)
  scores$Year <- data_clim$Year
  
  loadings <- res.pca$var$coord[, 1:nb_cp_keep, drop=FALSE]
  
  # Growth
  f_growth <- file.path(growth_dir, paste0("Chrono_", cyc, ".csv"))
  if(!file.exists(f_growth)) next
  
  df_growth <- read.csv(f_growth) %>% select(Year, std)
  df_pcr <- inner_join(df_growth, scores, by = "Year") %>% drop_na()
  if(nrow(df_pcr) < 20) next
  
  # Bootstrap
  vars_cp <- setdiff(colnames(scores), "Year")
  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))
  
  set.seed(123)
  for(i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(std ~ ., data = df_b %>% select(all_of(vars_cp), std))
    boot_coefs_cp[i, ] <- coef(mod)[-1]
  }
  
  # Back-transformation
  mean_coefs_cp <- colMeans(boot_coefs_cp, na.rm = TRUE)
  coefs_monthly <- mean_coefs_cp %*% t(loadings)
  
  boot_monthly <- boot_coefs_cp %*% t(loadings)
  inf95 <- apply(boot_monthly, 2, quantile, 0.025)
  sup95 <- apply(boot_monthly, 2, quantile, 0.975)
  
  month_names <- paste0(var_code, "_", rownames(loadings))
  
  res_df <- data.frame(
    Cycle = cyc,
    Term = month_names,
    Coef = as.vector(coefs_monthly),
    Inf95 = inf95,
    Sup95 = sup95
  )
  
  results_solo[[cyc]] <- res_df
}

# ---------------------------------------------------------
# 5. Merge all results
# ---------------------------------------------------------

all_results <- bind_rows(c(final_results_list, results_solo))

final_df <- all_results %>%
  mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", "")) %>%
  separate(Term, into = c("Variable", "Month"), sep = "_", fill = "right")

cat("\nPCR analysis completed successfully.\n")
