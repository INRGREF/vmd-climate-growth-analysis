# =========================================================
# 04 - PCR REGRESSION 
# =========================================================

# -----------------------------
# 1. PATHS AND FILE LISTING
# -----------------------------


fichiers_c <- list.files(output_dir, pattern = "\\.txt$", full.names = TRUE)

# -----------------------------
# 2. INITIALIZATION
# -----------------------------
final_results_list <- list()
noms_acp_dispo <- names(liste_scores_acp_all)
start_year <- 1902
end_year   <- 2001
n_boot     <- 2000

# -----------------------------
# 3. PCR PROCESSING LOOP
# -----------------------------
for (f in fichiers_c) {
  
  # A. Extract exact cycle ID (e.g., "28")
  id_croiss <- gsub("cycle_|_mean\\.txt|temperature_|precipitation_", "", basename(f))
  if (grepl("trend", id_croiss)) id_croiss <- "trend"
  
  # B. Exact matching with PCA keys
  nom_acp_match <- noms_acp_dispo[grepl(paste0("_", id_croiss, "$"), noms_acp_dispo)]
  
  if (id_croiss == "trend") {
    nom_acp_match <- noms_acp_dispo[grepl("trend", noms_acp_dispo)]
  }
  
  if (length(nom_acp_match) == 0) next
  
  cat("✅ Processing:", id_croiss, "with PCA key:", nom_acp_match[1], "\n")
  nom_key <- nom_acp_match[1]
  
  # -----------------------------
  # 4. DATA LOADING AND SYNC
  # -----------------------------
  df_growth <- read.table(f, header = TRUE) %>% 
    select(Year, std) %>% filter(Year >= start_year & Year <= end_year)
  
  acp_data  <- liste_scores_acp_all[[nom_key]]
  df_scores <- acp_data$scores %>% filter(Year >= start_year & Year <= end_year)
  loadings  <- as.matrix(acp_data$loadings)
  
  df_pcr <- inner_join(df_growth, df_scores, by = "Year")
  if (nrow(df_pcr) < 15) next
  
  # -----------------------------
  # 5. BOOTSTRAP PCR
  # -----------------------------
  vars_cp <- setdiff(colnames(df_scores), "Year")
  formule <- as.formula(paste("std ~", paste(vars_cp, collapse = " + ")))
  boot_coefs_cp <- matrix(NA, nrow = n_boot, ncol = length(vars_cp))
  
  set.seed(123)
  for (i in 1:n_boot) {
    df_b <- df_pcr[sample(nrow(df_pcr), replace = TRUE), ]
    mod <- lm(formule, data = df_b)
    cf <- coef(mod)[-1]
    if (length(cf) == length(vars_cp)) boot_coefs_cp[i, ] <- cf
  }
  
  # -----------------------------
  # 6. BACK-TRANSFORMATION
  # -----------------------------
  boot_coefs_cp[is.na(boot_coefs_cp)] <- 0
  mean_coefs_cp <- colMeans(boot_coefs_cp)
  coefs_mensuels <- as.vector(mean_coefs_cp %*% t(loadings))
  boot_mensuels <- boot_coefs_cp %*% t(loadings)
  inf95 <- apply(boot_mensuels, 2, quantile, 0.025, na.rm = TRUE)
  sup95 <- apply(boot_mensuels, 2, quantile, 0.975, na.rm = TRUE)
  
  final_results_list[[id_croiss]] <- data.frame(
    Cycle = id_croiss, Terme = rownames(loadings),
    Coef = coefs_mensuels, Inf95 = as.vector(inf95), Sup95 = as.vector(sup95)
  )
}

# -----------------------------
# 7. FINAL DATA RECONSTRUCTION
# -----------------------------
final_df_pcr <- bind_rows(final_results_list) %>%
  mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", "")) %>%
  mutate(
    Variable = ifelse(grepl("_P$|_P_|^P_|precip", Terme, ignore.case=T), "P", "T"),
    Mois = stringr::str_extract(Terme, "Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|Fev|Avr|Mai|Aou")
  )

# -----------------------------
# 8. OUTPUT SUMMARY
# -----------------------------
cat("\nDone! Processed cycles:", paste(unique(final_df_pcr$Cycle), collapse = ", "), "\n")
