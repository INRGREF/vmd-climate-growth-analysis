library(FactoMineR)
library(dplyr)
library(tidyr)

# -----------------------------
# 1. PATHS
# -----------------------------
doss_precip <- "results/VMD_precipitation_modes/Cycle_analysis/Merged_cycles/Biological_year"
doss_temp   <- "results/VMD_temperature_modes/Cycle_analysis/Merged_cycles/Biological_year"

# -----------------------------
# 2. FILE LISTING
# -----------------------------
files_p <- list.files(doss_precip, pattern = "_bio\\.txt$", full.names = TRUE)
files_t <- list.files(doss_temp, pattern = "_bio\\.txt$", full.names = TRUE)

# -----------------------------
# 3. GLOBAL STORAGE INITIALIZATION
# -----------------------------
liste_scores_acp_all <- list()

# -----------------------------
# 4. ROBUST PCA FUNCTION (PVP Criterion & Variance Display)
# -----------------------------
run_acp <- function(data_final, nom_cyc) {
  
  if (is.null(data_final) || nrow(data_final) < 20) return(NULL)
  
  years <- data_final$Year
  X <- data_final %>% select(-Year)
  X[is.na(X)] <- 0
  
  if (ncol(X) == 1) {
    cat("   ->", nom_cyc, ": 1 PC retained (100% variance)\n")
    score_val <- as.data.frame(scale(X[]))
    colnames(score_val) <- paste0(nom_cyc, "_CP1")
    score_val$Year <- years
    
    loadings <- matrix(1, nrow = 1, ncol = 1)
    rownames(loadings) <- colnames(X)
    colnames(loadings) <- paste0(nom_cyc, "_CP1")
    
    return(list(scores = score_val, loadings = loadings))
    
  } else {
    res <- PCA(X, scale.unit = TRUE, ncp = ncol(X), graph = FALSE)
    
    # 1. Eigenvalue extraction and Variance %
    eig <- res$eig[, 1]
    pct_var <- res$eig[, 2] # Variance percentage per PC
    
    # 2. PVP Criterion (Cumulative Product)
    cum_prod_eig <- cumprod(eig)
    ncp <- max(which(cum_prod_eig >= 1))
    ncp <- max(1, min(ncp, length(eig)))
    
    # 3. Calculation of cumulative variance for retained PCs
    cum_var_retained <- sum(pct_var[1:ncp])
    
    cat("   ->", nom_cyc, ":", ncp, "PC(s) retained | Cum. Variance:", round(cum_var_retained, 2), "%\n")
    
    scores <- as.data.frame(res$ind$coord[, 1:ncp, drop = FALSE])
    colnames(scores) <- paste0(nom_cyc, "_CP", 1:ncp)
    scores$Year <- years
    
    loadings <- res$var$coord[, 1:ncp, drop = FALSE]
    rownames(loadings) <- colnames(X)
    
    return(list(scores = scores, loadings = loadings))
  }
}

# -----------------------------
# 5. COMBINED PCA (P + T)
# -----------------------------
cat("\n--- STARTING COMBINED ANALYSES ---\n")
for (f_p in files_p) {
  # nom_cyc contient déjà "cycle_X" d'après le nom du fichier
  nom_cyc <- gsub("precipitation_|_bio\\.txt", "", basename(f_p))
  f_t <- file.path(doss_temp, paste0("temperature_", nom_cyc, "_bio.txt"))
  
  if (!file.exists(f_t)) next
  
  data_p <- read.table(f_p, header = TRUE, check.names = FALSE)
  data_t <- read.table(f_t, header = TRUE, check.names = FALSE)
  
  colnames(data_p)[colnames(data_p) != "Year"] <- paste0(colnames(data_p)[colnames(data_p) != "Year"], "_P")
  colnames(data_t)[colnames(data_t) != "Year"] <- paste0(colnames(data_t)[colnames(data_t) != "Year"], "_T")
  
  data_comb <- inner_join(data_p, data_t, by = "Year")
  
  res <- run_acp(data_comb, nom_cyc)
  if (!is.null(res)) {
    liste_scores_acp_all[[nom_cyc]] <- res
  }
}

# -----------------------------
# 6. SOLO PRECIPITATION PCA
# -----------------------------
cat("\n--- STARTING SOLO PRECIP ANALYSES ---\n")
for (f_p in files_p) {
  nom_cyc <- gsub("precipitation_|_bio\\.txt", "", basename(f_p))
  if (nom_cyc %in% names(liste_scores_acp_all)) next
  
  data_p <- read.table(f_p, header = TRUE, check.names = FALSE)
  colnames(data_p)[colnames(data_p) != "Year"] <- paste0(colnames(data_p)[colnames(data_p) != "Year"], "_P")
  
  res <- run_acp(data_p, nom_cyc)
  if (!is.null(res)) {
    liste_scores_acp_all[[nom_cyc]] <- res
  }
}

# -----------------------------
# 7. SOLO TEMPERATURE PCA
# -----------------------------
cat("\n--- STARTING SOLO TEMP ANALYSES ---\n")
for (f_t in files_t) {
  nom_cyc <- gsub("temperature_|_bio\\.txt", "", basename(f_t))
  if (nom_cyc %in% names(liste_scores_acp_all)) next
  
  data_t <- read.table(f_t, header = TRUE, check.names = FALSE)
  colnames(data_t)[colnames(data_t) != "Year"] <- paste0(colnames(data_t)[colnames(data_t) != "Year"], "_T")
  
  res <- run_acp(data_t, nom_cyc)
  if (!is.null(res)) {
    liste_scores_acp_all[[nom_cyc]] <- res
  }
}

# -----------------------------
# 8. OUTPUT SUMMARY
# -----------------------------
cat("\n✔ PCA processing complete.\n")
print(names(liste_scores_acp_all))
