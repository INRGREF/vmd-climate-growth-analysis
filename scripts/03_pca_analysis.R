# =========================================================
# 03 - PCA FINAL
# =========================================================

library(FactoMineR)
library(dplyr)

results_dir <- "results"

temp_dir <- file.path(results_dir, "VMD_climate_temperature", "Cycle_analysis")
prec_dir <- file.path(results_dir, "VMD_climate_precipitation", "Cycle_analysis")

process <- function(dir) {

  files <- list.files(dir, pattern = "_Cycle_analysis.txt$", full.names = TRUE)

  pca_list <- list()

  for (f in files) {

    df <- read.table(f, header = TRUE)
    cycle <- gsub("_Cycle_analysis.txt", "", basename(f))

    x <- df[, -1]
    x[is.na(x)] <- 0

    res <- PCA(x, scale.unit = TRUE, graph = FALSE)

    eig <- res$eig[,1]
    ncp <- max(1, sum(eig > mean(eig)))

    pca_list[[cycle]] <- list(
      scores = as.data.frame(res$ind$coord[,1:ncp]),
      loadings = res$var$coord[,1:ncp]
    )
  }

  return(pca_list)
}

pca_temp <- process(temp_dir)
pca_prec <- process(prec_dir)

pca_results <- c(pca_temp, pca_prec)

cat("\n✔ PCA OK\n")
