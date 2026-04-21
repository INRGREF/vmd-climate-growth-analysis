
# =========================================================
# 03 - PCA ANALYSIS ON VMD-BASED CLIMATE CYCLES
# Temperature + Precipitation (bio-climatic year)
# Pinus halepensis study
# =========================================================

library(FactoMineR)
library(dplyr)


# =========================================================
# Paths aligned with VMD workflow
# =========================================================

base_temp_folder <- "results/VMD_temperature_modes"
base_precip_folder <- "results/VMD_precipitation_modes"

temp_dir <- file.path(base_temp_folder, "Cycle_analysis")
precip_dir <- file.path(base_precip_folder, "Cycle_analysis")

# -----------------------------
# 2. Function: transform to biological year
# NOTE: To apply a different definition of the biological year, only the
# df_prev (end-of-year months) and df_curr (start-of-year months) blocks
# need to be modified; all subsequent processing steps remain unchanged.
# -----------------------------

process_bioclimate_year <- function(input_dir) {
  
  files <- list.files(input_dir, pattern = "^Cycle_.*\\.txt$", full.names = TRUE)
  files <- files[!grepl("\\.transf\\.txt$", files)]
  
  for (f in files) {
    
    cat("\nProcessing file:", basename(f), "\n")
    
    df <- read.table(f, header = TRUE, check.names = FALSE)
    
    # -----------------------------
    # Month detection (robust naming)
    # -----------------------------
    
    Oct <- grep("Oct", names(df), value = TRUE)
    Nov <- grep("Nov", names(df), value = TRUE)
    Dec <- grep("Dec", names(df), value = TRUE)
    Jan <- grep("Jan", names(df), value = TRUE)
    Feb <- grep("Fev|Feb", names(df), value = TRUE)
    Mar <- grep("Mar", names(df), value = TRUE)
    Apr <- grep("Avr|Apr", names(df), value = TRUE)
    May <- grep("Mai|May", names(df), value = TRUE)
    Jun <- grep("Jun", names(df), value = TRUE)
    Jul <- grep("Jul", names(df), value = TRUE)
    Aug <- grep("Aou|Aug", names(df), value = TRUE)
    Sep <- grep("Sep", names(df), value = TRUE)
    
    # -----------------------------
    # Previous year (Oct–Dec)
    # -----------------------------
    
    prev <- df[, c("Year", Oct, Nov, Dec)]
    prev$Year <- prev$Year + 1
    colnames(prev) <- c("Year", "Oct_prev", "Nov_prev", "Dec_prev")
    
    # -----------------------------
    # Current year (Jan–Sep)
    # -----------------------------
    
    curr <- df[, c("Year", Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep)]
    colnames(curr) <- c(
      "Year",
      "Jan", "Feb", "Mar", "Apr", "May",
      "Jun", "Jul", "Aug", "Sep"
    )
    
    # -----------------------------
    # Merge biological year
    # -----------------------------
    
    bio_year <- merge(prev, curr, by = "Year")
    
    # -----------------------------
    # Export
    # -----------------------------
    
    out_file <- gsub("\\.txt$", ".transf.txt", f)
    
    write.table(
      bio_year,
      file = out_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}

# -----------------------------
# 3. Run transformation
# -----------------------------

process_bioclimate_year(precip_dir)
process_bioclimate_year(temp_dir)

cat("\n✔ Biological year transformation completed\n")

# =========================================================
# 4. PCA ON EACH VMD CLIMATE CYCLE
# =========================================================

files_precip <- list.files(precip_dir, pattern = "\\.transf\\.txt$", full.names = TRUE)

pca_results <- list()

for (f_p in files_precip) {
  
  cycle_name <- gsub("_Precipitation\\.transf\\.txt", "", basename(f_p))
  cat("\n--- PCA for cycle:", cycle_name, "---\n")
  
  # -----------------------------
  # Load precipitation
  # -----------------------------
  
  data_p <- read.table(f_p, header = TRUE, check.names = FALSE)
  
  # -----------------------------
  # Match temperature file
  # -----------------------------
  
  f_t <- file.path(temp_dir, paste0(cycle_name, "_Temperature.transf.txt"))
  
  if (!file.exists(f_t)) {
    cat("Missing temperature file for:", cycle_name, "\n")
    next
  }
  
  data_t <- read.table(f_t, header = TRUE, check.names = FALSE)
  
  # -----------------------------
  # Merge climate variables
  # -----------------------------
  
  climate_data <- inner_join(data_p, data_t, by = "Year", suffix = c("_P", "_T"))
  years <- climate_data$Year
  
  # -----------------------------
  # PCA input matrix
  # -----------------------------
  
  pca_input <- climate_data %>% select(-Year)
  pca_input[is.na(pca_input)] <- 0
  
  # -----------------------------
  # PCA computation
  # -----------------------------
  
  res_pca <- PCA(
    pca_input,
    scale.unit = TRUE,
    ncp = ncol(pca_input),
    graph = FALSE
  )
  
  # -----------------------------
  # Component selection (Guiot criterion)
  # -----------------------------
  
  eigenvalues <- res_pca$eig[, 1]
  threshold <- mean(eigenvalues)
  n_pc <- sum(eigenvalues > threshold)
  
  cat("Selected PCs:", n_pc, "/", length(eigenvalues), "\n")
  
  # -----------------------------
  # Extract scores
  # -----------------------------
  
  scores <- as.data.frame(res_pca$ind$coord[, 1:n_pc])
  colnames(scores) <- paste0(cycle_name, "_PC", 1:n_pc)
  scores$Year <- years
  
  # -----------------------------
  # Store results
  # -----------------------------
  
  pca_results[[cycle_name]] <- list(
    scores = scores,
    loadings = res_pca$var$coord[, 1:n_pc],
    variance = res_pca$eig[n_pc, 3]
  )
}

cat("\n✔ PCA analysis completed for all VMD cycles\n")
