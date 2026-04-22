# =========================================================
# 06 - VARIANCE CONTRIBUTION ANALYSIS (FINAL CORRIGÉ)
# VMD modes + growth cycles
# Pinus halepensis (Tunisia)
# =========================================================

# ---------------------------------------------------------
# 1. PATHS (HARMONIZED PIPELINE)
# ---------------------------------------------------------

results_dir <- "results"

growth_dir <- file.path(results_dir, "VMD_growth")

temp_dir   <- file.path(results_dir, "VMD_climate_temperature")
prec_dir   <- file.path(results_dir, "VMD_climate_precipitation")

cycles_dir <- file.path(growth_dir, "SYNTHESE_GROWTH_CYCLES")

output_dir <- file.path(growth_dir, "variance_analysis")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------
# 2. VARIANCE OF VMD MODES (CLIMATE + GROWTH)
# ---------------------------------------------------------

compute_variance_modes <- function(dir_path) {
  
  files <- list.files(dir_path,
                      pattern = "VMD_.*\\.csv$|VMD_.*\\.txt$",
                      full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  var_list <- list()
  mode_names <- NULL
  
  for (f in files) {
    
    data <- tryCatch(
      read.table(f, header = TRUE, check.names = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(data)) next
    
    # remove year column safely
    data <- data[, !names(data) %in% c("Year", "Annee", "Year.")]
    
    if (ncol(data) < 1) next
    
    vars <- apply(data, 2, var, na.rm = TRUE)
    
    var_list[[length(var_list) + 1]] <- vars
    mode_names <- union(mode_names, names(vars))
  }
  
  if (length(var_list) == 0) return(NULL)
  
  var_matrix <- matrix(NA,
                       nrow = length(var_list),
                       ncol = length(mode_names))
  
  colnames(var_matrix) <- mode_names
  
  for (i in seq_along(var_list)) {
    var_matrix[i, names(var_list[[i]])] <- var_list[[i]]
  }
  
  mean_var <- colMeans(var_matrix, na.rm = TRUE)
  total_var <- sum(mean_var, na.rm = TRUE)
  
  data.frame(
    Mode = names(mean_var),
    MeanVariance = mean_var,
    Contribution_percent = (mean_var / total_var) * 100
  )
}

# ---------------------------------------------------------
# 3. VARIANCE CONTRIBUTION (CLIMATE MODES)
# ---------------------------------------------------------

cat("\nComputing climate variance contribution...\n")

climate_temp_var <- compute_variance_modes(temp_dir)
climate_prec_var <- compute_variance_modes(prec_dir)

if (!is.null(climate_temp_var)) {
  write.csv(climate_temp_var,
            file.path(output_dir, "variance_temperature_modes.csv"),
            row.names = FALSE)
}

if (!is.null(climate_prec_var)) {
  write.csv(climate_prec_var,
            file.path(output_dir, "variance_precipitation_modes.csv"),
            row.names = FALSE)
}

# ---------------------------------------------------------
# 4. VARIANCE CONTRIBUTION (GROWTH MODES)
# ---------------------------------------------------------

cat("\nComputing growth variance contribution...\n")

growth_modes <- list.files(growth_dir,
                           pattern = "VMD_.*\\.txt$",
                           full.names = TRUE)

if (length(growth_modes) > 0) {
  
  var_list <- list()
  mode_names <- NULL
  
  for (f in growth_modes) {
    
    data <- read.table(f, header = TRUE, check.names = FALSE)
    
    data <- data[, !names(data) %in% c("Year", "Annee")]
    
    vars <- apply(data, 2, var, na.rm = TRUE)
    
    var_list[[length(var_list) + 1]] <- vars
    mode_names <- union(mode_names, names(vars))
  }
  
  var_matrix <- matrix(NA,
                       nrow = length(var_list),
                       ncol = length(mode_names))
  
  colnames(var_matrix) <- mode_names
  
  for (i in seq_along(var_list)) {
    var_matrix[i, names(var_list[[i]])] <- var_list[[i]]
  }
  
  mean_var <- colMeans(var_matrix, na.rm = TRUE)
  total_var <- sum(mean_var, na.rm = TRUE)
  
  results_modes <- data.frame(
    Mode = names(mean_var),
    MeanVariance = mean_var,
    Contribution_percent = (mean_var / total_var) * 100
  )
  
  write.csv(results_modes,
            file.path(output_dir, "variance_contribution_growth_modes.csv"),
            row.names = FALSE)
}

# ---------------------------------------------------------
# 5. VARIANCE CONTRIBUTION (GROWTH CYCLES)
# ---------------------------------------------------------

cat("\nComputing growth cycles variance...\n")

files_cycles <- list.files(cycles_dir,
                          pattern = "Synthese_.*\\.txt$",
                          full.names = TRUE)

file_trend <- file.path(cycles_dir,
                        "Synthese_Tous_Arbres_Tendance.txt")

files_all <- c(files_cycles, file_trend)
cycle_names <- c(
  gsub(".*Cycle_|\\.txt", "", basename(files_cycles)),
  "Trend"
)

mean_variances <- numeric(length(files_all))

for (i in seq_along(files_all)) {
  
  if (!file.exists(files_all[i])) next
  
  data <- read.table(files_all[i], header = TRUE)
  data <- data[, !names(data) %in% c("Year", "Annee")]
  
  vars <- apply(data, 2, var, na.rm = TRUE)
  mean_variances[i] <- mean(vars, na.rm = TRUE)
}

total_var <- sum(mean_variances, na.rm = TRUE)

results_cycles <- data.frame(
  Cycle = cycle_names,
  MeanVariance = mean_variances,
  Contribution_percent = (mean_variances / total_var) * 100
)

# sort cycles
cycles_only <- results_cycles[results_cycles$Cycle != "Trend", ]
cycles_only <- cycles_only[order(as.numeric(cycles_only$Cycle)), ]

results_cycles <- rbind(
  cycles_only,
  results_cycles[results_cycles$Cycle == "Trend", ]
)

write.csv(results_cycles,
          file.path(output_dir, "variance_contribution_cycles.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# 6. FREQUENCY TABLE (TREE DISTRIBUTION)
# ---------------------------------------------------------

period_file <- file.path(growth_dir,
                         "tree_vmd_periods_rounded.txt")

if (file.exists(period_file)) {
  
  period_table <- read.table(period_file,
                             header = TRUE,
                             check.names = FALSE)
  
  classify_cycle <- function(x) {
    if (is.na(x)) return(NA)
    if (x == "tendance" || x == "trend") return("Trend")
    x <- as.numeric(x)
    if (x >= 2  & x <= 5)  return("2–5")
    if (x >= 6  & x <= 9)  return("6–9")
    if (x >= 10 & x <= 15) return("10–15")
    if (x >= 16 & x <= 25) return("16–25")
    return(NA)
  }
  
  class_matrix <- apply(period_table, c(1,2), classify_cycle)
  
  classes <- c("2–5","6–9","10–15","16–25","Trend")
  cycles <- colnames(period_table)
  
  freq_table <- matrix(0,
                       nrow = length(classes),
                       ncol = length(cycles),
                       dimnames = list(classes, cycles))
  
  for (i in seq_along(cycles)) {
    for (j in seq_along(classes)) {
      freq_table[j,i] <- sum(class_matrix[,i] == classes[j], na.rm = TRUE)
    }
  }
  
  write.csv(as.data.frame(freq_table),
            file.path(output_dir,
                      "frequency_classes_per_cycle.csv"),
            row.names = TRUE)
}

cat("\n✔ Variance contribution analysis completed\n")
