# =========================================================
# 06_variance_contribution.R
# Variance contribution analysis of VMD modes and cycles
# Pinus halepensis growth dataset (Tunisia)
# =========================================================

# ---------------------------------------------------------
# 1. Define paths (aligned with project structure)
# ---------------------------------------------------------

growth_dir <- "results/VMD_growth"
modes_dir  <- growth_dir                    # fichiers VMD arbres
cycles_dir <- file.path(growth_dir, "SYNTHESE_GROWTH_CYCLES")

output_dir <- file.path(growth_dir, "variance_analysis")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------
# 2. Variance contribution of VMD modes (tree-level)
# ---------------------------------------------------------

files <- list.files(modes_dir,
                    pattern = "VMD_.*\\.txt$",
                    full.names = TRUE)

var_list <- list()
mode_names <- NULL

for (i in seq_along(files)) {
  
  data <- read.table(files[i], header = TRUE)
  
  # Remove year column
  data <- data[, !(names(data) %in% "Annee")]
  
  # Variance per mode
  vars <- apply(data, 2, var, na.rm = TRUE)
  
  var_list[[i]] <- vars
  mode_names <- union(mode_names, names(vars))
}

# Align modes across trees
var_matrix <- matrix(NA,
                     nrow = length(var_list),
                     ncol = length(mode_names))

colnames(var_matrix) <- mode_names

for(i in seq_along(var_list)){
  var_matrix[i, names(var_list[[i]])] <- var_list[[i]]
}

# Mean variance and contribution
mean_var_modes <- colMeans(var_matrix, na.rm = TRUE)
total_var <- sum(mean_var_modes)

contribution_modes <- (mean_var_modes / total_var) * 100

results_modes <- data.frame(
  Mode = names(mean_var_modes),
  MeanVariance = mean_var_modes,
  Contribution_percent = contribution_modes
)

print(results_modes)

# Save
write.csv(results_modes,
          file.path(output_dir, "variance_contribution_modes.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# 3. Variance contribution by growth cycle
# ---------------------------------------------------------

files_cycles <- list.files(cycles_dir,
                          pattern = "Synthese_Tous_Arbres_Cycle_.*\\.txt$",
                          full.names = TRUE)

file_trend <- file.path(cycles_dir, "Synthese_Tous_Arbres_Tendance.txt")

files_all <- c(files_cycles, file_trend)

cycle_names <- c(
  gsub(".*Cycle_|\\.txt", "", basename(files_cycles)),
  "Trend"
)

mean_variances <- numeric(length(files_all))

for(i in seq_along(files_all)){
  
  data <- read.table(files_all[i], header = TRUE)
  
  data <- data[, !(names(data) %in% "Annee")]
  
  vars <- apply(data, 2, var, na.rm = TRUE)
  
  mean_variances[i] <- mean(vars, na.rm = TRUE)
}

total_var <- sum(mean_variances, na.rm = TRUE)
contributions <- (mean_variances / total_var) * 100

results_cycles <- data.frame(
  Cycle = cycle_names,
  MeanVariance = mean_variances,
  Contribution_percent = contributions
)

# Sort cycles numerically
cycles_only <- results_cycles[results_cycles$Cycle != "Trend", ]
cycles_only <- cycles_only[order(as.numeric(cycles_only$Cycle)), ]

results_cycles <- rbind(
  cycles_only,
  results_cycles[results_cycles$Cycle == "Trend", ]
)

print(results_cycles)

# Save
write.csv(results_cycles,
          file.path(output_dir, "variance_contribution_cycles.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# 4. Frequency of trees per cycle class
# ---------------------------------------------------------

period_table <- read.table(
  file.path(growth_dir, "tree_vmd_periods_rounded.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# Classification function
classify_cycle <- function(x){
  if(is.na(x)) return(NA)
  if(x == "tendance" || x == "trend") return("Trend")
  x <- as.numeric(x)
  if(x >= 2  & x <= 5)  return("2–5")
  if(x >= 6  & x <= 9)  return("6–9")
  if(x >= 10 & x <= 15) return("10–15")
  if(x >= 16 & x <= 25) return("16–25")
  return(NA)
}

class_matrix <- apply(period_table, c(1,2), classify_cycle)

classes <- c("2–5","6–9","10–15","16–25","Trend")
cycles <- colnames(period_table)

freq_table <- matrix(0,
                     nrow = length(classes),
                     ncol = length(cycles),
                     dimnames = list(classes, cycles))

for(i in seq_along(cycles)){
  for(j in seq_along(classes)){
    freq_table[j,i] <- sum(class_matrix[,i] == classes[j], na.rm = TRUE)
  }
}

freq_table <- as.data.frame(freq_table)

print(freq_table)

# Save
write.csv(freq_table,
          file.path(output_dir, "frequency_classes_per_cycle.csv"),
          row.names = TRUE)

cat("\nVariance contribution analysis completed.\n")
