# =========================================================
# 02 - FULL VMD DECOMPOSITION PIPELINE 
# Climate (Temperature + Precipitation)
# Growth (Tree-ring chronologies)
# Pinus halepensis study
# =========================================================

library(VMDecomp)
library(dplyr)
library(dplR)

# =========================================================
# 1. PROJECT STRUCTURE
# =========================================================

results_dir <- "results"

climate_temp_dir <- file.path(results_dir, "VMD_temperature_modes")
climate_prec_dir <- file.path(results_dir, "VMD_precipitation_modes")
growth_dir       <- file.path(results_dir, "VMD_growth")

cycle_temp_dir <- file.path(climate_temp_dir, "Cycle_analysis")
cycle_prec_dir <- file.path(climate_prec_dir, "Cycle_analysis")
cycle_growth_dir <- file.path(growth_dir, "Cycle_analysis")
if (!dir.exists(cycle_growth_dir)) {
  dir.create(cycle_growth_dir, recursive = TRUE)
}

dirs <- c(
  climate_temp_dir, climate_prec_dir,
  growth_dir, cycle_temp_dir, cycle_prec_dir,cycle_growth_dir
)

lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. OPTIMAL K FUNCTION 
# =========================================================

find_k_optimal <- function(series, max_k = 8, threshold_years = 1) {
  
  for (k in 2:max_k) {
    
    res <- vmd(
      series,
      alpha = 2000,
      tau = 0,
      K = k,
      DC = FALSE,
      init = 1,
      tol = 1e-7
    )
    
    # Final frequencies
    omegas <- sort(res$omega[nrow(res$omega), ])
    
    # Avoid division by zero
    omegas <- omegas[omegas > 0]
    
    # Conversion to periods (years)
    periods <- 1 / omegas
    
    # Ascending sort (important)
    periods <- sort(periods)
    
    # Difference between time scales
    min_diff <- min(diff(periods))
    
    # Physical criterion in years
    
    if (min_diff < threshold_years) {
      return(k - 1)
    }
  }
  
  return(max_k)
}

# =========================================================
# 3. CLIMATE VMD FUNCTION (TEMPERATURE + PRECIPITATION)
# =========================================================

run_vmd_climate <- function(file_path, output_dir, label) {
  
  data <- read.table(file_path, header = TRUE)
  
  monthly <- data[, 2:13]
  colnames(monthly) <- c(
    "Jan","Feb","Mar","Apr","May","Jun",
    "Jul","Aug","Sep","Oct","Nov","Dec"
  )
  
  vmd_results <- list()
  
  for (m in colnames(monthly)) {
    
    series <- monthly[[m]]
    k_opt <- find_k_optimal(series)
    
    res <- vmd(series, alpha = 2000, tau = 0,
               K = k_opt, DC = FALSE, init = 1, tol = 1e-7)
    
    vmd_results[[m]] <- list(
      modes = res$u,
      frequencies = res$omega[nrow(res$omega), ],
      K = k_opt
    )
  }
  
  # -----------------------------
  # PERIOD TABLE
  # -----------------------------
  
  maxK <- max(sapply(vmd_results, function(x) x$K))
  
  period_mat <- matrix(NA, 12, maxK)
  rownames(period_mat) <- names(vmd_results)
  colnames(period_mat) <- paste0("IMF_", 1:maxK)
  
  for (m in names(vmd_results)) {
    
    p <- sort(1 / vmd_results[[m]]$frequencies, decreasing = TRUE)
    
    # Numerical stability and data cleaning
    p[!is.finite(p)] <- NA
    p[p > 200] <- NA
    
    period_mat[m, 1:length(p)] <- p
  }
  
  # Rounding to nearest integer years
  period_mat <- round(period_mat, 0)
  
  
  # =========================================================
  # TREND DETECTION (IMF_1 → "trend")
  # =========================================================
  
  # IMPORTANT: convert to character to support "trend" as a label
  
  period_mat <- as.data.frame(period_mat, stringsAsFactors = FALSE)
  
  for (m in rownames(period_mat)) {
    
    periods <- suppressWarnings(as.numeric(period_mat[m, ]))
    periods <- periods[!is.na(periods)]
    
    if (length(periods) < 2) next
    
    # Comparison between IMF1 and remaining cycles
    
    max_cycle <- max(periods[-1], na.rm = TRUE)
    imf1 <- periods[1]
    
    # check for trend condition
    
    if (!is.na(imf1) && imf1 > max_cycle) {
      period_mat[m, 1] <- "trend"
    }
  }
  
  # =========================================================
  # RESTORE NUMERIC TYPES (IMF2+ remain numeric)
    # =========================================================
  
  if (ncol(period_mat) > 1) {
    period_mat[, -1] <- apply(period_mat[, -1], 2, function(x) as.numeric(x))
  }
  
  # =========================================================
  # EXPORT
  # =========================================================
  
  write.table(
    period_mat,
    file = file.path(output_dir, paste0(label, "_vmd_periods_rounded.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
  # -----------------------------
  # EXPORT MODES
  # -----------------------------
  
  modes_dir <- file.path(output_dir, "modes")
  dir.create(modes_dir, showWarnings = FALSE)
  
  for (m in names(vmd_results)) {
    
    modes <- vmd_results[[m]]$modes
    if (nrow(modes) != nrow(data)) modes <- t(modes)
    
    df <- cbind(Year = data[,1], as.data.frame(modes))
    colnames(df) <- c("Year", paste0("IMF_", 1:ncol(modes)))
    
    write.csv(df,
              file.path(modes_dir, paste0("VMD_", m, ".csv")),
              row.names = FALSE)
  }
  
  return(vmd_results)
}

# =========================================================
# 4. RUN CLIMATE VMD
# =========================================================

temp_results <- run_vmd_climate(
  "raw_data/temperature.txt",
  climate_temp_dir,
  "temperature"
)

prec_results <- run_vmd_climate(
  "raw_data/precipitation.txt",
  climate_prec_dir,
  "precipitation"
)

cat("\n✔ Climate VMD completed\n")

# =========================================================
# 5. GROWTH VMD (TREE-RINGS)
# =========================================================

tree_data <- rw_index

vmd_tree_results <- list()

for (i in 2:ncol(tree_data)) {
  
  tree_id <- colnames(tree_data)[i]
  series <- tree_data[, i]
  years <- tree_data[, 1]
  
  valid <- which(!is.na(series))
  clean <- series[valid]
  
  if (length(clean) > 30) {
    
    k_opt <- find_k_optimal(clean)
    
    res <- vmd(clean, alpha = 2000, tau = 0,
               K = k_opt, DC = FALSE, init = 1, tol = 1e-7)
    
    aligned <- matrix(NA, length(years), k_opt)
    modes <- res$u
    if (nrow(modes) != length(clean)) modes <- t(modes)
    
    aligned[valid, ] <- modes
    
    df <- cbind(Year = years, as.data.frame(aligned))
    colnames(df) <- c("Year", paste0("IMF_", 1:k_opt))
    
    write.table(df,
                file.path(growth_dir, paste0("VMD_", tree_id, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    vmd_tree_results[[tree_id]] <- list(
      frequencies = res$omega[nrow(res$omega), ],
      K = k_opt
    )
  }
}

# -----------------------------
# TREE PERIODS
# -----------------------------

maxK <- max(sapply(vmd_tree_results, function(x) x$K))

growth_periods <- matrix(NA, length(vmd_tree_results), maxK)
rownames(growth_periods) <- names(vmd_tree_results)
colnames(growth_periods) <- paste0("IMF_", 1:maxK)

for (t in names(vmd_tree_results)) {
  
  p <- sort(1 / vmd_tree_results[[t]]$frequencies, decreasing = TRUE)
  
  # Numerical stability
  p[!is.finite(p)] <- NA
  p[p > 200] <- NA
  
  growth_periods[t, 1:length(p)] <- p
}

# rounding to nearest full years

growth_periods <- round(growth_periods, 0)

# =========================================================
# TREND DETECTION (TREE-RINGS) 
# =========================================================

growth_periods <- as.data.frame(growth_periods, stringsAsFactors = FALSE)

for (t in rownames(growth_periods)) {
  
  periods <- suppressWarnings(as.numeric(growth_periods[t, ]))
  periods_clean <- periods[!is.na(periods)]
  
  if (length(periods_clean) < 1) next
  
  imf1 <- periods[1]
  max_cycle <- max(periods[-1], na.rm = TRUE)
  
  # ROBUST CORRECTION
  
  if (is.na(imf1) || imf1 > max_cycle) {
    growth_periods[t, 1] <- "trend"
  }
}

# -----------------------------
# EXPORT
# -----------------------------

write.table(
  growth_periods,
  file = file.path(growth_dir, "tree_vmd_periods_rounded.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

cat("\n✔ Growth VMD completed\n")

# =========================================================
# ORGANIZATION OF VMD MODES INTO CLIMATE CYCLES
# (TREND INTEGRATED WITHIN IMF1)
# =========================================================

library(dplyr)

run_cycle_pipeline <- function(base_folder, variable_name, years_ref = 1901:2023) {
  
  library(dplyr)
  
  # -----------------------------
  # PATHS
  # -----------------------------
  output_folder <- file.path(base_folder, "Cycle_analysis")
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  period_file <- file.path(base_folder, paste0(variable_name, "_vmd_periods_rounded.txt"))
  
  if (!file.exists(period_file)) {
    stop("Missing file: ", period_file)
  }
  
  period_table <- read.table(
    period_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  months <- rownames(period_table)
  
  # =========================================================
  # STEP 1: IMF → cycles (per month)
  # =========================================================
  
  for (month in months) {
    
    month_folder <- file.path(output_folder, paste0("Analysis_", month))
    dir.create(month_folder, recursive = TRUE, showWarnings = FALSE)
    
    vmd_file <- file.path(base_folder, "modes", paste0("VMD_", month, ".csv"))
    
    if (!file.exists(vmd_file)) next
    
    vmd_data <- read.csv(vmd_file, check.names = FALSE)
    names(vmd_data) <- gsub('"', "", names(vmd_data))
    
    cycle_values <- na.omit(as.character(period_table[month, ]))
    
    imf_index <- 1
    
    for (cycle in cycle_values) {
      
      col_name <- if (cycle == "trend") "IMF_1" else paste0("IMF_", imf_index)
      
      if (!col_name %in% names(vmd_data)) {
        imf_index <- imf_index + 1
        next
      }
      
      values <- vmd_data[[col_name]]
      
      if (all(is.na(values))) {
        imf_index <- imf_index + 1
        next
      }
      
      out_file <- file.path(
        month_folder,
        if (cycle == "trend") {
          paste0(month, "_trend.txt")
        } else {
          paste0(month, "_cycle_", cycle, "yr.txt")
        }
      )
      
      write.table(
        data.frame(Year = vmd_data$Year, Value = values),
        file = out_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      
      imf_index <- imf_index + 1
    }
  }
  
  # =========================================================
  # STEP 2: GROUP BY CYCLE
  # =========================================================
  
  base_cycle <- output_folder
  grouped_dir <- file.path(base_cycle, "Grouped_cycles")
  dir.create(grouped_dir, showWarnings = FALSE, recursive = TRUE)
  
  month_folders <- list.dirs(base_cycle, recursive = FALSE, full.names = TRUE)
  month_folders <- month_folders[grepl("Analysis_", month_folders)]
  
  cycle_storage <- list()
  
  for (mf in month_folders) {
    
    month_name <- gsub("Analysis_", "", basename(mf))
    files <- list.files(mf, pattern = "\\.txt$", full.names = TRUE)
    
    for (f in files) {
      
      data <- read.table(f, header = TRUE)
      
      if (grepl("trend", f)) {
        cycle_name <- "trend"
      } else {
        cycle_name <- paste0("cycle_", gsub(".*cycle_|yr\\.txt", "", f))
      }
      
      if (!cycle_name %in% names(cycle_storage)) {
        cycle_storage[[cycle_name]] <- list()
      }
      
      cycle_storage[[cycle_name]][[month_name]] <- data
    }
  }
  
  # export grouped
  for (cycle in names(cycle_storage)) {
    
    cycle_dir <- file.path(grouped_dir, paste0("Cycle_", cycle))
    dir.create(cycle_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (month in names(cycle_storage[[cycle]])) {
      
      write.table(
        cycle_storage[[cycle]][[month]],
        file = file.path(cycle_dir, paste0(month, "_", cycle, ".txt")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }
  
  # =========================================================
  # STEP 3: MERGE MATRICES
  # =========================================================
  
  merged_dir <- file.path(base_cycle, "Merged_cycles")
  dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)
  
  months_ref <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")
  
  analysis_folders <- list.dirs(base_cycle, recursive = FALSE, full.names = TRUE)
  analysis_folders <- analysis_folders[grepl("Analysis_", analysis_folders)]
  
  cycle_storage_mat <- list()
  
  for (folder in analysis_folders) {
    
    month <- gsub("Analysis_", "", basename(folder))
    files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
    
    for (f in files) {
      
      data <- read.table(f, header = TRUE)
      
      if (grepl("trend", f)) {
        cycle_name <- "trend"
      } else {
        cycle_name <- paste0("cycle_", gsub(".*cycle_|yr\\.txt", "", f))
      }
      
      if (!cycle_name %in% names(cycle_storage_mat)) {
        
        mat <- matrix(NA, nrow = length(years_ref), ncol = 12)
        rownames(mat) <- years_ref
        colnames(mat) <- months_ref
        
        cycle_storage_mat[[cycle_name]] <- mat
      }
      
      idx_y <- match(data$Year, years_ref)
      idx_m <- match(month, months_ref)
      
      cycle_storage_mat[[cycle_name]][idx_y, idx_m] <- data$Value
    }
  }
  
  for (cycle in names(cycle_storage_mat)) {
    
    write.table(
      cycle_storage_mat[[cycle]],
      file = file.path(merged_dir, paste0(variable_name, "_", cycle, ".txt")),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }
  
  cat("\n✔ Completed:", variable_name, "\n")
}

run_cycle_pipeline("results/VMD_temperature_modes", "temperature")
run_cycle_pipeline("results/VMD_precipitation_modes", "precipitation")

# =========================================================
# GROUPING TREE-RING VMD MODES BY CYCLE
# =========================================================

library(dplyr)

# -----------------------------
# 1. PATHS
# -----------------------------
base_folder <- "results/VMD_growth"
period_file <- file.path(base_folder, "tree_vmd_periods_rounded.txt")
vmd_folder <- base_folder

output_dir <- file.path(base_folder, "Grouped_cycles")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. TIME REFERENCE 
# -----------------------------
years_ref <- 1851:2001   # <<< MUST MATCH YOUR DATA

cat("\nTime range:", min(years_ref), "→", max(years_ref), "\n")

# -----------------------------
# 3. LOAD PERIOD TABLE
# -----------------------------
period_table <- read.table(
  period_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tree_names <- rownames(period_table)

# -----------------------------
# 4. STORAGE (cycle → matrix)
# -----------------------------
cycle_storage <- list()

# =========================================================
# 5. MAIN LOOP (TREES)
# =========================================================
for (tree in tree_names) {
  
  cat("\nProcessing tree:", tree, "\n")
  
  vmd_file <- file.path(vmd_folder, paste0("VMD_", tree, ".txt"))
  
  if (!file.exists(vmd_file)) {
    cat("⚠ Missing:", vmd_file, "\n")
    next
  }
  
  vmd_data <- read.table(vmd_file, header = TRUE, check.names = FALSE)
  
  # security
  if (!"Year" %in% names(vmd_data)) {
    cat("⚠ No Year column:", tree, "\n")
    next
  }
  
  # cycles for this tree-ring chronology
  
  cycles <- as.character(period_table[tree, ])
  cycles <- cycles[!is.na(cycles)]
  
  imf_index <- 1
  
  for (cyc in cycles) {
    
    col_name <- paste0("IMF_", imf_index)
    
    if (!(col_name %in% names(vmd_data))) {
      imf_index <- imf_index + 1
      next
    }
    
    values <- as.numeric(vmd_data[[col_name]])
    
    if (all(is.na(values))) {
      imf_index <- imf_index + 1
      next
    }
    
    # -----------------------------
    # INIT MATRIX FOR THIS CYCLE
    # -----------------------------
    if (!cyc %in% names(cycle_storage)) {
      
      mat <- matrix(NA,
                    nrow = length(years_ref),
                    ncol = length(tree_names))
      
      rownames(mat) <- years_ref
      colnames(mat) <- tree_names
      
      cycle_storage[[cyc]] <- mat
    }
    
    # -----------------------------
    # FILL MATRIX 
    # -----------------------------
    idx_years <- match(vmd_data$Year, years_ref)
    
    valid <- !is.na(idx_years) & !is.na(values)
    
    if (sum(valid) > 0) {
      cycle_storage[[cyc]][idx_years[valid], tree] <- values[valid]
    }
    
    imf_index <- imf_index + 1
  }
}

# =========================================================
# 6. EXPORT
# =========================================================
for (cyc in names(cycle_storage)) {
  
  out_file <- file.path(output_dir, paste0("cycle_", cyc, ".txt"))
  
  write.table(
    cycle_storage[[cyc]],
    file = out_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

cat("\n✔ Tree cycles grouped successfully\n")


# -----------------------------
# ROBUST MEAN CYCLE CREATION
# -----------------------------

# Using the dplR package function (standard in dendrochronology)
# or MASS for robust mean calculation
library(dplR)

# -----------------------------
# 1. Paths
# -----------------------------
input_dir  <- "results/VMD_growth/Grouped_cycles"
output_dir <- file.path(input_dir, "Mean_chronologies")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. List cycle files
# -----------------------------
files <- list.files(input_dir,
                    pattern = "^cycle_.*\\.txt$",
                    full.names = TRUE)

# -----------------------------
# 3. Loop over cycles
# -----------------------------
for (f in files) {
  
  cat("\nProcessing:", basename(f), "\n")
  
  # -----------------------------
  # Load data
  # -----------------------------
  df <- read.table(f,
                   header = TRUE,
                   check.names = FALSE)
  
  # -----------------------------
  # Convert to proper format for chron()
  # -----------------------------
  years <- as.numeric(rownames(df))
  
  # Ensure the "Year" column exists for safety
  
  if ("Year" %in% colnames(df)) {
    years <- df$Year
    df <- df[, -which(names(df) == "Year")]
  }
  
  rownames(df) <- years
  
  # -----------------------------
  # Convert to numeric
  # -----------------------------
  df[] <- lapply(df, as.numeric)
  
  # -----------------------------
  # Remove empty columns (trees lacking this cycle)
    # -----------------------------
  df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
  
  if (ncol(df) == 0) {
    cat("⚠ No valid series in:", basename(f), "\n")
    next
  }
  
  # -----------------------------
  # Compute chronology (robust mean)
  # -----------------------------
  crn <- chron(df,
               biweight = TRUE,
               prewhiten = FALSE)
  
  # -----------------------------
  # Extract outputs
  # -----------------------------
  result <- data.frame(
    Year = as.numeric(rownames(crn)),
    std  = crn[, "std"],   # moyenne robuste
    n    = crn[, "samp.depth"]  # nombre de séries
  )
  
  # -----------------------------
  # Output file name
  # -----------------------------
  cycle_name <- tools::file_path_sans_ext(basename(f))
  
  out_file <- file.path(
    output_dir,
    paste0(cycle_name, "_mean.txt")
  )
  
  # -----------------------------
  # Save
  # -----------------------------
  write.table(
    result,
    file = out_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("✔ Saved:", out_file, "\n")
}

cat("\n✔ All cycle chronologies computed successfully\n")
