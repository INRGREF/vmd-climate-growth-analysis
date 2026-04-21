
# =========================================================
# VMD DECOMPOSITION OF MONTHLY CLIMATE VARIABLES (TEMPERATURE)
# Pinus halepensis climate-growth study
# =========================================================
# The same methodological framework is applied to precipitation data
# to ensure a consistent multi-variable climate decomposition approach
# across hydroclimatic drivers.
# =========================================================


library(VMDecomp)
library(dplyr)

# -----------------------------
# 1. Load data
# -----------------------------

data <- read.table(
  "raw_data/temperature.txt",
  header = TRUE
)

# Assume:
# Column 1 = Year
# Columns 2–13 = monthly values

temp_monthly <- data[, 2:13]

colnames(temp_monthly) <- c(
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

# -----------------------------
# 2. Function: optimal K selection
# -----------------------------

find_k_optimal <- function(series, max_k = 8, threshold = 0.06) {
  
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
    
    omegas <- sort(res$omega[nrow(res$omega), ])
    min_diff <- min(diff(omegas))
    
    if (min_diff < threshold) {
      return(k - 1)
    }
  }
  
  return(max_k)
}

# -----------------------------
# 3. VMD decomposition (month-by-month)
# -----------------------------

vmd_results <- list()

for (month in colnames(temp_monthly)) {
  
  cat("\nProcessing month:", month, "\n")
  
  series <- temp_monthly[[month]]
  
  # Optimal number of modes
  k_opt <- find_k_optimal(series)
  cat("Optimal K:", k_opt, "\n")
  
  # VMD decomposition
  vmd_out <- vmd(
    series,
    alpha = 2000,
    tau = 0,
    K = k_opt,
    DC = FALSE,
    init = 1,
    tol = 1e-7
  )
  
  # Store results
  vmd_results[[month]] <- list(
    modes = vmd_out$u,
    frequencies = vmd_out$omega[nrow(vmd_out$omega), ],
    K = k_opt
  )
}

# -----------------------------
# 4. Extract periods (1/f)
# -----------------------------

max_K <- max(sapply(vmd_results, function(x) x$K))

period_matrix <- matrix(
  NA,
  nrow = 12,
  ncol = max_K
)

rownames(period_matrix) <- names(vmd_results)
colnames(period_matrix) <- paste0("IMF_", 1:max_K)

for (month in names(vmd_results)) {
  
  freqs <- vmd_results[[month]]$frequencies
  periods <- 1 / freqs
  
  period_matrix[month, 1:length(periods)] <- sort(periods, decreasing = TRUE)
}

# Save summary table
write.table(
  round(period_matrix, 2),
  file = "results/temperature_vmd_periods_rounded.txt",
  sep = "\t",
  quote = FALSE
)

# -----------------------------
# 5. Export VMD modes
# -----------------------------

output_folder <- "results/VMD_temperature_modes"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

for (month in names(vmd_results)) {
  
  modes <- vmd_results[[month]]$modes
  
  # transpose check
  if (nrow(modes) != nrow(data)) {
    modes <- t(modes)
  }
  
  df_modes <- as.data.frame(modes)
  colnames(df_modes) <- paste0("IMF_", 1:ncol(df_modes))
  
  df_out <- cbind(Year = data[, 1], df_modes)
  
  write.csv(
    df_out,
    file = file.path(output_folder, paste0("VMD_", month, ".csv")),
    row.names = FALSE
  )
}

# -----------------------------
# 6. Quick visualization example
# -----------------------------

jan_modes <- vmd_results[["Jan"]]$modes

plot(
  data[,1],
  jan_modes[1, ],
  type = "l",
  main = "IMF1 - January Temperature",
  xlab = "Year",
  ylab = "Amplitude"
)


# =========================================================
# ORGANIZATION OF VMD MODES INTO CLIMATE CYCLES
# USING ROUNDED PERIOD TABLE
# Temperature dataset (Pinus halepensis study)
# =========================================================

library(dplyr)

# -----------------------------
# 1. Define paths
# -----------------------------

base_folder <- "results/VMD_temperature_modes"
output_folder <- file.path(base_folder, "Cycle_analysis")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# -----------------------------
# 2. Load rounded period table
# -----------------------------

period_table <- read.table(
  file.path(base_folder, "temperature_vmd_periods_rounded.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

months <- rownames(period_table)

# -----------------------------
# 3. Loop over months
# -----------------------------

for (month in months) {
  
  cat("\nProcessing month:", month, "\n")
  
  # Create output folder per month
  month_folder <- file.path(output_folder, paste0("Analysis_", month))
  if (!dir.exists(month_folder)) {
    dir.create(month_folder, recursive = TRUE)
  }
  
  # -----------------------------
  # 4. Load VMD file
  # -----------------------------
  
  vmd_file <- file.path(base_folder, paste0("VMD_", month, ".csv"))
  
  if (!file.exists(vmd_file)) {
    message("Missing VMD file for month: ", month)
    next
  }
  
  vmd_data <- read.csv(vmd_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Clean column names if needed
  names(vmd_data) <- gsub('"', "", names(vmd_data))
  
  # -----------------------------
  # 5. Extract cycles from rounded table
  # -----------------------------
  
  cycle_values <- as.character(period_table[month, ])
  cycle_values <- cycle_values[!is.na(cycle_values)]
  
  if (length(cycle_values) == 0) {
    message("No cycles for month: ", month)
    next
  }
  
  # -----------------------------
  # 6. Export each cycle
  # -----------------------------
  
  imf_index <- 1
  
  for (cycle in cycle_values) {
    
    col_name <- if (cycle == "trend") {
      "IMF_1"
    } else {
      paste0("IMF_", imf_index)
    }
    
    if (!(col_name %in% names(vmd_data))) {
      message("Missing column: ", col_name)
      imf_index <- imf_index + 1
      next
    }
    
    values <- vmd_data[[col_name]]
    
    if (all(is.na(values))) {
      imf_index <- imf_index + 1
      next
    }
    
    # -----------------------------
    # 7. Save cycle series
    # -----------------------------
    
    out_file <- file.path(
      month_folder,
      if (cycle == "trend") {
        paste0(month, "_trend.txt")
      } else {
        paste0(month, "_cycle_", cycle, "yr.txt")
      }
    )
    
    out_df <- data.frame(
      Year = vmd_data$Year,
      Value = values
    )
    
    write.table(
      out_df,
      file = out_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    imf_index <- imf_index + 1
  }
  
  message("Completed month: ", month)
}

cat("\n✔ Cycle organization completed using rounded periods.\n")


# =========================================================
# TREE-RING VMD DECOMPOSITION (MULTI-SERIES ALIGNMENT)
# Pinus halepensis growth analysis
# =========================================================

library(VMDecomp)

# -----------------------------
# 1. Load tree-ring index data
# -----------------------------

tree_data <- rw_index

# column 1 = Year
# other columns = individual tree chronologies

# -----------------------------
# 2. Output directory
# -----------------------------

output_dir <- "results/VMD_growth"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# 3. Optimal K selection function
# -----------------------------

find_k_optimal <- function(series, max_k = 7, threshold = 0.04) {
  
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
    
    omegas <- sort(res$omega[nrow(res$omega), ])
    min_diff <- min(diff(omegas))
    
    if (min_diff < threshold) {
      return(k - 1)
    }
  }
  
  return(max_k)
}

# -----------------------------
# 4. Storage object
# -----------------------------

vmd_tree_results <- list()

# -----------------------------
# 5. VMD decomposition per tree
# -----------------------------

for (i in 2:ncol(tree_data)) {
  
  tree_id <- colnames(tree_data)[i]
  cat("\nProcessing tree:", tree_id, "\n")
  
  raw_series <- tree_data[, i]
  years <- tree_data[, 1]
  
  # remove missing values
  valid_idx <- which(!is.na(raw_series))
  clean_series <- as.numeric(raw_series[valid_idx])
  
  if (length(clean_series) > 30) {
    
    # optimal number of modes
    k_optimal <- find_k_optimal(clean_series)
    cat("Optimal K:", k_optimal, "\n")
    
    # VMD decomposition
    vmd_out <- vmd(
      clean_series,
      alpha = 2000,
      tau = 0,
      K = k_optimal,
      DC = FALSE,
      init = 1,
      tol = 1e-7
    )
    
    # store results
    vmd_tree_results[[tree_id]] <- list(
      modes = vmd_out$u,
      frequencies = vmd_out$omega[nrow(vmd_out$omega), ],
      K = k_optimal
    )
    
    # -----------------------------
    # ALIGNMENT TO FULL TIME SERIES
    # -----------------------------
    
    aligned_matrix <- matrix(
      NA,
      nrow = length(years),
      ncol = k_optimal
    )
    
    modes <- vmd_out$u
    
    # transpose if needed
    if (nrow(modes) != length(clean_series)) {
      modes <- t(modes)
    }
    
    aligned_matrix[valid_idx, ] <- modes
    
    # final output table
    df_out <- cbind(
      Year = years,
      as.data.frame(aligned_matrix)
    )
    
    colnames(df_out) <- c("Year", paste0("IMF_", 1:k_optimal))
    
    # save per tree
    write.table(
      df_out,
      file = file.path(output_dir, paste0("VMD_", tree_id, ".txt")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    cat("Export completed for:", tree_id, "\n")
  }
}

# -----------------------------
# 6. Periods summary table
# -----------------------------

if (length(vmd_tree_results) > 0) {
  
  max_k <- max(sapply(vmd_tree_results, function(x) x$K))
  
  period_table <- matrix(
    NA,
    nrow = length(vmd_tree_results),
    ncol = max_k
  )
  
  rownames(period_table) <- names(vmd_tree_results)
  colnames(period_table) <- paste0("IMF_", 1:max_k)
  
  for (tree in names(vmd_tree_results)) {
    
    periods <- sort(1 / vmd_tree_results[[tree]]$frequencies, decreasing = TRUE)
    period_table[tree, 1:length(periods)] <- periods
  }
  
  write.table(
    round(period_table, 2),
    file = file.path(output_dir, "tree_vmd_periods_rounded.txt"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  cat("\n✔ Tree-ring VMD decomposition completed successfully\n")
}         


# =========================================================
# ORGANIZATION OF TREE-RING VMD MODES INTO GROWTH CYCLES
# USING ROUNDED PERIOD TABLE
# Pinus halepensis growth dataset
# =========================================================

library(dplR)

# -----------------------------
# 1. Paths (aligned with VMD workflow)
# -----------------------------

growth_dir <- "results/VMD_growth"

synthesis_dir <- file.path(growth_dir, "SYNTHESE_GROWTH_CYCLES")
chronology_dir <- file.path(growth_dir, "CHRONOLOGIES_GROWTH_CYCLES")

if (!dir.exists(synthesis_dir)) dir.create(synthesis_dir, recursive = TRUE)
if (!dir.exists(chronology_dir)) dir.create(chronology_dir, recursive = TRUE)

# -----------------------------
# 2. Load rounded period table
# -----------------------------

period_table <- read.table(
  file.path(growth_dir, "tree_vmd_periods_rounded.txt"),
  header = TRUE,
  row.names = 1,
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# -----------------------------
# 3. Identify cycles
# -----------------------------

raw_values <- as.vector(as.matrix(period_table))

numeric_cycles <- sort(
  unique(as.numeric(raw_values[raw_values != "trend"])),
  decreasing = TRUE
)

numeric_cycles <- na.omit(numeric_cycles)

all_targets <- c("trend", numeric_cycles)

# -----------------------------
# 4. Main loop (cycle-by-cycle grouping)
# -----------------------------

for (target in all_targets) {
  
  is_trend <- (target == "trend")
  
  file_label <- if (is_trend) {
    "Trend"
  } else {
    paste0("Cycle_", target, "_years")
  }
  
  synthesis_file <- paste0("Growth_Synthesis_", file_label, ".txt")
  chrono_file <- paste0("Growth_Chronology_", file_label, ".csv")
  
  cat("\n--- Processing:", target, "---\n")
  
  synthesis_df <- NULL
  trees_found <- 0
  
  # -----------------------------
  # 5. Loop over trees
  # -----------------------------
  
  for (tree_id in rownames(period_table)) {
    
    row_values <- as.character(period_table[tree_id, ])
    target_col <- which(row_values == target)
    
    if (length(target_col) > 0) {
      
      vmd_file <- file.path(growth_dir, paste0("VMD_", tree_id, ".txt"))
      
      if (file.exists(vmd_file)) {
        
        vmd_data <- read.table(vmd_file, header = TRUE, check.names = FALSE)
        
        col_index <- target_col[1] + 1  # shift for Year column
        
        if (col_index <= ncol(vmd_data)) {
          
          temp_df <- data.frame(
            Year = as.integer(vmd_data[, 1]),
            Value = as.numeric(vmd_data[, col_index])
          )
          
          colnames(temp_df) <- c("Year", tree_id)
          
          if (is.null(synthesis_df)) {
            synthesis_df <- temp_df
          } else {
            synthesis_df <- merge(synthesis_df, temp_df, by = "Year", all = TRUE)
          }
          
          trees_found <- trees_found + 1
        }
      }
    }
  }
  
  # -----------------------------
  # 6. Chronology computation (dplR)
  # -----------------------------
  
  if (!is.null(synthesis_df) && ncol(synthesis_df) > 1) {
    
    synthesis_df <- synthesis_df[order(synthesis_df$Year), ]
    
    rwl_data <- synthesis_df
    rownames(rwl_data) <- rwl_data$Year
    rwl_data$Year <- NULL
    
    # Robust chronology
    chrono <- chron(rwl_data, biweight = TRUE)
    
    # -----------------------------
    # A. Full synthesis output
    # -----------------------------
    
    synthesis_df$CHRONO_STD <- chrono$std
    
    write.table(
      synthesis_df,
      file = file.path(synthesis_dir, synthesis_file),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      na = "NA"
    )
    
    # -----------------------------
    # B. Chronology output
    # -----------------------------
    
    chrono_df <- data.frame(
      Year = as.numeric(rownames(chrono)),
      std = chrono$std,
      sample_depth = chrono$samp.depth
    )
    
    write.csv(
      chrono_df,
      file = file.path(chronology_dir, chrono_file),
      row.names = FALSE
    )
    
    cat(" OK - Trees processed:", trees_found, "\n")
    
  } else {
    cat(" NO DATA FOUND FOR:", target, "\n")
  }
}                      
