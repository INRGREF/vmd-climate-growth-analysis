# =========================================
# 01 - DATA PREPROCESSING & PROJECT SETUP
# VMD Climate–Growth Analysis (Pinus halepensis)
# =========================================

# -------------------------
# Project structure
# -------------------------
# Standard structure:
#
# results/
# ├── VMD_growth/
# ├── VMD_temperature_modes/
# ├── VMD_precipitation_modes/
# ├── figures/
# └── data/

# =========================================
# 01 - DATA PREPROCESSING FIXED 
# =========================================

library(dplR)

# -------------------------
# Load raw data
# -------------------------

# Replace the filename below with your current dataset name
rw_raw <- read.table(
  "raw_data/YOUR_NEW_FILENAME.txt",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "")
)

# -------------------------
# Extract Year
# -------------------------

years <- as.numeric(rw_raw[, 1])

# -------------------------
# Extract tree data
# -------------------------

rw_values <- rw_raw[, -1]

# -------------------------
# CLEAN STEP 
# Remove trees with too many NA
# -------------------------

valid_cols <- colSums(!is.na(rw_values)) > 10   # adjustable threshold

rw_values <- rw_values[, valid_cols]

# -------------------------
# OPTIONAL: remove years with ALL NA
# -------------------------

valid_rows <- rowSums(!is.na(rw_values)) > 0

years <- years[valid_rows]
rw_values <- rw_values[valid_rows, ]

# -------------------------
# Detrending
# -------------------------

detrend_out <- detrend(
  rw_values,
  method = "ModNegExp",
  pos.slope = FALSE,
  constrain.nls = "when.fail",
  return.info = TRUE
)

rw_index <- as.data.frame(detrend_out$series)

# -------------------------
# ADD YEAR BACK 
# -------------------------

rw_index <- cbind(Year = years[1:nrow(rw_index)], rw_index)

# -------------------------
# Save processed data
# -------------------------

write.table(
  rw_index,
  file = "results/data/PH_TUN_ring_width_index.txt",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# -------------------------
# Regional chronology
# -------------------------

regional_chronology <- chron(rw_index[, -1])

write.table(
  regional_chronology,
  file = "results/data/PH_TUN_regional_chronology.txt",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

cat("\n✔ Preprocessing fixed: stable Year + cleaned series\n")

