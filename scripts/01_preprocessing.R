# =========================================
# 01 - DATA PREPROCESSING & PROJECT SETUP
# VMD Climate–Growth Analysis (Pinus halepensis)
# =========================================

# -------------------------
# Global project path
# -------------------------
base_dir <- "VMD_Climate_Growth_Project"

raw_dir     <- file.path(base_dir, "raw_data")
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")

# -------------------------
# Load libraries
# -------------------------
library(dplR)

# -------------------------
# Load raw tree-ring data
# -------------------------
rw_raw <- read.table(
  file.path(raw_dir, "PH_TUN_raw_ring_width.txt"),
  na.strings = "NA"
)

# -------------------------
# Detrending (Negative Exponential)
# -------------------------
detrend_out <- detrend(
  rw_raw,
  method = "ModNegExp",
  pos.slope = FALSE,
  constrain.nls = "when.fail",
  return.info = TRUE
)

rw_index <- as.data.frame(detrend_out$series)

# -------------------------
# Save processed data
# -------------------------
write.table(
  rw_index,
  file = file.path(data_dir, "PH_TUN_ring_width_index.txt"),
  row.names = FALSE
)

# -------------------------
# Regional chronology
# -------------------------
regional_chronology <- chron(rw_index)

write.table(
  regional_chronology,
  file = file.path(results_dir, "PH_TUN_regional_chronology.txt"),
  row.names = FALSE,
  sep = "\t"
)

# End of preprocessing script
