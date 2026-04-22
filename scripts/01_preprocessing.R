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

# -------------------------
# Load libraries
# -------------------------

library(dplR)

# -------------------------
# Load raw tree-ring data
# -------------------------

rw_raw <- read.table(
  "raw_data/PH_TUN_raw_ring_width.txt",
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
  file = "data/PH_TUN_ring_width_index.txt",
  row.names = FALSE
)

# -------------------------
# Regional chronology
# -------------------------

regional_chronology <- chron(rw_index)

write.table(
  regional_chronology,
  file = "results/PH_TUN_regional_chronology.txt",
  row.names = FALSE,
  sep = "\t"
)

# End of preprocessing script
