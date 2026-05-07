# =========================================================
# VARIANCE CONTRIBUTION PAR CYCLE (ARBRES)
# =========================================================

library(dplyr)

folder <- "results/VMD_growth/Grouped_cycles"

files <- list.files(folder,
                    pattern="^cycle_.*\\.txt$",
                    full.names=TRUE)

mean_variances <- c()
cycle_names <- c()

for(f in files){
  
  data <- read.table(f, header=TRUE, check.names=FALSE)
  
  # enlever colonne Year si présente
  if("Year" %in% colnames(data)){
    data <- data[, -which(names(data)=="Year")]
  }
  
  vars <- apply(data, 2, var, na.rm=TRUE)
  
  mean_variances <- c(mean_variances, mean(vars, na.rm=TRUE))
  
  cyc <- gsub("cycle_|\\.txt","", basename(f))
  cycle_names <- c(cycle_names, cyc)
}

# dataframe
results <- data.frame(
  Cycle = cycle_names,
  MeanVariance = mean_variances
)

# =========================
# TRI CORRECT (NUM + TREND)
# =========================
results$Cycle_Num <- suppressWarnings(as.numeric(results$Cycle))

results <- results %>%
  arrange(is.na(Cycle_Num), Cycle_Num)

results$Cycle[is.na(results$Cycle_Num)] <- "Trend"

# =========================
# CONTRIBUTION
# =========================
total_var <- sum(results$MeanVariance, na.rm=TRUE)

results$Contribution_percent <- (results$MeanVariance / total_var) * 100

print(results)

# sauvegarde
write.csv(results,
          file.path(folder, "Contribution_cycles_arbres.csv"),
          row.names=FALSE)

# =========================================================
# FREQUENCE DES CYCLES PAR CLASSE (ARBRES) - VERSION CORRIGÉE
# =========================================================

df <- read.table(
  "results/VMD_growth/tree_vmd_periods_rounded.txt",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# =========================
# FONCTION CLASSIFICATION
# =========================
classify_cycle <- function(x){
  
  if(is.na(x)) return(NA)
  
  # trend (texte)
  if(tolower(x) == "trend") return("Trend")
  
  x <- suppressWarnings(as.numeric(x))
  
  if(is.na(x)) return(NA)
  
  if(x >= 2 & x <= 5) return("2-5")
  if(x >= 6 & x <= 9) return("6-9")
  if(x >= 10 & x <= 15) return("10-15")
  if(x >= 16 & x <= 25) return("16-25")
  if(x > 25) return(">25")
  
  return(NA)
}

# =========================
# CLASSIFICATION MATRIX
# =========================
class_matrix <- apply(df, c(1,2), classify_cycle)

# =========================
# CLASSES FINAL ORDER
# =========================
classes <- c("2-5","6-9","10-15","16-25",">25","Trend")
cycles <- colnames(df)

# =========================
# FREQUENCY TABLE
# =========================
freq_table <- matrix(
  0,
  nrow = length(classes),
  ncol = length(cycles),
  dimnames = list(classes, cycles)
)

for(i in seq_along(cycles)){
  for(j in seq_along(classes)){
    freq_table[j,i] <- sum(class_matrix[,i] == classes[j], na.rm = TRUE)
  }
}

freq_table <- as.data.frame(freq_table)

print(freq_table)

# =========================
# EXPORT
# =========================
write.csv(
  freq_table,
  "results/VMD_growth/Frequence_classes_par_cycle.csv",
  row.names = TRUE
)

# =========================================================
# VARIANCE CONTRIBUTION CLIMAT PAR CYCLE
# =========================================================

library(dplyr)

calc_variance_cycles <- function(folder, varname){
  
  files <- list.files(folder,
                      pattern = paste0("^", varname, "_cycle_.*\\.txt$"),
                      full.names = TRUE)
  
  # ajouter trend
  trend_file <- file.path(folder, paste0(varname, "_trend.txt"))
  
  if(file.exists(trend_file)){
    files <- c(files, trend_file)
  }
  
  results <- data.frame()
  
  for(f in files){
    
    data <- read.table(f, header=TRUE, check.names=FALSE)
    
    # enlever colonne Year si présente
    if("Year" %in% colnames(data)){
      data <- data[, -which(names(data)=="Year")]
    }
    
    # variance globale du cycle (toutes colonnes/mois)
    vars <- apply(data, 2, var, na.rm=TRUE)
    mean_var <- mean(vars, na.rm=TRUE)
    
    cyc <- gsub(paste0(varname,"_|\\.txt"), "", basename(f))
    
    results <- rbind(results,
                     data.frame(Cycle = cyc,
                                MeanVariance = mean_var))
  }
  
  # tri
  results$Cycle_Num <- suppressWarnings(as.numeric(gsub("cycle_","",results$Cycle)))
  
  results <- results %>%
    arrange(is.na(Cycle_Num), Cycle_Num)
  
  results$Cycle[is.na(results$Cycle_Num)] <- "Trend"
  
  # contribution
  total <- sum(results$MeanVariance, na.rm=TRUE)
  
  results$Contribution_percent <- (results$MeanVariance / total) * 100
  
  return(results)
}

# =========================
# RUN
# =========================

temp_results <- calc_variance_cycles(
  "results/VMD_temperature_modes/Cycle_analysis/Merged_cycles",
  "temperature"
)

prec_results <- calc_variance_cycles(
  "results/VMD_precipitation_modes/Cycle_analysis/Merged_cycles",
  "precipitation"
)

print(temp_results)
print(prec_results)

# sauvegarde
write.csv(temp_results,
          "results/VMD_temperature_modes/Cycle_analysis/Merged_cycles/Contribution_temperature_cycles.csv",
          row.names=FALSE)

write.csv(prec_results,
          "results/VMD_precipitation_modes/Cycle_analysis/Merged_cycles/Contribution_precipitation_cycles.csv",
          row.names=FALSE)
