# =========================================================
# 05_figure_climate_growth_heatmap.R
# Visualization of multi-scale climate–growth relationships
# PCR bootstrap results (Pinus halepensis, Tunisia)
# =========================================================

# ---------------------------------------------------------
# 1. Load libraries
# ---------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------------------------
# 2. Define output path (project structure)
# ---------------------------------------------------------

figure_dir <- "results/figures"

if(!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# ---------------------------------------------------------
# 3. Merge PCR results (multi + single variable cycles)
# ---------------------------------------------------------

all_list <- c(final_results_list, results_solo)
valid_list <- Filter(is.data.frame, all_list)

if(length(valid_list) > 0) {
  
  final_df_all <- bind_rows(valid_list)
  
  # Significance based on bootstrap confidence intervals
  if(all(c("Inf95", "Sup95") %in% colnames(final_df_all))) {
    final_df_all <- final_df_all %>%
      mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", ""))
  } else {
    final_df_all$Signif <- ""
    warning("Missing Inf95 or Sup95 in some results.")
  }
  
} else {
  stop("No valid PCR results found.")
}

# ---------------------------------------------------------
# 4. Data preparation for heatmap visualization
# ---------------------------------------------------------

df_graph <- final_df_all %>%
  
  # Climate variable identification (Temperature vs Precipitation)
  mutate(Variable_Type = case_when(
    grepl("(^P_|_P$|_P_)", Term) ~ "Precipitation",
    grepl("(^T_|_T$|_T_)", Term) ~ "Mean Temperature",
    TRUE ~ "Climate Signal"
  )) %>%
  
  # Extract month (robust to naming variations)
  mutate(Month = case_when(
    grepl("Oct", Term, ignore.case = TRUE) ~ "Oct", 
    grepl("Nov", Term, ignore.case = TRUE) ~ "Nov", 
    grepl("Dec", Term, ignore.case = TRUE) ~ "Dec",
    grepl("Jan", Term, ignore.case = TRUE) ~ "Jan", 
    grepl("Feb|Fev", Term, ignore.case = TRUE) ~ "Feb", 
    grepl("Mar", Term, ignore.case = TRUE) ~ "Mar",
    grepl("Apr|Avr", Term, ignore.case = TRUE) ~ "Apr", 
    grepl("May|Mai", Term, ignore.case = TRUE) ~ "May", 
    grepl("Jun", Term, ignore.case = TRUE) ~ "Jun",
    grepl("Jul", Term, ignore.case = TRUE) ~ "Jul", 
    grepl("Aug|Aou", Term, ignore.case = TRUE) ~ "Aug", 
    grepl("Sep", Term, ignore.case = TRUE) ~ "Sep"
  )) %>%
  
  # Biological year ordering (Oct → Sep)
  mutate(Month = factor(Month, levels = c(
    "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", 
    "Apr", "May", "Jun", "Jul", "Aug", "Sep"
  ))) %>%
  #######################################################################
# NOTE: If a different definition of the biological year is required,
# modify the ordering of the 'Month' factor accordingly. This step controls
# the seasonal alignment of climate variables but does not affect the
# underlying PCR or VMD analyses, ensuring methodological flexibility.
##############################################################################
  # Extract numeric cycle (for ordering)
  mutate(Cycle_Num = as.numeric(gsub("[^0-9]", "", Cycle))) %>%
  mutate(Cycle_Sort = ifelse(is.na(Cycle_Num), 999, Cycle_Num)) %>%
  mutate(Cycle_Label = ifelse(Cycle_Sort == 999, "Trend", as.character(Cycle_Num))) %>%
  mutate(Cycle_Label = reorder(Cycle_Label, Cycle_Sort)) %>%
  
  drop_na(Month)

# ---------------------------------------------------------
# 5. Heatmap visualization
# ---------------------------------------------------------

limit_val <- max(abs(df_graph$Coef), na.rm = TRUE) * 0.6 

gg <- ggplot(df_graph, aes(x = Month, y = Cycle_Label, fill = Coef)) +
  
  geom_tile(color = "white", linewidth = 0.3) +
  
  # Significance markers
  geom_text(aes(label = ifelse(Signif == "*" & abs(Coef) > 0.001, "*", "")), 
            color = "black", size = 6, vjust = 0.75) +
  
  # Diverging color scale
  scale_fill_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = 0,
    limit = c(-limit_val, limit_val),
    oob = scales::squish,
    name = "Reg. Coeff."
  ) +
  
  # Facets by climate variable
  facet_wrap(~Variable_Type, ncol = 1, scales = "free_y") +
  
  # Clean theme
  theme_minimal() +
  labs(
    x = "Month (Biological Year)",
    y = "Periodicity (years)"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", color = "black", size = 13),
    axis.title = element_text(face = "bold", size = 13),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display
print(gg)

# ---------------------------------------------------------
# 6. Export figure
# ---------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "Figure_climate_growth_heatmap.png"),
  plot = gg,
  width = 14,
  height = 12,
  dpi = 300
)

cat("\nFigure successfully saved in:", figure_dir, "\n")
