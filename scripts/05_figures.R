# =========================================================
# 05 - CLIMATE–GROWTH HEATMAP (FINAL CORRIGÉ)
# Visualization of PCR bootstrap results
# Pinus halepensis (Tunisia)
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ---------------------------------------------------------

figure_dir <- "results/figures"
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# ---------------------------------------------------------
# 2. CHECK INPUT DATA (from Script 4)
# ---------------------------------------------------------

if (!exists("final_results_list") || !exists("results_solo")) {
  stop("Run Script 4 first (PCR results missing).")
}

all_list <- c(final_results_list, results_solo)
valid_list <- Filter(is.data.frame, all_list)

if (length(valid_list) == 0) {
  stop("No PCR results available.")
}

final_df_all <- bind_rows(valid_list)

# ---------------------------------------------------------
# 3. SIGNIFICANCE
# ---------------------------------------------------------

final_df_all <- final_df_all %>%
  mutate(Signif = ifelse(Inf95 * Sup95 > 0, "*", ""))

# ---------------------------------------------------------
# 4. PREPARE DATA
# ---------------------------------------------------------

df_graph <- final_df_all %>%
  
  mutate(
    Variable_Type = case_when(
      grepl("(^P_|_P_)", Term) ~ "Precipitation",
      grepl("(^T_|_T_)", Term) ~ "Temperature",
      TRUE ~ "Climate Signal"
    )
  ) %>%
  
  mutate(
    Month = case_when(
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
      grepl("Sep", Term, ignore.case = TRUE) ~ "Sep",
      TRUE ~ NA_character_
    )
  ) %>%
  
  mutate(Month = factor(Month,
                        levels = c("Oct","Nov","Dec","Jan","Feb","Mar",
                                   "Apr","May","Jun","Jul","Aug","Sep"))) %>%
  
  mutate(
    Cycle_Num = as.numeric(gsub("[^0-9]", "", Cycle)),
    Cycle_Sort = ifelse(is.na(Cycle_Num), 999, Cycle_Num),
    Cycle_Label = ifelse(Cycle_Sort == 999, "Trend", as.character(Cycle_Num))
  ) %>%
  
  mutate(Cycle_Label = reorder(Cycle_Label, Cycle_Sort)) %>%
  drop_na(Month)

# ---------------------------------------------------------
# 5. SAFETY CHECK
# ---------------------------------------------------------

if (nrow(df_graph) == 0) {
  stop("No valid data for heatmap.")
}

# ---------------------------------------------------------
# 6. HEATMAP
# ---------------------------------------------------------

limit_val <- max(abs(df_graph$Coef), na.rm = TRUE) * 0.6

gg <- ggplot(df_graph, aes(x = Month, y = Cycle_Label, fill = Coef)) +
  
  geom_tile(color = "white", linewidth = 0.3) +
  
  geom_text(aes(label = ifelse(Signif == "*" & abs(Coef) > 0.001, "*", "")),
            size = 5, vjust = 0.7) +
  
  scale_fill_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = 0,
    limit = c(-limit_val, limit_val),
    oob = scales::squish,
    name = "Coefficient"
  ) +
  
  facet_wrap(~Variable_Type, ncol = 1, scales = "free_y") +
  
  theme_minimal() +
  
  labs(
    x = "Biological Year",
    y = "Cycle Periodicity (years)"
  ) +
  
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

print(gg)

# ---------------------------------------------------------
# 7. SAVE FIGURE
# ---------------------------------------------------------

ggsave(
  filename = file.path(figure_dir,
                       "Figure_climate_growth_heatmap.png"),
  plot = gg,
  width = 14,
  height = 12,
  dpi = 300
)

cat("\n✔ Heatmap saved successfully\n")
