#=========================================================
# VMD-PCR CLIMATE-GROWTH HEATMAP 
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ─────────────────────────────────────────────────────────
# 1. DATA PREPARATION
# ─────────────────────────────────────────────────────────
df_graph <- final_df_pcr %>%
  mutate(Cycle_Num = as.numeric(gsub("[^0-9]", "", Cycle))) %>%
  mutate(Variable_Type = case_when(
    grepl("_P$|_P_", Terme) ~ "Precipitation",
    grepl("_T$|_T_", Terme) ~ "Mean Temperature",
    TRUE ~ "Climate Signal"
  )) %>%
  mutate(Variable_Type = factor(Variable_Type,
                                levels = c("Climate Signal",
                                           "Mean Temperature",
                                           "Precipitation")))

# ─────────────────────────────────────────────────────────
# 2. MONTH MANAGEMENT (BIOLOGICAL YEAR: OCT → SEP)
# NOTE: The biological year window is flexible and can be adjusted 
# to suit specific research objectives or species-specific phenology 
# (e.g., shifting the start/end months to match regional growth seasons).
# ─────────────────────────────────────────────────────────
month_levels <- c("Oct","Nov","Dec","Jan","Feb","Mar",
                  "Apr","May","Jun","Jul","Aug","Sep")

# Full month labels for X-axis 
month_labels <- c("Oct","Nov","Dec","Jan","Feb","Mar",
                  "Apr","May","Jun","Jul","Aug","Sep")

df_graph <- df_graph %>%
  mutate(Mois_Short = case_when(
    grepl("Oct",     Terme, ignore.case = TRUE) ~ "Oct",
    grepl("Nov",     Terme, ignore.case = TRUE) ~ "Nov",
    grepl("Dec",     Terme, ignore.case = TRUE) ~ "Dec",
    grepl("Jan",     Terme, ignore.case = TRUE) ~ "Jan",
    grepl("Feb|Fev", Terme, ignore.case = TRUE) ~ "Feb",
    grepl("Mar",     Terme, ignore.case = TRUE) ~ "Mar",
    grepl("Apr|Avr", Terme, ignore.case = TRUE) ~ "Apr",
    grepl("May|Mai", Terme, ignore.case = TRUE) ~ "May",
    grepl("Jun",     Terme, ignore.case = TRUE) ~ "Jun",
    grepl("Jul",     Terme, ignore.case = TRUE) ~ "Jul",
    grepl("Aug|Aou", Terme, ignore.case = TRUE) ~ "Aug",
    grepl("Sep",     Terme, ignore.case = TRUE) ~ "Sep"
  )) %>%
  mutate(Mois_Short = factor(Mois_Short, levels = month_levels)) %>%
  drop_na(Mois_Short)

# ─────────────────────────────────────────────────────────
# 3. Y-AXIS ORDERING (ascending periodicity, Trend on top)
# ─────────────────────────────────────────────────────────
df_graph <- df_graph %>%
  mutate(
    Cycle_Sort  = ifelse(is.na(Cycle_Num), 999, Cycle_Num),
    Cycle_Label = ifelse(Cycle_Sort == 999, "Trend", as.character(Cycle_Sort))
  )

ordres_cycles <- df_graph %>%
  select(Cycle_Label, Cycle_Sort) %>%
  distinct() %>%
  arrange(Cycle_Sort) %>%
  pull(Cycle_Label)

df_graph$Cycle_Label <- factor(df_graph$Cycle_Label, levels = ordres_cycles)

# ─────────────────────────────────────────────────────────
# 4. COLOR SCALE & SIGNIFICANCE
# ─────────────────────────────────────────────────────────

# Symmetric limit at 95th percentile 
# while preserving gradient contrast across the bulk of the data
limit_val <- quantile(abs(df_graph$Coef), 0.95, na.rm = TRUE)

# Palette: RdBu-derived, slightly muted for print legibility
col_neg  <- "#b2182b"   # warm red  — negative coefficients
col_zero <- "#f5f5f5"   # off-white — near-zero / neutral
col_pos  <- "#2166ac"   # cool blue — positive coefficients

# Significant tiles: we add a small filled circle 

df_graph <- df_graph %>%
  mutate(Signif_bool = (Signif == "*" & abs(Coef) > 0.001))

# ─────────────────────────────────────────────────────────
# 5. FACET STRIP LABELS

# ─────────────────────────────────────────────────────────
strip_labels <- c(
  "Climate Signal"    = "Climate signal",
  "Mean Temperature"  = "Mean temperature",
  "Precipitation"     = "Precipitation"
)

# ─────────────────────────────────────────────────────────
# 6. PLOT
# ─────────────────────────────────────────────────────────
gg <- ggplot(df_graph,
             aes(x    = Mois_Short,
                 y    = Cycle_Label,
                 fill = Coef)) +
  
  # — Tiles —
  geom_tile(color     = "white",
            linewidth  = 0.25) +       # thinner separator 
  
  # — Significance marker 
  geom_point(
    data  = filter(df_graph, Signif_bool),
    aes(x = Mois_Short, y = Cycle_Label),
    inherit.aes = FALSE,
    shape  = 16,           
    size   = 1.4,
    colour = "black"
  ) +
  
  # — Colour scale —
  scale_fill_gradient2(
    low      = col_neg,
    mid      = col_zero,
    high     = col_pos,
    midpoint = 0,
    limits   = c(-limit_val, limit_val),
    oob      = squish,
    name     = "Regression\ncoefficient",
    guide    = guide_colorbar(
      title.position  = "top",
      title.hjust     = 0.5,
      barwidth        = unit(0.4, "cm"),
      barheight       = unit(4.5, "cm"),
      ticks.linewidth = 0.5,
      frame.colour    = "grey50",
      frame.linewidth = 0.4
    )
  ) +
  
  # — Facets —
  facet_wrap(
    ~ Variable_Type,
    ncol    = 1,
    scales  = "free_y",
    labeller = as_labeller(strip_labels)
  ) +
  
  # — Axis scales —
  # PAR :
  scale_x_discrete(
    labels = c("Oct" = "Oct(t-1)", "Nov" = "Nov(t-1)", "Dec" = "Dec(t-1)",
               "Jan" = "Jan",      "Feb" = "Feb",      "Mar" = "Mar",
               "Apr" = "Apr",      "May" = "May",      "Jun" = "Jun",
               "Jul" = "Jul",      "Aug" = "Aug",      "Sep" = "Sep"),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  
  # — Labels —
  labs(
    x       = "Month (biological year)",
    y       = "Periodicity (years)",
    caption = ""
  ) +
  
  # — Theme —
  theme_minimal(base_size = 10, base_family = "serif") +
  
  theme(
    # Facet strips
    strip.background   = element_rect(fill = "grey94", colour = "grey70", linewidth = 0.4),
    strip.text         = element_text(face = "bold", size = 12,
                                      margin = margin(3, 0, 3, 0)),
    
    # Grid
    panel.grid         = element_blank(),
    panel.border       = element_rect(fill = NA, colour = "grey70", linewidth = 0.4),
    
    # Axes
    axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x       = element_text(size = 12, face = "bold",
                                      margin = margin(t = 5)),
    axis.title.y       = element_text(size = 12, face = "bold",
                                      margin = margin(r = 5)),
    axis.ticks         = element_line(colour = "grey60", linewidth = 0.3),
    axis.ticks.length  = unit(2, "pt"),
    
    # Legend
    legend.position    = "right",
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 8),
    legend.margin      = margin(0, 0, 0, 4),
    
    # Panel spacing
    panel.spacing      = unit(0.5, "lines"),
    
    # Caption
    plot.caption       = element_text(size = 6.5, colour = "grey40",
                                      hjust = 0, margin = margin(t = 6)),
    
    # Overall margins
    plot.margin        = margin(6, 6, 6, 6)
  )

print(gg)

# ─────────────────────────────────────────────────────────
# 7. EXPORT
#
# Dimensions follow standard journal column widths:
#   Single column : 8.5 cm
#   Double column : 17.5 cm  
#
# Resolution:
#   PDF (vector) 
#   PNG 600 dpi  
#   TIFF 600 dpi 
# ─────────────────────────────────────────────────────────
out_dir   <- "results/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
base_name <- file.path(out_dir, "PCR_heatmap_Q1")

# Vector PDF 
ggsave(
  filename = paste0(base_name, ".pdf"),
  plot     = gg,
  width    = 17.5,
  height   = 16,
  units    = "cm",
  device   = cairo_pdf
)

# PNG 600 dpi 
ggsave(
  filename = paste0(base_name, ".png"),
  plot     = gg,
  width    = 22,
  height   = 15,
  units    = "cm",
  dpi      = 600
)

# TIFF 600 dpi 
# ggsave(
#   filename    = paste0(base_name, ".tiff"),
#   plot        = gg,
#   width       = 17.5,
#   height      = 16,
#   units       = "cm",
#   dpi         = 600,
#   compression = "lzw"
# )

cat("\n✔ Files exported to:", file.path(getwd(), out_dir), "\n")
