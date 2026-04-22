# =========================================================
# 05 - HEATMAP CLIMATE–GROWTH (FINAL VERSION)
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# =========================================================
# LOAD FINAL RESULTS (IMPORTANT)
# =========================================================

final_df <- read.csv("results/pcr_final_results.csv")

# =========================================================
# PREP DATA
# =========================================================

df <- final_df %>%
  mutate(Month = gsub(".*_", "", Term),
         Variable = gsub("_.*", "", Term))

df$Month <- factor(df$Month,
                   levels = c("Oct","Nov","Dec","Jan","Feb","Mar",
                              "Apr","May","Jun","Jul","Aug","Sep"))

# =========================================================
# HEATMAP
# =========================================================

lim <- max(abs(df$Coef), na.rm = TRUE)

p <- ggplot(df, aes(x = Month, y = Cycle, fill = Coef)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       midpoint = 0,
                       limit = c(-lim, lim)) +
  facet_wrap(~Variable, ncol = 1) +
  theme_minimal()

ggsave("results/figures/heatmap_final.png", p, width = 12, height = 10)

cat("\n✔ Figure saved\n")
