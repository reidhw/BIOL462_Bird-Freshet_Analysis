library(tidyverse)

dat <- read.csv("gbif_with_freshet_updated.csv")

# ── Panel A: Annual freshet magnitude ─────────────────────────────────────────

panel_a_data <- dat %>%
  group_by(year) %>%
  summarise(annual_discharge = mean(avg_yearly_max_discharge, na.rm = TRUE),
            .groups = "drop")

p_annual <- ggplot(panel_a_data, aes(x = factor(year), y = annual_discharge)) +
  geom_col(fill = "#2c7bb6", width = 0.65) +
  labs(
    x     = "Year",
    y     = expression("Mean annual maximum discharge (m"^3*"/s)"),
    title = "A"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank()
  )

# ── Panel B: Mean monthly max discharge by year ───────────────────────────────

month_labels <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")

panel_b_data <- dat %>%
  group_by(year, month) %>%
  summarise(avg_monthly_max = mean(avg_monthly_max_discharge, na.rm = TRUE),  # column confirmed in updated CSV
            .groups = "drop") %>%
  mutate(
    month = factor(month, levels = 1:12, labels = month_labels),
    year  = factor(year)
  )

# Color palette — 6 muted, colorblind-friendly colors
year_colors <- c("#2c7bb6", "#abd9e9", "#74c476", "#fdae61", "#d7191c", "#8856a7")

p_monthly <- ggplot(panel_b_data,
                    aes(x = month, y = avg_monthly_max,
                        fill = year)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = year_colors, name = "Year") +
  labs(
    x     = "Month",
    y     = expression("Mean monthly maximum discharge (m"^3*"/s)"),
    title = "B"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

# ── Combine and save ──────────────────────────────────────────────────────────

library(patchwork)

combined <- p_annual / p_monthly +
  plot_annotation(
    caption = paste(
      "Fig. 1. Freshet characterization across the study period (2018–2023).",
      "(A) Mean annual maximum discharge per year.",
      "(B) Mean monthly maximum discharge by year.",
      "Discharge data from [your source]."
    )
  ) &
  theme(plot.caption = element_text(hjust = 0, size = 9))

ggsave("Fig1_freshet_characterization.png",
       plot   = combined,
       width  = 8,
       height = 7,
       dpi    = 300)

print("Saved: Fig1_freshet_characterization.png")
