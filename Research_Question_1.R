library(tidyverse)
library(MASS) 

data <- read.csv("gbif_with_monthly_discharge_updated.csv")

# --- Filter for Bank Swallow ---
bank_swallow <- data %>%
  filter(species == "Riparia riparia")

# --- Monthly abundance with correct discharge column ---
monthly_abundance <- bank_swallow %>%
  group_by(year, month) %>%
  summarise(
    abundance  = n(),
    discharge  = mean(avg_monthly_max_discharge, na.rm = TRUE),  # FIX 1
    .groups = "drop"
  )

# --- Breeding season filter ---
breeding_data <- monthly_abundance %>%
  filter(month %in% c(4, 5, 6, 7, 8, 9))

# --- Negative binomial GLM ---
glm_model <- glm.nb(abundance ~ discharge, data = breeding_data)
summary(glm_model)

# --- check overdispersion ---
deviance(glm_model) / df.residual(glm_model)

# --- QQ plot ---
qqnorm(resid(glm_model))
qqline(resid(glm_model))

# --- Plot 1: Single line ---
grey_abundance <- ggplot(breeding_data, aes(x = discharge, y = abundance)) +
  geom_point(size = 3) +
  geom_smooth(method = "glm.nb",                    
              se = TRUE, color = "grey40") +
  labs(
    x = "Mean Monthly Max Discharge (m³/s)",
    y = "Bank Swallow Abundance",
    title = "Bank Swallow Abundance vs Peak Freshet Discharge",
    subtitle = "Breeding season (Apr–Sept), 2018–2023"
  ) +
  theme_classic()

# --- Plot 2: Coloured by month ---
colour_abundance <- ggplot(breeding_data, aes(x = discharge, y = abundance, color = factor(month))) +
  geom_point(size = 3) +
  geom_smooth(method = "glm.nb", se = FALSE) +
  scale_color_brewer(palette = "Dark2",
                     name = "Month",
                     labels = c("4" = "Apr", "5" = "May", "6" = "Jun",
                                "7" = "Jul", "8" = "Aug", "9" = "Sept")) +
  labs(
    x = "Mean Monthly Max Discharge (m³/s)",
    y = "Bank Swallow Abundance",
    title = "Bank Swallow Abundance vs Peak Freshet Discharge",
    subtitle = "Breeding season (Apr–Sept), 2018–2023"
  ) +
  theme_classic()

grey_abundance + colour_abundance +
  plot_annotation(
    title = "Bank Swallow Abundance vs. Fraser River Freshet Magnitude",
    subtitle = "Breeding season (Apr–Sep), 2018–2023",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )
