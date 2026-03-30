library(tidyverse)
library(birdsize)
library(stringr)
library(patchwork)

# --- Load & clean ---
df_biomass <- read.csv("gbif_with_freshet_updated.csv")

data("known_species")
known_species <- known_species %>%
  mutate(scientific_name = paste(genus, species))
valid_species <- known_species$scientific_name

df_biomass$species <- trimws(df_biomass$species)
df_biomass$species <- word(df_biomass$species, 1, 2)
df_biomass$species <- paste(
  str_to_title(word(df_biomass$species, 1)),
  str_to_lower(word(df_biomass$species, 2))
)
df_biomass <- df_biomass %>% filter(species %in% valid_species)

# --- Abundance per species per year+month ---
annual_species_abundance <- df_biomass %>%
  group_by(year, month, species) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  left_join(
    df_biomass %>%
      group_by(year, month) %>%
      summarise(
        freshet_magnitude = mean(avg_monthly_max_discharge, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("year", "month")
  )

# --- Biomass simulation ---
annual_biomass <- annual_species_abundance %>%
  group_split(year, month) %>%
  map_df(function(month_df) {
    comm <- month_df %>%
      rename(scientific_name = species) %>%
      dplyr::select(scientific_name, abundance)
    sims <- community_generate(comm, abundance_column_name = "abundance")
    tibble(
      year              = unique(month_df$year),
      month             = unique(month_df$month),
      freshet_magnitude = unique(month_df$freshet_magnitude),
      total_abundance   = sum(month_df$abundance),
      total_biomass_g   = sum(sims$individual_mass)
    )
  })

annual_biomass_season <- annual_biomass %>%
  filter(month %in% c(4, 5, 6, 7, 8, 9))

model_biomass <- lm(log(total_biomass_g) ~ log(freshet_magnitude),
                    data = annual_biomass_season)

summary(model_biomass)
qqnorm(resid(model_biomass)); qqline(resid(model_biomass))
shapiro.test(resid(model_biomass))

# --- Annotation ---
r2    <- round(summary(model_biomass)$r.squared, 3)
pval  <- formatC(summary(model_biomass)$coefficients[2, 4], format = "e", digits = 2)
annotation <- paste0("R² = ", r2, "\np = ", pval)

year_colors <- scale_color_brewer(palette = "Dark2", name = "Year")

# --- Plot 1: Log-log ---
p1 <- ggplot(annual_biomass_season,
             aes(x = log(freshet_magnitude),
                 y = log(total_biomass_g),
                 color = factor(year))) +
  geom_smooth(method = "lm", se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_point(size = 3) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.4,
           label = annotation, size = 3.5) +
  year_colors +
  labs(
    x = "log(Mean Monthly Max Discharge)",
    y = "log(Total Monthly Biomass (g))",
    title = "Log-Log Scale"
  ) +
  theme_classic() +
  theme(legend.position = "none")

# --- Plot 2: Raw scale ---
pred_data <- data.frame(
  freshet_magnitude = seq(min(annual_biomass_season$freshet_magnitude),
                          max(annual_biomass_season$freshet_magnitude),
                          length.out = 200)
) %>%
  mutate(
    fit   = exp(predict(model_biomass, newdata = ., interval = "confidence")[, "fit"]),
    lower = exp(predict(model_biomass, newdata = ., interval = "confidence")[, "lwr"]),
    upper = exp(predict(model_biomass, newdata = ., interval = "confidence")[, "upr"])
  )

p2 <- ggplot() +
  geom_ribbon(data = pred_data,
              aes(x = freshet_magnitude, ymin = lower, ymax = upper),
              fill = "grey85") +
  geom_line(data = pred_data,
            aes(x = freshet_magnitude, y = fit),
            color = "grey40", linewidth = 0.8) +
  geom_point(data = annual_biomass_season,
             aes(x = freshet_magnitude,
                 y = total_biomass_g,
                 color = factor(year)),
             size = 3) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.4,
           label = annotation, size = 3.5) +
  year_colors +
  labs(
    x = "Mean Monthly Max Discharge (m³/s)",
    y = "Total Monthly Biomass (g)",
    title = "Raw Scale"
  ) +
  theme_classic() +
  theme(legend.position = "right")

# --- Combine ---
p1 + p2 +
  plot_annotation(
    title    = "Aerial Insectivore Biomass vs. Fraser River Freshet Magnitude",
    subtitle = "Breeding season (Apr–Sep), 2018–2023",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("freshet_biomass_plot.png", width = 12, height = 5.5, dpi = 300)