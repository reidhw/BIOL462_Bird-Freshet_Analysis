# ============================================================
# nMDS: Aerial Insectivore Community Composition x Freshet
# ============================================================
# Data:         gbif_with_freshet_updated.csv
# Response:     species abundance community matrix (year-level)
# Predictor:    avg_yearly_max_discharge (freshet peak, m3/s)
# Aggregation:  all observations collapsed to year-level vectors
#               (6 points, one per year 2018-2023)
# Tests:        PERMANOVA (adonis2) + envfit + Pearson correlation
# ============================================================

# ── Libraries ────────────────────────────────────────────────
library(tidyverse)
library(vegan)
library(conflicted)
library(ggrepel)

# Resolve common masking conflicts up front
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag",    "dplyr")

# ── 1. Load & tidy data ──────────────────────────────────────
dat <- read.csv("gbif_with_freshet_updated.csv")

# Remove rows with no species
dat <- dat %>%
  filter(!is.na(species), species != "")

# ── 2. Build year-level community matrix ─────────────────────
# Each row = one year
# Each column = total abundance of that species across all sites that year
comm_long <- dat %>%
  group_by(year, species) %>%
  summarise(abundance = n(), .groups = "drop")

comm_wide <- comm_long %>%
  pivot_wider(names_from  = species,
              values_from = abundance,
              values_fill = 0)

# Pure species matrix (rows = years, cols = species)
comm <- comm_wide %>%
  select(-year) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
rownames(comm) <- comm_wide$year

cat("Community matrix dimensions:", nrow(comm), "years,",
    ncol(comm), "species\n")

# ── 3. Build metadata table ───────────────────────────────────
# One freshet value per year (avg_yearly_max_discharge)
freshet_lut <- dat %>%
  group_by(year) %>%
  summarise(freshet_max = first(avg_yearly_max_discharge), .groups = "drop") %>%
  arrange(year)

cat("\nFreshet values by year:\n")
print(freshet_lut)

# ── 4. Hellinger transformation & Bray-Curtis distance ───────
# Hellinger transformation down-weights dominant species
# and is appropriate for zero-inflated abundance data
comm_hell <- decostand(comm, method = "hellinger")
bc_dist   <- vegdist(comm_hell, method = "bray")

cat("\nBray-Curtis dissimilarity matrix:\n")
print(round(as.matrix(bc_dist), 3))

# ── 5. nMDS ordination ───────────────────────────────────────
# With 6 points the solution is fast and stress is typically very low
set.seed(42)
nmds <- metaMDS(bc_dist,
                k             = 2,
                trymax        = 100,
                autotransform = FALSE)

cat("\nnMDS stress:", round(nmds$stress, 4), "\n")
# Stress guide: < 0.05 excellent, < 0.10 good, < 0.20 acceptable

# Stressplot to visually confirm ordination fit
stressplot(nmds, main = "Stressplot")

# ── 6. Extract scores ─────────────────────────────────────────
# Site (year) scores
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores <- site_scores %>%
  mutate(year        = as.numeric(rownames(site_scores)),
         freshet_max = freshet_lut$freshet_max[match(year, freshet_lut$year)])

# Species scores via weighted averages
# (required when passing a dist object to metaMDS)
species_scores <- as.data.frame(wascores(scores(nmds, display = "sites"),
                                         comm_hell))
colnames(species_scores) <- c("NMDS1", "NMDS2")
species_scores$species <- rownames(species_scores)

# ── 7. PERMANOVA (adonis2) ───────────────────────────────────
# Tests whether community composition differs by freshet and/or year
# NOTE: with only 6 observations there are only 6! = 720 possible
# permutations, so power is limited. A significant result is meaningful
# but non-significance does not rule out a real effect.

perm_meta <- data.frame(
  year        = freshet_lut$year,
  freshet_max = freshet_lut$freshet_max
)

set.seed(42)
perm_freshet <- adonis2(bc_dist ~ freshet_max,
                        data         = perm_meta,
                        permutations = 719,   # all possible permutations
                        method       = "bray")

set.seed(42)
perm_year <- adonis2(bc_dist ~ year,
                     data         = perm_meta,
                     permutations = 719,
                     method       = "bray")

cat("\n── PERMANOVA: Freshet Maximum ───────────────────────────\n")
print(perm_freshet)

cat("\n── PERMANOVA: Year ──────────────────────────────────────\n")
print(perm_year)

# ── 8. envfit: freshet vector ────────────────────────────────
# Tests whether the freshet gradient aligns with ordination axes
# More appropriate than PERMANOVA for continuous predictors with few points
ef <- envfit(nmds,
             env          = perm_meta["freshet_max"],
             permutations = 719)

cat("\n── envfit: freshet vector ───────────────────────────────\n")
print(ef)

# ── 9. Pearson correlation: freshet vs nMDS axes ──────────────
cat("\n── Pearson correlation: freshet vs nMDS axes ────────────\n")
cat("NMDS1:", round(cor(site_scores$NMDS1, site_scores$freshet_max), 3), "\n")
cat("NMDS2:", round(cor(site_scores$NMDS2, site_scores$freshet_max), 3), "\n")

# ── 10. Plot: nMDS without freshet arrow (main figure) ────────
year_labels <- setNames(
  paste0(freshet_lut$year, " (", freshet_lut$freshet_max, " m3/s)"),
  freshet_lut$year
)

p <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = factor(year)),
             size  = 5,
             alpha = 0.9) +
  geom_text_repel(aes(label  = year,
                      colour = factor(year)),
                  size        = 4,
                  fontface    = "bold",
                  show.legend = FALSE) +
  geom_text_repel(data        = species_scores,
                  aes(x       = NMDS1, y = NMDS2, label = species),
                  size        = 2.8,
                  colour      = "grey40",
                  fontface    = "italic",
                  inherit.aes = FALSE) +
  scale_colour_brewer(palette = "Dark2",
                      name    = "Year (freshet max)",
                      labels  = year_labels) +
  annotate("text",
           x     = min(site_scores$NMDS1),
           y     = max(site_scores$NMDS2),
           label = paste("Stress =", round(nmds$stress, 3)),
           hjust = 0, size = 3.5) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "right",
        legend.key.size  = unit(0.6, "cm"),
        axis.title       = element_text(size = 11),
        plot.title       = element_text(face = "bold")) +
  labs(title    = "nMDS - Aerial Insectivore Communities by Year",
       subtitle = NULL,
       x = "NMDS1", y = "NMDS2") +
  guides(size = "none")

ggsave("nmds_freshet_plot.png", plot = p,
       width = 9, height = 6.5, dpi = 300)

print(p)

# ── 11. Plot: nMDS with freshet arrow (supplementary) ─────────
# Arrow scaled by r² to honestly represent weak relationship
arrow_scale <- ef$vectors$r
arrow_df    <- as.data.frame(ef$vectors$arrows * arrow_scale * 0.3)
label_df    <- as.data.frame(ef$vectors$arrows * arrow_scale * 0.35)

p2 <- p +
  geom_segment(data        = arrow_df,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow       = arrow(length = unit(0.3, "cm")),
               colour      = "black",
               linewidth   = 0.9,
               inherit.aes = FALSE) +
  geom_text(data        = label_df,
            aes(x       = NMDS1, y = NMDS2, label = "Freshet max"),
            colour      = "black",
            size        = 3.5,
            inherit.aes = FALSE) +
  labs(title    = "nMDS - Aerial Insectivore Communities by Year",
       subtitle = "Freshet arrow scaled by r² (non-significant, p = 0.857)",
       x = "NMDS1", y = "NMDS2")

ggsave("nmds_freshet_vector_plot.png", plot = p2,
       width = 9, height = 6.5, dpi = 300)

print(p2)

