library(tidyverse)
library(vegan)

dat <- read.csv("gbif_with_discharge.csv")

dat <- dat %>%
  separate(geometry, into = c("x","y"), sep = ",", remove = FALSE) %>%
  mutate(x = as.numeric(gsub("c\\(|\\)", "", x)),
         y = as.numeric(gsub("\\)", "", y)),
         site = paste0(round(x, -3), "_", round(y, -3)))

#removing empty species rows
dat <- dat %>%
  filter(!is.na(species), species != "")
#Building the community matrix
comm_long <- dat %>%
  group_by(site, year, species) %>%
  summarise(abundance = n(), .groups = "drop")

comm_matrix <- comm_long %>%
  unite(site_year, site, year, remove = FALSE) %>%
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

comm_matrix <- comm_long %>%
  unite(site_year, site, year, remove = FALSE) %>%
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

# Drop the ID column explicitly
comm <- comm_matrix %>%
  select(-site_year) %>%   # only species
  as.matrix()

# Forcing numeric
comm <- apply(comm, 2, as.numeric)
rownames(comm) <- comm_matrix$site_year

#Helping to clean the data
bad_cols <- sapply(comm, function(x) any(is.na(as.numeric(x))))
bad_cols

names(bad_cols[bad_cols == TRUE])

str(comm)
summary(comm)
range(comm)

# Drop non-species columns
comm <- comm_matrix %>%
  select(-site_year, -site, -year) %>%   # remove metadata
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Put site-year back as rownames
rownames(comm) <- comm_matrix$site_year


meta <- comm_matrix %>%
  separate(site_year, into = c("site","year"), sep = "_")

comm <- as.matrix(comm_matrix[,-1])
rownames(comm) <- comm_matrix$site_year

str(comm)

#Standardizing
comm_hell <- decostand(comm, method = "hellinger")

#freshet variable
freshet <- dat %>%
  group_by(year) %>%
  summarise(freshet = mean(mean_discharge_max_year))

meta$year <- as.numeric(meta$year)
freshet$year <- as.numeric(freshet$year)

meta2 <- meta %>%
  left_join(freshet, by = "year")

bc_dist <- vegdist(comm_hell, method = "bray")

nmds <- metaMDS(comm_hell,
                distance = "bray",
                k = 2,
                trymax = 20)

nmds$stress

site_scores <- scores(nmds, display="sites")  # coordinates
species_scores <- scores(nmds, display="species")


# Extract years from rownames
years <- sub(".*_", "", rownames(site_scores))


#nmds plot
year_factor <- as.factor(years)

plot(site_scores[,1], site_scores[,2],
     col = year_factor,
     pch = 19,
     xlab = "NMDS1",
     ylab = "NMDS2")

legend("topright",
       legend = levels(year_factor),
       col = 1:length(levels(year_factor)),
       pch = 19)

# Add stress value
legend("topleft",
       legend = paste("Stress =", round(nmds$stress, 3)),
       bty = "n")
