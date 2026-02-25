# Leveraging citizen-science to estimate the population size of non-breeding Surfbirds (Calidris virgata) in Chile
# Authors: Erik M. Sandvig, Gabriela Contreras, Franco Villalobos, Marion Díaz, Giannira Álvarez, Carlos Silva & Fernando Medrano
# Date: February 2026

# Load required packages
library(tidyverse)
library(lubridate)
library(zoo)
library(scales)
library(sf)
library(brms)
library(glmmTMB)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)

# DATA IMPORT AND PREPARATION
# Import standardized survey data
CCAP <- read_csv("eBirdDataAbr25.csv")

# Import survey transect metadata
listas <- read.csv("ListasCCAP2025.csv", sep = ";")
names(listas)[names(listas) == "Submission.ID"] <- "Submission ID"

# Filter for shorebirds (Charadriiformes subset)
playeras <- CCAP %>%
  filter(`Taxonomic Order` >= 5713, `Taxonomic Order` <= 6171)

# Filter for 2025 surveys
playeras25 <- playeras %>%
  filter(Date >= as.Date("2025-01-01"))

# Filter for Surfbird only
plr <- playeras25 %>%
  filter(`Common Name` == "Playero de las rompientes" | `Taxonomic Order` == 6112)

# JOIN SURVEY DATA WITH EBIRD RECORDS
# Join transect metadata with Surfbird observations
# Also pull Date/Time from full CCAP for checklists without Surfbird detections
df <- left_join(listas, plr, by = "Submission ID") %>%
  left_join(
    CCAP %>%
      distinct(`Submission ID`, .keep_all = TRUE) %>%
      dplyr::select(`Submission ID`, Date_ccap = Date, Time_ccap = Time,
                    Duration_ccap = `Duration (Min)`, N_obs_ccap = `Number of Observers`),
    by = "Submission ID"
  ) %>%
  mutate(
    Date = coalesce(Date, Date_ccap),
    Time = coalesce(Time, Time_ccap),
    `Duration (Min)` = coalesce(`Duration (Min)`, Duration_ccap),
    `Number of Observers` = coalesce(`Number of Observers`, N_obs_ccap),
    dt = ymd_hms(paste(Date, Time))
  ) %>%
  dplyr::select(-Date_ccap, -Time_ccap, -Duration_ccap, -N_obs_ccap)

# ADD TIDE HEIGHT (INTERPOLATED)
# Import tide data
mareas <- read.csv2("TablaMareasRocosos.csv", colClasses = c("hora" = "character"))
mareas <- mareas %>%
  mutate(fecha_hora = dmy_hms(paste0(fecha, " ", str_sub(hora, 1, 2), ":",
                                      str_sub(hora, 3, 4), ":00"))) %>%
  arrange(lista_id, fecha_hora)

# Interpolate tide height at survey times
new <- dplyr::select(df, fecha_hora = dt)
new$altura <- NA
tides_sub <- dplyr::select(mareas, fecha_hora, altura)
new <- rbind(new, tides_sub) %>% arrange(fecha_hora)
new$pred_altura_marea <- na.approx(new$altura, x = new$fecha_hora, na.rm = FALSE)

# Merge interpolated tide height back to main data
df <- merge(df, dplyr::select(new, fecha_hora, pred_altura_marea),
            by.x = "dt", by.y = "fecha_hora", all.x = TRUE)

# CREATE ZERO-FILLED DATASET FOR MODELING
surfbird <- c("Playero de las rompientes")
all_listas <- unique(listas$`Submission ID`)

# Checklist-level metadata
checklist_info <- df %>%
  distinct(`Submission ID`, .keep_all = TRUE) %>%
  dplyr::select(`Submission ID`, Transecto, Region, kms, Time,
                `Duration (Min)`, pred_altura_marea, Lat, `Number of Observers`)

# Expand grid for all checklist × species combinations
combinaciones <- expand.grid(
  `Submission ID` = all_listas,
  `Common Name` = surfbird,
  stringsAsFactors = FALSE
)

# Join to create zero-filled dataset
df_modelo <- combinaciones %>%
  left_join(checklist_info, by = "Submission ID") %>%
  left_join(
    df %>% dplyr::select(`Submission ID`, `Common Name`, Count),
    by = c("Submission ID", "Common Name")
  ) %>%
  mutate(Count = replace_na(Count, 0))

# Replace NAs with mean value from data
df_modelo <- df_modelo %>%
  mutate(`Duration (Min)` = replace_na(`Duration (Min)`, 23.91))

# PREPARE VARIABLES FOR MODELING
# Rename variables to avoid special characters
df_modelo <- df_modelo %>%
  dplyr::rename(
    Duration_min = `Duration (Min)`,
    N_obs = `Number of Observers`,
    pred_marea = pred_altura_marea
  ) %>%
  mutate(
    # Replace NAs in N_obs with 1
    N_obs = replace_na(N_obs, 1),
    # Cap duration at 60 minutes
    Duration_capped = pmin(Duration_min, 60, na.rm = TRUE),
    # Log duration
    Duration_log = log(Duration_min),
    # Scale latitude
    Zlat = rescale(Lat)
  )

# FIT BAYESIAN ZERO-INFLATED MODE
# Bayesian zero-inflated negative binomial with regional random effects
mod_bayes <- brm(
  Count ~ Duration_log + (1|Region),
  family = zero_inflated_negbinomial(),
  data = df_modelo,
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(2), class = "sd")
  ),
  chains = 4,
  iter = 10000,
  warmup = 5000,
  cores = 4,
  control = list(adapt_delta = 0.99)
)

# Check model summary
summary(mod_bayes)

# Posterior predictive check
pp_check(mod_bayes, ndraws = 100)

# PREPARE SPATIAL PREDICTION GRID (ROCKY COASTLINE)
# Read rocky coastline shapefile
playas_rocosas <- st_read("PLayasRocas.kml")

# Assign regions based on transect names
playas_rocosas <- playas_rocosas %>%
  mutate(Region = case_when(
    Name %in% paste0(1:15, "A") ~ "Arica",
    Name %in% paste0(16:116, "A") ~ "Tarapaca",
    Name %in% paste0(117:247, "A") ~ "Antofagasta",
    Name %in% paste0(248:344, "A") ~ "Atacama",
    Name %in% paste0(345:413, "A") ~ "Coquimbo",
    Name %in% paste0(414:486, "A") ~ "Valparaiso",
    Name %in% paste0(487:521, "A") ~ "Ohiggins",
    Name %in% paste0(522:553, "A") ~ "Maule",
    Name %in% paste0(554:563, "A") ~ "Ñuble",
    Name %in% paste0(564:623, "A") ~ "Biobío",
    Name %in% paste0(624:631, "A") ~ "Araucania",
    Name %in% paste0(632:667, "A") ~ "Los Rios",
    Name %in% paste0(668:758, "A") ~ "Los Lagos",
    Name %in% paste0(414:522, "Au") ~ "Valparaiso",
    Name %in% paste0(498:499, "a") ~ "Valparaiso",
    TRUE ~ NA_character_
  ))

# Transform to UTM for distance calculations
playas_rocosas_utm <- st_transform(playas_rocosas, 32719)

# Generate points every 500m along rocky coastline
puntos_lista <- lapply(1:nrow(playas_rocosas_utm), function(i) {
  linea <- st_geometry(playas_rocosas_utm[i, ])
  largo <- st_length(linea)

  if (as.numeric(largo) >= 500) {
    puntos <- st_line_sample(linea, density = 1/500, type = "regular")
  } else {
    puntos <- st_line_sample(linea, sample = 0.5)
  }

  puntos <- st_cast(puntos, "POINT")
  if (length(puntos) > 0) {
    st_sf(
      Transecto = paste0("T", i),
      Region = playas_rocosas_utm$Region[i],
      geometry = puntos,
      crs = st_crs(playas_rocosas_utm)
    )
  } else {
    NULL
  }
})

puntos_sf <- do.call(rbind, puntos_lista)
puntos_latlon <- st_transform(puntos_sf, 4326)

# Extract coordinates
puntos_coords <- puntos_latlon %>%
  mutate(
    Lat = st_coordinates(.)[,2],
    Long = st_coordinates(.)[,1]
  ) %>%
  st_drop_geometry()

puntos_coords$Zlat <- rescale(puntos_coords$Lat)

# Create prediction grid with standardized survey conditions
pred_grid_clean <- puntos_coords %>%
  filter(!is.na(Region)) %>%
  mutate(
    Duration_log = log(30)  # Log of standard 30-minute survey
  ) %>%
  as.data.frame()

# PREDICT ABUNDANCE ACROSS ROCKY COASTLINE
# Get posterior draws of expected values
pred_fitted_draws <- fitted(mod_bayes, newdata = pred_grid_clean,
                            type = "response", summary = FALSE)

# Sum by region for each posterior draw
regions <- pred_grid_clean$Region
unique_regions <- unique(regions[!is.na(regions)])

regional_totals_bayes <- do.call(rbind, lapply(unique_regions, function(r) {
  idx <- which(regions == r)
  region_sims <- rowSums(pred_fitted_draws[, idx])
  data.frame(
    Region = r,
    km_rocky = length(idx) * 0.5,
    Mean = mean(region_sims),
    Median = median(region_sims),
    CI_low = quantile(region_sims, 0.025),
    CI_high = quantile(region_sims, 0.975),
    row.names = NULL
  )
}))

# National total
national_sims <- rowSums(pred_fitted_draws)
national_total <- data.frame(
  Region = "TOTAL",
  km_rocky = sum(regional_totals_bayes$km_rocky),
  Mean = mean(national_sims),
  Median = median(national_sims),
  CI_low = quantile(national_sims, 0.025),
  CI_high = quantile(national_sims, 0.975)
)

# CREATE SUMMARY TABLE WITH OBSERVED DATA
# Observed data by region
observed_by_region <- df_modelo %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    observed = sum(Count),
    n_detections = sum(Count > 0)
  )

# Survey effort by region
effort_by_region <- df_modelo %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_transects = n(),
    total_duration_hrs = sum(Duration_min, na.rm = TRUE) / 60,
    km_surveyed = sum(kms, na.rm = TRUE)
  )

# Final summary table
summary_table <- regional_totals_bayes %>%
  left_join(observed_by_region, by = "Region") %>%
  left_join(effort_by_region, by = "Region") %>%
  mutate(
    observed_density = round(observed / km_surveyed, 2),
    predicted_density = round(Mean / km_rocky, 2),
    coverage_pct = round(km_surveyed / km_rocky * 100, 1),
    Estimate = paste0(round(Mean, 0), " (", round(CI_low, 0), " - ", round(CI_high, 0), ")")
  )

# Order regions from north to south
region_order <- c("Arica", "Tarapaca", "Antofagasta", "Atacama", "Coquimbo",
                  "Valparaiso", "Ohiggins", "Maule", "Ñuble", "Biobío",
                  "Araucania", "Los Rios", "Los Lagos")

summary_table <- summary_table %>%
  mutate(Region = factor(Region, levels = region_order)) %>%
  arrange(Region) %>%
  dplyr::select(
    Region,
    km_rocky,
    n_transects,
    km_surveyed,
    coverage_pct,
    n_detections,
    observed,
    observed_density,
    Estimate,
    predicted_density
  )

# Print results
print(summary_table)
print(national_total)

# Add predictions to the grid
pred_grid_clean$pred_mean <- fitted(mod_bayes, newdata = pred_grid_clean,
                                    type = "response")[, "Estimate"]

# Figure 1: Map of Surfbird survey locations along Chile's rocky coastline
# Get Chile country outline
chile <- ne_countries(scale = "medium", country = "Chile", returnclass = "sf")
south_america <- ne_countries(scale = "medium", continent = "South America", returnclass = "sf")
chile_regions <- ne_states(country = "Chile", returnclass = "sf")
# Create survey data with coordinates (using your df_modelo data)

# Simulated data
survey_points <- listas %>%
  left_join(
    df_modelo %>%
      dplyr::group_by(Transecto) %>%
      dplyr::summarise(total_count = sum(Count, na.rm = TRUE)),
    by = "Transecto"
  ) %>%
  dplyr::mutate(
    total_count = replace_na(total_count, 0),
    detected = total_count > 0
  )

# Convert to sf object
survey_sf <- st_as_sf(survey_points, coords = c("Long", "Lat"), crs = 4326)

# Regional summaries for labeling
region_labels <- df_modelo %>%
  left_join(listas %>% dplyr::select(Transecto, Long), by = "Transecto") %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Lat = mean(Lat, na.rm = TRUE),
    Long = mean(Long, na.rm = TRUE) - 2,
    total_birds = sum(Count, na.rm = TRUE),
    n_sites = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(Lat))


# Panel A: Observed Counts
map_observed <- ggplot() +
  # Regional boundaries
  geom_sf(data = chile_regions, fill = "gray95", color = "gray60", size = 0.3) +

  # Survey points - no detection (gray)
  geom_point(data = survey_points %>% filter(!detected),
             aes(x = Long, y = Lat),
             color = "gray50", size = 1.5, alpha = 0.6) +

  # Survey points - detection (orange, size by count)
  geom_point(data = survey_points %>% filter(detected),
             aes(x = Long, y = Lat, size = total_count),
             color = "#D55E00", alpha = 0.8) +

  scale_size_continuous(name = "Count",
                        range = c(2, 10),
                        breaks = c(50, 100, 150, 200)) +

  coord_sf(xlim = c(-76, -66), ylim = c(-44, -17), expand = FALSE) +

  annotation_scale(location = "br", width_hint = 0.25,
                   pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
                   text_cex = 0.7) +

  labs(title = "A) Observed counts", x = NULL, y = "Latitude") +

  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_rect(fill = "aliceblue", color = NA),
    legend.position.inside = c(0.2, 0.2),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 8)
  )

# Panel B: Model predictions
# Map region names to match your data
chile_regions <- chile_regions %>%
  dplyr::mutate(Region = case_when(
    name == "Arica y Parinacota" ~ "Arica",
    name == "Tarapacá" ~ "Tarapaca",
    name == "Libertador General Bernardo O'Higgins" ~ "Ohiggins",
    name == "La Araucanía" ~ "Araucania",
    name == "Los Ríos" ~ "Los Rios",
    TRUE ~ name
  ))

# Filter to coastal regions
coastal_regions <- c("Arica", "Tarapaca", "Antofagasta", "Atacama",
                     "Coquimbo", "Valparaiso", "Ohiggins",
                     "Maule", "Ñuble", "Biobío", "Araucania", "Los Rios", "Los Lagos")

chile_coastal <- chile_regions %>%
  filter(Region %in% coastal_regions)

# Create regional centroids for point placement
regional_centroids <- regional_totals_bayes %>%
  left_join(
    survey_points %>%
      dplyr::group_by(Region) %>%
      dplyr::summarise(
        Long = mean(Long, na.rm = TRUE),
        Lat = mean(Lat, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Region"
  ) %>%
  dplyr::mutate(
    # Calculate coefficient of variation (CV) as measure of uncertainty
    CV = (CI_high - CI_low) / (2 * Mean),
    # Certainty = inverse of CV (for alpha mapping)
    Certainty = 1 / CV
  )

map_predicted <- ggplot() +
  geom_sf(data = chile_regions, fill = "gray95", color = "gray60", size = 0.3) +
  geom_point(data = regional_centroids,
             aes(x = Long, y = Lat, size = Mean, alpha = Certainty),
             color = "#0072B2") +
  scale_size_continuous(name = "Predicted\npopulation",
                        range = c(3, 15),
                        breaks = c(1000, 3000, 6000, 9000),
                        labels = scales::comma) +
  scale_alpha_continuous(name = "Certainty",
                         range = c(0.3, 1),
                         guide = "none") +  # Hide alpha legend for simplicity
  coord_sf(xlim = c(-76, -66), ylim = c(-44, -17), expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.25,
                   pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
                   text_cex = 0.7) +
  labs(title = "B) Predicted population",
       subtitle = "(transparency = certainty)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_rect(fill = "aliceblue", color = NA),
    legend.position.inside = c(0.2, 0.25),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8, face = "italic"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8)
  )


combined_maps <- plot_grid(map_observed, map_predicted,
                           ncol = 2, align = "h", rel_widths = c(1, 0.95))
