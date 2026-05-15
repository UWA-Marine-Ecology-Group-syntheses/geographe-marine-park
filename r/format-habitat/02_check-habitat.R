# Clear obejcts from the environment
rm(list = ls())

library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(dplyr)
library(tidyr)
library(stringr)
library(ggbeeswarm)
library(RColorBrewer)
library(leaflet)
library(leaflet.minicharts)
library(here)

name <- "2024-04_Geographe_BOSS"

metadata <- read_csv("data/raw/MEG_Labsheets - 2024-04_Geographe_BOSS_metadata.csv") %>%
  #dplyr::rename(sample = opcode) %>% #if BOSS put hash in front, if BRUV remove #
  dplyr::select(sample, longitude_dd, latitude_dd, date_time, location, site, depth_m, successful_count, successful_length, successful_habitat_panoramic) %>% #Check forward facing camera or habitat panoramic and make sure it matches what is in the lab sheet
  glimpse()

habitat <- read_csv("data/tidy/2024-04_Geographe_AMP_benthos-count.csv") %>%
  dplyr::mutate(sample = str_replace(sample, "2024-04_Geographe_BOSS_", "")) %>%
  dplyr::mutate(habitat = case_when(
    level_2 %in% "Macroalgae" ~ level_2, 
    level_2 %in% "Seagrasses" ~ level_2, 
    level_3 %in% "Consolidated (hard)" ~ level_3,
    level_3 %in% "Unconsolidated (soft)" ~ level_3, 
    level_2 %in% "Sponges" ~ "Sessile invertebrates", 
    level_2 %in% "Sessile invertebrates" ~ level_2, 
    level_2 %in% "Bryozoa" ~ "Sessile invertebrates", 
    level_2 %in% "Cnidaria" ~ "Sessile invertebrates")) %>%
  dplyr::select(sample, habitat, count) %>%
  group_by(sample, habitat) %>%
  dplyr::summarise(count = sum(count)) %>%
  dplyr::mutate(total_points_annotated = sum(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = "habitat", values_from = "count", values_fill = 0) %>%
  dplyr::mutate(reef = Macroalgae + Seagrasses + `Sessile invertebrates` + `Consolidated (hard)`) %>%
  pivot_longer(cols = c("Macroalgae", "Seagrasses", "Sessile invertebrates", "Consolidated (hard)", "Unconsolidated (soft)", "reef"), names_to = "habitat", values_to = "count") %>%
  glimpse()


habitat.missing.metadata <- anti_join(habitat, metadata, by = c("sample")) %>%
  glimpse()

metadata.missing.habitat <- anti_join(metadata, habitat, by = c("sample")) %>%
  glimpse()

tidy.relief <- read_csv("data/tidy/2024-04_Geographe_AMP_benthos-relief.csv") %>%
  dplyr::mutate(sample = str_replace(sample, "2024-04_Geographe_BOSS_", "")) %>%
  uncount(count) %>%
  filter(!is.na(level_5)) %>%
  group_by(sample) %>%
  dplyr::summarise(mean.relief = mean(as.numeric(level_5)), sd.relief = sd(as.numeric(level_5), na.rm = T)) %>%
  ungroup() %>%
  glimpse()

relief.missing.metadata <- anti_join(tidy.relief, metadata, by = c("sample")) %>%
  glimpse()

metadata.missing.relief <- anti_join(metadata, tidy.relief, by = c("sample")) %>%
  glimpse()

tidy.habitat <- metadata %>%
  left_join(habitat) %>%
  left_join(tidy.relief) %>%
  dplyr::mutate(longitude_dd = as.numeric(longitude_dd),
                latitude_dd = as.numeric(latitude_dd)) %>%
  clean_names() %>%
  glimpse()

plot.relief <- read_csv("data/tidy/2024-04_Geographe_AMP_benthos-relief.csv") %>%
  group_by(campaignid, sample, level_5) %>%
  dplyr::summarise(count = sum(count)) %>%
  ungroup() %>%
  dplyr::mutate(class.relief = as.factor(level_5)) %>%
  glimpse()

ggplot() +
  geom_quasirandom(data = tidy.habitat, aes(x = (count/total_points_annotated), y = habitat), groupOnX = F, method = "quasirandom", alpha = 0.25, size = 1.8, width = 0.2) +
  labs(x = "Number of points", y = "") +
  theme_classic()

ggplot() +
  geom_quasirandom(data = plot.relief, aes(x = count, y = class.relief), groupOnX = F, method = "quasirandom", alpha = 0.25, size = 1.8, width = 0.05) +
  labs(x = "Number of points", y = "Relief (0-5)") +
  theme_classic()

cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tidy.habitat$habitat)))

plot.habitat <- tidy.habitat %>%
  pivot_wider(names_from = "habitat", values_from = "count", names_prefix = "broad.") %>%
  glimpse()

plot.leaflet <- tidy.habitat %>%
  filter(!habitat %in% "reef") %>%
  pivot_wider(names_from = "habitat", values_from = "count", names_prefix = "broad.") %>%
  glimpse()

leaflet() %>%
  addTiles(group = "Open Street Map") %>%
  addProviderTiles('Esri.WorldImagery', group = "World Imagery") %>%
  addLayersControl(baseGroups = c("World Imagery", "Open Street Map"), options = layersControlOptions(collapsed = FALSE)) %>%
  addMinicharts(plot.leaflet$longitude_dd, plot.leaflet$latitude_dd, 
                type = "pie", colorPalette = cols, 
                chartdata = plot.leaflet[grep("broad", names(plot.leaflet))], 
                width = 40, 
                transitionTime = 0) %>%
  setView(mean(as.numeric(plot.leaflet$longitude_dd)),
          mean(as.numeric(plot.leaflet$latitude_dd)), zoom = 12)

hab.name <- 'Sessile invertebrates'

overzero <-  tidy.habitat %>%
  filter(habitat %in% hab.name & count > 0)

equalzero <- tidy.habitat %>%
  filter(habitat %in% hab.name & count == 0)

bubble.plot <- leaflet(data = tidy.habitat) %>%
  addTiles() %>%
  addProviderTiles('Esri.WorldImagery', group = "World Imagery") %>%
  addLayersControl(baseGroups = c("Open Street Map", "World Imagery"), options = layersControlOptions(collapsed = FALSE))

if (nrow(overzero)) {
  bubble.plot <- bubble.plot %>%
    addCircleMarkers(data = overzero, lat = ~ latitude_dd, lng = ~ longitude_dd, radius = ~ count + 3, fillOpacity = 0.5, stroke = FALSE, label = ~ as.character(sample))}

if (nrow(equalzero)) {
  bubble.plot <- bubble.plot %>%
    addCircleMarkers(data = equalzero, lat = ~ latitude_dd, lng = ~ longitude_dd, radius = 2, fillOpacity = 0.5, color = "white", stroke = FALSE, label = ~ as.character(sample))}
bubble.plot

library(tidyverse)

# ---- 1. Select and reshape habitat columns ----
df_long <- plot.habitat %>%
  dplyr::select (depth_m, starts_with("broad.")) %>%  # adjust if needed
  pivot_longer(
    cols = starts_with("broad."),
    names_to = "habitat",
    values_to = "count"
  ) %>%
  mutate(habitat = str_remove(habitat, "broad\\.")) %>%
  filter(count > 0)

# ---- 2. Plot ----
plot_alldata <- ggplot(df_long, aes(x = depth_m, weight = count)) +
  geom_histogram(binwidth = 5, colour = "black") +  # adjust binwidth as needed
  facet_wrap(~ habitat, scales = "free_y") +
  labs(
    x = "Depth (m)",
    y = "Count",
    title = "Depth distribution by habitat"
  ) +
  theme_minimal()

plot_alldata

# 3. 0-100m in depth

df_long100 <- df_long %>%
  filter(depth_m >= 0, depth_m <= 100)

plot_100 <- ggplot(df_long100, aes(x = depth_m, weight = count)) +
  geom_histogram(binwidth = 5, colour = "black") +  # adjust binwidth as needed
  facet_wrap(~ habitat, scales = "free_y") +
  labs(
    x = "Depth (m)",
    y = "Count",
    title = "Depth distribution by habitat"
  ) +
  theme_minimal()

plot_100

# check for deep macroalgae or seagrasses
deep_plants <- plot.habitat %>%
  dplyr::select (sample, depth_m, starts_with("broad.")) %>%  # adjust if needed
  pivot_longer(
    cols = starts_with("broad."),
    names_to = "habitat",
    values_to = "count"
  ) %>%
  mutate(habitat = str_remove(habitat, "broad\\.")) %>%
  filter(count > 0) %>%
  filter(depth_m > 70, habitat %in% c("Macroalgae", "Seagrasses"))


#histogram of all samples
ggplot(metadata, aes(x = depth_m)) +
  geom_histogram(binwidth = 5, colour = "black") +  # adjust binwidth as needed
  labs(
    x = "Depth (m)",
    y = "Count",
  ) +
  theme_minimal()


