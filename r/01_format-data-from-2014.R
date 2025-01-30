library(tidyverse)
library(CheckEM)

caab_codes <- CheckEM::australia_life_history %>%
  dplyr::distinct(caab, family, genus, species) %>%
  rename(caab_code = caab)

GeographeAMP_complete_count <- readRDS("H:/south-west-network/data/geographe/raw/GeographeAMP_complete_count.RDS")
GeographeAMP_complete_length <- readRDS("H:/south-west-network/data/geographe/raw/GeographeAMP_complete_length.RDS")

test_count <- GeographeAMP_complete_count %>%
  dplyr::filter(count > 0)

test_length <- GeographeAMP_complete_length %>%
  dplyr::filter(number > 0)

# Read in metadata ----
metadata <- read_csv("data/geographe/raw/temp/2007-2014-Geographe-stereo-BRUVs.checked.metadata.csv") %>%
  dplyr::filter(campaignid %in% "2014-12_Geographe.Bay_stereoBRUVs") %>%
  CheckEM::clean_names() %>%
  dplyr::rename(latitude_dd = latitude,
                longitude_dd = longitude,
                depth_m = depth,
                observer_count = observer) %>%
  dplyr::select(-c(raw_hdd_number, con_hdd_number, state_zone, commonwealth_zone)) %>%
  dplyr::mutate(observer_length = "Unknown") %>%
  dplyr::mutate(date = if_else(sample %in% "68", "2014-12-10", date)) %>%
  dplyr::mutate(time = if_else(sample %in% "68", "14:07:00", time)) %>%
  dplyr::mutate(time = if_else(str_length(time) == 5, str_c(time, ":00"), time)) %>%
  dplyr::mutate(date_time = paste0(date, "T", time, "+08:00")) %>%
  dplyr::mutate(depth_m = if_else(depth_m %in% "N/A", "0", depth_m)) %>%
  dplyr::mutate(observer_count = if_else(successful_count %in% "Yes" & is.na(observer_count), "Unknown", observer_count)) %>%
  dplyr::mutate(depth_m = as.numeric(depth_m)) %>%
  dplyr::select(-c(date, time)) %>%
  dplyr::glimpse()

names(metadata)
unique(metadata$campaignid)
unique(metadata$date) %>% sort()
unique(metadata$time) %>% sort()
unique(metadata$latitude_dd) %>% sort()
unique(metadata$longitude_dd) %>% sort()
unique(metadata$location) %>% sort()
unique(metadata$site) %>% sort()
unique(metadata$depth_m) %>% sort()
unique(metadata$successful_count) %>% sort()
unique(metadata$successful_length) %>% sort()
unique(metadata$successful_length) %>% sort()
unique(metadata$status) %>% sort()

# read in Count ----
count_raw <- read_csv("data/geographe/raw/temp/2007-2014-Geographe-stereo-BRUVs.complete.maxn.csv") %>%
  dplyr::filter(campaignid %in% "2014-12_Geographe.Bay_stereoBRUVs") %>%
  dplyr::filter(maxn > 0) %>%
  CheckEM::clean_names() %>%
  dplyr::select(campaignid, sample, family, genus, species, maxn) %>%
  dplyr::rename(count = maxn) %>%
  dplyr::glimpse()

# run count data with synonyms list ----
count <- dplyr::left_join(count_raw, CheckEM::aus_synonyms) %>%
  dplyr::mutate(genus = ifelse(!genus_correct%in%c(NA), genus_correct, genus)) %>%
  dplyr::mutate(species = ifelse(!is.na(species_correct), species_correct, species)) %>%
  dplyr::mutate(family = ifelse(!is.na(family_correct), family_correct, family)) %>%
  dplyr::select(-c(family_correct, genus_correct, species_correct)) %>%
  dplyr::group_by(sample, family, genus, species) %>%
  dplyr::slice(which.max(count)) %>%
  left_join(caab_codes) %>%
  dplyr::filter(!is.na(caab_code)) %>% # only removes 3 unknowns
  ungroup()

# read in length ----
length_raw <- read_csv("data/geographe/raw/temp/2007-2014-Geographe-stereo-BRUVs.expanded.length.csv") %>%
  dplyr::filter(campaignid %in% "2014-12_Geographe.Bay_stereoBRUVs") %>%
  dplyr::select(campaignid, sample, family, genus, species, length, range) %>%
  dplyr::rename(length_mm = length,
                range_mm = range) %>%
  dplyr::mutate(precision_mm = 1, rms_mm = 1)

names(length_raw)

# TODO check that the number column should be 1?

# fix synonyms ----
length <- length_raw %>%
  dplyr::left_join(., CheckEM::aus_synonyms) %>%
  dplyr::mutate(genus = ifelse(!genus_correct%in%c(NA), genus_correct, genus)) %>%
  dplyr::mutate(species = ifelse(!is.na(species_correct), species_correct, species)) %>%
  dplyr::mutate(family = ifelse(!is.na(family_correct), family_correct, family)) %>%
  dplyr::select(-c(family_correct, genus_correct, species_correct)) %>%
  dplyr::mutate(number = 1) %>%
  left_join(caab_codes) %>%

  dplyr::filter(!(sample %in% "492" & species %in% "vittiger")) %>%

  glimpse()

names(length)

# Read in habitat ----
habitat <- read_csv("data/geographe/raw/temp/2014-12_Geographe.Bay_stereoBRUVs_habitat.csv") %>%
  semi_join(metadata)



# Find length that doesn't have count ----



# write data to be uploaded ----
write.csv(metadata, "data/geographe/uploads/2014-12_Geographe.Bay_stereoBRUVs_metadata.csv", row.names = F)
write.csv(length, "data/geographe/uploads/2014-12_Geographe.Bay_stereoBRUVs_length.csv", row.names = F)
write.csv(count, "data/geographe/uploads/2014-12_Geographe.Bay_stereoBRUVs_count.csv", row.names = F)
write.csv(habitat, "data/geographe/uploads/2014-12_Geographe.Bay_stereoBRUVs_habitat.csv", row.names = F)
