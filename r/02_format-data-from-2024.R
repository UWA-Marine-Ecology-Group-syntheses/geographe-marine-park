library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(tidyverse)
library(googlesheets4)
library(sf)
library(terra)
library(here)

name <- "2024-04_Geographe_stereo-BRUVs"


# Read in metadata
metadata <- read_metadata(here::here("data/raw/em export/"), method = "BRUVs") %>% # Change here to "DOVs"
  dplyr::select(campaignid, sample, 
                status, 
                longitude_dd, latitude_dd, 
                observer_count, observer_length,
                date_time, 
                location, site, 
                depth_m, 
                successful_count, 
                successful_length, 
                successful_habitat_forwards, 
                successful_habitat_backwards) %>%
  glimpse()

write_csv(metadata, paste0("data/uploads/", name, "_metadata.csv"))

saveRDS(metadata, file = here::here(paste0("data/tidy/",
                                           name, "_Metadata.rds")))

metadata_sf <- st_as_sf(metadata, coords = c("longitude_dd", "latitude_dd"), crs = 4326)
regions <- st_as_sf(CheckEM::aus_regions, crs = st_crs(4326))
regions <- st_transform(regions, 4326) %>%
  dplyr::select(REGION)

metadata <- st_join(metadata_sf, regions, join = st_nearest_feature) %>%
  dplyr::rename(marine_region = REGION) %>%
  dplyr::mutate(sample = as.character(sample)) %>%
  as.data.frame() %>%
  dplyr::select(-c(geometry)) %>%
  glimpse()


# Read in the count and length data
points <- read_points(here::here("data/raw/em export/")) %>%
  glimpse()

maxn_points <- points %>% 
  dplyr::mutate(species = ifelse(species %in% c("sp", "sp1", "sp2"), "spp", as.character(species))) %>%
  dplyr::group_by(campaignid, sample, filename, periodtime, frame, family, genus, species, stage) %>% # If you have MaxN'd by stage (e.g. Adult, Juvenile) add stage here
  dplyr::mutate(number = as.numeric(number)) %>%
  dplyr::summarise(maxn = sum(number)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(campaignid, sample, family, genus, species, stage) %>%
  dplyr::slice(which.max(maxn)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(maxn)) %>%
  dplyr::select(-frame) %>%
  tidyr::replace_na(list(maxn = 0)) %>%
  dplyr::mutate(maxn = as.numeric(maxn)) %>%
  dplyr::filter(maxn > 0) %>%
  dplyr::inner_join(metadata, by = join_by(campaignid, sample)) %>%
  dplyr::filter(successful_count %in% c("Yes")) %>% 
  dplyr::filter(maxn > 0) %>%
  dplyr::select(campaignid, sample, family, genus, species, maxn, stage) %>%
  dplyr::glimpse()

maxn <- bind_rows(get0("maxn_points"), get0("maxn_counts")) # this works even if you only have one type of data

em_length3dpoints <- read_em_length(here::here("data/raw/em export/")) %>%
  dplyr::mutate(species = ifelse(species %in% c("sp", "sp1", "sp2"), "spp", as.character(species))) %>%
  dplyr::select(-c(comment))%>% # there is a comment column in metadata, so you will need to remove this column from EM data
  dplyr::inner_join(metadata, by = join_by(sample, campaignid)) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  # dplyr::rename(length_mm = length) %>%
  glimpse() 

# If only EventMeasure data then length only includes Length and 3D points data
# If only Generic data then length only includes generic length data
# If both exist, then length includes both Length and 3D points and generic length data
length <- bind_rows(get0("em_length3dpoints"), get0("gen_length")) # this works even if you only have one type of data


# Format data for checks
count <- maxn %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  dplyr::select(campaignid, sample, family, genus, species, maxn) %>%
  tidyr::complete(nesting(campaignid, sample), nesting(family, genus, species)) %>%
  tidyr::replace_na(list(maxn = 0)) %>%
  group_by(campaignid, sample, family, genus, species) %>%
  dplyr::summarise(count = sum(maxn)) %>%
  ungroup() %>%
  mutate(scientific = paste(family, genus, species, sep = " "))%>%
  dplyr::select(campaignid, sample, scientific, count)%>%
  spread(scientific, count, fill = 0)

count_families <- maxn %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(scientific = paste(family, genus, species, sep = " ")) %>%
  filter(!(family %in% "Unknown")) %>%
  dplyr::select(c(family, genus, species, scientific)) %>%
  distinct()

complete_count <- count %>%
  pivot_longer(names_to = "scientific", values_to = "count",
               cols = 3:ncol(.)) %>%
  inner_join(count_families, by = c("scientific")) %>%
  full_join(metadata)%>%
  filter(successful_count %in% "Yes") %>%
  glimpse()

complete_length <- length %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  dplyr::filter(!family %in% "Unknown")%>%
  # First make one row for every length measurement
  dplyr::mutate(number = as.numeric(number)) %>%
  tidyr::uncount(number) %>%
  dplyr::mutate(number = 1) %>% 
  # Add in missing samples
  dplyr::right_join(metadata) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  # Complete the data (add in zeros for every species)
  dplyr::select(campaignid, sample, family, genus, species, length_mm, number, any_of(c("range", "rms", "precision"))) %>% # this will keep EM only columns
  tidyr::complete(nesting(campaignid, sample), nesting(family, genus, species)) %>%
  replace_na(list(number = 0)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(number)) %>%
  dplyr::mutate(length_mm = as.numeric(length_mm)) %>%
  left_join(., metadata) %>%
  glimpse()


# Check data
number_of_samples <- metadata %>%
  dplyr::distinct(campaignid, sample)

message(paste(nrow(number_of_samples), "unique samples in the metadata"))

duplicate_samples <- metadata %>%
  dplyr::group_by(campaignid, sample) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 1)

message(paste(nrow(duplicate_samples), "samples duplicated in the metadata"))

metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")), 
                successful_count, successful_length) %>%
  distinct()

samples <- maxn %>%
  distinct(campaignid, sample)

missing_count <- anti_join(metadata_samples, samples, by = join_by(campaignid, sample))

message(paste(nrow(missing_count), "samples in the metadata missing count data"))

missing_metadata <- anti_join(samples, metadata_samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_metadata), "samples in count data missing metadata"))

metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")), 
                successful_count, successful_length) %>%
  distinct()

samples <- length %>%
  distinct(campaignid, sample)

missing_length <- anti_join(metadata_samples, samples, by = join_by(campaignid, sample))

message(paste(nrow(missing_length), "samples in metadata missing length data"))

missing_metadata <- anti_join(samples, metadata_samples, by = join_by(campaignid, sample))

message(paste(nrow(missing_metadata), "samples in length data missing metadata"))

periods <- read_periods(here::here("data/raw/em export/")) %>%
  glimpse()

periods_without_end <- periods %>%
  dplyr::filter(has_end == 0)

message(paste(nrow(periods_without_end), "periods without an end"))
glimpse(periods_without_end)

metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")), successful_count, successful_length) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sample = as.factor(sample))

periods_samples <- periods %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period"))) %>%
  distinct()

missing_periods <- anti_join(metadata_samples, periods_samples) %>%
  dplyr::select(!sample)

message(paste(nrow(missing_periods), "samples missing period"))
glimpse(missing_periods)

points_outside_periods <- points %>%
  dplyr::filter(period %in% c("NA", NA, NULL, "")) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, number, frame)

message(paste(nrow(points_outside_periods), "points outside a period"))
glimpse(points_outside_periods)

lengths_outside_periods <- em_length3dpoints %>%
  dplyr::filter(period %in% c("NA", NA, NULL, "")) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, number)

message(paste(nrow(lengths_outside_periods), "lengths/3D points outside period"))
glimpse(lengths_outside_periods)

period_length <- 60 # in minutes

periods_wrong <- periods %>%
        dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), time_start, time_end, has_end) %>%
        dplyr::distinct() %>%
        dplyr::mutate(period_time = round(time_end - time_start)) %>%
        dplyr::filter(!period_time %in% period_length)

message(paste(nrow(periods_wrong), "periods not", period_length, "minutes long"))
glimpse(periods_wrong)

total_count <- sum(complete_count$count)
message(paste(total_count, "fish counted in the count data"))

total_length <- sum(complete_length$number)
message(paste(total_length, "fish counted in the length data"))

length.vs.maxn <- complete_length %>%
  dplyr::group_by(campaignid, sample, family, genus, species) %>%
  dplyr::summarise(length_maxn = sum(number)) %>%
  dplyr::ungroup() %>%
  dplyr::full_join(complete_count) %>%
  tidyr::replace_na(list(length_maxn = 0)) %>%
  filter(!length_maxn %in% count) %>%
  glimpse()

message(paste(nrow(length.vs.maxn), "MaxN fish not measured for length"))

points_without_number <- points %>%
  filter(number %in% c("NA", NA, 0, NULL, "", " "))

message(paste(nrow(points_without_number), "points in the _Points.txt file that do not have a number"))
glimpse(points_without_number)

lengths_without_number <- em_length3dpoints %>%
  filter(number %in% c("NA", NA, 0, NULL, "", " "))

message(paste(nrow(lengths_without_number), "lengths or 3D points in the EMObs that do not have a number"))
glimpse(lengths_without_number)

synonyms_in_count <- dplyr::left_join(complete_count, CheckEM::aus_synonyms) %>%
      dplyr::filter(!is.na(genus_correct)) %>%
      dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
      dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
      dplyr::select('old name', 'new name') %>%
      dplyr::distinct()

message(paste(nrow(synonyms_in_count), "synonyms used in the count data"))
glimpse(synonyms_in_count)

synonyms_in_length <- dplyr::left_join(complete_length, CheckEM::aus_synonyms) %>%
      dplyr::filter(!is.na(genus_correct)) %>%
      dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
      dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
      dplyr::select('old name', 'new name') %>%
      dplyr::distinct()

message(paste(nrow(synonyms_in_length), "synonyms used in the length data"))
glimpse(synonyms_in_length)

complete_count <- dplyr::left_join(complete_count, CheckEM::aus_synonyms) %>%
  dplyr::mutate(genus = ifelse(!genus_correct%in%c(NA), genus_correct, genus)) %>%
  dplyr::mutate(species = ifelse(!is.na(species_correct), species_correct, species)) %>%
  dplyr::mutate(family = ifelse(!is.na(family_correct), family_correct, family)) %>%
  dplyr::select(-c(family_correct, genus_correct, species_correct)) %>%
  dplyr::mutate(scientific = paste(family, genus, species))

complete_length <- dplyr::left_join(complete_length, CheckEM::aus_synonyms) %>%
  dplyr::mutate(genus = ifelse(!genus_correct%in%c(NA), genus_correct, genus)) %>%
  dplyr::mutate(species = ifelse(!is.na(species_correct), species_correct, species)) %>%
  dplyr::mutate(family = ifelse(!is.na(family_correct), family_correct, family)) %>%
  dplyr::select(-c(family_correct, genus_correct, species_correct)) %>%
  dplyr::mutate(scientific = paste(family, genus, species))

count_species_not_observed_region <- complete_count %>%
  dplyr::distinct(campaignid, sample, family, genus, species, marine_region, count) %>%
  dplyr::anti_join(., expand_life_history(CheckEM::australia_life_history), by = c("family", "genus", "species", "marine_region")) %>%
  dplyr::filter(count > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, marine_region) %>%
  dplyr::distinct() %>%
  dplyr::rename('marine region not observed in' = marine_region) %>%
  dplyr::semi_join(., CheckEM::australia_life_history, by = c("family", "genus", "species"))

message(paste(nrow(count_species_not_observed_region), "species not observed in the region before"))

length_species_not_observed_region <- complete_length %>%
  dplyr::distinct(campaignid, sample, family, genus, species, marine_region, number) %>%
  dplyr::anti_join(., expand_life_history(CheckEM::australia_life_history), by = c("family", "genus", "species", "marine_region")) %>%
  dplyr::filter(number > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, marine_region) %>%
  dplyr::distinct() %>%
  dplyr::rename('marine region not observed in' = marine_region) %>%
  dplyr::semi_join(., CheckEM::australia_life_history, by = c("family", "genus", "species"))

message(paste(nrow(length_species_not_observed_region), "species not observed in the region before"))
glimpse(length_species_not_observed_region)

count_species_not_in_list <- complete_count %>%
  dplyr::anti_join(., CheckEM::australia_life_history, by = c("family", "genus", "species")) %>%
  dplyr::filter(count > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species) %>%
  dplyr::distinct() 

message(paste(nrow(count_species_not_in_list), "species not in chosen life history list"))
glimpse(count_species_not_in_list)

length_species_not_in_list <- complete_length %>%
  dplyr::anti_join(., CheckEM::australia_life_history, by = c("family", "genus", "species")) %>%
  dplyr::filter(number > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species) %>%
  dplyr::distinct() 

message(paste(nrow(length_species_not_in_list), "species not in chosen life history list"))
glimpse(length_species_not_in_list)

incorrect_lengths <- left_join(complete_length, create_min_max(CheckEM::australia_life_history, minimum = 0.15, maximum = 0.85)) %>%
  dplyr::filter(length_mm < min_length_mm | length_mm > max_length_mm) %>%
  mutate(reason = ifelse(length_mm < min_length_mm, "too small", "too big")) %>%
  dplyr::select(campaignid, sample, family, genus, species, length_mm, min_length_mm, max_length_mm, length_max_mm, reason, any_of(c("em_comment", "frame_left")), length_max_mm) %>%
  mutate(difference = ifelse(reason %in% c("too small"), (min_length_mm - length_mm), (length_mm - max_length_mm))) %>%
  dplyr::mutate(percent_of_fb_max = (length_mm/length_max_mm)*100) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, length_mm, min_length_mm, max_length_mm, length_max_mm, reason, any_of(c("em_comment", "frame_left")), difference, percent_of_fb_max)

too_small <- incorrect_lengths %>%
  dplyr::filter(reason %in% "too small")

too_big <- incorrect_lengths %>%
  dplyr::filter(reason %in% "too big")

message(paste(nrow(too_small), "lengths are too small"))
glimpse(too_small)

message(paste(nrow(too_big), "lengths are too big"))
glimpse(too_big)

rms_limit <- 20 # in mm

over_rms <- complete_length %>%
  dplyr::filter(as.numeric(rms) > rms_limit)

message(paste(nrow(over_rms), "lengths over RMS limit"))
glimpse(over_rms)

precision_limit <- 10 # in %

over_precision <- complete_length %>%
  dplyr::mutate(precision_percent = (precision/length_mm) *100) %>%
  dplyr::filter(precision_percent > precision_limit)

message(paste(nrow(over_precision), "lengths over precision limit"))
glimpse(over_precision)

range_limit <- 8 # in metres

over_range <- complete_length %>%
  dplyr::filter(as.numeric(range) > (range_limit* 1000))

message(paste(nrow(over_range), "lengths over range limit"))
glimpse(over_range)

saveRDS(complete_count,
        file = here::here(paste0("data/staging/",
                                 name, "_complete-count.rds")))
saveRDS(complete_length,
        file = here::here(paste0("data/staging/",
                                 name, "_complete-length.rds")))

# Format data for upload to GA and final checks
codes <- australia_life_history %>%
  dplyr::select(family, genus, species, caab_code)

count_upload <- maxn %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  {message(paste(length(which(.$family %in% "Unknown")), "rows removed because family is 'Unknown'"));
    .} %>%
  dplyr::filter(!family %in% "Unknown")%>%
  dplyr::rename(count = maxn) %>%
  left_join(codes)

length_upload <- em_length3dpoints %>% # This includes only EM data, not generic length
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  {message(paste(length(which(.$family %in% "Unknown")), "rows removed because family is 'Unknown'"));
    .} %>%
  dplyr::filter(!family %in% "Unknown")%>%
  left_join(codes) %>%
  dplyr::rename(rms_mm = rms, range_mm = range, precision_mm = precision, count = number)

# Check missing caab codes (okay if from sp, sp1 etc.)
n_no_caab_count <- length(count_upload$species[which(is.na(count_upload$caab_code))])
sp_no_caab_count <- unique(count_upload$species[which(is.na(count_upload$caab_code))])
n_no_caab_length <- length(length_upload$species[which(is.na(length_upload$caab_code))])
sp_no_caab_length <- unique(length_upload$species[which(is.na(length_upload$caab_code))])

message(paste(n_no_caab_count, "count rows without caab codes from species:", paste(sp_no_caab_count, collapse = ", ")))
message(paste(n_no_caab_length, "length rows without caab codes from species:", paste(sp_no_caab_length, collapse = ", ")))

# Save GA upload data
write_csv(count_upload, paste0("data/uploads/", name, "_count.csv"))
write_csv(length_upload, paste0("data/uploads/", name, "_length.csv"))
