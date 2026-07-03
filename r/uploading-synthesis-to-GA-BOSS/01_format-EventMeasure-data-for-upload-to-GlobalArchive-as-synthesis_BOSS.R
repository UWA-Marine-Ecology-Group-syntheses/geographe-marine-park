# Before you use this script please read the below instructions!!
# This workflow assumes all of your metadata is correctly formatted, and has been cleaned by using CheckEM. CheckEM available here: https://marine-ecology.shinyapps.io/CheckEM/
# Please see the CheckEM usermanual for the correct format(https://globalarchivemanual.github.io/CheckEM/articles/manuals/CheckEM_user_guide.html)

# Load libraries -----
library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM") # Run this if you do not have CheckEM installed.
library(CheckEM)
library(tidyverse)
library(here)
library(sf)
library(dplyr)

# Set name for synthesis, will be used as prefix for your files to upload
name <- "geographe-marine-park"

# All data needs to be saved in a folder structure "data/raw/" for the script to work, or change the direcory when reading in files.

# Read in metadata ----
# Note check if you have used opcode, opcode or period and change the below code accordingly.


metadata <- read_metadata(here::here("data/raw/BOSS/"), method = "BOSS") %>% # Change here to "DOVs"
  dplyr::select(campaignid, sample, 
                status, 
                longitude_dd, latitude_dd, 
                observer_count, observer_length,
                date_time, 
                location, site, 
                depth_m, 
                successful_count, 
                successful_length, 
                successful_habitat_panoramic, 
                successful_habitat_downward,
                observer_habitat_panoramic,
                observer_habitat_downward) %>%
  rename(period = sample) %>% # use this line if you need to rename opcode to opcode
  dplyr::filter(campaignid %in% c("2021-03_Geographe_BOSS", "2024-04_Geographe_BOSS")) %>%
  glimpse()

#checks for duplicates of campaign IDs and opcodes
duplicates <- metadata %>%
  dplyr::group_by(campaignid, period) %>%
  dplyr::filter(n() > 1) 

#checks for duplicates in coordinates
metadata %>%
  dplyr::group_by(latitude_dd, longitude_dd) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::arrange(latitude_dd, longitude_dd)

unique(metadata$campaignid)

metadata_number_of_samples <- metadata %>%
  group_by(campaignid) %>%
  dplyr::summarise(n=n())

metadata_number_of_samples_count <- metadata %>%
  filter(successful_count %in% "Yes") %>%
  group_by(campaignid) %>%
  dplyr::summarise(n=n())

metadata_number_of_samples_length <- metadata %>%
  filter(successful_length %in% "Yes") %>%
  group_by(campaignid) %>%
  dplyr::summarise(n=n())

sum(metadata_number_of_samples_length$n)
sum(metadata_number_of_samples_count$n)
sum(metadata_number_of_samples$n)

write_csv(metadata, paste0("data/uploads/BOSS/", name, "_metadata.csv"))

# # Read in the maxn and length data ----
maxn <- read_points(here::here("data/raw/BOSS//"), method = "BOSS") %>%
  dplyr::mutate(species = ifelse(species %in% c("sp", "sp1", "sp2", "sp10"), "spp", as.character(species))) %>%
  dplyr::group_by(campaignid, period, filename, periodtime, frame, family, genus, species, stage) %>% # If you have MaxN'd by stage (e.g. Adult, Juvenile) add stage here
  dplyr::mutate(number = as.numeric(number)) %>%
  dplyr::summarise(maxn = sum(number)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(stage = "AD") %>%
  dplyr::group_by(campaignid, period, family, genus, species, stage) %>%
  dplyr::slice(which.max(maxn)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(maxn)) %>%
  dplyr::select(-frame) %>%
  tidyr::replace_na(list(maxn = 0)) %>%
  dplyr::mutate(maxn = as.numeric(maxn)) %>%
  dplyr::filter(maxn > 0) %>%
  dplyr::inner_join(metadata, by = join_by(campaignid, period)) %>%
  dplyr::filter(successful_count %in% c("Y")) %>%
  dplyr::filter(maxn > 0) %>%
  dplyr::select(campaignid, period, family, genus, species, maxn, stage) %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%

  {message(paste(length(which(.$family %in% "Unknown")), "rows removed because family is 'Unknown'"));
    .} %>%
dplyr::filter(!family %in% "Unknown") # %>%
# dplyr::mutate(stage_join = case_when(
#     stage %in% c("M", "F", "J") ~ "AD",
#     TRUE ~ stage))

unique(maxn$stage)
#unique(maxn$stage_join)
metadata %>%
  count(campaignid, period) %>%
  filter(n > 1)

metadata %>%
  count(campaignid, period) %>%
  filter(n > 1)
# 
# length <- read_em_length(here::here("data/raw/")) %>%
#   dplyr::mutate(species = ifelse(species %in% c("sp", "sp1", "sp2", "sp10"), "spp", as.character(species))) %>%
#   dplyr::select(-c(comment))%>% # there is a comment column in metadata, so you will need to remove this column from EM data
#   dplyr::inner_join(metadata, by = join_by(opcode, campaignid)) %>%
#   dplyr::filter(successful_length %in% "Yes") %>%
#   dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
#   dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
#   dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species)))  %>%
#   
#   {message(paste(length(which(.$family %in% "Unknown")), "rows removed because family is 'Unknown'"));
#     .} %>%
#   dplyr::filter(!family %in% "Unknown")%>%
#   dplyr::mutate(stage = "AD")
# 
# unique(length$stage)
# 
# Format data for upload to GA and final checks -----
codes <- australia_life_history %>%
  dplyr::select(family, genus, species, caab_code)

count_upload <- maxn %>%
  dplyr::rename(count = maxn) %>%
  dplyr::mutate(across(c(family, genus, species), as.character)) %>%
  dplyr::left_join(., CheckEM::aus_synonyms) %>%
  dplyr::mutate(
    genus   = dplyr::coalesce(as.character(genus_correct),   genus),
    species = dplyr::coalesce(as.character(species_correct), species),
    family  = dplyr::coalesce(as.character(family_correct),  family)
  ) %>%
  dplyr::select(-c(family_correct, genus_correct, species_correct)) %>%
  dplyr::mutate(scientific = paste(family, genus, species)) %>%
  left_join(codes) %>%
  glimpse()
# 
# length(unique(count_upload$opcode))
# length(unique(length_upload$opcode))
# count_upload_check <- count_upload %>%
#   distinct(campaignid, opcode) %>%
#   group_by(campaignid) %>%
#   dplyr::summarise(n=n())
# 
# length_upload_check <- length_upload %>%
#   distinct(campaignid, opcode) %>%
#   group_by(campaignid) %>%
#   dplyr::summarise(n=n())
# 
synonyms_in_count <- dplyr::left_join(count_upload, CheckEM::aus_synonyms) %>%
  dplyr::filter(!is.na(genus_correct)) %>%
  dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
  dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
  dplyr::select('old name', 'new name') %>%
  dplyr::distinct()
#   
# 
# length_upload <- length %>% # This includes only EM data, not generic length
#   mutate(family = if_else(genus %in% "Anoplocapros", "Aracanidae", family)) %>%
#   mutate(family = if_else(genus %in% "Aracana", "Aracanidae", family)) %>% 
#   mutate(genus = if_else(genus %in% "Dasyatis", "Bathytoshia", genus)) %>%
#   #mutate(genus = if_else(genus %in% "Dasyatis", "Bathytoshia", genus)) %>% #turn this on if you need it 
#   dplyr::left_join(., CheckEM::aus_synonyms) %>%
#   dplyr::mutate(genus = ifelse(!genus_correct%in%c(NA), genus_correct, genus)) %>%
#   dplyr::mutate(species = ifelse(!is.na(species_correct), species_correct, species)) %>%
#   dplyr::mutate(family = ifelse(!is.na(family_correct), family_correct, family)) %>%
#   dplyr::select(-c(family_correct, genus_correct, species_correct)) %>%
#   dplyr::mutate(scientific = paste(family, genus, species)) %>% # turn this on if you need it 
#   left_join(codes) %>%
#   dplyr::rename(rms_mm = rms, range_mm = range, precision_mm = precision, count = number) %>%
#   dplyr::select(-period) %>%
#   # filter(!is.na(caab_code)) %>%
#   select(campaignid, opcode, sample, caab_code, count, family, genus, species, length_mm, precision_mm, range_mm, rms_mm, stage)
# 
# names(length_upload) %>% sort()
# 
# # Check missing caab codes (okay if from sp, sp1 etc.)
# n_no_caab_count <- length(count_upload$species[which(is.na(count_upload$caab_code))])
# sp_no_caab_count <- unique(count_upload$species[which(is.na(count_upload$caab_code))])
# n_no_caab_length <- length(length_upload$species[which(is.na(length_upload$caab_code))])
# sp_no_caab_length <- unique(length_upload$species[which(is.na(length_upload$caab_code))])
# 
# species_2_update <- count_upload %>%
#   filter(is.na(caab_code)) %>%
#   distinct(campaignid, opcode, family, genus, species, caab_code)
# 
# # Save GA upload data ----
write_csv(count_upload, paste0("data/uploads/BOSS/", name, "_count.csv"))

mismatches <- dplyr::anti_join(count_upload, metadata, by = c("campaignid", "period"))
mismatches %>% dplyr::distinct(campaignid, period)  
head(mismatches)

class(metadata$period)
class(count_upload$period)
# make sure neither has trailing decimals when printed
head(metadata$period, 20)
head(count_upload$period, 20)

metadata_upload <- metadata %>% dplyr::mutate(period = as.integer(period))
count_upload_final <- count_upload %>% dplyr::mutate(period = as.integer(period))

metadata %>%
  dplyr::filter(is.na(suppressWarnings(as.integer(period)))) %>%
  dplyr::select(campaignid, period)
metadata %>%
  dplyr::filter(period == "" | grepl("^\\s|\\s$", period) | is.na(period)) %>%
  dplyr::select(campaignid, period)

dplyr::anti_join(count_upload, metadata, by = c("campaignid", "period")) %>%
  dplyr::distinct(campaignid, period)

# Re-write fresh, right now, from the confirmed-clean objects
write_csv(metadata, "data/uploads/BOSS/geographe-marine-park_metadata.csv")
write_csv(count_upload, "data/uploads/BOSS/geographe-marine-park_count.csv")

# Then read it back in as a completely independent check - not from the R objects in memory
check_meta  <- read_csv("data/uploads/BOSS/geographe-marine-park_metadata.csv")
check_count <- read_csv("data/uploads/BOSS/geographe-marine-park_count.csv")

names(check_meta)   # confirm "period" is really the header, not "sample"



# shows the raw character codes so you can spot the difference
# write_csv(length_upload, paste0("data/uploads/", name, "_length.csv"))
# 
# 
# # Read in AusBathyTopo bathymetry
# 
# bathy <- terra::rast(
#   
#   "data/spatial/rasters/AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
#   
# )
# 
# names(bathy) <- "bathy_depth_m"
# 
# # Convert metadata to spatial points
# 
# metadata_sf <- st_as_sf(
#   
#   metadata,
#   
#   coords = c("longitude_dd", "latitude_dd"),
#   
#   crs = 4326,
#   
#   remove = FALSE
#   
# )
# 
# # Reproject points to bathymetry CRS if needed
# 
# metadata_vect <- terra::vect(metadata_sf)
# 
# if (!terra::same.crs(metadata_vect, bathy)) {
#   
#   metadata_vect <- terra::project(metadata_vect, terra::crs(bathy))
#   
# }
# 
# # Extract bathymetry value
# 
#   
#   dplyr::select(bathy_depth_m)
# 
# # Combine with metadata and replace missing depth_m
# 
# metadata_bathy <- bind_cols(metadata, bathy_values) %>%
#   
#   mutate(
#     
#     depth_m_original = depth_m,
#     
#     depth_m = if_else(
#       
#       is.na(depth_m) | depth_m == "",
#       
#       as.character(bathy_depth_m),
#       
#       as.character(depth_m)
#       
#     ))
# 
# glimpse(metadata_bathy)
# 
# maxn_check <- count_upload %>%
#   group_by(campaignid, opcode, family, genus, species, stage) %>%
#   dplyr::summarise(
#     maxn_total = sum(count, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Summarise number of length measurements
# length_check <- length_upload %>%
#   group_by(campaignid, opcode, family, genus, species, stage) %>%
#   dplyr::summarise(
#     n_lengths = sum(count),
#     length_count_total = sum(count, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Compare lengths against MaxN
# comparison <- maxn_check %>%
#   full_join(length_check,
#             by = c("campaignid", "opcode", "family", "genus", "species", "stage")) %>%
#   mutate(
#     maxn_total = replace_na(maxn_total, 0),
#     n_lengths = replace_na(n_lengths, 0),
#     length_count_total = replace_na(length_count_total, 0),
#     excess_lengths = n_lengths - maxn_total,
#     excess_count = length_count_total - maxn_total
#   )
# 
# # Rows where there are more lengths than MaxN
# problems <- comparison %>%
#   filter(n_lengths > maxn_total |
#            length_count_total > maxn_total |
#            n_lengths < maxn_total |
#            length_count_total < maxn_total) %>%
#   left_join(metadata) %>%
#   filter(successful_length %in% 'Yes')
# 
# problems
# 
# points %>%
#   count(longitude, latitude) %>%
#   filter(n > 1)
