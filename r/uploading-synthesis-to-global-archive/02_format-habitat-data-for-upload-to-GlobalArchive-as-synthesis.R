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

# Set name for synthesis, will be used as prefix for your files to upload
name <- "test-synthesis"

# All data needs to be saved in a folder structure "data/raw/" for the script to work, or change the direcory when reading in files.

# YOU NEED TO COPY THIS SCRIPT FOR EACH TM CAMPAIGN YOU WOULD LIKE TO INCLUDES

# Read in metadata ----
# Note check if you have used opcode, opcode or period and change the below code accordingly.
metadata <- read_csv(paste0("data/uploads/", name, "_metadata.csv")) %>%
  dplyr::filter(campaignid %in% "2022-11_Salisbury_stereo-BRUVs") # change campaign here

# Read in TM data ----
points <- read_TM(here::here("data/raw/"),
                  sample = "opcode") %>%
  dplyr::filter(campaignid %in% "2022-11_Salisbury_stereo-BRUVs") # change campaign here

unique(points$relief_annotated)

habitat <- points %>%
  dplyr::filter(relief_annotated %in% c("no", NA, "No")) %>%
  dplyr::select(campaignid, sample, starts_with("level"), scientific) %>%
  dplyr::rename(opcode = sample) %>%
  glimpse()

relief <- points %>%
  dplyr::filter(relief_annotated %in% c("yes", "Yes")) %>%
  dplyr::select(campaignid, sample, starts_with("level"), scientific) %>%
  dplyr::rename(opcode = sample) %>%
  glimpse()

num.points <- 20

wrong.points.habitat <- habitat %>%
  group_by(campaignid, opcode) %>%
  summarise(points.annotated = n()) %>%
  left_join(metadata) %>%
  dplyr::mutate(expected = case_when(
    successful_habitat_forward %in% "Yes" & successful_habitat_backward %in% "Yes" ~ num.points * 2, 
    successful_habitat_forward %in% "Yes" & successful_habitat_backward %in% "No" ~ num.points * 1, 
    successful_habitat_forward %in% "No" & successful_habitat_backward %in% "Yes" ~ num.points * 1, 
    successful_habitat_forward %in% "No" & successful_habitat_backward %in% "No" ~ num.points * 0)) %>%
  dplyr::filter(!points.annotated == expected) %>%
  glimpse()

wrong.points.relief <- relief %>%
  group_by(campaignid, opcode) %>%
  summarise(points.annotated = n()) %>%
  left_join(metadata) %>%
  dplyr::mutate(expected = case_when(successful_habitat_forward %in% "Yes" & successful_habitat_backward %in% "Yes" ~ num.points * 2,
                                     successful_habitat_forward %in% "Yes" & successful_habitat_backward %in% "No" ~ num.points * 1,
                                     successful_habitat_forward %in% "No" & successful_habitat_backward %in% "Yes" ~ num.points * 1,
                                     successful_habitat_forward %in% "No" & successful_habitat_backward %in% "No" ~ num.points * 0)) %>%
  dplyr::filter(!points.annotated == expected) %>%
  glimpse()

habitat.missing.metadata <- anti_join(habitat, metadata, by = c("campaignid", "opcode")) %>%
  glimpse()

metadata.missing.habitat <- anti_join(metadata, habitat, by = c("campaignid", "opcode")) %>%
  glimpse()

catami_cols <- c("level_1" = NA,
                 "level_2" = NA,
                 "level_3" = NA,
                 "level_4" = NA,
                 "level_5" = NA,
                 "level_6" = NA,
                 "level_7" = NA,
                 "level_8" = NA,
                 "family" = NA,
                 "genus" = NA,
                 "species" = NA)

tidy.habitat <- habitat %>%
  dplyr::mutate(count = 1) %>%
  add_column(!!!catami_cols[!names(catami_cols) %in% names(.)]) %>%
  left_join(catami) %>%
  dplyr::select(campaignid, opcode, count, starts_with("level"), family, genus, species) %>%
  dplyr::filter(!level_2 %in% c("","Unscorable", NA)) %>%
  group_by(campaignid, opcode, across(starts_with("level")), family, genus, species) %>%
  dplyr::summarise(count = sum(count)) %>%
  ungroup() %>%
  dplyr::select(campaignid, opcode, level_1, everything()) %>%
  glimpse()

tidy.relief <- relief %>%
  dplyr::select(campaignid, opcode, starts_with("level"), scientific) %>%
  dplyr::filter(!level_2 %in% c("","Unscorable", NA)) %>%
  dplyr::mutate(count = 1) %>%
  add_column(!!!catami_cols[!names(catami_cols) %in% names(.)]) %>%
  left_join(catami) %>%
  group_by(campaignid, opcode, across(starts_with("level"))) %>%
  dplyr::summarise(count = sum(count)) %>%
  ungroup() %>%
  glimpse()

campaignid <- unique(metadata$campaignid)

write_csv(tidy.habitat, paste0("data/uploads/", campaignid, "_benthos-count.csv"))
write_csv(tidy.relief, paste0("data/uploads/", campaignid, "_benthos-relief.csv"))