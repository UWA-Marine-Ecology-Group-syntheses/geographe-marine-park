###
# Project: Australia BRUV synthesis 2024/25
# Data:    Squidle + Benthic annotations
# Task:    Format Squidle + data for import into GA v2
# author:  Claude Spencer
# date:    July 2025
##

# Clear obejcts from the environment
rm(list = ls())
# install.packages('remotes')
# library('remotes')
# options(timeout=9999999)
# 
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
# Load libraries
library(tidyverse)
library(CheckEM)
# library(SQAPI)


# Load the schema file to get caab codes ----
## OPTION 1. Manually download the schema file ----

schema <- read.csv("data/raw/squidle-schema_20250325.csv") %>%
  dplyr::mutate(caab_code = as.character(vocab_registry.caab)) %>%
  dplyr::select(uuid, caab_code) %>%
  glimpse()

raw_benthos <- read.csv("data/raw/2024-04_Geographe_BOSS_habitat.csv")

benthos <- raw_benthos %>%
  clean_names() %>%
  dplyr::rename(campaignid = any_of(c("point_media_deployment_campaign_key",
                                      "point_media_deployment_campaign_name")),
                sample = any_of(c("point_media_deployment_key",
                                  "point_media_deployment_name")),
                uuid = label_uuid) %>%
  left_join(schema) %>%
  dplyr::mutate(caab_code = if_else(label_lineage_names %in% c("Unscorable", "3 Field of View > Open"), "00000001", caab_code)) %>%
  left_join(CheckEM::catami) 

# Check here that they are all open water - if they are not you need to talk to Brooke about how to fix

benthos_missing <- benthos %>%
  filter(caab_code %in% NA) %>%
  distinct(label_lineage_names)

benthos_missing <- benthos %>%
  filter(level_2 %in% NA) %>%
  distinct(label_lineage_names)

benthos_clean <- benthos %>%
  dplyr::mutate(count = 1) %>%
  group_by(campaignid, sample, caab_code, across(starts_with("level")), family, genus, species) %>%
  dplyr::summarise(count = sum(count)) %>%
  ungroup() %>%
  glimpse()

write.csv(benthos_clean, file = paste0("data/tidy/",
                                 unique(benthos_clean$campaignid),
                                 "_benthos-count.csv"), # Change here for relief
          row.names = F)





raw_relief <- read.csv("data/raw/2024-04_Geographe_BOSS_relief.csv")


relief <- raw_relief %>%
  clean_names() %>%
  dplyr::rename(campaignid = any_of(c("point_media_deployment_campaign_key",
                                      "point_media_deployment_campaign_name")),
                sample = any_of(c("point_media_deployment_key",
                                  "point_media_deployment_name")),
                uuid = label_uuid) %>%
  left_join(schema) %>%
  dplyr::mutate(caab_code = if_else(label_lineage_names %in% c("Unscorable", "3 Field of View > Open"), "00000001", caab_code)) %>%
  left_join(CheckEM::catami) 

# Check here that they are all open water - if they are not you need to talk to Brooke about how to fix

relief_missing <- relief %>%
  filter(caab_code %in% NA) %>%
  distinct(label_lineage_names)


relief_clean <- relief %>%
  dplyr::mutate(count = 1) %>%
  group_by(campaignid, sample, caab_code, across(starts_with("level"))) %>%
  dplyr::summarise(count = sum(count)) %>%
  ungroup() %>%
  glimpse()

write.csv(relief_clean, file = paste0("data/tidy/",
                                       unique(relief_clean$campaignid),
                                       "_benthos-relief.csv"), # Change here for relief
          row.names = F)


