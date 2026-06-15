library(readxl)
library(dplyr)

metadataold <- read.csv("data/raw/2007-2014-Geographe-stereo-BRUVs.checked.metadata.csv")
newmetadata <- read.csv("data/raw/2014-12_Geographe-bay_stereo-BRUVs_metadata.csv")

sp1_clean <- metadataold %>%
  filter(campaignid != "2007-03_Capes.MF_stereoBRUVs")

unique(sp1_clean$campaignid)

newmetadata2 <- newmetadata %>%
  rename(sample = opcode)

missing_in_new <- anti_join(
  sp1_clean,
  newmetadata2,
  by = "sample"
)

missing_in_new %>%
  select(sample)

missing_in_old <- anti_join(
  newmetadata2,
  sp1_clean,
  by = "sample"
)

missing_in_old %>%
  select(sample)
