## Project: Effects of grass–endophyte symbiosis and herbivory on population demography across climatic and geographic gradients
## Purpose: Create variables that most accurately reflect the climate across study area.
## Authors: Jacob Moutouama
## Date last modified: 2024-08-03


## Clear workspace and load packages
rm(list = ls())

## Load reauire package 
library(tidyverse)      # data wrangling and plotting
library(prism)          # download and process PRISM climate data
library(raster)         # raster processing
library(terra)          # raster processing alternative
library(stringr)        # string manipulation
library(magrittr)       # pipes and functional programming
library(readxl)         # read Excel files
library(ggsci)          # color palettes for ggplot2
library(lubridate)      # date manipulation

## Load and compile PRISM raster data
options(prism.path = "/Users/jacobmoutouama/Documents/prism/")
# Download PRISM data (commented out as already done)
# get_prism_monthlys(type="ppt", years=1990:2025, mon=1:12, keepZip = TRUE)
# get_prism_monthlys(type="tmean", years=1990:2025, mon=1:12, keepZip = TRUE)

climate_data <- prism_archive_ls() %>%
  pd_stack()

climate_crs <- climate_data@crs@projargs

## Load and format garden site data
read.csv(
  "https://www.dropbox.com/scl/fi/si346imz380lpdgo9yekr/Study_site.csv?rlkey=yiue42npkzzu9w8dggjr4fent&dl=1",
  stringsAsFactors = FALSE
) %>%
  arrange(latitude) -> garden

garden_sites <- as.data.frame(garden)
coordinates(garden_sites) <- c("longitude", "latitude")
proj4string(garden_sites) <- CRS(climate_crs)

##  Extract climate values at garden sites
climate_garden <- data.frame(
  coordinates(garden_sites),
  garden_sites$site_code,
  extract(climate_data, garden_sites)
)

## Reshape climate data
climate_garden <- climate_garden %>%
  gather(date, value, 4:ncol(climate_garden))

# Clean column headers
climate_garden$date <- gsub("PRISM_", "", climate_garden$date) %>%
  gsub("stable_4kmM3_", "", .) %>%
  gsub("provisional_4kmM3_", "", .) %>%
  gsub("_bil", "", .)

# Split header into climate type, year, month
climate_garden <- separate(climate_garden, "date", into = c("clim", "YearMonth"), sep = "_")
climate_garden <- separate(climate_garden, "YearMonth", into = c("year", "month"), sep = 4)

# Reshape wide for each climate variable
climate_garden_1995_2025 <- climate_garden %>%
  unique() %>%
  spread(clim, value) %>%
  rename(lon = longitude, lat = latitude, site = garden_sites.site_code) %>%
  mutate(year = as.numeric(year), month = as.numeric(month)) %>%
  filter(year > 1994)

## Inspect climate data
#summary(climate_garden_1995_2025)
# table(climate_garden_1995_2025$year, climate_garden_1995_2025$month)

## Subset climate for 2023–2025 period
climate_garden_2023_2025 <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date > as.Date("2023-05-01") & date < as.Date("2025-06-01"))

## Compute site-level summary
climate_garden_2023_2025 %>%
  group_by(site) %>%
  summarise(
    sum_ppt = sum(ppt, na.rm = TRUE),
    mean_temp = mean(tmean, na.rm = TRUE),
  ) -> prism_means

prism_means <- as.data.frame(prism_means)
# saveRDS(prism_means, "/Users/jacobmoutouama/Downloads/endo-range-limits/Data/prism_means.rds")

## Load census data
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1")

dat23 <- read.csv("https://www.dropbox.com/scl/fi/fkwm0dan6nx2eaeyxjrjw/census2023.csv?rlkey=hy9209t53j9n7vxhta7axl5jk&dl=1")
datini23 <- right_join(x = datini, y = dat23, by = "Tag_ID")
datini23 %>%
  dplyr::select(Site, Species, date_23) %>%
  distinct() -> date_sp_site23

dat24 <- read.csv("https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1")

# Combine initial and 2023 census
combined_data <- bind_rows(
  datini[,c("Site","Species","Plot","Tag_ID","Population","Endo")],
  datini23[,c("Site","Species","Plot","Tag_ID","Population","Endo")]
) %>%
  distinct(Tag_ID, .keep_all = TRUE)

dat24_sp_site_tag <- left_join(dat24, combined_data, by = "Tag_ID") %>%
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_24","brokenInf_24"))) %>%
  filter(!is.na(Species))

dat24_sp_site_tag %>%
  dplyr::select(Site, Species, date_24) %>%
  distinct() %>%
  na.omit() -> date_sp_site24

dat25 <- read.csv("https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1")
dat25 %>%
  dplyr::select(Site, Species, date_25) %>%
  distinct() %>%
  na.omit() -> date_sp_site25

## Prepare unique census dates for each period
date_sp_site23_unique <- date_sp_site23 %>%
  group_by(Site, Species) %>%
  summarise(date_23 = min(as.Date(date_23)), .groups = "drop")

date_sp_site24_unique <- date_sp_site24 %>%
  group_by(Site, Species) %>%
  summarise(date_24 = max(as.Date(date_24)), .groups = "drop")

census_dates_23_24 <- date_sp_site23_unique %>%
  left_join(date_sp_site24_unique, by = c("Site", "Species")) %>%
  rename(start_date = date_23, end_date = date_24) %>%
  mutate(
    end_date = if_else(is.na(end_date), start_date + years(1), end_date),
    census_year = paste(year(end_date)),
    start_year  = year(start_date),
    start_month = month(start_date),
    end_year    = year(end_date),
    end_month   = month(end_date)
  )
view(census_dates_23_24)
# head(census_dates_23_24)

## Compute cumulative precipitation and mean temperature per species × site for 2023–2024
climate_garden_2023_2024 <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date > as.Date("2023-05-01") & date < as.Date("2024-06-01"))

climate_2023_2024 <- climate_garden_2023_2024 %>%
  inner_join(census_dates_23_24, by = c("site" = "Site"), relationship = "many-to-many") %>%
  filter((year > start_year | (year == start_year & month >= start_month)) &
           (year < end_year | (year == end_year & month <= end_month))) %>%
  group_by(Species, site, census_year) %>%
  summarise(
    cum_ppt    = sum(ppt, na.rm = TRUE),
    mean_tmean = mean(tmean, na.rm = TRUE),
    .groups    = "drop"
  )

## Prepare unique census dates for 2024–2025
date_sp_site24_unique <- date_sp_site24 %>%
  group_by(Site, Species) %>%
  summarise(date_24 = min(as.Date(date_24)), .groups = "drop")

date_sp_site25_unique <- date_sp_site25 %>%
  group_by(Site, Species) %>%
  summarise(date_25 = max(as.Date(date_25)), .groups = "drop")

census_dates_24_25 <- date_sp_site24_unique %>%
  left_join(date_sp_site25_unique, by = c("Site", "Species")) %>%
  rename(start_date = date_24, end_date = date_25) %>%
  mutate(
    end_date = if_else(is.na(end_date), start_date + years(1), end_date),
    census_year = paste(year(end_date)),
    start_year  = year(start_date),
    start_month = month(start_date),
    end_year    = year(end_date),
    end_month   = month(end_date)
  )

head(census_dates_24_25)
## Compute cumulative climate per species × site for 2024–2025
climate_garden_2024_2025 <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date >= as.Date("2024-05-01") & date <= as.Date("2025-06-01"))

climate_2024_2025 <- climate_garden_2024_2025 %>%
  inner_join(census_dates_24_25, by = c("site" = "Site"), relationship = "many-to-many") %>%
  filter((year > start_year | (year == start_year & month >= start_month)) &
           (year < end_year | (year == end_year & month <= end_month))) %>%
  group_by(Species, site, census_year) %>%
  summarise(
    cum_ppt    = sum(ppt, na.rm = TRUE),
    mean_tmean = mean(tmean, na.rm = TRUE),
    .groups    = "drop"
  )

## Combine all census periods into a single dataset
climate_census_years <- bind_rows(climate_2023_2024, climate_2024_2025)
# view(climate_census_years)
# saveRDS(climate_census_years, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/climate_census_years.rds")
