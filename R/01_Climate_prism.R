## Project: Aridity and herbivory shape demographic benefits of fungal endophytes in cool-season grasses
## Purpose: Create variables that most accurately reflect the climate across study area.
## Authors: Jacob K. Moutouama
## Date last modified (Y-M-D): 2026-11-17

## Clear workspace and load packages
rm(list = ls())
source(here::here("R", "config.R"))   # shared paths + data_path() helper
library(tidyverse)
library(prism)
library(raster)
library(stringr)
library(magrittr)
library(readxl)
library(ggsci)
library(lubridate)


## Load PRISM climate data
prism_set_dl_dir(PRISM_DIR)   # default: Data/prism (override with env var PRISM_DIR)

## Download the PRISM monthly rasters this script needs (skipped if already present).
## What is needed: monthly precipitation ("ppt") and mean temperature ("tmean"),
## 800 m resolution, 1995-2025, all 12 months (~744 rasters; several GB, can take hours).
## PRISM throttles repeated downloads; if it stops part-way just re-run this block.
if (!exists("DOWNLOAD_PRISM")) DOWNLOAD_PRISM <- TRUE
if (DOWNLOAD_PRISM) {
  for (v in c("ppt", "tmean")) {
    get_prism_monthlys(type = v, years = 1995:2025, mon = 1:12,
                       resolution = "800m", keepZip = FALSE)
  }
}
prism_rasters <- prism_archive_ls() %>% pd_stack()

## Load and format garden site data
garden_sites_raw <- read.csv(
  data_path("Study_site.csv"),
  stringsAsFactors = FALSE
) %>% arrange(latitude)

garden_sites_sp <- as.data.frame(garden_sites_raw)
coordinates(garden_sites_sp) <- c("longitude", "latitude")
crs(garden_sites_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## Extract climate values at garden sites
climate_garden_long <- data.frame(
  coordinates(garden_sites_sp),
  garden_sites_sp$site_code,
  raster::extract(prism_rasters, garden_sites_sp)
) %>%
  gather(date, value, 4:ncol(.)) %>%
  mutate(
    layer = gsub("^prism_|_us_30s_", "", date),
    variable = str_extract(layer, "^[a-zA-Z]+"),
    ym = str_extract(layer, "[0-9]{6}"),
    year = as.integer(substr(ym, 1, 4)),
    month = as.integer(substr(ym, 5, 6))
  )

climate_garden_wide_1995_2025 <- climate_garden_long %>%
  distinct() %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  rename(
    lon = longitude,
    lat = latitude,
    site = garden_sites_sp.site_code
  ) %>%
  mutate(year = as.integer(year), month = as.integer(month)) %>%
  filter(year >= 1995 & year <= 2025)

climate_garden_filtered_2023_2025 <- climate_garden_wide_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date > as.Date("2023-05-01") & date < as.Date("2025-06-01"))

## Helper function for Tag_ID cleanup
clean_tag <- function(df) {
  df %>%
    mutate(
      Tag_ID = as.character(Tag_ID) %>% str_trim() %>% str_squish()
    )
}

##  Load census data
census_initial_plantings <- read.csv(
  data_path("Initialdata.csv")
) %>% clean_tag()

census_2023_raw <- read.csv(
  data_path("census2023.csv")
) %>% clean_tag()

# Fix D_LE4_1 stale planting date — correct to actual 2023 census date for LAF Plot 4
census_2023_raw <- census_2023_raw %>%
  mutate(date_23 = if_else(Tag_ID == "D_LE4_1", "2023-06-13", date_23))

census_2023_metadata <- right_join(
  x = census_initial_plantings,
  y = census_2023_raw,
  by = "Tag_ID"
)
census_2023_unique <- census_2023_metadata %>%
  dplyr::select(Site, Plot, Species, date_23) %>%
  distinct() %>%
  filter(!is.na(Site) & !is.na(Species) & !is.na(Plot))



census_2024_raw <- read.csv(
  data_path("census2024.csv")
) %>% clean_tag()

census_2024_initial <- read.csv(
  data_path("data_ini_2024.csv")
) %>% clean_tag() %>%
  mutate(date_24 = mdy(date_24))

metadata_all_previous <- bind_rows(
  census_initial_plantings %>% dplyr::select(Site, Species, Plot, Tag_ID, Population, Endo),
  census_2023_metadata %>% dplyr::select(Site, Species, Plot, Tag_ID, Population, Endo)
) %>% distinct(Tag_ID, .keep_all = TRUE)

census_2024_with_metadata <- census_2024_raw %>%
  left_join(metadata_all_previous %>% dplyr::select(Tag_ID, Site, Plot, Species, Population, Endo),
            by = "Tag_ID") %>%
  dplyr::select(-any_of(c("Spikelet_A","Spikelet_B","Spikelet_C","digit",
                   "attachedInf_24","brokenInf_24")))
# census_2024_with_metadata %>% 
#   filter(Site=="KER")

census_2024_full <- census_2024_raw %>%
  left_join(metadata_all_previous %>% dplyr::select(Tag_ID, Site, Plot, Species, Population, Endo),
            by = "Tag_ID") %>%
  left_join(census_2024_initial %>% dplyr::select(Tag_ID, Site, Species, Plot),
            by = "Tag_ID", suffix = c("", "_ini")) %>%
  mutate(
    Site = coalesce(Site, Site_ini),
    Species = coalesce(Species, Species_ini),
    Plot = coalesce(Plot, Plot_ini)
  ) %>%
  dplyr::select(-any_of(c("Spikelet_A","Spikelet_B","Spikelet_C","digit",
                   "attachedInf_24","brokenInf_24",
                   "Site_ini","Species_ini","Plot_ini")))

census_2025_raw <- read.csv(
  data_path("census_2025.csv")
) %>% clean_tag()

census_2025_with_metadata <- census_2025_raw %>%
  left_join(census_2024_full %>% dplyr::select(Tag_ID, Site, Species, Plot), by = "Tag_ID") %>%
  mutate(
    Site = coalesce(Site.x, Site.y),
    Species = coalesce(Species.x, Species.y),
    Plot = Plot
  ) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y"))

##  Census dates & windows
census_dates_2023 <- census_2023_unique %>%
  mutate(date_23 = as.Date(date_23)) %>%
  filter(lubridate::year(date_23) == 2023) %>%
  group_by(Site, Species) %>%
  summarise(date_23 = min(date_23), .groups = "drop")


census_dates_2024 <- census_2024_with_metadata %>%
  dplyr::select(Site, Plot, Species, date_24) %>%
  distinct() %>%
  na.omit() %>%
  group_by(Site, Species) %>%
  summarise(date_24 = max(as.Date(date_24)), .groups = "drop")

census_dates_2024_full <- census_2024_full %>%
  dplyr::select(Site, Plot, Species, date_24) %>%
  distinct() %>%
  na.omit() %>%
  group_by(Site, Species) %>%
  summarise(date_24 = min(as.Date(date_24)), .groups = "drop")

census_dates_2025 <- census_2025_with_metadata %>%
  dplyr::select(Site, Plot, Species, date_25) %>%
  distinct() %>%
  na.omit() %>%
  group_by(Site, Species) %>%
  summarise(date_25 = max(as.Date(date_25)), .groups = "drop")

census_windows_2023_2024 <- census_dates_2023 %>%
  left_join(census_dates_2024, by = c("Site","Species")) %>%
  rename(start_date = date_23, end_date = date_24) %>%
  mutate(
    start_year  = year(start_date),
    start_month = month(start_date),
    end_year    = year(end_date),
    end_month   = month(end_date)
  )

census_windows_2024_2025 <- census_dates_2024_full %>%
  full_join(census_dates_2025, by = c("Site","Species")) %>%
  rename(start_date = date_24, end_date = date_25) %>%
  mutate(
    census_year = year(end_date),
    start_year  = year(start_date),
    start_month = month(start_date),
    end_year    = year(end_date),
    end_month   = month(end_date)
  )

census_windows_by_site_species_23_24 <- census_windows_2023_2024 %>%
  group_by(Site, Species) %>%
  summarise(
    start_year  = min(start_year),
    start_month = min(start_month),
    end_year    = max(end_year),
    end_month   = max(end_month),
    census_year = max(end_year),
    .groups = "drop"
  )

# Copy SON ELVI values to KER AGHY and KER POAU manually
census_windows_by_site_species_23_24 <- census_windows_by_site_species_23_24 %>%
  mutate(
    end_year  = ifelse(Site == "KER" & is.na(end_year), 2024, end_year),
    end_month = ifelse(Site == "KER" & is.na(end_month), 6, end_month),
    census_year = ifelse(Site == "KER" & is.na(census_year), 2024, census_year)
  )

census_windows_by_site_species_24_25 <- census_windows_2024_2025 %>%
  group_by(Site, Species) %>%
  summarise(
    start_year  = min(start_year),
    start_month = min(start_month),
    end_year    = max(end_year),
    end_month   = max(end_month),
    census_year = max(end_year),
    .groups = "drop"
  )

##  Climate summaries by census window
climate_monthly_all <- climate_garden_wide_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

climate_census_window_2023_2024 <- climate_monthly_all %>%
  cross_join(census_windows_by_site_species_23_24 %>% distinct(Species)) %>%
  inner_join(census_windows_by_site_species_23_24, by = c("site" = "Site", "Species" = "Species")) %>%
  filter(
    (year > start_year | (year == start_year & month >= start_month)) &
      (year < end_year   | (year == end_year & month <= end_month))
  ) %>%
  group_by(Species, site, census_year) %>%
  summarise(
    cum_ppt = sum(ppt, na.rm = TRUE),
    mean_tmean = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

climate_census_window_2024_2025 <- climate_monthly_all %>%
  cross_join(census_windows_by_site_species_24_25 %>% distinct(Species)) %>%
  inner_join(census_windows_by_site_species_24_25, by = c("site" = "Site", "Species" = "Species")) %>%
  filter(
    (year > start_year | (year == start_year & month >= start_month)) &
      (year < end_year   | (year == end_year & month <= end_month))
  ) %>%
  group_by(Species, site, census_year) %>%
  summarise(
    cum_ppt = sum(ppt, na.rm = TRUE),
    mean_tmean = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

climate_census_year_summary <- bind_rows(
  climate_census_window_2023_2024,
  climate_census_window_2024_2025
)

# census_2023_unique %>%
#   mutate(date_23 = as.Date(date_23)) %>%
#   filter(lubridate::year(date_23) != 2023)
# 
# census_windows_by_site_species_23_24 %>%
#   filter(Site == "LAF") %>%
#   arrange(Species)
# 
# climate_census_year_summary %>%
#   filter(Species == "ELVI", census_year == 2024) %>%
#   arrange(desc(cum_ppt))
# 
# climate_census_year_summary %>%
#   group_by(Species, census_year) %>%
#   summarise(
#     min_ppt = min(cum_ppt),
#     max_ppt = max(cum_ppt),
#     n       = n(),
#     .groups = "drop"
#   ) %>%
#   print(n = Inf)

saveRDS(climate_census_year_summary, file.path(DATA_DIR, "climate_census_years.rds"))

##  Easternmost edge climate
species_edge_coords <- readRDS(
  file.path(DATA_DIR, "max_longitudes.rds")
)

species_edge_coords_sp <- species_edge_coords
coordinates(species_edge_coords_sp) <- c("longitude","latitude")

climate_species_edge_long <- data.frame(
  coordinates(species_edge_coords_sp),
  species = species_edge_coords_sp$species,
  raster::extract(prism_rasters, species_edge_coords_sp)
) %>%
  gather(date, value, 4:ncol(.)) %>%
  mutate(
    layer = gsub("^prism_|_us_30s_", "", date),
    variable = str_extract(layer, "^[a-zA-Z]+"),
    ym = str_extract(layer, "[0-9]{6}"),
    year = as.integer(substr(ym,1,4)),
    month = as.integer(substr(ym,5,6))
  ) %>%
  distinct() %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  rename(lon = longitude, lat = latitude) %>%
  mutate(year = as.integer(year), month = as.integer(month))

climate_species_edge_annual_means <- climate_species_edge_long %>%
  group_by(species, year) %>%
  summarise(annual_ppt = sum(ppt, na.rm = TRUE), .groups = "drop") %>%
  group_by(species) %>%
  summarise(mean_annual_ppt = mean(annual_ppt, na.rm = TRUE), .groups = "drop")

saveRDS(climate_species_edge_annual_means,
        file.path(DATA_DIR, "prism_edge_yr_means.rds"))
