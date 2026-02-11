## Project: Effects of grass–endophyte symbiosis and herbivory on population demography across climatic and geographic gradients
## Purpose: Create variables that most accurately reflect the climate across study area.
## Authors: Jacob Moutouama
## Date last modified: 2024-08-03

## Clear workspace and load packages
rm(list = ls())
## Load reauire package 
library(tidyverse)      # data wrangling and plotting
library(prism)          # download and process PRISM climate data
library(terra)          # raster processing alternative
library(raster)         # raster processing
library(stringr)        # string manipulation
library(magrittr)       # pipes and functional programming
library(readxl)         # read Excel files
library(ggsci)          # color palettes for ggplot2
library(lubridate)      # date manipulation

## Load and compile PRISM raster data
prism_set_dl_dir("/Users/jacobmoutouama/Documents/prism")
# Download PRISM monthly precipitation (ppt) and mean temperature (tmean) at 800 m (commented out as already done)
# years <- 1990:2025
# months <- 1:12
# 
# # Precipitation
# get_prism_monthlys(
#   type = "ppt",
#   years = years,
#   mon = months,
#   resolution = "800m",
#   keepZip = TRUE
# )
# 
# # Mean temperature
# get_prism_monthlys(
#   type = "tmean",
#   years = years,
#   mon = months,
#   resolution = "800m",
#   keepZip = TRUE
# )

climate_data <- prism_archive_ls() %>%
  pd_stack()

## Load and format garden site data
read.csv(
  "https://www.dropbox.com/scl/fi/si346imz380lpdgo9yekr/Study_site.csv?rlkey=yiue42npkzzu9w8dggjr4fent&dl=1",
  stringsAsFactors = FALSE
) %>%
  arrange(latitude) -> garden

garden_sites <- as.data.frame(garden)
coordinates(garden_sites) <- c("longitude", "latitude")
# Assign lon/lat CRS explicitly (matches PRISM)
crs(garden_sites) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

##  Extract climate values at garden sites
climate_garden <- data.frame(
  coordinates(garden_sites),
  garden_sites$site_code,
  raster::extract(climate_data, garden_sites)
)

## Reshape climate data
climate_garden <- climate_garden %>%
  gather(date, value, 4:ncol(climate_garden))

# Clean column headers
climate_garden$layer <- gsub("^prism_|_us_30s_", "", climate_garden$date)

# Split header into climate type, year, month
climate_garden <- climate_garden %>%
  mutate(
    variable = str_extract(layer, "^[a-zA-Z]+"),
    ym       = str_extract(layer, "[0-9]{6}"),
    year     = as.integer(substr(ym, 1, 4)),
    month    = as.integer(substr(ym, 5, 6))
  )
# Reshape wide for each climate variable
climate_garden_1995_2025 <- climate_garden %>%
  unique() %>%                                # remove duplicate rows
  pivot_wider(
    names_from = variable,                     # use the 'variable' column for new columns
    values_from = value                        # fill with 'value'
  ) %>%
  rename(
    lon = longitude,
    lat = latitude,
    site = garden_sites.site_code
  ) %>%
  mutate(
    year = as.integer(year),
    month = as.integer(month)
  ) %>%
  filter(year >= 1995, year <= 2025)

## Inspect climate data
#summary(climate_garden_1995_2025)
# table(climate_garden_1995_2025$year, climate_garden_1995_2025$month)

## Subset climate for 2023–2025 period
climate_garden_2023_2025 <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date > as.Date("2023-05-01") & date < as.Date("2025-06-01"))

# climate_garden_2023_2025 %>%
#   filter(site=="KER" & year=="2025")

## Compute site-level summary
climate_garden_2023_2025 %>%
  # create a year column if not already present
  mutate(year = year(date)) %>%
  group_by(site, year) %>%
  summarise(
    sum_ppt = sum(ppt, na.rm = TRUE),
    mean_temp = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  ) -> prism_means_per_year

prism_means <- as.data.frame(prism_means_per_year)
#saveRDS(prism_means, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/prism_means.rds")

# Helper function to clean Tag_ID
clean_tag <- function(df) {
  df %>%
    mutate(
      Tag_ID = as.character(Tag_ID),        # ensure character
      Tag_ID = str_trim(Tag_ID),            # trim leading/trailing spaces
      Tag_ID = str_squish(Tag_ID)           # remove extra spaces inside
    )
}

## Load census data
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1")
dat23 <- read.csv("https://www.dropbox.com/scl/fi/fkwm0dan6nx2eaeyxjrjw/census2023.csv?rlkey=hy9209t53j9n7vxhta7axl5jk&dl=1")

# Ensure Tag_ID is numeric in both tables
datini <-clean_tag(datini)
dat23 <-clean_tag(dat23) 

datini23 <- right_join(x = datini, y = dat23, by = "Tag_ID")
date_sp_site23 <- datini23 %>%
  dplyr::select(Site, Plot, Species, date_23) %>%
  distinct() %>%
  filter(!is.na(Site) & !is.na(Species) & !is.na(Plot))

dat24 <- read.csv("https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1")
dat24ini<-read.csv("https://www.dropbox.com/scl/fi/zncghunh1p9ull9j1jhkp/data_ini_2024.csv?rlkey=fkhbp3sdvb65va0gg4rwp4ra4&dl=1")
dat24 <- clean_tag(dat24)
#unique(dat24$date_24)
#which(dat24$date_24 == "CO+IR 2-24-1 E+")
#which(dat24$date_24 == "2026-06-24")
#which(dat24$date_24 == "2025-05-23")
dat24ini <-clean_tag(dat24ini)
dat24ini <- dat24ini %>%
  mutate(
    date_24 = mdy(date_24)  # converts "2-26-2024" to "2024-02-26"
  )
#unique(dat24ini$date_24)
# Combine initial and 2023 census
combined_data <- bind_rows(
  datini[,c("Site","Species","Plot","Tag_ID","Population","Endo")],
  datini23[,c("Site","Species","Plot","Tag_ID","Population","Endo")]
) %>%
  distinct(Tag_ID, .keep_all = TRUE)

dat24_sp_site_tag <- dat24 %>%
  # join combined_data for extra columns including Site, Plot, Species, Population, Endo
  left_join(
    combined_data %>% dplyr::select(Tag_ID, Site, Plot, Species, Population, Endo),
    by = "Tag_ID"
  ) %>%
  # join initial 2024 plantings to fill any missing metadata
  left_join(
    dat24ini %>% dplyr::select(Tag_ID, Site, Species, Plot),
    by = "Tag_ID",
    suffix = c("", "_ini")
  ) %>%
  # coalesce to fill missing Site, Species, Plot
  mutate(
    Site = coalesce(Site, Site_ini),
    Species = coalesce(Species, Species_ini),
    Plot = coalesce(Plot, Plot_ini)
  ) %>%
  # remove unwanted measurement columns and the temporary "_ini" columns
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit",
                          "attachedInf_24","brokenInf_24",
                          "Site_ini","Species_ini","Plot_ini")))
#unique(dat24_sp_site_tag$date_24)
dat24_sp_site_tag %>%
  dplyr::select(Site,Plot, Species, date_24) %>%
  distinct() %>%
  na.omit() -> date_sp_site24
#unique(date_sp_site24$date_24)
dat25 <- read.csv("https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1")
# Ensure Tag_ID is character and trim whitespace
dat25 <- clean_tag(dat25)
#unique(dat25$date_25)

dat25_plot <- dat25 %>%
  left_join(
    dat24_sp_site_tag %>% dplyr::select(Tag_ID, Site, Species, Plot),
    by = "Tag_ID"
  ) %>%
  mutate(
    Site    = coalesce(Site.x, Site.y),
    Species = coalesce(Species.x, Species.y),
    Plot    = Plot  # <-- this will fill Plot from dat24
  ) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y"))

# Check duplicates
dat25_plot %>% count(Tag_ID) %>% filter(n > 1)
#unique(dat25_plot$date_25)

dat25_plot %>%
  dplyr::select(Site,Plot,Species, date_25) %>%
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
    start_year  = lubridate::year(start_date),
    start_month = lubridate::month(start_date),
    end_year    = lubridate::year(end_date),
    end_month   = lubridate::month(end_date)
  )
#view(census_dates_23_24)
# head(census_dates_23_24)

## Compute cumulative precipitation and mean temperature per species × site for 2023–2024
#  Summarise main 2024 data by Site, Plot, Species
date_sp_site24_unique <- date_sp_site24 %>%
  group_by(Site, Plot, Species) %>%
  summarise(date_24_main = min(as.Date(date_24)), .groups = "drop")

# Summarise 2025 data
date_sp_site25_unique <- date_sp_site25 %>%
  group_by(Site, Plot, Species) %>%
  summarise(date_25 = max(as.Date(date_25)), .groups = "drop")

# Combine 2024 and 2025 for census dates
census_dates_24_25 <- date_sp_site24_unique %>%
  full_join(date_sp_site25_unique, by = c("Site", "Species")) %>%
  rename(start_date = date_24, end_date = date_25) %>%
  mutate(
    # If 2025 date missing, leave as NA (don't backfill)
    census_year = year(end_date),
    start_year  = year(start_date),
    start_month = month(start_date),
    end_year    = year(end_date),
    end_month   = month(end_date)
  )

# Check
census_dates_24_25
view(census_dates_24_25)
# COLLAPSE CENSUS WINDOWS TO SITE × SPECIES
# Ensure census_site_species tables assign census_year = end_year
census_site_species_23_24 <- census_dates_23_24 %>%
  group_by(Site, Species) %>%
  summarise(
    start_year  = min(start_year),
    start_month = min(start_month),
    end_year    = max(end_year),
    end_month   = max(end_month),
    census_year = max(end_year),  # assign census_year to the end year
    .groups = "drop"
  )

census_site_species_24_25 <- census_dates_24_25 %>%
  group_by(Site, Species) %>%
  summarise(
    start_year  = min(start_year),
    start_month = min(start_month),
    end_year    = max(end_year),
    end_month   = max(end_month),
    census_year = max(end_year),
    .groups = "drop"
  )

# ---- CLIMATE DATA ----
climate_monthly <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

# ---- 2023–2024 ----
climate_2023_2024 <- climate_monthly %>%
  inner_join(census_site_species_23_24, by = c("site" = "Site"), relationship = "many-to-many") %>%
  # Keep only months inside each site × species census window
  filter(
    (year > start_year | (year == start_year & month >= start_month)) &
      (year < end_year   | (year == end_year   & month <= end_month))
  ) %>%
  group_by(Species, site, census_year) %>%
  summarise(
    cum_ppt    = sum(ppt, na.rm = TRUE),
    mean_tmean = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 2024–2025 ----
climate_2024_2025 <- climate_monthly %>%
  inner_join(census_site_species_24_25, by = c("site" = "Site"), relationship = "many-to-many") %>%
  filter(
    (year > start_year | (year == start_year & month >= start_month)) &
      (year < end_year   | (year == end_year   & month <= end_month))
  ) %>%
  group_by(Species, site, census_year) %>%
  summarise(
    cum_ppt    = sum(ppt, na.rm = TRUE),
    mean_tmean = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

# ---- COMBINE ALL YEARS ----
climate_census_years <- bind_rows(climate_2023_2024, climate_2024_2025)

# Check unique census years
unique(climate_census_years$census_year)

# view(climate_census_years)
# saveRDS(climate_census_years, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/climate_census_years.rds")

# Identify the 30-year easternmost point
#Load species coordinates
max_long_table <- readRDS(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/max_longitudes.rds"
)

max_long_yr <- max_long_table 

coordinates(max_long_yr) <- c("longitude", "latitude")

climate_sites_max <- data.frame(
  coordinates(max_long_yr),
  species = max_long_yr$species,
  extract(climate_data, max_long_yr)
)

climate_sites_max <- climate_sites_max %>%
  gather(date, value, 4:ncol(.))

climate_sites_max$layer <- gsub("^prism_|_us_30s_", "", climate_sites_max$date)

climate_sites_max <- climate_sites_max %>%
  mutate(
    variable = str_extract(layer, "^[a-zA-Z]+"),
    ym       = str_extract(layer, "[0-9]{6}"),
    year     = as.integer(substr(ym, 1, 4)),
    month    = as.integer(substr(ym, 5, 6))
  )
climate_sites_max <- climate_sites_max %>%
  distinct() %>%
  pivot_wider(
    names_from  = variable,
    values_from = value
  ) %>%
  rename(
    lon = longitude,
    lat = latitude
  ) %>%
  mutate(
    year  = as.integer(year),
    month = as.integer(month)
  )

prism_edge_yr_means <- climate_sites_max %>%
  group_by(species, year) %>%
  summarise(
    annual_ppt = sum(ppt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  summarise(
    mean_annual_ppt = mean(annual_ppt, na.rm = TRUE),
    .groups = "drop"
  )


saveRDS(prism_edge_yr_means, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/prism_edge_yr_means.rds")

