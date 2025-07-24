# Project:
# Purpose: Create variables that most accurately reflect the climate across study area.
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 2024-08-03
rm(list = ls())
# load packages
library(tidyverse) # a suite of packages for wrangling and tidying data
library(prism) # package to access and download climate data
library(raster) # the climate data comes in raster files- this package helps process those
library(stringr)
library(magrittr)
library(readxl) # read excel data
library(ggsci) # package for color blind color in ggplot 2
library(corrplot) # visualize the correlation
library(terra)
library(zoo)
library(SPEI)
library(smplot2)
library(Evapotranspiration)
library(lubridate)
# PRISM data ----
# First, set a file path where prism data will be stored
options(prism.path = "/Users/jm200/Documents/Prism Range limit/")
# get_prism_monthlys(type="ppt",years=1990:2025,mon=1:12,keepZip = TRUE)
# get_prism_monthlys(type = "tmean", years = 1990:2025, mon = 1:12, keepZip = TRUE)
# get_prism_monthlys(type = "vpdmin", years = 1990:2025, mon = 1:12, keepZip = TRUE)
# get_prism_monthlys(type = "vpdmax", years = 1990:2025, mon = 1:12, keepZip = TRUE)
# Grab the prism data and compile the files
climate_data <- prism_archive_ls() %>%
  pd_stack(.)
climate_crs <- climate_data@crs@projargs
# Convert these locations to format that can be matched to Prism climate data
read.csv("https://www.dropbox.com/scl/fi/zm29qvqkc1bab8bpbc8wz/Study_site.csv?rlkey=ma4efsaa95tvvuca20pyc4d1f&dl=1", stringsAsFactors = F) %>%
  arrange(latitude) -> garden ## common garden populations
garden_sites <- as.data.frame(garden)
coordinates(garden_sites) <- c("longitude", "latitude")
proj4string(garden_sites) <- CRS(climate_crs)
# Extract climatic data from the raster stack for those sites
climate_garden <- data.frame(
  coordinates(garden_sites),
  garden_sites$site_code,
  extract(climate_data, garden_sites)
)
# Reshape data. Col 1:3 are lat, long, and site ID. Col 4:ncol are climate data
# Column headers include date and climate type info
climate_garden <- climate_garden %>%
  gather(date, value, 4:ncol(climate_garden))
# The column header includes the date and data type, but also some other metadata that we don't need
# Here, I remove the extra info from the column header
climate_garden$date <- gsub("PRISM_", "", climate_garden$date) %>%
  gsub("stable_4kmM3_", "", .) %>%
  gsub("provisional_4kmM3_", "", .) %>%
  gsub("_bil", "", .)

# Split header into type (precipitation or temperature), year, and month
climate_garden <- separate(climate_garden, "date",
  into = c("clim", "YearMonth"),
  sep = "_"
)
climate_garden <- separate(climate_garden, "YearMonth",
  into = c("year", "month"),
  sep = 4
)
# Reshape data-- make a separate column for temperature and precipitation and vpd
climate_garden <- unique(climate_garden)

climate_garden_1995_2025 <- climate_garden %>%
  spread(clim, value) %>%
  rename(lon = longitude, lat = latitude, site = garden_sites.site_code) %>%
  mutate(year = as.numeric(year), month = as.numeric(month)) %>%
  filter(year > 1994)

summary(climate_garden_1995_2025)

# Subset the data collection period
climate_garden_2023_2025 <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date > as.Date("2023-05-01") & date < as.Date("2025-06-01"))

## Plot the daily trend for temperature and soil moisture from start to end
climate_garden_2023_2025 %>%
  mutate(vpd = (vpdmax + vpdmin) / 2) %>%
  group_by(site) %>%
  summarise(
    sum_ppt = sum(ppt, na.rm = TRUE),
    mean_temp = mean(tmean, na.rm = TRUE),
    mean_vpd = mean(vpd, na.rm = TRUE)
  ) -> prism_means

# Export prism data for the data collection period
prism_means <- as.data.frame(prism_means)
saveRDS(prism_means, "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Data/prism_means.rds")
#  Time scale provide insight into long-term drought trends.
# SPEI > 0: Wet conditions.
# SPEI < 0: Dry conditions (drought).
# SPEI ≤ -1.5: Moderate to severe drought.

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/Climate_prism.pdf", width = 14, height = 10, useDingbats = F)
par(mar = c(5, 5, 2, 3), mfrow = c(2, 2))
barplot(prism_means[order(prism_means[, 2], decreasing = FALSE), ][, 2], names.arg = prism_means[order(prism_means[, 2], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "", ylim = c(0, 3500))
mtext("Precipitation", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("A", side = 3, adj = 0, cex = 1.2)
barplot(prism_means[order(prism_means[, 3], decreasing = FALSE), ][, 3], names.arg = prism_means[order(prism_means[, 3], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "", ylim = c(0, 25))
mtext("Temperature", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("B", side = 3, adj = 0, cex = 1.2)
barplot(prism_means[order(prism_means[, 4], decreasing = FALSE), ][, 4], names.arg = prism_means[order(prism_means[, 4], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "")
mtext("Vapor Pressure Deficit", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("C", side = 3, adj = 0, cex = 1.2)

dev.off()

# site_names <- c(
#   "LAF" = "Lafayette",
#   "HUN" = "Huntville",
#   "KER" = "Kerville",
#   "BAS" = "Bastrop",
#   "COL" = "College Station",
#   "BFL" = "Brackenridge",
#   "SON" = "Sonora"
# )
# climate_garden_2023_2025 %>%
#   ggplot(aes(x = as.Date(date), y = ppt)) +
#   geom_line(aes(colour = site)) +
#   # ggtitle("d")+
#   # scale_color_manual(values = cbp1)+
#   # scale_fill_manual(values = cbp1)+
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(size = 4.5, color = "black", angle = 0),
#     plot.title = element_text(size = 14, color = "black", angle = 0)
#   ) +
#   labs(y = "Daily precipitation  (°C)", x = "Month") +
#   # ylim=c(0,100)+
#   # facet_grid(~site,labeller = labeller(site=site_names))+
#   facet_grid(~ factor(site, levels = c("HUN", "LAF", "COL", "BAS", "BFL", "KER", "SON"))) +
#   geom_hline(data = prism_means, aes(yintercept = sum_ppt, colour = site)) -> figpptsite_prism
#

# Extract all the value for the each yera prior data collection
datini <- read.csv("https://www.dropbox.com/scl/fi/9lemc8zv5fdi07r7j1jqq/Initialdata.csv?rlkey=gmryfg74motuf7h49s1aul0tm&dl=1", stringsAsFactors = F)
dat23 <- read.csv("https://www.dropbox.com/scl/fi/o18136j5t9sd94irtd46i/census2023.csv?rlkey=7byyglfehzqo5e197fp1jsasw&dl=1", stringsAsFactors = F)
datini23 <- right_join(x = datini, y = dat23, by = c("Tag_ID"))
datini23 %>%
  dplyr::select(Site, Species, date_23) %>%
  distinct() -> date_sp_site23
dat24 <- read.csv("https://www.dropbox.com/scl/fi/w6qdatnkosrjsvvub1vsf/census2024.csv?rlkey=tt1qbdaratj2zvt4igkxz89n9&dl=1", stringsAsFactors = F)

datini24 <- right_join(x = datini23, y = dat24, by = c("Tag_ID"))
datini23_spike <- datini23 %>% filter(!is.na(Tag_ID))
datini24 %>%
  dplyr::select(Site, Species, date_24) %>%
  distinct() %>%
  na.omit() -> date_sp_site24

dat25 <- read.csv("https://www.dropbox.com/scl/fi/x6gttyivrogzug0w3iveo/census_2025.csv?rlkey=jidxofw54907a73gdnn7t9n3i&dl=1", stringsAsFactors = F)
dat25 %>%
  dplyr::select(Site, Species, date_25) %>%
  distinct() %>%
  na.omit() -> date_sp_site25

climate_garden_2022_2025 <- climate_garden_1995_2025 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  filter(date > as.Date("2022-05-01") & date < as.Date("2025-06-01"))

# For 2023
#  Convert date_23 column to Date format if needed
obs_table_23 <- date_sp_site23 %>%
  mutate(date_23 = mdy(date_23)) # If it's character like "5/10/25"
# Ensure climate date column is in Date format
climate <- climate_garden_2022_2025 %>%
  mutate(date = ymd(date),
         vpd = (vpdmax + vpdmin)/2) # Only if it's not already a Date object
# For each observation, calculate cumulative and average tmean
climate_23 <- obs_table_23 %>%
  rowwise() %>%
  mutate(
    cumulative_pptmean = sum(
      climate %>%
        filter(
          site == Site,
          date >= date_23 - years(1),
          date <= date_23
        ) %>%
        pull(ppt),
      na.rm = TRUE
    ),
    average_tmean = mean(
      climate %>%
        filter(
          site == Site,
          date >= date_23 - years(1),
          date <= date_23
        ) %>%
        pull(tmean),
      na.rm = TRUE
    ),
    average_vpdmean = mean(
      climate %>%
        filter(
          site == Site,
          date >= date_23 - years(1),
          date <= date_23
        ) %>%
        pull(vpd),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

view(climate_23)

# For 2024
#  Convert date_24 column to Date format if needed
obs_table_24 <- date_sp_site24 %>%
  mutate(date_24 = mdy(date_24)) # If it's character like "5/10/25"
# For each observation, calculate cumulative and average tmean
climate_24 <- obs_table_24 %>%
  rowwise() %>%
  mutate(
    cumulative_pptmean = sum(
      climate %>%
        filter(
          site == Site,
          date >= date_24 - years(1),
          date <= date_24
        ) %>%
        pull(ppt),
      na.rm = TRUE
    ),
    average_tmean = mean(
      climate %>%
        filter(
          site == Site,
          date >= date_24 - years(1),
          date <= date_24
        ) %>%
        pull(tmean),
      na.rm = TRUE
    ),
    average_vpdmean = mean(
      climate %>%
        filter(
          site == Site,
          date >= date_24 - years(1),
          date <= date_24
        ) %>%
        pull(vpd),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

view(climate_24)

# For 2025
#  Convert date_25 column to Date format if needed
obs_table_25 <- date_sp_site25 %>%
  mutate(date_25 = mdy(date_25)) # If it's character like "5/10/25"
# For each observation, calculate cumulative and average tmean
climate_25 <- obs_table_25 %>%
  rowwise() %>%
  mutate(
    cumulative_pptmean = sum(
      climate %>%
        filter(
          site == Site,
          date >= date_25 - years(1),
          date <= date_25
        ) %>%
        pull(ppt),
      na.rm = TRUE
    ),
    average_tmean = mean(
      climate %>%
        filter(
          site == Site,
          date >= date_25 - years(1),
          date <= date_25
        ) %>%
        pull(tmean),
      na.rm = TRUE
    ),
    average_vpdmean = mean(
      climate %>%
        filter(
          site == Site,
          date >= date_25 - years(1),
          date <= date_25
        ) %>%
        pull(vpd),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

view(climate_25)



