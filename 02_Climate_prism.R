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
# PRISM data ----
# First, set a file path where prism data will be stored
options(prism.path = "/Users/jm200/Documents/Prism Range limit/")
#get_prism_monthlys(type="ppt",years=1994:2024,mon=1:12,keepZip = TRUE)
#get_prism_monthlys(type = "tmean", years = 1994:2024, mon = 1:12, keepZip = TRUE)
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
# Reshape data-- make a separate column for temperature and precipitation
climate_garden <- unique(climate_garden)

climate_garden_1994_2024 <- climate_garden %>%
  spread(clim, value) %>%
  rename(lon = longitude, lat = latitude, site = garden_sites.site_code) %>% 
  mutate(year=as.numeric(year),month=as.numeric(month)) %>% 
  filter(year>=1994)

summary(climate_garden_1994_2024)
# Separate the data per site

climate_garden_1994_2024 %>%
  filter(site == "SON") -> climate_garden_SON
climate_garden_1994_2024 %>%
  filter(site == "KER") -> climate_garden_KER
climate_garden_1994_2024 %>%
  filter(site == "BAS") -> climate_garden_BAS
climate_garden_1994_2024 %>%
  filter(site == "BFL") -> climate_garden_BFL
climate_garden_1994_2024 %>%
  filter(site == "LAF") -> climate_garden_LAF
climate_garden_1994_2024 %>%
  filter(site == "COL") -> climate_garden_COL
climate_garden_1994_2024 %>%
  filter(site == "HUN") -> climate_garden_HUN

# Load the required libraries
# Function to calculate SPEI using the Thornthwaite method for PET
# Function to calculate SPEI using Thornthwaite PET over a two-year (24-month) scale
calculate_SPEI <- function(precipitation, temperature, latitude, site_names, scale, start_year) {
  # Ensure input vectors have the correct length
  if (length(precipitation) != length(temperature)) {
    stop("Precipitation and temperature vectors must have the same length.")
  }
  # Create a dataframe for the site data
  months <- rep(1:12, length(precipitation) / 12)
  data <- data.frame(
    site = rep(site_names, each = length(precipitation) / length(site_names)), # Assign site names
    month = months,
    latitude = rep(latitude, length(precipitation)),  # Repeat latitude for all months
    precipitation = precipitation,
    temperature = temperature
  )
  # Create a Date column based on start year and month
  data$year <- rep(seq(start_year, length.out = length(precipitation) / 12, by = 1), each = 12)[1:length(precipitation)]
  data$date <- as.Date(paste(data$year, data$month, "01", sep = "-"))
  # Calculate PET using the Thornthwaite method from the SPEI package
  data$PET <- thornthwaite(data$temperature, latitude)  # Thornthwaite PET calculation
  # Calculate water balance: Water Balance = Precipitation - PET
  data$water_balance <- data$precipitation - data$PET
  # Calculate the SPEI using the water balance and specified scale
  spei_result <- spei(data$water_balance, scale = scale)
  # Return the required columns
  result <- data.frame(
    date = data$date, 
    site = data$site,  
    month = data$month,
    tmean = data$temperature,  
    ppt = data$precipitation,
    PET = data$PET,  
    SPEI = spei_result$fitted  # Extract the fitted SPEI values
  )
  return(result)
}
# Apply function to each site
climate_garden_SON_SPEI <- calculate_SPEI(climate_garden_SON$ppt, climate_garden_SON$tmean, latitude=30.27347,site_names="SON", scale = 3, start_year = 2000)
climate_garden_KER_SPEI <- calculate_SPEI(climate_garden_KER$ppt, climate_garden_KER$tmean, latitude=30.02844, site_names="KER",scale = 3, start_year = 2000)
climate_garden_BAS_SPEI <- calculate_SPEI(climate_garden_BAS$ppt, climate_garden_BAS$tmean, latitude=30.09234, site_names="BAS",scale = 3, start_year = 2000)
climate_garden_BFL_SPEI <- calculate_SPEI(climate_garden_BFL$ppt, climate_garden_BFL$tmean, latitude=30.28487,site_names="BFL",scale = 3, start_year = 2000)
climate_garden_LAF_SPEI <- calculate_SPEI(climate_garden_LAF$ppt, climate_garden_LAF$tmean, latitude=30.30477,site_names="LAF", scale = 3, start_year = 2000)
climate_garden_COL_SPEI <- calculate_SPEI(climate_garden_COL$ppt, climate_garden_COL$tmean, latitude=30.56930,site_names="COL", scale = 3, start_year = 2000)
climate_garden_HUN_SPEI <- calculate_SPEI(climate_garden_HUN$ppt, climate_garden_HUN$tmean, latitude=30.74281,site_names="HUN", scale = 3, start_year = 2000)

# Put all the sites together
climate_garden_SPEI <- bind_rows(
  climate_garden_HUN_SPEI,
  climate_garden_COL_SPEI,
  climate_garden_LAF_SPEI,
  climate_garden_BAS_SPEI,
  climate_garden_BFL_SPEI,
  climate_garden_KER_SPEI,
  climate_garden_SON_SPEI
)

# Subset the data collection period
climate_garden_SPEI %>%
  filter(date > as.Date("2023-05-01") & date < as.Date("2024-06-01"))-> climate_garden_SPEI_2023_2024

## Plot the daily trend for temperature and soil moisture from start to end
climate_garden_SPEI_2023_2024 %>%
  group_by(site) %>%
  summarise(
    sum_ppt = sum(ppt),
    mean_temp = mean(tmean),
    mean_pet = mean(PET),
    mean_spei = mean(SPEI)
  ) -> prism_means

prism_means <- as.data.frame(prism_means)
saveRDS(prism_means,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Data/prism_means.rds")
#  Time scale provide insight into long-term drought trends.
# SPEI > 0: Wet conditions.
# SPEI < 0: Dry conditions (drought).
# SPEI ≤ -1.5: Moderate to severe drought.

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Climate_prism.pdf", width = 14, height = 10, useDingbats = F)
par(mar = c(5, 5, 2, 3), mfrow = c(2, 2))
barplot(prism_means[order(prism_means[, 2], decreasing = FALSE), ][, 2], names.arg = prism_means[order(prism_means[, 2], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "", ylim = c(0, 2000))
mtext("Precipitation", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("A", side = 3, adj = 0, cex = 1.2)
barplot(prism_means[order(prism_means[, 3], decreasing = FALSE), ][, 3], names.arg = prism_means[order(prism_means[, 3], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "", ylim = c(0, 25))
mtext("Temperature", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("B", side = 3, adj = 0, cex = 1.2)
barplot(prism_means[order(prism_means[, 4], decreasing = FALSE), ][, 4], names.arg = prism_means[order(prism_means[, 4], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "")
mtext("Precipitation Evapotranspiration ", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("C", side = 3, adj = 0, cex = 1.2)
barplot(prism_means[order(prism_means[, 5], decreasing = FALSE), ][, 5], names.arg = prism_means[order(prism_means[, 5], decreasing = FALSE), ][, 1], col = "#E69F00", xlab = "Sites", ylab = "Mean", main = "")
mtext("Standardized Precipitation Evapotranspiration Index", side = 3, adj = 0.5, cex = 1.2, line = 0.3)
mtext("D", side = 3, adj = 0, cex = 1.2)
dev.off()

site_names <- c(
  "LAF" = "Lafayette",
  "HUN" = "Huntville",
  "KER" = "Kerville",
  "BAS" = "Bastrop",
  "COL" = "College Station",
  "BFL" = "Brackenridge",
  "SON" = "Sonora"
)
climate_garden_SPEI_2023_2024 %>%
  ggplot(aes(x = as.Date(date), y = ppt)) +
  geom_line(aes(colour = site)) +
  # ggtitle("d")+
  # scale_color_manual(values = cbp1)+
  # scale_fill_manual(values = cbp1)+
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 4.5, color = "black", angle = 0),
    plot.title = element_text(size = 14, color = "black", angle = 0)
  ) +
  labs(y = "Daily precipitation  (°C)", x = "Month") +
  # ylim=c(0,100)+
  # facet_grid(~site,labeller = labeller(site=site_names))+
  facet_grid(~ factor(site, levels = c("HUN", "LAF", "COL", "BAS", "BFL", "KER", "SON"))) +
  geom_hline(data = prism_means, aes(yintercept = sum_ppt, colour = site)) -> figpptsite_prism

ggplot(data = climate_garden_SPEI_2023_2024, mapping = aes(x = tmean, y = ppt)) +
  geom_point(shape = 21, fill = "#0f993d", color = "white", size = 3) +
  sm_statCorr()
