# Project: Distance from niche centroid and distance from range center using species climatic niche modeling
# Purpose: Import and prepare climate and occurrence data for ellipsoid modeling and calculating distance from range edge
# Note: Raster files are stored locally due to size constraints
# Author: Jacob Moutouama
# Last modified: 2024-08-03
# Clean environment
rm(list = ls(all = TRUE)) # Remove all objects from the environment
# Load required packages
library(tidyverse) # Data manipulation and visualization
library(terra) # Raster data handling
library(geojsonio) # Reading GeoJSON files
library(sp) # Spatial data classes and functions
library(ggspatial) # Spatial plotting with ggplot2
library(ntbox) # Niche modeling tools including ellipsoid modeling
library(raster) # Raster data (older package used in ntbox)
library(stringr) # String manipulation
library(prism) # PRISM climate data access
library(RColorBrewer) # Color palettes for maps
library(dplyr)
library(geosphere)
library(mgcv)
options(rgl.useNULL = TRUE)
library(rgl)

# Set random seed for reproducibility
set.seed(13)
# Load PRISM climate data----
# making a folder to store prism data
options(prism.path = "/Users/jacobmoutouama/Documents/prism/")
# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
# get_prism_monthlys(type = "tmean", years = 1990:2024, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "ppt", years = 1990:2024, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "vpdmin", years = 1990:2025, mon = 1:12, keepZip = TRUE)
# get_prism_monthlys(type = "vpdmax", years = 1990:2025, mon = 1:12, keepZip = TRUE)
# pulling out values to get normals for time periods

# Load and compute monthly temperature
tmean_annual_norm <- terra::mean(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023
  )
)))
tmean_spring_norm <- terra::mean(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023,
    mon = 1:4
  )
)))
tmean_summer_norm <- terra::mean(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023,
    mon = 5:8
  )
)))
tmean_autumn_norm <- terra::mean(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023,
    mon = 9:12
  )
)))

# calculating standard deviation in temp
tmean_annual_sd <- terra::stdev(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023
  )
)))
tmean_spring_sd <- terra::stdev(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023,
    mon = 1:4
  )
)))
tmean_summer_sd <- terra::stdev(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023,
    mon = 5:8
  )
)))
tmean_autumn_sd <- terra::stdev(terra::rast(pd_stack(
  prism_archive_subset(
    type = "tmean",
    temp_period = "monthly",
    year = 1993:2023,
    mon = 9:12
  )
)))


# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual <- ppt_spring <- ppt_summer <- ppt_autumn <- ppt_winter <- list()
for (y in 1993:2023) {
  ppt_annual[[y]] <- sum(terra::rast(pd_stack(
    prism_archive_subset(
      type = "ppt",
      temp_period = "monthly",
      year = y
    )
  )))
  ppt_spring[[y]] <- sum(terra::rast(pd_stack(
    prism_archive_subset(
      type = "ppt",
      temp_period = "monthly",
      year = y,
      mon = 1:4
    )
  )))
  ppt_summer[[y]] <- sum(terra::rast(pd_stack(
    prism_archive_subset(
      type = "ppt",
      temp_period = "monthly",
      year = y,
      mon = 5:8
    )
  )))
  ppt_autumn[[y]] <- sum(terra::rast(pd_stack(
    prism_archive_subset(
      type = "ppt",
      temp_period = "monthly",
      year = y,
      mon = 9:12
    )
  )))
}

# Taking the mean of the cumulative precipitation values
ppt_annual_norm <- terra::mean(terra::rast(unlist(ppt_annual)))
ppt_spring_norm <- terra::mean(terra::rast(unlist(ppt_spring)))
ppt_summer_norm <- terra::mean(terra::rast(unlist(ppt_summer)))
ppt_autumn_norm <- terra::mean(terra::rast(unlist(ppt_autumn)))

# calculating the standard deviation in precip
ppt_annual_sd <- terra::stdev(terra::rast(unlist(ppt_annual)))
ppt_spring_sd <- terra::stdev(terra::rast(unlist(ppt_spring)))
ppt_summer_sd <- terra::stdev(terra::rast(unlist(ppt_summer)))
ppt_autumn_sd <- terra::stdev(terra::rast(unlist(ppt_autumn)))

# Study area shapefile
US <- terra::vect(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp"
)
US_land <- US[(
  !US$State_Name %in% c(
    "HAWAII",
    "ALASKA",
    "ARIZONA",
    "COLORADO",
    "UTAH",
    "NEVADA",
    "NEW MEXICO",
    "IDAHO",
    "MONTANA",
    "WYOMING",
    "CALIFORNIA",
    "WASHINGTON",
    "OREGON"
  )
), ]
US_land_reprojected <- terra::project(US_land, crs(tmean_annual_norm))
plot(US_land_reprojected)

# Crop the study area -----
tmean_spring_norm <- terra::crop(tmean_spring_norm, US_land_reprojected, mask = TRUE)
tmean_summer_norm <- terra::crop(tmean_summer_norm, US_land_reprojected, mask = TRUE)
tmean_autumn_norm <- terra::crop(tmean_autumn_norm, US_land_reprojected, mask = TRUE)
tmean_spring_sd <- terra::crop(tmean_spring_sd, US_land_reprojected, mask = TRUE)
tmean_summer_sd <- terra::crop(tmean_summer_sd, US_land_reprojected, mask = TRUE)
tmean_autumn_sd <- terra::crop(tmean_autumn_sd, US_land_reprojected, mask = TRUE)

ppt_spring_norm <- terra::crop(ppt_spring_norm, US_land_reprojected, mask = TRUE)
ppt_summer_norm <- terra::crop(ppt_summer_norm, US_land_reprojected, mask = TRUE)
ppt_autumn_norm <- terra::crop(ppt_autumn_norm, US_land_reprojected, mask = TRUE)
ppt_spring_sd <- terra::crop(ppt_spring_sd, US_land_reprojected, mask = TRUE)
ppt_summer_sd <- terra::crop(ppt_summer_sd, US_land_reprojected, mask = TRUE)
ppt_autumn_sd <- terra::crop(ppt_autumn_sd, US_land_reprojected, mask = TRUE)

# Now stacking all climatic variables including updated VPD min and max -----
US_land_clim <- terra::rast(
  list(
    tmean_spring_norm,
    tmean_summer_norm,
    tmean_autumn_norm,
    tmean_spring_sd,
    tmean_summer_sd,
    tmean_autumn_sd,
    ppt_spring_norm,
    ppt_summer_norm,
    ppt_autumn_norm,
    ppt_spring_sd,
    ppt_summer_sd,
    ppt_autumn_sd
  )
)

# Naming the layers -----
names(US_land_clim) <- c(
  "tmean_spring_norm",
  "tmean_summer_norm",
  "tmean_autumn_norm",
  "tmean_spring_sd",
  "tmean_summer_sd",
  "tmean_autumn_sd",
  "ppt_spring_norm",
  "ppt_summer_norm",
  "ppt_autumn_norm",
  "ppt_spring_sd",
  "ppt_summer_sd",
  "ppt_autumn_sd"
)

# Read rasters as a stack
# Visualize climate data (optional)
# plot(clim_stack)
US_land_clim_stack <- stack(US_land_clim)
plot(US_land_clim_stack)

# Function to thin occurrences - Keep only one per climatic data pixel
thin_occurrences <- function(occ_data, clim_raster) {
  occ_data <- occ_data %>%
    filter(!is.na(lat) & !is.na(lon)) %>%
    unique()
  
  # Assign raster cell ID to each occurrence point
  occ_data$cell <- terra::cellFromXY(clim_raster, occ_data[, c("lon", "lat")])
  
  # Keep only one occurrence per pixel (cell)
  occ_thinned <- occ_data %>%
    group_by(cell) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-cell) # Remove temporary column
  
  return(occ_thinned)
}

# Agrostis hyemalis-----
# Apply thinning function to each species dataset
aghy_occ_raw <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/ijl1i7964qxzcxvm7blz8/aghy_occ_raw.rds?rlkey=msbs3xzkjc719uld8yw7hb6cj&dl=1"
  )
)
# names(aghy_occ_raw)
aghy_occ_raw %>%
  filter(
    !is.na(lat),
    !is.na(lon),
    !is.na(year),
    # !is.na(coordinatePrecision),!is.na(coordinateUncertaintyInMeters),
    # coordinatePrecision <= 0.01,
    coordinateUncertaintyInMeters <= 1000,
    lon >= -102.6458,
    lat < 45,
    country == "United States"
  ) %>%
  distinct() %>%
  dplyr::select(country, lon, lat, year) %>%
  arrange(lat) -> aghy

# Thin occurrences
aghy_occ_thinned <- thin_occurrences(aghy, US_land_clim)

# Plot thinned occurrence points for Agrostis hyemalis
plot(US_land_reprojected)
points(aghy_occ_thinned[, c("lon", "lat")],
       pch = 20,
       cex = 0.5,
       col = "red")

# Model calibration selection using Minimum Volume Ellipsoids (MVEs)
train_index_aghy <- sample(1:nrow(aghy_occ_thinned), 0.80 * nrow(aghy_occ_thinned))
test_index_aghy <- setdiff(1:nrow(aghy_occ_thinned), train_index_aghy)
# Split occurrences into train and test
aghy_train <- aghy_occ_thinned[train_index_aghy, ]
aghy_test <- aghy_occ_thinned[test_index_aghy, ]

# Extract environmental information for both train and test data
aghy_etrain <- raster::extract(US_land_clim_stack, aghy_train[, c("lon", "lat")], df = TRUE)
#sum(na.omit(aghy_etrain))
aghy_etrain <- na.omit(aghy_etrain)[, -1]
#summary(aghy_etrain)

aghy_etest <- raster::extract(US_land_clim_stack, aghy_test[, c("lon", "lat")], df = TRUE)
aghy_etest <- na.omit(aghy_etest)[, -1]
#summary(aghy_etest)

# Find correlated environmental variables
env_varsL_aghy <- ntbox::correlation_finder(cor(aghy_etrain, method = "spearman"),
                                            threshold = 0.75,
                                            verbose = F)
env_vars_aghy <- env_varsL_aghy$descriptors
print(env_vars_aghy)

# Fit ellipsoid models
nvarstest <- c(3,4)
level <- 0.99
env_bg <- ntbox::sample_envbg(US_land_clim_stack, 10000)
omr_criteria <- 0.06
proc <- TRUE

e_select_aghy <- ntbox::ellipsoid_selection(
  env_train = aghy_etrain,
  env_test = aghy_etest,
  env_vars = env_vars_aghy,
  level = level,
  nvarstest = nvarstest,
  env_bg = env_bg,
  omr_criteria = omr_criteria,
  proc = proc
)

# Display the first 10 rows of the results
head(e_select_aghy, 10)

# Best ellipsoid model for "omr_criteria"
bestvarcomb_aghy <- stringr::str_split(e_select_aghy$fitted_vars, ",")[[1]]
best_mod_aghy <- ntbox::cov_center(
  aghy_etrain[, bestvarcomb_aghy],
  mve = TRUE,
  level = 0.99,
  vars = 1:length(bestvarcomb_aghy)
)

# Projection model in geographic space
mProj_aghy <- ntbox::ellipsoidfit(
  US_land_clim_stack[[bestvarcomb_aghy]],
  centroid = best_mod_aghy$centroid,
  covar = best_mod_aghy$covariance,
  level = 0.99,
  lw = 10,
  size = 2
)

if (length(bestvarcomb_aghy) == 3) {
  rgl::rglwidget(reuse = TRUE)
}

# Mahalanobis distance for common garden populations
aghy_clim <- terra::rast(list(tmean_spring_norm,ppt_spring_norm,
                              ppt_summer_norm,ppt_autumn_norm,ppt_summer_sd))

names(aghy_clim) <- c("tmean_spring_norm","ppt_spring_norm",
                      "ppt_summer_norm","ppt_autumn_norm","ppt_summer_sd")

aghy_clim_stack <- stack(aghy_clim)
plot(aghy_clim_stack)

garden <- read.csv(
  "https://www.dropbox.com/scl/fi/si346imz380lpdgo9yekr/Study_site.csv?rlkey=yiue42npkzzu9w8dggjr4fent&dl=1",
  stringsAsFactors = FALSE
) %>%
  unique() %>%
  arrange(latitude)

garden %>%
  filter(Species == "AGHY") -> garden_AGHY

garden_aghy_clim <- raster::extract(aghy_clim_stack, garden_AGHY[, c("longitude", "latitude")], df = TRUE)[, -1]
mhd_aghy <- stats::mahalanobis(garden_aghy_clim[, bestvarcomb_aghy],
                               center = best_mod_aghy$centroid,
                               cov = best_mod_aghy$covariance)
distance_aghy <- data.frame(garden_AGHY, distance = mhd_aghy)
plot(mhd_aghy ~ longitude, data = distance_aghy)
cor.test(distance_aghy$longitude, distance_aghy$distance)


# Extract ellipsoid suitability layer (values between 0 and 1, inside ellipsoid)
# suitability_map <- mProj_aghy$suitRaster
# plot(suitability_map)
# Extract the suitability value for each population
# suit_vals <- terra::extract(suitability_map, garden_AGHY[, c("longitude", "latitude")])

# Compute distance from edge as: distance = 1 - suitability score
# Inside niche: suitability close to 1, so distance is small
# Near the edge: suitability close to 0, so distance is large
# niche_edge_dist <- 1 - suit_vals

# Combine into final dataframe
# distance_aghy <- data.frame(garden_AGHY, distance = niche_edge_dist)

# Plot
# plot(niche_edge_dist ~ longitude, data = distance_aghy)
# cor.test(distance_aghy$longitude, distance_aghy$distance)

# Elymus virginicus----
# elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE)
# saveRDS(elvi_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/elvi_occ_raw.rds")
elvi_occ_raw <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/0ssa5gepxyz28b7ykw1x8/elvi_occ_raw.rds?rlkey=4dx0q4lw2112droh73hmh7xte&dl=1"
  )
)

elvi_occ_raw %>%
  filter(
    !is.na(lat),
    !is.na(lon),
    !is.na(year),
    # !is.na(coordinatePrecision),!is.na(coordinateUncertaintyInMeters),
    # coordinatePrecision <= 0.01,
    coordinateUncertaintyInMeters <= 1000,
    lon >= -102.6458,
    # lat < 45,
    country == "United States"
  ) %>%
  distinct() %>%
  dplyr::select(country, lon, lat, year) %>%
  arrange(lat) -> elvi

# Thin occurrences
elvi_occ_thinned <- thin_occurrences(elvi, US_land_clim)

# Plot thinned occurrences
plot(US_land_reprojected)
points(elvi_occ_thinned[, c("lon", "lat")],
       pch = 20,
       cex = 0.5,
       col = "red")

# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_elvi <- sample(1:nrow(elvi_occ_thinned), 0.80 * nrow(elvi_occ_thinned))
test_index_elvi <- setdiff(1:nrow(elvi_occ_thinned), train_index_elvi)

# Split occurrences into train and test
elvi_train <- elvi_occ_thinned[train_index_elvi, ]
elvi_test <- elvi_occ_thinned[test_index_elvi, ]

# Extracts the environmental information for both train and test data
elvi_etrain <- raster::extract(US_land_clim_stack, elvi_train[, c("lon", "lat")], df = TRUE)
elvi_etrain <- na.omit(elvi_etrain)
elvi_etrain <- elvi_etrain[, -1]

elvi_etest <- raster::extract(US_land_clim_stack, elvi_test[, c("lon", "lat")], df = TRUE)
elvi_etest <- na.omit(elvi_etest)
elvi_etest <- elvi_etest[, -1]

env_varsL_elvi <- ntbox::correlation_finder(cor(elvi_etrain, method = "spearman"),
                                            threshold = 0.75,
                                            verbose = F)
env_vars_elvi <- env_varsL_elvi$descriptors
print(env_vars_elvi)

# Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3 dimensions
nvarstest <- c(3,4)

# Now we use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_elvi <- ntbox::ellipsoid_selection(
  env_train = elvi_etrain,
  env_test = elvi_etest,
  env_vars = env_vars_elvi,
  level = level,
  nvarstest = nvarstest,
  env_bg = env_bg,
  omr_criteria = omr_criteria,
  proc = proc
)

# Let’s see the first 10 rows of the results
head(e_select_elvi, 10)

# Best ellipsoid model for "omr_criteria"
bestvarcomb_elvi <- stringr::str_split(e_select_elvi$fitted_vars, ",")[[1]]

# Ellipsoid model (environmental space)
best_mod_elvi <- ntbox::cov_center(
  elvi_etrain[, bestvarcomb_elvi],
  mve = T,
  level = 0.99,
  vars = 1:length(bestvarcomb_elvi)
)

# Projection model in geographic space
mProj_elvi <- ntbox::ellipsoidfit(
  US_land_clim_stack[[bestvarcomb_elvi]],
  centroid = best_mod_elvi$centroid,
  covar = best_mod_elvi$covariance,
  level = 0.99,
  size = 3
)
if (length(bestvarcomb_elvi) == 3) {
  rgl::rglwidget(reuse = TRUE)
}

elvi_clim <- terra::rast(list(tmean_spring_norm,ppt_spring_norm,
                              ppt_summer_norm,ppt_spring_sd,ppt_summer_sd))

names(elvi_clim) <- c("tmean_spring_norm","ppt_spring_norm",
                      "ppt_summer_norm","ppt_spring_sd","ppt_summer_sd")

elvi_clim_stack <- stack(elvi_clim)
plot(elvi_clim_stack)

garden %>%
  filter(Species == "ELVI") -> garden_ELVI
elvi_garden_clim <- raster::extract(elvi_clim_stack, garden_ELVI[, c("longitude", "latitude")], df = TRUE)[, -1]

# Mahalanobis distance for common garden populations
mhd_elvi <- stats::mahalanobis(elvi_garden_clim[, bestvarcomb_elvi],
                               center = best_mod_elvi$centroid,
                               cov = best_mod_elvi$covariance)
distance_elvi <- data.frame(garden_ELVI, distance = mhd_elvi)
plot(mhd_elvi ~ longitude, data = distance_elvi)
cor.test(distance_elvi$longitude, distance_elvi$distance)


# Poa autumnalis---
# poa_occ_raw <- gbif(genus="Poa", species="autumnalis", download=TRUE)
# saveRDS(poa_occ_raw, file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/ELVI Model output/occurence/poa_occ_raw.rds")
poa_occ_raw <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/oip7ndyf0d99rqxcqxb0q/poa_occ_raw.rds?rlkey=920uql1gd4gahnh8utw9fz96l&dl=1"
  )
)

poa_occ_raw %>%
  filter(
    !is.na(lat),
    !is.na(lon),
    !is.na(year),
    # !is.na(coordinatePrecision),!is.na(coordinateUncertaintyInMeters),
    # coordinatePrecision <= 0.01,
    # coordinateUncertaintyInMeters <= 1000,
    lon >= -102.6458,
    lat < 45,
    country == "United States"
  ) %>%
  distinct() %>%
  dplyr::select(country, lon, lat, year) %>%
  arrange(lat) -> poa

# Thin occurrences
poa_occ_thinned <- thin_occurrences(poa, US_land_clim)

# Plot thinned occurrences
plot(US_land_reprojected)
points(poa_occ_thinned[, c("lon", "lat")],
       pch = 20,
       cex = 0.5,
       col = "red")

# Model calibration selection using Minimum Volume Ellipsoids (MVEs).
# Random sample indexes
train_index_poa <- sample(1:nrow(poa_occ_thinned), 0.80 * nrow(poa_occ_thinned))
test_index_poa <- setdiff(1:nrow(poa_occ_thinned), train_index_poa)

# Split occurrences into train and test
poa_train <- poa_occ_thinned[train_index_poa, ]
poa_test <- poa_occ_thinned[test_index_poa, ]

# Extracts the environmental information for both train and test data
poa_etrain <- raster::extract(US_land_clim_stack, poa_train[, c("lon", "lat")], df = TRUE)
poa_etrain <- na.omit(poa_etrain)
poa_etrain <- poa_etrain[, -1]

poa_etest <- raster::extract(US_land_clim_stack, poa_test[, c("lon", "lat")], df = TRUE)
poa_etest <- na.omit(poa_etest)
poa_etest <- poa_etest[, -1]

env_varsL_poa <- ntbox::correlation_finder(cor(poa_etrain, method = "spearman"),
                                           threshold = 0.75,
                                           verbose = F)
env_vars_poa <- env_varsL_poa$descriptors
print(env_vars_poa)

# Now we specify the number of variables to fit the ellipsoid models; in the example, we will fit for 3 dimensions
nvarstest <- c(3,4)

# Now we use the function ellipsoid_selection to run the model calibration and selection protocol
e_select_poa <- ntbox::ellipsoid_selection(
  env_train = poa_etrain,
  env_test = poa_etest,
  env_vars = env_vars_poa,
  level = level,
  nvarstest = nvarstest,
  env_bg = env_bg,
  omr_criteria = omr_criteria,
  proc = proc
)

# Let’s see the first 10 rows of the results
head(e_select_poa, 10)

# Best ellipsoid model for "omr_criteria"
bestvarcomb_poa <- stringr::str_split(e_select_poa$fitted_vars, ",")[[1]]

# Ellipsoid model (environmental space)
best_mod_poa <- ntbox::cov_center(
  poa_etrain[, bestvarcomb_poa],
  mve = T,
  level = 0.99,
  vars = 1:length(bestvarcomb_poa)
)

# Projection model in geographic space
mProj_poa <- ntbox::ellipsoidfit(
  US_land_clim_stack[[bestvarcomb_poa]],
  centroid = best_mod_poa$centroid,
  covar = best_mod_poa$covariance,
  level = 0.99,
  size = 3
)
if (length(bestvarcomb_poa) == 3) {
  rgl::rglwidget(reuse = TRUE)
}

poa_clim <- terra::rast(list(tmean_spring_norm,ppt_spring_norm,
                             ppt_summer_norm,ppt_autumn_norm,  
                             ppt_spring_sd,ppt_summer_sd,ppt_autumn_sd))
names(poa_clim) <- c("tmean_spring_norm","ppt_spring_norm",
                     "ppt_summer_norm","ppt_autumn_norm",  
                     "ppt_spring_sd","ppt_summer_sd","ppt_autumn_sd")

poa_clim_stack <- stack(poa_clim)
plot(poa_clim_stack)

garden %>%
  filter(Species == "POAU") -> garden_POAU
poa_garden_clim <- raster::extract(poa_clim_stack, garden_POAU[, c("longitude", "latitude")], df = TRUE)[, -1]

# Mahalanobis distance for common garden populations
mhd_poa <- stats::mahalanobis(poa_garden_clim[, bestvarcomb_poa],
                              center = best_mod_poa$centroid,
                              cov = best_mod_poa$covariance)
distance_poa <- data.frame(garden_POAU, distance = mhd_poa)
plot(mhd_poa ~ longitude, data = distance_poa)
cor.test(distance_poa$longitude, distance_poa$distance)

# Combine distance data for all species
distance_species <- bind_rows(distance_aghy, distance_elvi, distance_poa)
Species.label <- c("AGHY", "ELVI", "POAU")
names(Species.label) <- c("A. hyemalis", "E. virginicus", "P. autumnalis")

# Function to calculate the distance to the centroid using Earth's curvature
calculate_distance_to_centroid <- function(occurrence_data, garden_data) {
  # Filter the occurrence data (remove NA values)
  occurrence_data <- occurrence_data %>%
    filter(!is.na(lat) & !is.na(lon))
  
  # Calculate the centroid (mean latitude and longitude)
  centroid_lon <- mean(occurrence_data$lon, na.rm = TRUE)
  centroid_lat <- mean(occurrence_data$lat, na.rm = TRUE)
  
  # Calculate the distance from each garden point to the centroid using the Vincenty formula (takes curvature into account)
  garden_data$distance_to_centroid <- distVincentySphere(
    cbind(garden_data$longitude, garden_data$latitude),
    c(centroid_lon, centroid_lat)
  )
  
  return(garden_data)
}
# Assuming you have the `elvi_occ_thinned` data and `garden_clim` data
garden_coordinate_aghy <- garden_AGHY[, 2:3]
aghy_with_distance <- calculate_distance_to_centroid(aghy_occ_thinned, garden_coordinate_aghy)
aghy_geo_distance <- data.frame(site_code = garden_AGHY$site_code, aghy_with_distance)
aghy_geo_distance$Species <- rep("AGHY", nrow(aghy_geo_distance))

elvi_with_distance <- calculate_distance_to_centroid(elvi_occ_thinned, garden_ELVI[, 2:3])
elvi_geo_distance <- data.frame(site_code = garden_ELVI$site_code, elvi_with_distance)
elvi_geo_distance$Species <- rep("ELVI", nrow(elvi_geo_distance))

poa_with_distance <- calculate_distance_to_centroid(poa_occ_thinned, garden_POAU[, 2:3])
poa_geo_distance <- data.frame(site_code = garden_POAU$site_code, poa_with_distance)
poa_geo_distance$Species <- rep("POAU", nrow(poa_geo_distance))
geo_distance <- bind_rows(aghy_geo_distance, elvi_geo_distance, poa_geo_distance)

distance_species <- bind_rows(distance_aghy, distance_elvi, distance_poa)
Species.label <- c("AGHY", "ELVI", "POAU")
names(Species.label) <- c("A. hyemalis", "E. virginicus", "P. autumnalis")
distance_species <- cbind(distance_species, geo_distance = geo_distance[, 4])
# saveRDS(
#   distance_species,
#   "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/distance_species.rds"
# )

pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/distance_plot.pdf",
  width = 9,
  height = 8
)
# Set up 2x2 layout
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))

# Custom colors and labels
cols <- c(
  "AGHY" = "#009E73",
  "ELVI" = "#D55E00",
  "POAU" = "#0072B2"
)
labels <- c(
  "AGHY" = expression(italic("A. hyemalis")),
  "ELVI" = expression(italic("E. virginicus")),
  "POAU" = expression(italic("P. autumnalis"))
)

# Panel labels
panel_labels <- c("A", "B", "C", "D")

# First 3 plots: Longitude vs. Mahalanobis distance by species
for (i in 1:length(cols)) {
  sp <- names(cols)[i]
  dat <- subset(distance_species, Species == sp)
  
  # Fit GAM
  gam_model <- gam(distance ~ s(longitude, k = 4), data = dat)
  pred_long <- seq(min(dat$longitude), max(dat$longitude), length.out = 200)
  pred <- predict(gam_model,
                  newdata = data.frame(longitude = pred_long),
                  se.fit = TRUE)
  
  # Plot points
  plot(
    dat$longitude,
    dat$distance,
    col = cols[sp],
    pch = 16,
    xlab = "Longitude",
    ylab = "Mahalanobis distance",
    main = labels[[sp]],
    cex.lab = 1.3,
    cex.axis = 1.1
  )
  
  # Add smooth line with CI from GAM
  polygon(
    c(pred_long, rev(pred_long)),
    c(pred$fit + 2 * pred$se.fit, rev(pred$fit - 2 * pred$se.fit)),
    col = adjustcolor(cols[sp], alpha.f = 0.2),
    border = NA
  )
  lines(pred_long, pred$fit, col = cols[sp], lwd = 2)
  
  # Add panel label (A, B, C, D)
  mtext(
    panel_labels[i],
    side = 3,
    line = 0.2,
    at = par("usr")[1],
    adj = 0,
    cex = 1.25,
    font = 1
  )
}

# 4th plot: log-log plot with GAM
# Filter to valid values
valid <- distance_species$geo_distance > 0 &
  distance_species$distance > 0
log_geo <- log(distance_species$geo_distance[valid])
log_dist <- log(distance_species$distance[valid])
species_colors <- cols[distance_species$Species[valid]]

# Fit GAM to log-log
log_data <- data.frame(log_geo = log_geo, log_dist = log_dist)
gam_log <- gam(log_dist ~ s(log_geo, k = 4), data = log_data)

# Predictions
new_logx <- seq(min(log_geo), max(log_geo), length.out = 200)
pred_log <- predict(gam_log,
                    newdata = data.frame(log_geo = new_logx),
                    se.fit = TRUE)

# Plot log-log data
plot(
  log_geo,
  log_dist,
  col = species_colors,
  pch = 16,
  xlab = expression(log * " (Distance from geographic center)"),
  ylab = expression(log * " (Mahalanobis distance)"),
  main = "",
  cex.lab = 1.3,
  cex.axis = 1.1
)

# Add confidence band and smooth from GAM
polygon(
  c(new_logx, rev(new_logx)),
  c(
    pred_log$fit + 2 * pred_log$se.fit,
    rev(pred_log$fit - 2 * pred_log$se.fit)
  ),
  col = adjustcolor("black", alpha.f = 0.2),
  border = NA
)
lines(new_logx, pred_log$fit, col = "black", lwd = 2)

# Add p-value from the GAM summary
gam_p <- summary(gam_log)$s.table[1, "p-value"]
p_label <- ifelse(gam_p < 0.001, "p < 0.001", paste0("p = ", signif(gam_p, 3)))

# Add p-value label (for GAM)
# text(x = min(log_geo) + 0.1,
#      y = max(log_dist) - 0.1,
#      labels = p_label,
#      cex = 1.5,
#      font = 1)

# Add panel label (D)
mtext(
  panel_labels[4],
  side = 3,
  line = 0.2,
  at = par("usr")[1],
  adj = 0,
  cex = 1.25,
  font = 1
)

dev.off()

# Plot the SDM map

pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/SDM.pdf",
  width = 12,
  height = 10,
  useDingbats = F
)
par(mar = c(5, 5, 2, 3), mfrow = c(2, 2))
raster::plot(
  mProj_aghy$suitRaster,
  main = "",
  xlab = "Longitude",
  ylab = "Latitude",
  cex.lab = 1.5,
  col = rev(terrain.colors(100)),
  legend = TRUE
)
# points(aghy[,c("lon","lat")],pch=23,cex=0.3,col="grey")
# plot(garden_map_aghy,add=T,pch = 3,col="black",cex =2)
# plot(source_map_aghy,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("A",
      side = 3,
      adj = 0,
      cex = 1.25)
mtext(
  ~ italic("A. hyemalis"),
  side = 3,
  adj = 0.5,
  cex = 1.2,
  line = 0.3
)
raster::plot(
  mProj_elvi$suitRaster,
  main = "",
  xlab = "Longitude",
  ylab = "",
  cex.lab = 1.5,
  col = rev(terrain.colors(100)),
  legend = TRUE
)
# points(elvi[,c("lon","lat")],pch=23,cex=0.3,col="grey")
# plot(garden_map_elvi,add=T,pch = 3,col="black",cex =2)
# plot(source_map_elvi,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("B",
      side = 3,
      adj = 0,
      cex = 1.25)
mtext(
  ~ italic("E. virginicus"),
  side = 3,
  adj = 0.5,
  cex = 1.2,
  line = 0.3
)
raster::plot(
  mProj_poa$suitRaster,
  xlab = "Longitude",
  ylab = "Latitude",
  cex.lab = 1.5,
  col = rev(terrain.colors(100))
)
# points(poa[,c("lon","lat")],pch=23,cex=0.3,col="grey")
# plot(garden_map_poau,add=T,pch = 3,col="black",cex =2)
# plot(source_map_poau,add=T,pch = 21,col="black",bg="red",cex =1)
mtext("C",
      side = 3,
      
      adj = 0,
      cex = 1.25)
mtext(
  ~ italic("P. autumnalis"),
  side = 3,
  adj = 0.5,
  cex = 1.2,
  line = 0.3
)
# legend(-119, 25.5,
#        legend=c( "GBIF occurences","Common garden sites"),
#        pch = c(23,3),
#        pt.cex=c(1.5,1.5),
#        col = c("grey50","black"),
#        pt.bg=c("grey","black"),
#        cex = 1,
#        bty = "n",
#        horiz = F ,
# )
dev.off()
