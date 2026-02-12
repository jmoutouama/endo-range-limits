# Project: 
# Purpose: Maps of Common garden experiment for Agrostis hyemalis, Elymus virginicus and Poa autumnalis across a climatic gradient. 
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 

# remove all objects and clear workspace
rm(list = ls(all=TRUE))
# load packages
library(sp)
library(raster)
library(terra)
library(RColorBrewer)
library(tidyverse)
library(dismo)
library(prism)
library(MESS)
library(mgcv)
library(maps)
library(dplyr)
library(tidyr)
library(reshape2)

# Climatic data----
## Data from PRISM---- 
# making a folder to store prism data
prism_set_dl_dir("/Users/jacobmoutouama/Documents/prism")
prism_archive_ls()

# Common garden, natural  and source populations locations ---- 
read.csv("https://www.dropbox.com/scl/fi/si346imz380lpdgo9yekr/Study_site.csv?rlkey=yiue42npkzzu9w8dggjr4fent&dl=1", stringsAsFactors = F) %>% 
  # dplyr::select(latitude,longitude) %>%
  unique() %>% 
  arrange(latitude)->garden ## common garden populations

read.csv("https://www.dropbox.com/scl/fi/a6cycehor6k6jp7xvbw3c/source_pop.csv?rlkey=chftk4ym5rsc7edjug90ie6o7&dl=1", stringsAsFactors = F) %>% 
  # dplyr::select(latitude,longitude) %>% 
  unique() %>% 
  arrange(latitude)->source ## source populations

# Agrostis hyemalis-----
#aghy_occ_raw <- gbif(genus="Agrostis",species="hyemalis",download=TRUE) 
aghy_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/ijl1i7964qxzcxvm7blz8/aghy_occ_raw.rds?rlkey=msbs3xzkjc719uld8yw7hb6cj&dl=1"))
head(aghy_occ_raw) 
aghy_occ_raw %>% 
  filter(!is.na(lat) & !is.na(lon) & basisOfRecord=="HUMAN_OBSERVATION")->aghy_occ
aghy_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-106.6458) &  as.numeric(longitude <=-94.02083) & as.numeric(latitude >=25.85417) &  as.numeric(latitude <=32.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->aghy1

aghy_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-94) &  as.numeric(longitude <=-92) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=32.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->aghy2

aghy_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-92) &  as.numeric(longitude <=-89.5) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=31) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->aghy3

aghy<-rbind(aghy1,aghy2,aghy3)

# Elymus virginicus----
#dir.create("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Data/occurence")
#elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE) 
elvi_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/0ssa5gepxyz28b7ykw1x8/elvi_occ_raw.rds?rlkey=4dx0q4lw2112droh73hmh7xte&dl=1"))
#head(elvi_occ_raw) 
elvi_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  filter(!is.na(lat) & !is.na(lon)) %>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-106.6458) &  as.numeric(longitude <=-94.02083) & as.numeric(latitude >=25.85417) &  as.numeric(latitude <=33.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->elvi1

elvi_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-94) &  as.numeric(longitude <=-92) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=32.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->elvi2

elvi_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-92) &  as.numeric(longitude <=-89.5) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=31) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->elvi3

elvi<-rbind(elvi1,elvi2,elvi3)

# Poa autumnalis ---
#poa_occ_raw <- gbif(genus="Poa",species="autumnalis",download=TRUE) 
poa_occ_raw<-readRDS(url("https://www.dropbox.com/scl/fi/oip7ndyf0d99rqxcqxb0q/poa_occ_raw.rds?rlkey=920uql1gd4gahnh8utw9fz96l&dl=1"))
#head(poa_occ_raw) 
poa_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  filter(!is.na(lat) & !is.na(lon))%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-106.6458) &  as.numeric(longitude <=-94.02083) & as.numeric(latitude >=25.85417) &  as.numeric(latitude <=33.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->poa1

poa_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-95) &  as.numeric(longitude <=-92) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=32.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->poa2

poa_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-96) &  as.numeric(longitude <=-90) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=31) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->poa3

poau<-rbind(poa1,poa2,poa3)
#class(poau)
# Define the 30-year window (relative to 2025)
#start_year <- 2025 - 29  # 1996

# Combine species into one table
dat_all_species <- bind_rows(
  mutate(aghy, species = "aghy"),
  mutate(elvi, species = "elvi"),
  mutate(poau, species = "poau")
)

# # Annual eastern range edge for each species
# annual_east_range_edge_1996_2025 <- dat_all_species %>%
#   filter(year >= start_year, year <= 2025) %>%
#   group_by(species, year) %>%
#   slice_max(longitude, n = 1, with_ties = FALSE) %>%
#   ungroup()
# Find westernmost (most negative longitude) per species
west_edge <- dat_all_species %>%
  group_by(species) %>%                           # handle each species separately
  slice_min(longitude, n = 1, with_ties = FALSE) %>%  # pick the minimum longitude
  ungroup() %>%
  dplyr::select(species, longitude, latitude)

# Check result
west_edge

# Save as RDS
#saveRDS(west_edge, file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/max_longitudes.rds")

# Georeferencing the occurences -----
garden %>% 
  filter(Species=="AGHY")->garden_aghy
garden %>% 
  filter(Species=="ELVI")->garden_elvi
garden %>% 
  filter(Species=="POAU")->garden_poau

# Order sites west-to-east
ordered_sites <- garden_aghy$site[order(garden_aghy$lon)]

source %>% 
  filter(Species=="AGHY")->source_aghy
source %>% 
  filter(Species=="ELVI")->source_elvi
source %>% 
  filter(Species=="POAU")->source_poau

sp::coordinates(aghy) <- ~ longitude + latitude
sp::coordinates(elvi) <- ~ longitude + latitude
sp::coordinates(poau) <- ~ longitude + latitude

sp::coordinates(garden_aghy) <- ~ longitude + latitude
sp::coordinates(garden_elvi) <- ~ longitude + latitude
sp::coordinates(garden_poau) <- ~ longitude + latitude

sp::coordinates(source_aghy) <- ~ longitude + latitude
sp::coordinates(source_elvi) <- ~ longitude + latitude
sp::coordinates(source_poau) <- ~ longitude + latitude

CRS1 <- CRS("+init=epsg:4326") # WGS 84

crs(aghy) <- CRS1
crs(elvi) <- CRS1
crs(poau) <- CRS1

crs(garden_aghy) <- CRS1
crs(garden_elvi) <- CRS1
crs(garden_poau) <- CRS1

crs(source_aghy) <- CRS1
crs(source_elvi) <- CRS1
crs(source_poau) <- CRS1

# Climatic data----
prism_summary <- readRDS(url("https://www.dropbox.com/scl/fi/x231nlm6rtm96uqffsv03/prism_means.rds?rlkey=w29usaquhd1u1w3u7guf1rgbp&dl=1"))
# Study area shapefile ----
study_area<-terra::vect("/Users/jacobmoutouama/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
study_area <- study_area[(study_area$State_Name %in% c("TEXAS","LOUISIANA")), ]
#plot(study_area)
# Clip the climatic rasters
tmean_annual <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1996:2025,resolution  = "800m"))))
crs(tmean_annual)<-CRS1
crop_tmean_annual <- terra::crop(tmean_annual, study_area,mask=TRUE)
# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual <- list()
for(y in 1996:2025){
  ppt_annual[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y,resolution  = "800m"))))
}
# Taking the mean of the cumulative precipitation values
ppt_annual_norm <- terra::mean(terra::rast(unlist(ppt_annual)))
crs(ppt_annual_norm)<-CRS1
crop_ppt_annual <- terra::crop(ppt_annual_norm, study_area,mask=TRUE)
col_precip <- terrain.colors(30)
col_precip_rev <- rev(col_precip)

# Fenced plots change herbivory that plants experience
## Data
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

dat24 <- read.csv("https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1")
dat24ini<-read.csv("https://www.dropbox.com/scl/fi/zncghunh1p9ull9j1jhkp/data_ini_2024.csv?rlkey=fkhbp3sdvb65va0gg4rwp4ra4&dl=1")
dat24 <- clean_tag(dat24)
dat24ini <-clean_tag(dat24ini)

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
# datini23 %>%
#   mutate(Tag_ID_num = as.numeric(Tag_ID)) %>%
#   filter(Tag_ID_num >= 75 & Tag_ID_num <= 85)
#dat24_sp_site_tag %>% filter(Tag_ID %in% missing_tags$Tag_ID)

dat25 <- read.csv("https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1")
# Ensure Tag_ID is character and trim whitespace
dat25 <- clean_tag(dat25)

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
datherbivory <- read.csv("https://www.dropbox.com/scl/fi/2gnlfozxpd2u9gprzp9oi/herbivory.csv?rlkey=sz2cloqxtbc6ou29j97l3t10f&dl=1", stringsAsFactors = F)

## Merge the demographic data with the herbivory data
dat23_herb <- left_join(x = datini23, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
dat24_herb <- left_join(x = dat24_sp_site_tag, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
dat25_herb <- left_join(x = dat25_plot, y = datherbivory, by = c("Site","Plot", "Species")) # Merge the demographic data with the herbivory data

dat23_herb_plot <- dat23_herb %>%
  group_by(Site, Plot, Herbivory) %>%
  summarise(
    prop_plants_damaged = mean(tiller_Herb_23 > 0, na.rm = TRUE),
    .groups = "drop"
  )

#  Summarize per Site × Herbivory
summary_23 <- dat23_herb_plot %>%
  filter(!is.na(Site) & !is.na(Herbivory)) %>%  # remove rows with missing grouping variables
  group_by(Site, Herbivory) %>%
  summarise(
    Mean = mean(prop_plants_damaged, na.rm = TRUE),
    SE   = sd(prop_plants_damaged, na.rm = TRUE) / sqrt(sum(!is.na(prop_plants_damaged))),
    .groups = "drop"
  ) %>%
  filter(!is.nan(Mean))  # remove groups with all NA


# Define colors
# Extract site codes and longitudes
site_order_df <- data.frame(
  site = garden_aghy$site_code,
  lon  = coordinates(garden_aghy)[,1]
)

# Order sites west → east
site_order_df <- site_order_df[order(site_order_df$lon), ]
ordered_sites <- site_order_df$site
# Reorder columns of mean_matrix and se_matrix

herb_levels <- c(0, 1)
colors <- c("lightgreen", "salmon")
names(colors) <- herb_levels

mean_matrix <- sapply(ordered_sites, function(s){
  sapply(herb_levels, function(h){
    v <- summary_23$Mean[summary_23$Site == s & summary_23$Herbivory == h]
    if (length(v) == 0) NA else v
  })
})

se_matrix <- sapply(ordered_sites, function(s){
  sapply(herb_levels, function(h){
    v <- summary_23$SE[summary_23$Site == s & summary_23$Herbivory == h]
    if (length(v) == 0) NA else v
  })
})

# Calculate mean annual precipitation per site
# mean_ppt_per_site <- prism_summary %>%
#   group_by(site) %>%
#   summarise(mean_ppt = mean(sum_ppt, na.rm = TRUE))

# Join with garden site coordinates
# garden_sites_ppt <- data.frame(
#   site = garden_aghy$site_code,
#   lon  = coordinates(garden_aghy)[,1],
#   mean_ppt  = mean_ppt_per_site$mean_ppt[match(garden_aghy$site_code, mean_ppt_per_site$site)]
# )
# 
# # Order sites west → east
# garden_sites_ppt_west_east <- garden_sites_ppt[order(garden_sites_ppt$lon), ]

# # Prepare barplot values and labels
# mean_ppt_west_east <- garden_sites_ppt_west_east$mean_ppt
# names(mean_ppt_west_east) <- garden_sites_ppt_west_east$site

# Maps (Figure 1) ----
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/clim_map1.pdf",
    width=9, height=10.5)

# Define layout as a 3x2 grid
layout_matrix <- matrix(c(1,2,
                          3,4,
                          5,6), nrow=3, ncol=2, byrow=TRUE)

layout(layout_matrix, heights=c(1,1,0.9))  # adjust heights if you want barplot bigger

# Set default margins
par(mar=c(3,0,4,1), oma=c(0,2,1,0))

### Panel A
plot(crop_ppt_annual, xlab="Longitude", ylab="Latitude", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=TRUE)
plot(aghy, add=TRUE, pch=23, col="grey50", bg="grey", cex=0.55)
plot(garden_aghy, add=TRUE, pch=3, col="black", cex=2)
plot(source_aghy, add=TRUE, pch=21, col="black", bg="red", cex=1)
mtext(~italic("Agrostis hyemalis"), side=3, adj=0.5, cex=1.2, line=0.2)
mtext("A", side=3, adj=0, cex=1.25, line=0.2)
mtext("ppt (mm)", side=3, adj=1.21, cex=0.6, line=-1.2)
map.scale(
  x = -95,       # longitude position of scale bar
  y = 28,        # latitude position of scale bar
  relwidth = 0.2,  # relative width of the scale bar
  metric = TRUE,   # use metric units (km)
  cex = 0.8 ,      # size of text
  ratio = FALSE   # removes the 1:16 ratio label
)
### Panel B
plot(crop_ppt_annual, xlab="Longitude", ylab="", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=TRUE)
plot(elvi, add=TRUE, pch=23, col="grey50", bg="grey", cex=0.55)
plot(garden_elvi, add=TRUE, pch=3, col="black", cex=2)
plot(source_elvi, add=TRUE, pch=21, col="black", bg="red", cex=1)
mtext(~italic("Elymus virginicus"), side=3, adj=0.5, cex=1.2, line=0.2)
mtext("B", side=3, adj=0, cex=1.25, line=0.2)
mtext("ppt (mm)", side=3, adj=1.21, cex=0.6, line=-1.2)
map.scale(
  x = -95,       # longitude position of scale bar
  y = 28,        # latitude position of scale bar
  relwidth = 0.2,  # relative width of the scale bar
  metric = TRUE,   # use metric units (km)
  cex = 0.8 ,      # size of text
  ratio = FALSE   # removes the 1:16 ratio label
)
### Panel C
par(mar=c(0,3,3.75,1))
plot(crop_ppt_annual, xlab="Longitude", ylab="Latitude", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=TRUE)
plot(poau, add=TRUE, pch=23, col="grey50", bg="grey", cex=0.55)
plot(garden_poau, add=TRUE, pch=3, col="black", cex=2)
plot(source_poau, add=TRUE, pch=21, col="black", bg="red", cex=1)
mtext(~italic("Poa autumnalis"), side=3, adj=0.5, cex=1.2, line=0.7)
mtext("C", side=3, adj=0, cex=1.25, line=0.3)
mtext("ppt (mm)", side=3, adj=1.21, cex=0.6, line=-1.2)
map.scale(
  x = -95,       # longitude position of scale bar
  y = 28,        # latitude position of scale bar
  relwidth = 0.2,  # relative width of the scale bar
  metric = TRUE,   # use metric units (km)
  cex = 0.8 ,      # size of text
  ratio = FALSE   # removes the 1:16 ratio label
)
legend(
  -105, 28,
  legend = c("GBIF occurrence", "Experimental site", "Source"),
  pch    = c(23, 3, 21),
  pt.bg  = c("grey", NA, "red"),
  col    = c("grey50", "black", "black"),
  pt.cex = c(0.55, 2, 1),
  bty    = "n",
  cex    = 0.9
)

### Panel D (barplot per site per year)
### Panel D (barplot per site per year)
par(mar=c(6,4,4,1))

# Convert prism_summary to a wide matrix: rows = sites, columns = years
prism_summary_census <- readRDS(url("https://www.dropbox.com/scl/fi/fj6aqhej58k9fjt9v02ap/climate_census_years.rds?rlkey=28dsn6b9o3x06wd1c2ehs5ama&dl=1"))
unique(prism_summary_census$census_year)
# Example: create ppt_matrix for one species, e.g., AGHY
ppt_matrix <- prism_summary_census %>%
  dplyr::group_by(site, census_year) %>%
  dplyr::summarise(cum_ppt = mean(cum_ppt), .groups = "drop") %>%  # average duplicates
  tidyr::pivot_wider(
    id_cols = site,
    names_from = census_year,
    values_from = cum_ppt
  ) %>%
  column_to_rownames("site") %>%
  as.matrix()

# Reorder sites west → east
ppt_matrix <- ppt_matrix[ordered_sites, ]

# Colors for years
#year_colors <- c("2023" = "lightgreen", "2024" = "salmon", "2025" = "goldenrod")
year_colors <- c( "2024" = "salmon",       # pink-orange
                 "2025" = "lightgreen")  # blue


# Create barplot
bp <- barplot(
  t(ppt_matrix),              # transpose so bars are grouped by site
  beside = TRUE,              # grouped bars
  col = year_colors[colnames(ppt_matrix)],
  names.arg = rownames(ppt_matrix),
  ylim = c(0, max(ppt_matrix, na.rm=TRUE) * 1.1),
  las = 2,
  cex.lab = 1.2,
  cex.names = 0.75,
  xlab = "Sites",
  ylab = "Annual Precipitation (mm)"
)

mtext("D", side=3, adj=-0.06, cex=1.25, line=0.5)
box()

# Add legend
legend(
  "topleft",
  legend = colnames(ppt_matrix),
  fill = year_colors[colnames(ppt_matrix)],
  bty = "n",
  cex = 1
)

### Panel E
#par(mar=c(0, 0, 0, 0)) 
par(mar=c(4,2,2,1))  
#plot(0, 0, type="n", xlim=c(0,3.6), ylim=c(0,1.5), axes=FALSE, xlab="", ylab="", main="", asp=1)
plot(0, 0, type="n", xlim=c(0,3.6), ylim=c(0,1.8), axes=FALSE, xlab="", ylab="", main="", asp=1)

mtext("E", side=3, adj=0.07, cex=1.25, line=-1.75)

# fenced
#rect(0.1, 0.1, 1.6, 1.6, border="black", lwd=2)
rect(0.1, 0.2, 1.6, 1.5, border="black", lwd=2)
plants1_x <- rep(seq(0.375, 1.125, length.out=4), times=4)[-1]
plants1_y <- rep(seq(0.375, 1.125, length.out=4), each=4)[-1]
# points(plants1_x+0.1, plants1_y+0.1, pch=22, col="black", cex=0.75,bg = "black")
set.seed(13)  
signs1 <- sample(c("+", "−"), length(plants1_x), replace = TRUE)
text(plants1_x + 0.1, plants1_y + 0.1, labels = signs1, cex = 1.1, font = 2)

mtext("Fenced", side=3, at=0.85, line=-4, cex=1.2)

# unfenced
# rect(1.9, 0.1, 3.4, 1.6, border=NA)
# segments(1.9,0.1,3.4,0.1, col="black", lwd=2)
# segments(1.9,1.6,3.4,1.6, col="black", lwd=2)
rect(1.9, 0.2, 3.4, 1.5, border=NA)
segments(1.9,0.2,3.4,0.2, col="black", lwd=2)
segments(1.9,1.5,3.4,1.5, col="black", lwd=2)
plants2_x <- rep(seq(2.375, 3.125, length.out=4), times=4)[-1]
plants2_y <- rep(seq(0.375, 1.125, length.out=4), each=4)[-1]
# points(plants2_x-2+1.9, plants2_y+0.1, pch=22, col="black", cex=0.75,bg = "black")
set.seed(13)
signs2 <- sample(c("+", "−"), length(plants2_x), replace = TRUE)
text(plants2_x - 2 + 1.9, plants2_y + 0.1, labels = signs2, cex = 1.1, font = 2)

mtext("Unfenced", side=3, at=2.65, line=-4, cex=1.2)

### Panel F
par(mar=c(4,2,2,1))  # same for both E and F
mtext("F", side = 3, adj = 1.1, cex = 1.25, line = 0.5)

bp <- barplot(
  height = mean_matrix,
  beside = TRUE,
  names.arg = ordered_sites,
  col = c("grey80","grey15"),  # use herbivory colors
  ylim = c(0, max(mean_matrix + se_matrix, na.rm = TRUE) * 1.25),
  xlab = "Sites",
  ylab = "Proportion of damaged plants",
  las = 2,
  cex.lab = 1.2,
  cex.names = 0.75
)

# Add error bars
arrows(
  x0 = bp,
  y0 = mean_matrix - se_matrix,
  x1 = bp,
  y1 = mean_matrix + se_matrix,
  angle = 90, code = 3, length = 0.07
)
box()
# Add clear legend
legend(
  "topright",
  legend = c("Unfenced", "Fenced"),
  fill = c("grey80","grey15"),
  bty = "n",
  cex = 1.2
)

dev.off()

## Herbivory for 2024 and 2025
# Define herbivory levels and colors
herb_levels <- c(1,0)  # 1 = Fenced, 0 = Unfenced
colors <- c("grey15","grey80")
names(colors) <- herb_levels

# Function to summarize data per year
summarize_herb <- function(dat, year_col, tiller_col){
  dat %>%
    group_by(Site, Plot, Herbivory) %>%
    summarise(
      prop_plants_damaged = mean(.data[[tiller_col]] > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Site, Herbivory) %>%
    summarise(
      Mean = mean(prop_plants_damaged, na.rm = TRUE),
      SE   = sd(prop_plants_damaged, na.rm = TRUE)/sqrt(n()),
      .groups = "drop"
    ) %>%
    filter(!is.na(Site) & !is.na(Herbivory))
}

# Summarize 2024
summary_24 <- summarize_herb(dat24_herb, "date_24", "tiller_herb_24")
# Summarize 2025
summary_25 <- summarize_herb(dat25_herb, "date_25", "tiller_herb_25")  # assuming dat25_herb exists

# Combine into wide matrices, ordered by ordered_sites
make_matrix <- function(summary_df){
  mat <- sapply(ordered_sites, function(s){
    sapply(herb_levels, function(h){
      v <- summary_df$Mean[summary_df$Site==s & summary_df$Herbivory==h]
      if(length(v)==0) NA else v
    })
  })
  rownames(mat) <- herb_levels
  mat
}

make_se_matrix <- function(summary_df){
  mat <- sapply(ordered_sites, function(s){
    sapply(herb_levels, function(h){
      v <- summary_df$SE[summary_df$Site==s & summary_df$Herbivory==h]
      if(length(v)==0) NA else v
    })
  })
  rownames(mat) <- herb_levels
  mat
}

mean_matrix_24 <- make_matrix(summary_24)
se_matrix_24 <- make_se_matrix(summary_24)

mean_matrix_25 <- make_matrix(summary_25)
se_matrix_25 <- make_se_matrix(summary_25)

# Export combined plot
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/prop_damage_barplot_24_25.pdf",
    width = 8, height = 4)

par(mfrow=c(1,2), mar=c(5,4,3,1))  # side-by-side panels

# --- Panel 2024 ---
bp1 <- barplot(
  height = mean_matrix_24,
  beside = TRUE,
  names.arg = ordered_sites,
  col = colors,
  ylim = c(0, max(c(mean_matrix_24 + se_matrix_24,
                    mean_matrix_25 + se_matrix_25), na.rm = TRUE) * 1.2),
  ylab = "Proportion of damaged plants",
  xlab = "Sites",
  main = "2024",
  las = 2,
  cex.names = 0.75
)
arrows(
  x0 = bp1,
  y0 = mean_matrix_24 - se_matrix_24,
  x1 = bp1,
  y1 = mean_matrix_24 + se_matrix_24,
  angle = 90, code = 3, length = 0.05
)
box()
# --- Panel 2025 ---
bp2 <- barplot(
  height = mean_matrix_25,
  beside = TRUE,
  names.arg = ordered_sites,
  col = colors,
  ylim = c(0, max(c(mean_matrix_24 + se_matrix_24,
                    mean_matrix_25 + se_matrix_25), na.rm = TRUE) * 1.2),
  ylab = "",
  xlab = "Sites",
  main = "2025",
  las = 2,
  cex.names = 0.75
)
arrows(
  x0 = bp2,
  y0 = mean_matrix_25 - se_matrix_25,
  x1 = bp2,
  y1 = mean_matrix_25 + se_matrix_25,
  angle = 90, code = 3, length = 0.05
)
box()
# Shared legend
legend("topleft", legend = c("Fenced","Unfenced"), fill = colors, bty = "n", cex = 0.9)

dev.off()

## Herbivory per species  and per year
# Define herbivory levels and colors
herb_levels <- c(0,1)  # 0 = Fenced, 1 = Unfenced
colors <- c("grey80","grey15")
names(colors) <- herb_levels

# Function to summarize herbivory per dataset
summarize_herb <- function(dat, tiller_col){
  # Remove rows with NA in key columns
  dat <- dat %>% filter(!is.na(Species), !is.na(Site), !is.na(Herbivory))
  
  if(!tiller_col %in% colnames(dat)){
    stop(paste("Column", tiller_col, "not found in dataset"))
  }
  
  # Compute proportion damaged per Plot
  summary_plot <- dat %>%
    group_by(Species, Site, Plot, Herbivory) %>%
    summarise(
      prop_plants_damaged = mean(.data[[tiller_col]] > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Summarize per Species × Site × Herbivory
  summary_site <- summary_plot %>%
    group_by(Species, Site, Herbivory) %>%
    summarise(
      Mean = mean(prop_plants_damaged, na.rm = TRUE),
      SE   = if(n() > 1) sd(prop_plants_damaged, na.rm = TRUE)/sqrt(n()) else 0,
      .groups = "drop"
    )
  
  return(summary_site)
}

# Summarize all years
summary_23 <- summarize_herb(dat23_herb, "tiller_Herb_23")
summary_24 <- summarize_herb(dat24_herb, "tiller_herb_24")
summary_25 <- summarize_herb(dat25_herb, "tiller_herb_25")

# List of summaries and years
summary_list <- list("2023"=summary_23, "2024"=summary_24, "2025"=summary_25)
years <- names(summary_list)

# Function to create mean & SE matrices per species
make_matrices <- function(summary_df, species){
  df <- summary_df %>% filter(Species == species)
  
  mean_matrix <- sapply(ordered_sites, function(s){
    sapply(herb_levels, function(h){
      v <- df$Mean[df$Site==s & df$Herbivory==h]
      if(length(v)==0) NA else v
    })
  })
  rownames(mean_matrix) <- herb_levels
  
  se_matrix <- sapply(ordered_sites, function(s){
    sapply(herb_levels, function(h){
      v <- df$SE[df$Site==s & df$Herbivory==h]
      if(length(v)==0) NA else v
    })
  })
  rownames(se_matrix) <- herb_levels
  
  list(mean=mean_matrix, se=se_matrix)
}

# PDF output
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/prop_damage_species_years.pdf",
    width = 12, height = 8)

# Layout: rows = species, cols = years
species_list <- c("AGHY","ELVI","POAU")
par(mfrow=c(length(species_list), length(years)), mar=c(5,4,3,1))

# Loop over species and years
for(species in species_list){
  for(yr in years){
    mats <- make_matrices(summary_list[[yr]], species)
    mean_matrix <- mats$mean
    se_matrix <- mats$se
    
    bp <- barplot(
      height = mean_matrix,
      beside = TRUE,
      names.arg = ordered_sites,
      col = colors,
      ylim = c(0, max(unlist(lapply(summary_list, function(s){ 
        max(s$Mean + s$SE, na.rm=TRUE) 
      })), na.rm=TRUE) * 1.2),
      ylab = ifelse(yr=="2023", "Proportion of damaged plants",""),
      xlab = "Sites",
      main = paste(species, "-", yr),
      las = 2,
      cex.names = 0.75
    )
    
    arrows(
      x0 = bp,
      y0 = mean_matrix - se_matrix,
      x1 = bp,
      y1 = mean_matrix + se_matrix,
      angle = 90, code = 3, length = 0.05
    )
    box()
    # Add legend only on first panel
    if(species == species_list[1] & yr==years[1]){
      legend("topright", legend=c("Unfenced","Fenced"), fill=colors, bty="n", cex=1)
    }
  }
}

dev.off()

# Load climate census data
prism_summary_census <- readRDS(url(
  "https://www.dropbox.com/scl/fi/fj6aqhej58k9fjt9v02ap/climate_census_years.rds?rlkey=28dsn6b9o3x06wd1c2ehs5ama&dl=1"
))

# Ensure site is a factor with desired west-to-east order
ordered_sites <- c("LAF", "HUN", "COL", "BAS", "BFL", "KER", "SON")
prism_summary_census$site <- factor(prism_summary_census$site, levels = ordered_sites)

# Aggregate mean temperature by site and year
temp_data <- aggregate(mean_tmean ~ site + census_year, data = prism_summary_census, mean)

# Reshape into matrix: rows = years, columns = sites
temp_matrix <- dcast(temp_data, census_year ~ site, value.var = "mean_tmean")
rownames(temp_matrix) <- temp_matrix$census_year
temp_matrix <- as.matrix(temp_matrix[, -1, drop=FALSE])  # remove year column

# Years (rows)
years <- as.character(temp_matrix[,1])  # optional if you want to reference

# Colors for years (color-blind friendly)
cols <- c("#0072B2", "#D55E00", "#009E73")  # adjust depending on number of years
num_years <- nrow(temp_matrix)
if(length(cols) < num_years) cols <- rainbow(num_years)

# Export temperature-only barplot
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/temp_census_years.pdf",
    width = 8, height = 4)

par(mar=c(5,6,3,2), cex.lab=1.5, cex.axis=1.3)

# Transpose so years are grouped side-by-side
barplot(
  temp_matrix,        # transpose to get years beside each site
  beside = TRUE,
  col = cols[1:nrow(temp_matrix)],
  names.arg = ordered_sites,
  ylim = c(0, max(temp_matrix, na.rm = TRUE) * 1.8),
  xlab = "Site",
  ylab = "Mean Temperature (°C)",
  las = 1
)
legend("topleft", legend = rownames(temp_matrix), fill = cols[1:nrow(temp_matrix)],
       title = "Year", cex = 1.2)
box()
#mtext("Temperature", side = 3, adj = 0, line = 0.25, cex = 1.5)

dev.off()
