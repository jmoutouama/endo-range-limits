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
head(elvi_occ_raw) 
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
head(poa_occ_raw) 
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
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1", stringsAsFactors = F)
dat23 <- read.csv("https://www.dropbox.com/scl/fi/fkwm0dan6nx2eaeyxjrjw/census2023.csv?rlkey=hy9209t53j9n7vxhta7axl5jk&dl=1", stringsAsFactors = F)
dat24 <- read.csv("https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1", stringsAsFactors = F)
dat25<-read.csv("https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1", stringsAsFactors = F)
datherbivory <- read.csv("https://www.dropbox.com/scl/fi/2gnlfozxpd2u9gprzp9oi/herbivory.csv?rlkey=sz2cloqxtbc6ou29j97l3t10f&dl=1", stringsAsFactors = F)
## Match the data with initial data 
dat23_sp <- dat23 %>%
  left_join(
    datini %>%
      dplyr::select(Site, Species, Plot, Position, Tag_ID, Population, 
                    GreenhouseID, Clone, Endo),
    by = "Tag_ID"
  ) %>%
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C")))
#Combine datini and dat23
combined_data <- bind_rows(datini[,c("Site","Species","Plot","Tag_ID","Population","Endo")], dat23_sp[,c("Site","Species","Plot","Tag_ID","Population","Endo" )])%>% 
  distinct(Tag_ID, .keep_all = TRUE)

dat24_sp<-left_join(x = dat24, y = combined_data, by = c("Tag_ID")) %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit"))) %>% 
  filter(!is.na(Species))
## Merge the demographic data with the herbivory data
dat23_herb <- left_join(x = dat23_sp, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
dat24_herb <- left_join(x = dat24_sp, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
#dat25_herb <- left_join(x = dat25, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data

dat23_herb_plot <- dat23_herb %>%
  group_by(Site, Plot, Herbivory) %>%
  summarise(
    prop_plants_damaged = mean(tiller_Herb_23 > 0, na.rm = TRUE),
    .groups = "drop"
  )

#  Summarize per Site × Herbivory
summary_23 <- dat23_herb_plot %>%
  group_by(Site, Herbivory) %>%
  summarise(
    Mean = mean(prop_plants_damaged, na.rm = TRUE),
    SE   = sd(prop_plants_damaged, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

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
mean_ppt_per_site <- prism_summary %>%
  group_by(site) %>%
  summarise(mean_ppt = mean(sum_ppt, na.rm = TRUE))

# Join with garden site coordinates
garden_sites_ppt <- data.frame(
  site = garden_aghy$site_code,
  lon  = coordinates(garden_aghy)[,1],
  mean_ppt  = mean_ppt_per_site$mean_ppt[match(garden_aghy$site_code, mean_ppt_per_site$site)]
)

# Order sites west → east
garden_sites_ppt_west_east <- garden_sites_ppt[order(garden_sites_ppt$lon), ]

# Prepare barplot values and labels
mean_ppt_west_east <- garden_sites_ppt_west_east$mean_ppt
names(mean_ppt_west_east) <- garden_sites_ppt_west_east$site

# Maps (Figure 1) ----
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/clim_map.pdf",
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

### Panel D (barplot ordered largest → smallest)
par(mar=c(6,4,4,1))
barplot(
  mean_ppt_west_east,
  col="#E69F00",
  xlab="Sites",
  ylab="Cumulative Precipitation (mm)",
  ylim=c(0, max(mean_ppt_west_east)*1.1),
  las=2,
  cex.lab=1.2,
  cex.names=0.75
)
mtext("D", side=3, adj=-0.06, cex=1.25, line=0.5)
box()
### Panel E
par(mar=c(0, 0, 0, 0))  
plot(0, 0, type="n", xlim=c(0,3.6), ylim=c(0,1.5), axes=FALSE, xlab="", ylab="", main="", asp=1)
mtext("E", side=3, adj=0.07, cex=1.25, line=-1.75)

# fenced
rect(0.1, 0.1, 1.6, 1.6, border="black", lwd=2)
# plants1_x <- rep(seq(0.375, 1.125, length.out=4), times=4)[-1]
# plants1_y <- rep(seq(0.375, 1.125, length.out=4), each=4)[-1]
# points(plants1_x+0.1, plants1_y+0.1, pch=22, col="black", cex=0.75,bg = "black")
set.seed(13)  
signs1 <- sample(c("+", "−"), length(plants1_x), replace = TRUE)
text(plants1_x + 0.1, plants1_y + 0.1, labels = signs1, cex = 1.1, font = 2)

mtext("Fenced", side=3, at=0.85, line=-4, cex=0.9)

# unfenced
rect(1.9, 0.1, 3.4, 1.6, border=NA)
segments(1.9,0.1,3.4,0.1, col="black", lwd=2)
segments(1.9,1.6,3.4,1.6, col="black", lwd=2)
# plants2_x <- rep(seq(2.375, 3.125, length.out=4), times=4)[-1]
# plants2_y <- rep(seq(0.375, 1.125, length.out=4), each=4)[-1]
# points(plants2_x-2+1.9, plants2_y+0.1, pch=22, col="black", cex=0.75,bg = "black")
set.seed(13)
signs2 <- sample(c("+", "−"), length(plants2_x), replace = TRUE)
text(plants2_x - 2 + 1.9, plants2_y + 0.1, labels = signs2, cex = 1.1, font = 2)

mtext("Unfenced", side=3, at=2.65, line=-4, cex=0.9)

### Panel F
par(mar = c(8, 4, 2, 1))  # bottom, left, top, right
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

# Add clear legend
legend(
  "topright",
  legend = c("Unfenced", "Fenced"),
  fill = c("grey80","grey15"),
  bty = "n",
  cex = 1.2
)

dev.off()

## 2024 
dat24_herb_plot <- dat24_herb %>%
  group_by(Site, Plot, Herbivory) %>%
  summarise(
    prop_plants_damaged = mean(tiller_herb_24 > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Summarize per Site × Herbivory
summary_24 <- dat24_herb_plot %>%
  group_by(Site, Herbivory) %>%
  summarise(
    Mean = mean(prop_plants_damaged, na.rm = TRUE),
    SE   = sd(prop_plants_damaged, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

# Define colors
herb_levels <- c(0,1)  # 0 = Fenced, 1 = Unfenced
colors <- c("lightgreen", "salmon")
names(colors) <- herb_levels

# Create matrix for barplot (rows = Herbivory, columns = Site)
site_ids <- unique(summary_24$Site)
mean_matrix <- sapply(site_ids, function(s){
  sapply(herb_levels, function(h){
    v <- summary_24$Mean[summary_24$Site==s & summary_24$Herbivory==h]
    if(length(v)==0) NA else v
  })
})

# Create SE matrix
se_matrix <- sapply(site_ids, function(s){
  sapply(herb_levels, function(h){
    v <- summary_24$SE[summary_24$Site==s & summary_24$Herbivory==h]
    if(length(v)==0) NA else v
  })
})

# Export barplot as PDF
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/prop_damage_barplot24.pdf",
    width = 5, height = 4)

bp <- barplot(
  height = mean_matrix,
  beside = TRUE,
  names.arg = site_ids,
  col = colors,
  ylim = c(0, max(mean_matrix + se_matrix, na.rm=TRUE) * 1.2),
  ylab = "Proportion of damaged plants",
  xlab="Sites",
  main = "",
  las = 2,
  cex.names = 0.6
)
# Add error bars
arrows(
  x0 = bp,
  y0 = mean_matrix - se_matrix,
  x1 = bp,
  y1 = mean_matrix + se_matrix,
  angle = 90, code = 3, length = 0.05
)

# Add legend
legend(0,0.95, legend = c("Unfenced", "Fenced"), fill = colors, bty = "n", cex = 0.8)

dev.off()


# Summarize data per plot and herbivory
dat24_herb_plot_by_species <- dat24_herb %>%
  group_by(Species, Site, Plot, Herbivory) %>%
  summarise(
    prop_plants_damaged = mean(tiller_herb_24 > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Summarize per Site × Herbivory × Species
summary_24 <- dat24_herb_plot_by_species %>%
  group_by(Species, Site, Herbivory) %>%
  summarise(
    Mean = mean(prop_plants_damaged, na.rm = TRUE),
    SE   = sd(prop_plants_damaged, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

# Define colors and herbivory levels (1 = Fenced, 0 = Unfenced)
herb_levels <- c(1, 0)  # reordered
colors <- c("salmon","lightgreen")
names(colors) <- herb_levels

# Get species list
species_list <- unique(summary_24$Species)

# Export barplots as a multi-page PDF
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/prop_damage_barplot24_by_species.pdf",
    width = 5, height = 4)

for (sp in species_list) {
  
  sp_data <- summary_24 %>% filter(Species == sp)
  site_ids <- unique(sp_data$Site)
  
  # Create numeric mean matrix (rows = Herbivory, columns = Site)
  mean_matrix <- matrix(
    as.numeric(sapply(site_ids, function(s){
      sapply(herb_levels, function(h){
        v <- sp_data$Mean[sp_data$Site==s & sp_data$Herbivory==h]
        if(length(v)==0) NA else v
      })
    })),
    nrow = length(herb_levels), byrow = FALSE
  )
  
  # SE matrix
  se_matrix <- matrix(
    as.numeric(sapply(site_ids, function(s){
      sapply(herb_levels, function(h){
        v <- sp_data$SE[sp_data$Site==s & sp_data$Herbivory==h]
        if(length(v)==0) NA else v
      })
    })),
    nrow = length(herb_levels), byrow = FALSE
  )
  
  # Create barplot with italic species name
  bp <- barplot(
    height = mean_matrix,
    beside = TRUE,
    names.arg = site_ids,
    col = colors,
    ylim = c(0, max(mean_matrix + se_matrix, na.rm=TRUE) * 1.2),
    ylab = "Proportion of damaged plants",
    xlab = "Sites",
    main = substitute(italic(sp_name), list(sp_name = sp)),
    las = 2,
    cex.names = 0.6
  )
  
  # Add error bars
  arrows(
    x0 = bp,
    y0 = mean_matrix - se_matrix,
    x1 = bp,
    y1 = mean_matrix + se_matrix,
    angle = 90, code = 3, length = 0.05
  )
  
  # Add legend (Fenced = 1, Unfenced = 0)
  legend("topright", legend = c("Fenced", "Unfenced"), fill = colors, bty = "n", cex = 0.8)
}

dev.off()

# Climatic census ----
prism_summary_census <- readRDS(url("https://www.dropbox.com/scl/fi/fj6aqhej58k9fjt9v02ap/climate_census_years.rds?rlkey=28dsn6b9o3x06wd1c2ehs5ama&dl=1"))

# Aggregate cum_ppt by site and year
# Aggregate cum_ppt by site and year
site_data <- aggregate(cum_ppt ~ site + census_year, data = prism_summary_census, sum)

# Ensure site is a factor with the desired order
site_data$site <- factor(site_data$site, levels = c("LAF", "HUN", "COL", "BAS", "BFL", "KER", "SON"))

# Reshape into matrix: rows = sites, columns = years
library(reshape2)
bar_matrix <- dcast(site_data, site ~ census_year, value.var = "cum_ppt", fill = 0)
rownames(bar_matrix) <- bar_matrix$site
bar_matrix <- as.matrix(bar_matrix[, c("2024", "2025")])

# Side-by-side bar plot
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/climate_census_years.pdf",
    width = 5, height = 4)
barplot(
  t(bar_matrix),             # transpose so 2024 and 2025 are beside each other
  beside = TRUE,
  col = c("skyblue", "orange"),
  legend.text = c("2024", "2025"),
  args.legend = list(x = "topright"),
  main = "",
  ylim=c(0,7000),
  xlab = "Sites",
  ylab = "Cumulative Precipitation (mm)",
  las = 1
)
dev.off()

# Order sites west-to-east
ordered_sites <- garden_sites_ppt$site[order(garden_sites_ppt$lon)]

# Years
years <- sort(unique(prism_summary$year))

# Create matrices for barplot (years as rows, sites as columns) in west-east order
ppt_mat <- sapply(ordered_sites, function(s) prism_summary$sum_ppt[prism_summary$site==s])
temp_mat <- sapply(ordered_sites, function(s) prism_summary$mean_temp[prism_summary$site==s])
rownames(ppt_mat) <- years
rownames(temp_mat) <- years
colnames(ppt_mat) <- ordered_sites
colnames(temp_mat) <- ordered_sites

# Colors for years (color-blind friendly)
cols <- c("#0072B2", "#D55E00", "#009E73")
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/temp_pppt_site_year.pdf",
    width=10, height=12)

# Layout with custom heights (ppt taller, temp shorter)
layout(matrix(1:2, ncol=1), heights=c(2.5, 1))  # 2.5:1 ratio

# Increase label sizes
par(mar=c(5,6,2,2), cex.lab=1.5, cex.axis=1.3)  # margins, axis and label sizes

# --- Precipitation ---
barplot(ppt_mat, beside=TRUE, col=cols,
        ylim=c(0,max(ppt_mat)*1.2),
        names.arg=ordered_sites,
        xlab="Site", ylab="Annual Precipitation (mm)")
legend("topleft", legend=years, fill=cols, title="Year", cex=1.2)
box()
mtext("A", side=3, adj=0, line=0.25, cex=2)

# --- Temperature ---
barplot(temp_mat, beside=TRUE, col=cols,
        ylim=c(0,max(temp_mat)*1.2),
        names.arg=ordered_sites,
        xlab="Site", ylab="Mean Temperature (°C)")
box()
mtext("B", side=3, adj=0, line=0.25, cex=2)

dev.off()
