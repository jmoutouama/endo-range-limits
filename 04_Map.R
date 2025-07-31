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
# Climatic data----
## Data from PRISM---- 
# making a folder to store prism data
options(prism.path = "/Users/jacobmoutouama/Documents/prism/")
# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
# get_prism_monthlys(type = "tmean", years = 1994:2024, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "ppt", years = 1994:2024, mon = 1:12, keepZip = FALSE)
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
  filter(year %in% (1901:2024) & as.numeric(longitude >=-94) &  as.numeric(longitude <=-92) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=32.5) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->poa2

poa_occ_raw %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(longitude=lon,latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(longitude >=-92) &  as.numeric(longitude <=-89.5) & as.numeric(latitude >=29.5) &  as.numeric(latitude <=31) & country=="United States") %>% 
  unique() %>% 
  arrange(latitude)->poa3

poau<-rbind(poa1,poa2,poa3)

# Georeferencing the occurences -----
garden %>% 
  filter(Species=="AGHY")->garden_aghy
garden %>% 
  filter(Species=="ELVI")->garden_elvi
garden %>% 
  filter(Species=="POAU")->garden_poau

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

# Climatic and distance data----
distance_summary <- readRDS(url("https://www.dropbox.com/scl/fi/6rnf3ahwave2p7gaf9cqd/distance_species.rds?rlkey=cs4rtiee8h611brv9z1jv4prn&dl=1"))
distance_summary$geo_distance<-distance_summary$geo_distance/1000
distance_summary_ordered <- distance_summary[order(distance_summary$longitude), ]
# Study area shapefile ----
study_area<-terra::vect("/Users/jacobmoutouama/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
study_area <- study_area[(study_area$State_Name %in% c("TEXAS","LOUISIANA")), ]
#plot(study_area)
# Clip the climatic rasters
tmean_annual <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1994:2024))))
crs(tmean_annual)<-CRS1
crop_tmean_annual <- terra::crop(tmean_annual, study_area,mask=TRUE)
# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual <- list()
for(y in 1993:2023){
  ppt_annual[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
}
# Taking the mean of the cumulative precipitation values
ppt_annual_norm <- terra::mean(terra::rast(unlist(ppt_annual)))
crs(ppt_annual_norm)<-CRS1
crop_ppt_annual <- terra::crop(ppt_annual_norm, study_area,mask=TRUE)
col_precip <- terrain.colors(30)
col_precip_rev <- rev(col_precip)


# Maps (Figure 1) ----
pdf("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/clim_map.pdf",width=9,height=8)
op <- par(mfrow = c(2,2), mar=c(0,1,3.75,1), oma = c(0, 2, 1, 0)) 

# First plot (A)
plot(crop_ppt_annual, xlab="Longitude", ylab="Latitude", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=T)
plot(aghy, add=T, pch = 23, col="grey50", bg="grey", cex=0.55)
plot(garden_aghy, add=T, pch = 3, col="black", cex=2)
plot(source_aghy, add=T, pch = 21, col="black", bg="red", cex=1)
mtext(~ italic("Agrostis hyemalis"), side = 3, adj = 0.5, cex=1.2, line=0.2)
mtext("A", side = 3, adj = 0, cex=1.25, line=0.2)
legend(-106, 28, 
       legend = c("GBIF occurences", "Common garden sites", "Source populations"),
       pch = c(23, 3, 21),
       pt.cex = c(0.55, 1, 1),
       col = c("grey50", "black", "black"),
       pt.bg = c("grey", "black", "red"),
       cex = 0.7, 
       bty = "n", 
       horiz = F)

# Second plot (B)
plot(crop_ppt_annual, xlab="Longitude", ylab="", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=T)
plot(elvi, add=T, pch = 23, col="grey50", bg="grey", cex=0.55)
plot(garden_elvi, add=T, pch = 3, col="black", cex=2)
plot(source_elvi, add=T, pch = 21, col="black", bg="red", cex=1)
mtext(~ italic ("Elymus virginicus"), side = 3, adj = 0.5, cex=1.2, line=0.2)
mtext("B", side = 3, adj = 0, cex=1.25, line=0.2)
legend(-106, 28, 
       legend = c("GBIF occurences", "Common garden sites", "Source populations"),
       pch = c(23, 3, 21),
       pt.cex = c(0.55, 1, 1),
       col = c("grey50", "black", "black"),
       pt.bg = c("grey", "black", "red"),
       cex = 0.7, 
       bty = "n", 
       horiz = F)

# Third plot (C)
par(mar=c(0,3,3.75,1))
plot(crop_ppt_annual, xlab="Longitude", ylab="Latitude", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=T)
plot(poau, add=T, pch = 23, col="grey50", bg="grey", cex=0.55)
plot(garden_poau, add=T, pch = 3, col="black", cex=2)
plot(source_poau, add=T, pch = 21, col="black", bg="red", cex=1)
mtext( ~ italic("Poa autumnalis"), side = 3, adj = 0.5, cex=1.2, line=0.2)
mtext("C", side = 3, adj = 0, cex=1.25, line=0.2)
legend(-106, 28, 
       legend = c("GBIF occurences", "Common garden sites", "Source populations"),
       pch = c(23, 3, 21),
       pt.cex = c(0.55, 1, 1),
       col = c("grey50", "black", "black"),
       pt.bg = c("grey", "black", "red"),
       cex = 0.7, 
       bty = "n", 
       horiz = F)

# Fourth plot (D)
par(mar=c(5,4,3.75,1))  
plot(distance_summary_ordered$longitude[distance_summary_ordered$Species=="AGHY"], 
     distance_summary_ordered$geo_distance[distance_summary_ordered$Species=="AGHY"], 
     type = "l", lty = 1, xlab="Longitude", ylab="Distance from geographic center (Km)", 
     cex.lab=1.2, col="#009E73", cex.axis=0.8, ylim=c(500,1500))
points(distance_summary_ordered$longitude[distance_summary_ordered$Species=="AGHY"], 
       distance_summary_ordered$geo_distance[distance_summary_ordered$Species=="AGHY"], 
       pch = 16, cex = 2, col="#009E73")
lines(distance_summary_ordered$longitude[distance_summary_ordered$Species=="ELVI"], 
      distance_summary_ordered$geo_distance[distance_summary_ordered$Species=="ELVI"], 
      col="#D55E00")
points(distance_summary_ordered$longitude[distance_summary_ordered$Species=="ELVI"], 
       distance_summary_ordered$geo_distance[distance_summary_ordered$Species=="ELVI"], 
       pch = 16, cex = 2, col="#D55E00")
lines(distance_summary_ordered$longitude[distance_summary_ordered$Species=="POAU"], 
      distance_summary_ordered$geo_distance[distance_summary_ordered$Species=="POAU"], 
      col="#0072B2")
points(jitter(distance_summary_ordered$longitude[distance_summary_ordered$Species=="POAU"],amount = 0.1), 
       distance_summary_ordered$geo_distance[distance_summary_ordered$Species=="POAU"], 
       pch = 16, cex = 2, col="#0072B2")

legend("topright",
       legend = c("AGHY", "ELVI", "POAU"),
       col = c("#009E73", "#D55E00", "#0072B2"),
       lwd = 2,            # Line width for the curves
       lty = 1,            # Line type for the curves (solid)
       pch = 16,           # Point symbol (circle)
       pt.cex = 2,         # Size of the points in the legend
       cex = 1)

mtext("D", side = 3, adj = 0, cex = 1.25)

# Inset plot inside Panel D (move to bottom left)
# usr <- par("usr")
# x_inset_min <- usr[1] + 0.05 * (usr[2] - usr[1])
# x_inset_max <- usr[1] + 0.45 * (usr[2] - usr[1])
# y_inset_min <- usr[3] + 0.05 * (usr[4] - usr[3])
# y_inset_max <- usr[3] + 0.57 * (usr[4] - usr[3])
# 
# par(xpd = NA)
# par(fig = c(grconvertX(c(x_inset_min, x_inset_max), from="user", to="ndc"),
#             grconvertY(c(y_inset_min, y_inset_max), from="user", to="ndc")),
#     new = TRUE, mar = c(1, 1, 0.5, 0.5))
# 
# par(mgp = c(1.3, 0.5, 0))
# plot(distance_summary_ordered$longitude[distance_summary_ordered$Species=="AGHY"], 
#      distance_summary_ordered$distance[distance_summary_ordered$Species=="AGHY"], 
#      type = "l", col = "#000000", axes = TRUE, xlab = "", ylab = "DNC", ylim=c(0,40), cex.axis=0.5)
# points(distance_summary_ordered$longitude[distance_summary_ordered$Species=="AGHY"], 
#        distance_summary_ordered$distance[distance_summary_ordered$Species=="AGHY"], 
#        pch=16, col = "#000000")
# lines(distance_summary_ordered$longitude[distance_summary_ordered$Species=="ELVI"], 
#       distance_summary_ordered$distance[distance_summary_ordered$Species=="ELVI"], 
#       col = "#E69F00")
# points(distance_summary_ordered$longitude[distance_summary_ordered$Species=="ELVI"], 
#        distance_summary_ordered$distance[distance_summary_ordered$Species=="ELVI"], 
#        col = "#E69F00", pch=16)
# lines(distance_summary_ordered$longitude[distance_summary_ordered$Species=="POAU"], 
#       distance_summary_ordered$distance[distance_summary_ordered$Species=="POAU"], 
#       col = "#56B4E9")
# points(distance_summary_ordered$longitude[distance_summary_ordered$Species=="POAU"], 
#        distance_summary_ordered$distance[distance_summary_ordered$Species=="POAU"], 
#        col = "#56B4E9", pch=16)

par(op)
dev.off()





