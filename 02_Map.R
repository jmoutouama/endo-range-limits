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
prism_summary <- readRDS(url("https://www.dropbox.com/scl/fi/rjwgk98idmdatk025g6st/prism_means.rds?rlkey=dfooohwn3j2d5pew0ryjpultl&dl=1"))

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
dat25_herb <- left_join(x = dat25, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
demography_climate<-readRDS("/Users/jacobmoutouama/Desktop/Range/demography_climate.RDS")

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
herb_levels <- c(0, 1)  # 0 = Unfenced, 1 = Fenced
colors <- c("lightgreen", "salmon")
names(colors) <- herb_levels

# Create matrix for barplot (rows = Herbivory, columns = Site)
site_ids <- unique(summary_23$Site)
mean_matrix <- sapply(site_ids, function(s){
  sapply(herb_levels, function(h){
    v <- summary_23$Mean[summary_23$Site == s & summary_23$Herbivory == h]
    if (length(v) == 0) NA else v
  })
})

# Create SE matrix
se_matrix <- sapply(site_ids, function(s){
  sapply(herb_levels, function(h){
    v <- summary_23$SE[summary_23$Site == s & summary_23$Herbivory == h]
    if (length(v) == 0) NA else v
  })
})



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
mtext("P(mm)", side=3, adj=1.17, cex=0.6, line=-1.2)
### Panel B
plot(crop_ppt_annual, xlab="Longitude", ylab="", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=TRUE)
plot(elvi, add=TRUE, pch=23, col="grey50", bg="grey", cex=0.55)
plot(garden_elvi, add=TRUE, pch=3, col="black", cex=2)
plot(source_elvi, add=TRUE, pch=21, col="black", bg="red", cex=1)
mtext(~italic("Elymus virginicus"), side=3, adj=0.5, cex=1.2, line=0.2)
mtext("B", side=3, adj=0, cex=1.25, line=0.2)
mtext("P(mm)", side=3, adj=1.17, cex=0.6, line=-1.2)
### Panel C
par(mar=c(0,3,3.75,1))
plot(crop_ppt_annual, xlab="Longitude", ylab="Latitude", col=col_precip_rev, cex.lab=1.2)
plot(study_area, add=TRUE)
plot(poau, add=TRUE, pch=23, col="grey50", bg="grey", cex=0.55)
plot(garden_poau, add=TRUE, pch=3, col="black", cex=2)
plot(source_poau, add=TRUE, pch=21, col="black", bg="red", cex=1)
mtext(~italic("Poa autumnalis"), side=3, adj=0.5, cex=1.2, line=0.7)
mtext("C", side=3, adj=0, cex=1.25, line=0.3)
mtext("P(mm)", side=3, adj=1.17, cex=0.6, line=-1.2)
### Panel D (barplot ordered largest → smallest)
par(mar=c(6,4,4,1))
ordered_data <- prism_summary[order(prism_summary[,2], decreasing=TRUE), ]
bar_vals <- as.numeric(ordered_data[,2])
names(bar_vals) <- ordered_data[,1]
barplot(bar_vals,
        col="#E69F00",
        xlab="Sites", ylab="Cumulative Precipitation (mm)",
        ylim=c(0,3500),
        las=2,
        cex.lab=1.2,
        cex.names=0.75)
mtext("D", side=3, adj=-0.06, cex=1.25, line=0.5)

### Panel E
par(mar=c(0, 0, 0, 0))  
plot(0, 0, type="n", xlim=c(0,3.6), ylim=c(0,1.5), axes=FALSE, xlab="", ylab="", main="", asp=1)
mtext("E", side=3, adj=0.07, cex=1.25, line=-1.75)

# fenced
rect(0.1, 0.1, 1.6, 1.6, border="black", lwd=2)
plants1_x <- rep(seq(0.375, 1.125, length.out=4), times=4)[-1]
plants1_y <- rep(seq(0.375, 1.125, length.out=4), each=4)[-1]
points(plants1_x+0.1, plants1_y+0.1, pch=19, col="salmon", cex=1.5)
mtext("Fenced", side=3, at=0.85, line=-4, cex=0.9)

# unfenced
rect(1.9, 0.1, 3.4, 1.6, border=NA)
segments(1.9,0.1,3.4,0.1, col="black", lwd=2)
segments(1.9,1.6,3.4,1.6, col="black", lwd=2)
plants2_x <- rep(seq(2.375, 3.125, length.out=4), times=4)[-1]
plants2_y <- rep(seq(0.375, 1.125, length.out=4), each=4)[-1]
points(plants2_x-2+1.9, plants2_y+0.1, pch=19, col="lightgreen", cex=1.5)
mtext("Unfenced", side=3, at=2.65, line=-4, cex=0.9)
# Match site order in Panel F to that in Panel D
ordered_sites <- ordered_data[,1]
# Reorder the data used for Panel F
match_idx <- match(ordered_sites, site_ids)
mean_matrix <- mean_matrix[, match_idx, drop=FALSE]
se_matrix <- se_matrix[, match_idx, drop=FALSE]
site_ids <- ordered_sites

### Panel F
par(mar = c(8, 4, 2, 1))  # bottom, left, top, right
mtext("F", side = 3, adj = 1.1, cex = 1.25, line = 0.5)

bp <- barplot(
  height = mean_matrix,
  beside = TRUE,
  names.arg = site_ids,
  col = colors,
  ylim = c(0, max(mean_matrix + se_matrix, na.rm = TRUE) * 1.25),
  xlab="Sites",
  ylab = "Proportion of damaged plants",
  las = 2,
  cex.lab=1.2,
  cex.names=0.75
)

# Add error bars
arrows(
  x0 = bp,
  y0 = mean_matrix - se_matrix,
  x1 = bp,
  y1 = mean_matrix + se_matrix,
  angle = 90, code = 3, length = 0.07
)

# Add legend — shifted inside the plot for cleaner layout
legend(0,0.98,"topleft", legend = c("Unfenced", "Fenced"), fill = colors, bty = "n", cex = 1.2)

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
    width = 4, height = 5)

bp <- barplot(
  height = mean_matrix,
  beside = TRUE,
  names.arg = site_ids,
  col = colors,
  ylim = c(0, max(mean_matrix + se_matrix, na.rm=TRUE) * 1.2),
  ylab = "Proportion of plants with any damage",
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
legend(0,0.7, legend = c("Unfenced", "Fenced"), fill = colors, bty = "n", cex = 0.8)

dev.off()




