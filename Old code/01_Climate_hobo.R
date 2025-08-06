# Project: 
# Purpose: 
# Authors: Jacob Moutouama
# Date last modified (Y-M-D): 2024-08-03
rm(list = ls())
# load packages
set.seed(13)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(rmutil)
library(actuar)
#library(SPEI)
library(LaplacesDemon)
library(ggpubr)
library(raster)
library(rgdal)
library(readxl)
library(ggsci)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("scater")
# library(scater)
library(BiocManager)
library(swfscMisc)
jacob_path<-"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density"
# tom_path<-"C:/Users/tm9/Dropbox/github/ELVI-endophyte-density" 
choose_path<-jacob_path
## format date and separate year-month-day
list.files(path = paste0(choose_path,"/Data/HOBO data/"),  
           pattern = "*.xlsx", full.names = TRUE) %>% # Identify all excel files
  lapply(read_excel) %>%                              # Store all files in list
  bind_rows ->hobo_data_raw # get HOBO data

tidyr::separate(hobo_data_raw, "date",
                into = c('longdate', 'time'),
                sep= ' ') %>%
  tidyr::separate('longdate', # Separate the ‘longdate’ column into separate columns for month, day and year using the separate() function.
                  into = c('year','month', 'day'),
                  sep= '-',
                  remove = FALSE)->hobo_data_full 

## double check if the starting and ending dates are correct
hobo_data_full %>% 
  group_by(site) %>% 
  summarise(start=range(longdate)[1],
            end=range(longdate)[2],
            duration=as.Date(end)-as.Date(start))->hobo_dates 

## average over days to look at overall trend across sites
hobo_data_full %>% 
  group_by(longdate,site,day) %>% 
  summarise(daily_mean_moist=mean(water),daily_mean_temp=mean(temperature))->HOBO_daily

## Plot the daily trend for temperature and soil moisture from start to end
hobo_means<-HOBO_daily %>% 
  filter(longdate>as.Date("2023-05-01") & longdate<as.Date("2024-06-25")) %>% 
  group_by(site) %>% 
  summarise(mean_temp=mean(daily_mean_temp),
            mean_moisture=mean(daily_mean_moist))

hobo_means<-as.data.frame(hobo_means)

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Climate_hobo.pdf",width=14,height=5,useDingbats = F)
par(mar=c(5,5,2,3),mfrow=c(1,2))
barplot(hobo_means[order(hobo_means[,2],decreasing=FALSE),][,2],names.arg=hobo_means[order(hobo_means[,2],decreasing=FALSE),][,1],col="#E69F00",xlab="Sites", ylab="Mean",main="",ylim=c(0,30))
mtext("Soil temperature",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "A",side = 3, adj = 0,cex=1.2)
barplot(hobo_means[order(hobo_means[,3],decreasing=FALSE),][,3],names.arg=hobo_means[order(hobo_means[,3],decreasing=FALSE),][,1],col="#E69F00",xlab="Sites", ylab="Mean",main="",ylim=c(0,0.3))
mtext("Soil moisture",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "B",side = 3, adj = 0,cex=1.2)
dev.off()

data_plotclim<-data.frame(site=c(HOBO_daily$site,HOBO_daily$site),daily_mean_clim=c(HOBO_daily$daily_mean_temp,HOBO_daily$daily_mean_moist),date=c(HOBO_daily$longdate,HOBO_daily$longdate),clim=c(rep("temp",nrow(HOBO_daily)),rep("water",nrow(HOBO_daily))))

site_names <- c("LAF"="Lafayette",
                 "HUN"="Huntville",
                 "BAS"="Bastrop",
                "COL"="College Station",
                "KER" ="Kerville",
                "BFL" ="Brackenridge",
                "SON"="Sonora")

figtempsite<-ggplot(HOBO_daily, aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_temp))+
  geom_line(aes(colour=site))+
  ggtitle("a")+
  scale_fill_jco()+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.5,color="black", angle=0))+
  labs( y="Daily temperature  (°C)", x="")+
  facet_grid(~factor(site,levels=c("LAF","HUN","COL","BAS","BFL","KER","SON")))+
  geom_hline(data=hobo_means,aes(yintercept = mean_temp,colour=site))

figmoistsite<-ggplot(HOBO_daily, aes(x=as.Date(longdate, format= "%Y - %m - %d"), y=daily_mean_moist))+
  geom_line(aes(colour=site))+
  ggtitle("b")+
  scale_fill_jco()+
  theme_bw()+ 
  theme(legend.position = "none",
        axis.text.x = element_text(size=4.55, color="black",angle=0))+
  labs( y="Daily soil moisture (wfv)", x="Month")+
  facet_grid(~factor(site,levels=c("LAF","HUN","COL","BAS","BFL","KER","SON")))+
  geom_hline(data=hobo_means,aes(yintercept = mean_moisture,colour=site))

# pdf("paste0(choose_path/Figure/climatesite.pdf",height =5,width =12,useDingbats = F)
# (Figclimatesite<-ggpubr::ggarrange(figtempsite,figmoistsite,common.legend = FALSE,ncol = 1, nrow = 2))
# dev.off()