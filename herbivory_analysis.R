## This script reads in range limits census data and analyzes herbivory

# remove all objects and clear workspace
rm(list = ls(all=TRUE))

##tom's local wd
setwd("C:/Users/tm9/Dropbox/github/endo-range-limits-prevalence")

##load packages
library(tidyverse)
library(googlesheets4)
library(rstan)
library(googledrive)
library(readxl)
library(bayesplot)
library(tidybayes)
library(posterior)
library(magrittr)

##read in endo range limits database from google drive
##read in data as text so that nothing gets coerced
file_id<-"11m_-0KYtLtp5J0EqGWZ8RzNM-JfsDB3B"
drive_download(as_id(file_id), path = "Data/endo_range_limits_experiment.xlsx", overwrite = TRUE)
##plot-level treatment assignments
plots<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="Initial Plot Data", col_types = "text")
##initial plant data (including source and endo status)
init23<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="Initial Plant Data", col_types = "text")
##2024 replant
init24<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2024 Re-plant inital data", col_types = "text")
##2023-25 census
census23<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2023 Data Census", col_types = "text")
census24<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2024 Data Census", col_types = "text")
census25<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2025 Data Census", col_types = "text")


# data prep ---------------------------------------------------------------


##there is one unfortunated underscored tag id: 731_2, which was the second use of 731
##give it a new unique tag not used by anything else
##it only appears once in the initial data, and once in the 23 census
init23[init23$Tag_ID=="731_2","Tag_ID"]<-"10001"
census23[census23$Tag_ID=="731_2","Tag_ID"]<-"10001"

##drop all the initial plant tags starting with D (these were dead on first census)
init23 %<>% filter(substr(Tag_ID,1,1)!="D") %>% 
  ##coerce remaining tags, plot, and endo to integer
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Plot = as.integer(as.numeric(Plot)),
         Endo = as.integer(as.numeric(Endo)))
##warning OK: there is 1 NA tag ID -- dead ELVI in LAF

plots %<>% mutate(Plot = as.integer(as.numeric(Plot)),
                  Herbivory = case_when(Herbivory == "0.0" ~ "Access",
                                        Herbivory == "1.0" ~ "Exclusion"))

##merge initial plant and plot treatment data
init_plot <- left_join(init23,plots %>% select(Site,Plot,Species,Herbivory),
                       by=c("Site","Species","Plot"))
##2024 is more complicated because of replant
init24 %<>% filter(substr(Tag_ID,1,1)!="D") %>% 
  ##coerce remaining tags, plot, and endo to integer
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Plot = as.integer(as.numeric(Plot)),
         Endo = as.integer(as.numeric(Endo)),
         Herbivory = case_when(Herbivory_treatment == "0.0" ~ "Access",
                               Herbivory_treatment == "1.0" ~ "Exclusion"))

init_combo<-bind_rows(init_plot,init24 %>% select(names(init_plot))) %>% 
  filter(!is.na(Tag_ID))

census23 %<>% 
  ##drop dead plants
  filter(substr(Tag_ID,1,1)!="D") %>% 
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Tiller = as.integer(as.numeric(Tiller_23)),
         Inflor = as.integer(as.numeric(Inf_23)),
         TillerHerb = as.integer(as.numeric(tiller_Herb)))
#warning OK: one row of NA data, not sure why
##are all the 2023 census tags in the initial plant data? -- YES
census23[which(!(census23$Tag_ID %in% init_plot$Tag_ID)),]

##merge 2023 censusu and initial plant-plot data
dat23<-left_join(census23,init_plot,by="Tag_ID")

census24 %<>% 
  ##drop dead plants
  filter(substr(Tag_ID,1,1)!="D") %>% 
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Tiller = as.integer(as.numeric(Tiller_24)),
         AttachedInf = as.integer(as.numeric(attachedInf_24)),
         BrokenInf = as.integer(as.numeric(brokenInf_24)),
         TillerHerb = as.integer(as.numeric(tiller_Herb)))
  #warning OK
##are all the 2024 census tags in either init23 or init24?
sum(census24$Tag_ID %in% init23$Tag_ID)+sum(census24$Tag_ID %in% init24$Tag_ID)
length(census24$Tag_ID)
##it looks like there were 3 tagIDs that matched into both init23 and init24:
census24[(census24$Tag_ID %in% init23$Tag_ID) & (census24$Tag_ID %in% init24$Tag_ID),]
##tags 2097, 2098, and 2099 -- these are actually the same data in both places, 
##because they are the sole 3 survivors at SON

dat24<-left_join(census24,init_combo,by="Tag_ID",multiple="first")
##multiple=first deals with two init_combo rows for tag 2097

census25 %<>% 
  ##drop dead plants
  filter(substr(Tag_ID,1,1)!="D") %>% 
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Tiller = as.integer(as.numeric(Tiller_25)),
         AttachedInf = as.integer(as.numeric(attachedInf_25)),
         BrokenInf = as.integer(as.numeric(brokenInf_25)),
         TillerHerb = as.integer(as.numeric(tiller_Herb)))

dat25<-left_join(census25,init_combo,by=c("Site","Species","Tag_ID"))

##combine years
bind_rows(dat23 %>% 
            mutate(Year=2023,
                   AttachedInf=NA,
                   BrokenInf=NA) %>% 
            select(Year,Site,Species,Plot,Herbivory,Tag_ID,Population,Endo,Tiller,AttachedInf,BrokenInf,TillerHerb),
          dat24 %>% 
            mutate(Year=2024) %>% 
            select(Year,Site,Species,Plot,Herbivory,Tag_ID,Population,Endo,Tiller,AttachedInf,BrokenInf,TillerHerb),
          dat25 %>% 
            mutate(Year=2025) %>% 
            select(Year,Site,Species,Plot,Herbivory,Tag_ID,Population,Endo,Tiller,AttachedInf,BrokenInf,TillerHerb)
)->herb_dat_tidy

#write.csv(herb_dat_tidy,"Data/herbivory_tidy_data.csv")


# analysis ----------------------------------------------------------------


