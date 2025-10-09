# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis on  vital rate models (survival, growth, flowering and spikelet).
# Authors: Jacob Moutouama
# Date last modified (Y-M-D):
rm(list = ls())
# load packages
# remove.packages(c("StanHeaders", "rstan"))
# install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(rstan)
# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
# install.packages("countreg",repos = "http://R-Forge.R-project.org")
# library(countreg)
library(rmutil)
library(actuar)
# library(SPEI)
library(LaplacesDemon)
library(ggpubr)
library(raster)
# library(rgdal)
library(readxl)
library(ggsci)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("scater")
# library(scater)
library(BiocManager)
library(swfscMisc)
library(bayesplot)
library(extraDistr)
library(lubridate)
library(patchwork)
# Install cmdstanr from R
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# # Install CmdStan backend (once)
# cmdstanr::install_cmdstan()
#library(cmdstanr)
# Define some basic functions that we'll use later
quote_bare <- function(...) {
  substitute(alist(...)) %>%
    eval() %>%
    sapply(deparse)
}
set.seed(13)
# Demographic data -----
# Merge the demographic census
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1", stringsAsFactors = F)
dat23 <- read.csv("https://www.dropbox.com/scl/fi/fkwm0dan6nx2eaeyxjrjw/census2023.csv?rlkey=hy9209t53j9n7vxhta7axl5jk&dl=1", stringsAsFactors = F)
dat24 <- read.csv("https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1", stringsAsFactors = F)
dat25<-read.csv("https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1", stringsAsFactors = F)
datherbivory <- read.csv("https://www.dropbox.com/scl/fi/2gnlfozxpd2u9gprzp9oi/herbivory.csv?rlkey=sz2cloqxtbc6ou29j97l3t10f&dl=1", stringsAsFactors = F)

# Climatic data ----
climate_site <- readRDS(url("https://www.dropbox.com/scl/fi/arl44h1v0xoaymz0s4b5m/site_climate_summary.rds?rlkey=w7wuh62on061cm3cdhwl6gu9e&dl=1"))
# head(climate_site)

# calculate the average spikelet and inflorescence number for each census
dat23 %>%
  mutate(spikelet_23 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T)), digit = 0) -> dat23_spike
dat24 %>%
  mutate(spikelet_24 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_24 = round(rowSums(across(attachedInf_24:brokenInf_24), na.rm = T), digit = 0)) -> dat24_spike
dat25 %>%
  mutate(spikelet_25 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_25 = round(rowSums(across(attachedInf_25:brokenInf_25), na.rm = T), digit = 0)) -> dat25_spike

# Check for duplicate Tag_IDs
datini$Tag_ID <- as.character(datini$Tag_ID)
dat23_spike$Tag_ID <- as.character(dat23_spike$Tag_ID)
dat23_spike %>% count(Tag_ID) %>% filter(n > 1)
dat24_spike$Tag_ID <- as.character(dat24_spike$Tag_ID)
dat25_spike$Tag_ID <- as.character(dat25_spike$Tag_ID)
#view(dat24_spike)
#dat24_spike %>% count(Tag_ID) %>% filter(n > 1)
#dat25_spike %>% count(Tag_ID) %>% filter(n > 1)

## Merge the initial data to the 23 to get the Tag ID for each elements -----
datini23_spike <- dat23_spike %>%
  left_join(
    datini %>%
      dplyr::select(Site, Species, Plot, Position, Tag_ID, Population, 
                    GreenhouseID, Clone, Endo),
    by = "Tag_ID"
  ) %>%
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit")))

# view(datini23_spike)

#Combine datini and dat23
combined_data <- bind_rows(datini[,c("Site","Species","Plot","Tag_ID","Population","Endo")], datini23_spike[,c("Site","Species","Plot","Tag_ID","Population","Endo" )])%>% 
  distinct(Tag_ID, .keep_all = TRUE)

dat24_spike_sp_site_tag<-left_join(x = dat24_spike, y = combined_data, by = c("Tag_ID")) %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_24","brokenInf_24"))) %>% 
  filter(!is.na(Species))
#dat24_spike_sp_site_tag %>% count(Tag_ID) %>% filter(n > 1) # No duplicate 
# view(dat24_spike_sp_site_tag)

dat25_spike_sp_site_tag<- dat25_spike %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_25","brokenInf_25"))) 
#dat25_spike_sp_site_tag %>% count(Tag_ID) %>% filter(n > 1) # No duplicate 

dat2324 <- datini23_spike %>%
  left_join(
    dat24_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Inf_24, Tiller_24, tiller_herb_24, date_24, stroma_24, spikelet_24),
    by = "Tag_ID"
  ) 

#dat2324 %>% count(Tag_ID) %>% filter(n > 1) # No duplicate 
# dat2324 %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)

dat2425 <- dat24_spike_sp_site_tag %>%
  left_join(
    dat25_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Inf_25, Tiller_25, tiller_herb_25, date_25, stroma_25, spikelet_25),
    by = "Tag_ID"
  )

# Change variable names
dat2324%>%
  mutate(
    tiller_t = Tiller_23,
    tiller_t1 = Tiller_24,
    inf_t = Inf_23,
    inf_t1 = Inf_24,
    spikelet_t = spikelet_23,
    spikelet_t1 = spikelet_24,
    tiller_Herb_t = tiller_Herb_23,
    tiller_Herb_t1 = tiller_herb_24,
    date_t=date_23,
    date_t1=date_24,
    census_year=rep(2024,nrow(datini23_spike))
  ) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Tag_ID,
    Population,
    Endo,
    tiller_t,
    tiller_t1,
    inf_t,
    inf_t1,
    spikelet_t,
    spikelet_t1,
    tiller_Herb_t,
    tiller_Herb_t1,
    census_year
  ) -> dat2324_t_t1

# dat2324_t_t1 %>% count(Tag_ID) %>% filter(n > 1) 
# dat2324_t_t1 %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)
dat2425 %>%
  mutate(
    tiller_t = Tiller_24,
    tiller_t1 = Tiller_25,
    inf_t = Inf_24,
    inf_t1 = Inf_25,
    spikelet_t = spikelet_24,
    spikelet_t1 = spikelet_25,
    tiller_Herb_t = tiller_herb_24,
    tiller_Herb_t1 = tiller_herb_25,
    date_t = date_24,
    date_t1 = date_25,
    census_year=rep(2025,nrow(dat24_spike_sp_site_tag))
  ) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Tag_ID,
    Population,
    Endo,
    tiller_t,
    tiller_t1,
    inf_t,
    inf_t1,
    spikelet_t,
    spikelet_t1,
    tiller_Herb_t,
    tiller_Herb_t1,
    census_year
  ) -> dat2425_t_t1

# dat2425_t_t1 %>% count(Tag_ID) %>% filter(n > 1)
# dat2425_t_t1 %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)
dat_t_t1 <- rbind(dat2324_t_t1, dat2425_t_t1)
# Find duplicates by Tag ID within each census year
# dup_tags_per_year <- dat_t_t1 %>%
#   count(Tag_ID, census_year) %>%
#   filter(n > 1)
# dup_tags_per_year

## Merge the demographic data with the herbivory data -----
dat_t_t1_herb <- left_join(x = dat_t_t1, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
# head(dat_t_t1_herb)
# unique(dat_t_t1_herb$Species)
# view(dat_t_t1_herb)

## Merge the demographic data with the climatic data -----
climate_site_unique <- climate_site %>%
  distinct(site, Species, census_year, .keep_all = TRUE)

#view(climate_site_unique)
dat_t_t1_herb$census_year<-as.character(dat_t_t1_herb$census_year)
climate_site_unique$census_year<-as.character(climate_site_unique$census_year)

dat_t_t1_herb_clim <- dat_t_t1_herb %>%
  left_join(
    climate_site_unique,
    by = c("Site" = "site", "Species", "census_year")
  )

# dat_t_t1_herb_clim %>%
# group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)

dat_t_t1_herb_clim %>%
  mutate(
    site_year=interaction(Site,census_year),
    ##if the plant was dead or size was NA at the start of the transition year, survival is NA
    ##if the plant was alive at the start of the transition year, it survived if tillers_t1>0
    surv1 = ifelse(tiller_t>0,tiller_t1>0,NA),
    site_species_plot = interaction(Site, Species, Plot),
    ##if the plant was alive at the start of the transition year and it survived, growth is the log ratio of tiller counts, else NA
    ##note there that growth is conditional on survival, which I think is how it should be
    grow = ifelse(tiller_t>0 & tiller_t1>0,log(tiller_t1/tiller_t),NA_real_)
  ) -> demography_climate

## check for TagIDs  duplicated within years
# demography_climate %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)

## Explore the surival rate per species
demography_climate %>%
  group_by(Species) %>%
  summarise(
    n_total = n(),
    n_survived = sum(surv1, na.rm = TRUE),
    survival_rate = n_survived / n_total
  )
# names(demography_climate)
# view(demography_climate)
# summary(demography_climate)

## Running the stan model
sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  thin = 2,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4
)

# Survival----
## Read and format survival data to build the model
demography_climate_surv <- demography_climate %>%
  filter(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, site_species_plot, site_year, Endo, Herbivory,
    tiller_t, surv1, cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    # Grouping indices for random effects
    Site              = as.integer(factor(Site)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_year         = as.integer(factor(site_year)),
    site_species_plot = as.integer(factor(site_species_plot)),
    
    # Survival and predictors
    log_size_t0 = log(tiller_t),
    surv_t1     = as.integer(surv1),
    ppt         = log(cum_ppt)
  )

# summary(demography_climate_surv$cum_ppt)
# any(demography_climate_surv$cum_ppt <= 0)
# Convert into list for Stan
demography_surv_ppt <- list(
  nSpp       = n_distinct(demography_climate_surv$Species),
  nSite      = n_distinct(demography_climate_surv$Site),
  nsite_year = n_distinct(demography_climate_surv$site_year),
  nPop       = n_distinct(demography_climate_surv$Population),
  nPlot      = n_distinct(demography_climate_surv$site_species_plot),
  Spp        = demography_climate_surv$Species,
  site       = demography_climate_surv$Site,
  site_year  = demography_climate_surv$site_year,
  pop        = demography_climate_surv$Population,
  plot       = demography_climate_surv$site_species_plot,
  clim       = as.vector(demography_climate_surv$ppt),
  endo       = demography_climate_surv$Endo,      # already 0/1
  herb       = demography_climate_surv$Herbivory, # already 0/1
  size       = demography_climate_surv$log_size_t0,
  y          = demography_climate_surv$surv_t1,
  N          = nrow(demography_climate_surv)
)


fit_surv_abio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/surv_abio_endo.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)

fit_sur_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/surv_bio_endo.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)


fit_surv_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/survival.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)

# rstan::check_hmc_diagnostics(fit_surv_ppt)
#summary(fit_surv_ppt)$summary[, c("Rhat", "n_eff")]
posterior_surv_abio_bio_endo <- as.array(fit_surv_abio_bio_endo) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_abio_bio_endo,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()

## Save RDS file for further use
# saveRDS(fit_surv_abio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_abio_endo.rds')
# saveRDS(fit_sur_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_sur_bio_endo.rds')
# saveRDS(fit_surv_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_abio_bio_endo.rds')

# Growth----
## Read and format survival data to build the model
demography_climate%>%
  filter(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot,site_year, site_species_plot, Endo, Herbivory,
    tiller_t, grow,cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    site_year = as.integer(factor(site_year)),
    Endo = as.integer(Endo),
    Herbivory = as.integer(Herbivory)
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    grow = grow,
    ppt = log(cum_ppt),
  ) -> demography_climate_grow

## Separate each variable to use the same model stan
### Precipitation
demography_grow_ppt <- list(
  nSpp = demography_climate_grow$Species %>% n_distinct(),
  nSite = demography_climate_grow$Site %>% n_distinct(),
  nsite_year = demography_climate_grow$site_year %>% n_distinct(),
  nPop = demography_climate_grow$Population %>% n_distinct(),
  nPlot = demography_climate_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_grow$Species,
  site = demography_climate_grow$Site,
  site_year = demography_climate_grow$site_year,
  pop = demography_climate_grow$Population,
  plot = demography_climate_grow$site_species_plot,
  clim = as.vector(demography_climate_grow$ppt),
  endo = demography_climate_grow$Endo,
  herb = demography_climate_grow$Herbivory,
  size = demography_climate_grow$log_size_t0,
  y = demography_climate_grow$grow,
  N = nrow(demography_climate_grow)
)


fit_grow_abio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/grow_abio_endo.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

fit_grow_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/grow_bio_endo.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)


fit_grow_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/growth.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

# summary(fit_grow_ppt)$summary[, c("Rhat", "n_eff")]
posterior_grow_abio_bio_endo <- as.array(fit_grow_abio_bio_endo) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_abio_bio_endo,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()



## Save RDS file for further use
# saveRDS(fit_grow_abio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_abio_endo.rds')
# saveRDS(fit_grow_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_bio_endo.rds')
# saveRDS(fit_grow_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_abio_bio_endo.rds')

# Flowering----
demography_climate %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    site_year=as.integer(factor(site_year)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(Endo),
    Herbivory = as.integer(Herbivory)
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    flow_t1 = inf_t1,
    ppt = log(cum_ppt)
  ) -> demography_climate_flow

## Separate each variable to use the same model stan
### Precipitation 
demography_flow_ppt <- list(
  nSpp = demography_climate_flow$Species %>% n_distinct(),
  nSite = demography_climate_flow$Site %>% n_distinct(),
  nsite_year = demography_climate_flow$site_year %>% n_distinct(),
  nPop = demography_climate_flow$Population %>% n_distinct(),
  nPlot = demography_climate_flow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_flow$Species,
  site = demography_climate_flow$Site,
  site_year = demography_climate_flow$site_year,
  pop = demography_climate_flow$Population,
  plot = demography_climate_flow$site_species_plot,
  clim = as.vector(demography_climate_flow$ppt),
  endo = demography_climate_flow$Endo,
  herb = demography_climate_flow$Herbivory,
  size = demography_climate_flow$log_size_t0,
  y = demography_climate_flow$flow_t1,
  N = nrow(demography_climate_flow)
)

fit_flow_ppt_intercept <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering_ppt_intercept.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

fit_flow_ppt_abiotic <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering_ppt_abiotic.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

fit_flow_ppt_biotic <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/fowering_ppt_biotic.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)


fit_flow_ppt <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

#summary(fit_flow_ppt)$summary[, c("Rhat", "n_eff")]
posterior_flow_ppt <- as.array(fit_flow_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_ppt,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3],
                        bendoclim2[1], bendoclim2[2], bendoclim2[3]
                      )
) + theme_bw()



## Save RDS file for further use
# saveRDS(fit_flow_ppt_intercept, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_ppt_intercept.rds')
# saveRDS(fit_flow_ppt_abiotic, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_ppt_abiotic.rds')
# saveRDS(fit_flow_ppt_biotic, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_ppt_biotic.rds')
# saveRDS(fit_flow_ppt, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_ppt.rds')

# Spikelet----
demography_climate %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    site_year=as.integer(factor(site_year)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(Endo),
    Herbivory = as.integer(Herbivory)
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    spi_t1 = spikelet_t1,
    ppt = log(cum_ppt)
  ) -> demography_climate_spik

### Precipitation
demography_spik_ppt <- list(
  nSpp = demography_climate_spik$Species %>% n_distinct(),
  nSite = demography_climate_spik$Site %>% n_distinct(),
  nsite_year =demography_climate_spik$site_year %>% n_distinct(),
  nPop = demography_climate_spik$Population %>% n_distinct(),
  nPlot = demography_climate_spik$site_species_plot %>% n_distinct(),
  Spp = demography_climate_spik$Species,
  site = demography_climate_spik$Site,
  site_year = demography_climate_spik$site_year,
  pop = demography_climate_spik$Population,
  plot = demography_climate_spik$site_species_plot,
  clim = as.vector(demography_climate_spik$ppt),
  endo = demography_climate_spik$Endo,
  herb = demography_climate_spik$Herbivory,
  size = demography_climate_spik$log_size_t0,
  y = demography_climate_spik$spikelet_t1,
  N = nrow(demography_climate_spik)
)

fit_spik_ppt_intercept <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spik_ppt_intercept.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control =sim_pars$control)

fit_spik_ppt_abiotic <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spik_ppt_abiotic.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control =sim_pars$control)

fit_spik_ppt_biotic <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spik_ppt_biotic.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control =sim_pars$control)

fit_spik_ppt <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control =sim_pars$control)

# summary(fit_spik_ppt)$summary[, c("Rhat", "n_eff")]
posterior_spik_ppt <- as.array(fit_spik_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_ppt,
                      pars = quote_bare(
                        b0[1], b0[2],
                        bendo[1], bendo[2],
                        bherb[1], bherb[2],
                        bclim[1], bclim[2],
                        bendoclim[1], bendoclim[2],
                        bendoherb[1], bendoherb[2],
                        bclim2[1], bclim2[2],
                        bendoclim2[1], bendoclim2[2]
                      )
) + theme_bw()



## Save RDS file for further use
# saveRDS(fit_spik_ppt_intercept, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_ppt_intercept.rds')
# saveRDS(fit_spik_ppt_abiotic, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_ppt_abiotic.rds')
# saveRDS(fit_spik_ppt_biotic, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_ppt_biotic.rds')
# saveRDS(fit_spik_ppt, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_ppt.rds')

# Posterior predictive check----
## Helper functions
inv_logit <- function(x) 1 / (1 + exp(-x))

# Function to simulate posterior predictive samples for each vital rate
simulate_ppc <- function(pred_matrix, phi = NULL, family = c("bernoulli", "normal", "negbinomial")) {
  family <- match.arg(family)
  n_iter <- nrow(pred_matrix)
  N <- ncol(pred_matrix)
  y_rep <- matrix(NA, nrow = n_iter, ncol = N)
  
  if (family == "bernoulli") {
    probs <- inv_logit(pred_matrix)
    for (i in 1:n_iter) {
      y_rep[i, ] <- rbinom(n = N, size = 1, prob = probs[i, ])
    }
  } else if (family == "normal") {
    # phi is sd here
    for (i in 1:n_iter) {
      y_rep[i, ] <- rnorm(N, mean = pred_matrix[i, ], sd = phi[i])
    }
  } else if (family == "negbinomial") {
    # phi is size parameter
    for (i in 1:n_iter) {
      mu_i <- exp(pred_matrix[i, ])
      size_i <- phi[i]
      y_rep[i, ] <- rnbinom(n = N, size = size_i, mu = mu_i)
    }
  }
  return(y_rep)
}

bayesplot::color_scheme_set("blue")
## Survival Model full model----
fit_survival <- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_surv_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))
post_surv <- rstan::extract(fit_survival)
pred_surv <- post_surv$predS
y_surv <- demography_surv_ppt$y
y_rep_surv <- simulate_ppc(pred_surv, family = "bernoulli")
p_surv <- ppc_dens_overlay(y_surv, y_rep_surv[1:500, ]) + ggtitle("Survival")

## Growth Model full model----
fit_growth <- readRDS(url("https://www.dropbox.com/scl/fi/5oduhrkn3l0cu5b9soju5/fit_grow_ppt.rds?rlkey=mpxxl4aejowhdm29pij9jv8kk&dl=1"))
post_grow <- rstan::extract(fit_growth)
pred_grow <- post_grow$predG
sigma_grow <- post_grow$sigma
y_grow <- demography_grow_ppt$y
y_rep_grow <- simulate_ppc(pred_grow, phi = sigma_grow, family = "normal")
p_grow <- ppc_dens_overlay(y_grow, y_rep_grow[1:500, ]) + ggtitle("Growth")

## Flowering Model full model ----
fit_flowering <- readRDS(url("https://www.dropbox.com/scl/fi/1j3ln3jxk94s56c9j193q/fit_flow_ppt.rds?rlkey=ag5bdlhngtg2gsfbfbx15purr&dl=1"))
post_flow <- rstan::extract(fit_flowering)
pred_flow <- post_flow$predF
phi_flow <- post_flow$phi
y_flow <- demography_flow_ppt$y
y_rep_flow <- simulate_ppc(pred_flow, phi = phi_flow, family = "negbinomial")
p_flow <- ppc_dens_overlay(y_flow, y_rep_flow[1:500, ]) + ggtitle("Flowering")

## Spikelet Model full model ----
fit_spikelet <- readRDS(url("https://www.dropbox.com/scl/fi/fg562lkkl077j8qoegb8q/fit_spik_ppt.rds?rlkey=cbfyy7tq7tja3e9mxujv0scm6&dl=1"))
post_spik <- rstan::extract(fit_spikelet)
pred_spik <- post_spik$predF
phi_spik <- post_spik$phi
y_spik <- demography_spik_ppt$y
y_rep_spik <- simulate_ppc(pred_spik, phi = phi_spik, family = "negbinomial")
p_spik <- ppc_dens_overlay(y_spik, y_rep_spik[1:500, ]) + ggtitle("Spikelet")

# Combine all PPC plots full model ----
combined_plot <- (p_surv | p_grow) / (p_flow | p_spik) +
  plot_annotation(title = "")
print(combined_plot)

combined_plot <- (p_surv | p_grow) / (p_flow | p_spik) +
  plot_annotation(title = "") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# Save to PDF
ggsave(filename = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/combined_ppc_plots.pdf", plot = combined_plot, width = 7, height = 6)

## Survival Model abiotic----
fit_survival_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/ty2yavh70dopswv2wbrij/fit_surv_ppt_abiotic.rds?rlkey=sqd5etosja0pw41pocx6x9pyu&dl=1"))
post_surv_abiotic <- rstan::extract(fit_survival_abiotic)
pred_surv_abiotic <- post_surv_abiotic$predS
y_surv_abiotic <- demography_surv_ppt$y
y_rep_surv_abiotic <- simulate_ppc(pred_surv_abiotic, family = "bernoulli")
p_surv_abiotic <- ppc_dens_overlay(y_surv_abiotic, y_rep_surv_abiotic[1:500, ]) + ggtitle("Survival")

## Growth Model abiotic----
fit_growth_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/wwy6y5z2xpu5ukausel6m/fit_grow_ppt_abiotic.rds?rlkey=teyk9e49bjp8ojsyi860z1241&dl=1"))
post_grow_abiotic <- rstan::extract(fit_growth_abiotic)
pred_grow_abiotic <- post_grow_abiotic$pred
sigma_grow_abiotic <- post_grow_abiotic$sigma
y_grow_abiotic <- demography_grow_ppt$y
y_rep_grow_abiotic <- simulate_ppc(pred_grow_abiotic, phi = sigma_grow_abiotic, family = "normal")
p_grow_abiotic <- ppc_dens_overlay(y_grow_abiotic, y_rep_grow_abiotic[1:500, ]) + ggtitle("Growth")

## Flowering Model abiotic ----
fit_flowering_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/8amcccm6vtn9jxsoa97yk/fit_flow_ppt_abiotic.rds?rlkey=l87f6ipeyn5si41wozsz2c57k&dl=1"))
post_flow_abiotic <- rstan::extract(fit_flowering_abiotic)
pred_flow_abiotic <- post_flow_abiotic$pred
phi_flow_abiotic <- post_flow_abiotic$phi
y_flow_abiotic <- demography_flow_ppt$y
y_rep_flow_abiotic <- simulate_ppc(pred_flow_abiotic, phi = phi_flow_abiotic, family = "negbinomial")
p_flow_abiotic <- ppc_dens_overlay(y_flow_abiotic, y_rep_flow_abiotic[1:500, ]) + ggtitle("Flowering")

## Spikelet Model abiotic----
fit_spikelet_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/idg2f35yu055pv38ugp71/fit_spik_ppt_abiotic.rds?rlkey=c860tyv2b187p8w9a8ctrl3x9&dl=1"))
post_spik_abiotic <- rstan::extract(fit_spikelet_abiotic)
pred_spik_abiotic <- post_spik_abiotic$pred
phi_spik_abiotic <- post_spik_abiotic$phi
y_spik_abiotioc <- demography_spik_ppt$y
y_rep_spik_abiotic <- simulate_ppc(pred_spik_abiotic, phi = phi_spik_abiotic, family = "negbinomial")
p_spik_abiotic <- ppc_dens_overlay(y_spik_abiotioc, y_rep_spik_abiotic[1:500, ]) + ggtitle("Spikelet")

# Combine all PPC plots abiotic ----
combined_plot_abiotic <- (p_surv_abiotic | p_grow_abiotic) / (p_flow_abiotic | p_spik_abiotic) +
  plot_annotation(title = "")+
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
print(combined_plot_abiotic)
ggsave(filename = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/combined_ppc_plots_abiotic.pdf", plot = combined_plot_abiotic, width = 7, height = 6)


## Survival Model biotic----
fit_survival_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/ty2yavh70dopswv2wbrij/fit_surv_ppt_biotic.rds?rlkey=sqd5etosja0pw41pocx6x9pyu&dl=1"))
post_surv_biotic <- rstan::extract(fit_survival_biotic)
pred_surv_biotic <- post_surv_biotic$predS
y_surv_biotic <- demography_surv_ppt$y
y_rep_surv_biotic <- simulate_ppc(pred_surv_biotic, family = "bernoulli")
p_surv_biotic <- ppc_dens_overlay(y_surv_biotic, y_rep_surv_biotic[1:500, ]) + ggtitle("Survival")

## Growth Model biotic----
fit_growth_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/wwy6y5z2xpu5ukausel6m/fit_grow_ppt_biotic.rds?rlkey=teyk9e49bjp8ojsyi860z1241&dl=1"))
post_grow_biotic <- rstan::extract(fit_growth_biotic)
pred_grow_biotic <- post_grow_biotic$pred
sigma_grow_biotic <- post_grow_biotic$sigma
y_grow_biotic <- demography_grow_ppt$y
y_rep_grow_biotic <- simulate_ppc(pred_grow_biotic, phi = sigma_grow_biotic, family = "normal")
p_grow_biotic <- ppc_dens_overlay(y_grow_biotic, y_rep_grow_biotic[1:500, ]) + ggtitle("Growth")

## Flowering Model biotic ----
fit_flowering_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/8amcccm6vtn9jxsoa97yk/fit_flow_ppt_biotic.rds?rlkey=l87f6ipeyn5si41wozsz2c57k&dl=1"))
post_flow_biotic <- rstan::extract(fit_flowering_biotic)
pred_flow_biotic <- post_flow_biotic$pred
phi_flow_biotic <- post_flow_biotic$phi
y_flow_biotic <- demography_flow_ppt$y
y_rep_flow_biotic <- simulate_ppc(pred_flow_biotic, phi = phi_flow_biotic, family = "negbinomial")
p_flow_biotic <- ppc_dens_overlay(y_flow_biotic, y_rep_flow_biotic[1:500, ]) + ggtitle("Flowering")

## Spikelet Model biotic----
fit_spikelet_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/idg2f35yu055pv38ugp71/fit_spik_ppt_biotic.rds?rlkey=c860tyv2b187p8w9a8ctrl3x9&dl=1"))
post_spik_biotic <- rstan::extract(fit_spikelet_biotic)
pred_spik_biotic <- post_spik_biotic$pred
phi_spik_biotic <- post_spik_biotic$phi
y_spik_abiotioc <- demography_spik_ppt$y
y_rep_spik_biotic <- simulate_ppc(pred_spik_biotic, phi = phi_spik_biotic, family = "negbinomial")
p_spik_biotic <- ppc_dens_overlay(y_spik_abiotioc, y_rep_spik_biotic[1:500, ]) + ggtitle("Spikelet")

# Combine all PPC plots biotic ----
combined_plot_biotic <- (p_surv_biotic | p_grow_biotic) / (p_flow_biotic | p_spik_biotic) +
  plot_annotation(title = "")+
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
print(combined_plot_biotic)
ggsave(filename = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/combined_ppc_plots_biotic.pdf", plot = combined_plot_biotic, width = 7, height = 6)


