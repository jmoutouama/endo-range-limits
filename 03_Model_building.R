# Model building to compare linear vs quadratic effects and identify drivers of grass-endophyte symbiosis outcomes
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
climate_site_scaled <- climate_site %>%
  mutate(
    # Log-transform precipitation first to reduce skew
    ppt_log = log(cumulative_pptmean),
    # Standardize across all rows: mean = 0, SD = 1
    ppt_scaled = (ppt_log - mean(ppt_log, na.rm = TRUE)) / sd(ppt_log, na.rm = TRUE)
  )
#saveRDS(climate_site_scaled,"/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/climate_site_scaled.rds")

# calculate the average spikelet and inflorescence number for each census
dat23_spike <- dat23 %>%
  mutate(Flowered_23 = ifelse(Inf_23 > 0, 1, 0), 
    spikelet_23 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = TRUE), digits = 0))
dat24_spike <- dat24 %>%
  mutate(
    spikelet_24 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = TRUE), digits = 0),
    Inf_24 = round(rowSums(across(attachedInf_24:brokenInf_24), na.rm = TRUE), digits = 0),
    Flowered_24 = ifelse(Inf_24 > 0, 1, 0)
  )
dat25_spike <- dat25 %>%
  mutate(
    spikelet_25 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = TRUE), digits = 0),
    Inf_25 = round(rowSums(across(attachedInf_25:brokenInf_25), na.rm = TRUE), digits = 0),
    Flowered_25 = ifelse(Inf_25 > 0, 1, 0)
  )
# calculate the total spikelet and inflorescence number for each census
# dat23_spike <- dat23 %>%
#   mutate(Flowered_23 = ifelse(Inf_23 > 0, 1, 0), 
#     spikelet_23 = rowSums(across(Spikelet_A:Spikelet_C), na.rm = TRUE))
# 
# dat24_spike <- dat24 %>%
#   mutate(
#     spikelet_24 = rowSums(across(Spikelet_A:Spikelet_C), na.rm = TRUE),
#     Inf_24 = rowSums(across(attachedInf_24:brokenInf_24), na.rm = TRUE),
#     Flowered_24 = ifelse(Inf_24 > 0, 1, 0)  # 1 if any inflorescence, 0 otherwise
#   )
# 
# dat25_spike <- dat25 %>%
#   mutate(
#     spikelet_25 = rowSums(across(Spikelet_A:Spikelet_C), na.rm = TRUE),
#     Inf_25 = rowSums(across(attachedInf_25:brokenInf_25), na.rm = TRUE),
#     Flowered_25 = ifelse(Inf_25 > 0, 1, 0)  # 1 if any inflorescence, 0 otherwise
#   )


# Check for duplicate Tag_IDs
datini$Tag_ID <- as.character(datini$Tag_ID)
dat23_spike$Tag_ID <- as.character(dat23_spike$Tag_ID)
dat23_spike %>% count(Tag_ID) %>% filter(n > 1)
dat24_spike$Tag_ID <- as.character(dat24_spike$Tag_ID)
dat25_spike$Tag_ID <- as.character(dat25_spike$Tag_ID)
#view(dat24_spike)
# dat24_spike %>% count(Tag_ID) %>% filter(n > 1)
# dat25_spike %>% count(Tag_ID) %>% filter(n > 1)

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
      dplyr::select(Tag_ID, Flowered_24,Inf_24, Tiller_24, tiller_herb_24, date_24, stroma_24, spikelet_24),
    by = "Tag_ID"
  ) 

dat2425 <- dat24_spike_sp_site_tag %>%
  left_join(
    dat25_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Flowered_25,Inf_25, Tiller_25, tiller_herb_25, date_25, stroma_25, spikelet_25),
    by = "Tag_ID"
  )

# Change variable names
dat2324%>%
  mutate(
    tiller_t = Tiller_23,
    tiller_t1 = Tiller_24,
    Flowered_t=Flowered_23,
    Flowered_t1=Flowered_24,
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
    Flowered_t,
    Flowered_t1,
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
    Flowered_t=Flowered_24,
    Flowered_t1=Flowered_25,
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
    Flowered_t,
    Flowered_t1,
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
climate_site_unique <- climate_site_scaled %>%
  rename(census_year=year,cum_ppt=cumulative_pptmean) %>% 
  distinct(Site, Species, census_year, .keep_all = TRUE)

#view(climate_site_unique)
dat_t_t1_herb$census_year<-as.character(dat_t_t1_herb$census_year)
climate_site_unique$census_year<-as.character(climate_site_unique$census_year)

dat_t_t1_herb_clim <- dat_t_t1_herb %>%
  left_join(
    climate_site_unique,
    by = c("Site", "Species", "census_year")
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

#saveRDS(demography_climate,"/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/demography_climate.rds")

## check for TagIDs  duplicated within years
# demography_climate %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)

## Explore the surival rate per species
# demography_climate %>%
#   group_by(Species) %>%
#   summarise(
#     n_total = n(),
#     n_survived = sum(surv1, na.rm = TRUE),
#     survival_rate = n_survived / n_total
#   )
# names(demography_climate)
# view(demography_climate)
# summary(demography_climate)

## Running the stan model
sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4
)

# Survival----
## Read and format survival data to build the model
demography_climate_surv <- demography_climate %>%
  filter(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, site_species_plot, site_year, Endo, Herbivory,
    tiller_t, surv1, cum_ppt,ppt_scaled
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
    log_size_t0 = ppt_scaled,
    surv_t1     = as.integer(surv1),
    ppt         = ppt_scaled
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


fit_surv_abio_bio_endo_linear <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/survival_l.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)

posterior_surv_abio_bio_endo_linear <- as.array(fit_surv_abio_bio_endo_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_abio_bio_endo_linear,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bendoherbclim[1], bendoherbclim[2], bendoherbclim[3]
                      )
) + theme_bw()

fit_surv_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/survival.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)

#rstan::check_hmc_diagnostics(fit_surv_abio_bio_endo)
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
                        bendoherbclim[1], bendoherbclim[2], bendoherbclim[3],
                        bclim2[1], bclim2[2], bclim2[3]
                      )
) + theme_bw()

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

fit_surv_endo_clim <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/survival_endo_clim.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)

fit_surv_endo_herb <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/survival_endo_herb.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  control = sim_pars$control,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  seed = 13
)
## Save RDS file for further use
# saveRDS(fit_surv_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_abio_bio_endo.rds')
# saveRDS(fit_surv_abio_bio_endo_linear, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_abio_bio_endo_linear.rds')
# saveRDS(fit_surv_endo_clim, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_endo_clim.rds')
# saveRDS(fit_surv_endo_herb, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_endo_herb.rds')

# Growth----
## Read and format survival data to build the model
demography_climate%>%
  filter(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot,site_year, site_species_plot, Endo, Herbivory,
    tiller_t, grow,cum_ppt,ppt_scaled
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
    ppt = ppt_scaled,
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

fit_grow_abio_bio_endo_linear <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/growth_l.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

posterior_grow_abio_bio_endo_linear <- as.array(fit_grow_abio_bio_endo_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_abio_bio_endo_linear,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()


fit_grow_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/growth.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

posterior_grow_abio_bio_endo <- as.array(fit_grow_abio_bio_endo) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_abio_bio_endo,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3]
                      )
) + theme_bw()


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
                        bclim2[1], bclim2[2], bclim2[3]
                      )
) + theme_bw()

fit_grow_endo_clim <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/growth_endo_clim.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

fit_grow_endo_herb <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/growth_endo_herb.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

## Save RDS file for further use
# saveRDS(fit_grow_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_abio_bio_endo.rds')
# saveRDS(fit_grow_abio_bio_endo_linear, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_abio_bio_endo_linear.rds')
# saveRDS(fit_grow_endo_clim, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_endo_clim.rds')
# saveRDS(fit_grow_endo_herb, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_endo_herb.rds')

# Flowering----
demography_climate %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cum_ppt,ppt_scaled,Flowered_t1
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
    flow_t1 = Flowered_t1,
    ppt = ppt_scaled
  ) -> demography_climate_flow
#sum(demography_climate_flow$flow_t1 == 0)
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

fit_flow_abio_bio_endo_linear <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering_l.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

posterior_flow_abio_bio_endo_linear <- as.array(fit_flow_abio_bio_endo_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_abio_bio_endo_linear,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                        
                      )
) + theme_bw()


fit_flow_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

#summary(fit_flow_ppt)$summary[, c("Rhat", "n_eff")]
posterior_flow_abio_bio_endo <- as.array(fit_flow_abio_bio_endo) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_abio_bio_endo,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3]
                      )
) + theme_bw()

## Save RDS file for further use
# saveRDS(fit_flow_abio_bio_endo_linear, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_abio_bio_endo_linear.rds')
# saveRDS(fit_flow_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_abio_bio_endo.rds')

# Inflorescence----
demography_climate %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cum_ppt,ppt_scaled
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
    inf_t1 = inf_t1,
    ppt = ppt_scaled
  ) -> demography_climate_inf
#sum(demography_climate_flow$inf_t1 == 0)
## Separate each variable to use the same model stan
### Precipitation 
demography_inf_ppt <- list(
  nSpp = demography_climate_inf$Species %>% n_distinct(),
  nSite = demography_climate_inf$Site %>% n_distinct(),
  nsite_year = demography_climate_inf$site_year %>% n_distinct(),
  nPop = demography_climate_inf$Population %>% n_distinct(),
  nPlot = demography_climate_inf$site_species_plot %>% n_distinct(),
  Spp = demography_climate_inf$Species,
  site = demography_climate_inf$Site,
  site_year = demography_climate_inf$site_year,
  pop = demography_climate_inf$Population,
  plot = demography_climate_inf$site_species_plot,
  clim = as.vector(demography_climate_inf$ppt),
  endo = demography_climate_inf$Endo,
  herb = demography_climate_inf$Herbivory,
  size = demography_climate_inf$log_size_t0,
  y = demography_climate_inf$inf_t1,
  N = nrow(demography_climate_inf)
)

fit_inf_abio_bio_endo_linear <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/inflorescence_l.stan",
  data = demography_inf_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

posterior_inf_abio_bio_endo_linear <- as.array(fit_inf_abio_bio_endo_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_inf_abio_bio_endo_linear,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                       
                      )
) + theme_bw()


fit_inf_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/inflorescence.stan",
  data = demography_inf_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

#summary(fit_flow_ppt)$summary[, c("Rhat", "n_eff")]
posterior_flow_abio_bio_endo <- as.array(fit_flow_abio_bio_endo) # Converts to an array
bayesplot::mcmc_trace(posterior_flow_abio_bio_endo,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3],
                        bclim2[1], bclim2[2], bclim2[3]
                      )
) + theme_bw()

fit_inf_endo_clim <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/inflorescence_endo_clim.stan",
  data = demography_inf_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

fit_inf_endo_herb <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/inflorescence_endo_herb.stan",
  data = demography_inf_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control = sim_pars$control,
  seed = 13)

## Save RDS file for further use
# saveRDS(fit_inf_abio_bio_endo_linear, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_inf_abio_bio_endo_linear.rds')
# saveRDS(fit_inf_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_inf_abio_bio_endo.rds')
saveRDS(fit_inf_endo_clim, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_inf_endo_clim.rds')
saveRDS(fit_inf_endo_herb, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_inf_endo_herb.rds')

# Spikelet----
demography_climate %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  filter(tiller_t1 > 0,inf_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cum_ppt,ppt_scaled,inf_t1
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
    ppt = ppt_scaled
  ) -> demography_climate_spik

#sum(demography_climate_spik$spikelet_t1 == 0)
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

fit_spik_abio_bio_endo_linear <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet_l.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control =sim_pars$control,
  seed = 13)

posterior_spik_abio_bio_endo_linear <- as.array(fit_spik_abio_bio_endo_linear) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_abio_bio_endo_linear,
                      pars = quote_bare(
                        b0[1], b0[2],
                        bendo[1], bendo[2],
                        bherb[1], bherb[2],
                        bclim[1], bclim[2],
                        bendoclim[1], bendoclim[2],
                        bendoherb[1], bendoherb[2],
                        bendoherbclim[1],bendoherbclim[2]
                      )
) + theme_bw()

fit_spik_abio_bio_endo <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control =sim_pars$control,
  seed = 13)

# summary(fit_spik_ppt)$summary[, c("Rhat", "n_eff")]
posterior_spik_abio_bio_endo <- as.array(fit_spik_abio_bio_endo) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_abio_bio_endo,
                      pars = quote_bare(
                        b0[1], b0[2],
                        bendo[1], bendo[2],
                        bherb[1], bherb[2],
                        bclim[1], bclim[2],
                        bendoclim[1], bendoclim[2],
                        bendoherb[1], bendoherb[2],
                        bclim2[1], bclim2[2],
                        bendoherbclim[1],bendoherbclim[2]
                      )
) + theme_bw()

fit_spik_endo_clim <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet_endo_clim.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control =sim_pars$control,
  seed = 13)
fit_spik_endo_herb <- stan(
  file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet_endo_herb.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  chains = sim_pars$chains,
  control =sim_pars$control,
  seed = 13)
## Save RDS file for further use
# saveRDS(fit_spik_abio_bio_endo_linear, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_abio_bio_endo_linear.rds')
# saveRDS(fit_spik_abio_bio_endo, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_abio_bio_endo.rds')
saveRDS(fit_spik_endo_clim, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_endo_clim.rds')
saveRDS(fit_spik_endo_herb, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_endo_herb.rds')

# Posterior predictive check----
# Quadratic models
## Helper functions
inv_logit <- function(x) 1 / (1 + exp(-x))

# Function to simulate posterior predictive samples for each vital rate
# Function to simulate posterior predictive samples for each vital rate
simulate_ppc <- function(pred_matrix, phi = NULL, zi = NULL,
                         family = c("bernoulli", "normal", "negbinomial", "zinb")) {
  family <- match.arg(family)
  n_iter <- nrow(pred_matrix)
  N <- ncol(pred_matrix)
  y_rep <- matrix(NA, nrow = n_iter, ncol = N)
  
  if (family == "bernoulli") {
    probs <- inv_logit(pred_matrix)
    for (i in 1:n_iter) {
      y_rep[i, ] <- rbinom(N, size = 1, prob = probs[i, ])
    }
    
  } else if (family == "normal") {
    for (i in 1:n_iter) {
      y_rep[i, ] <- rnorm(N, mean = pred_matrix[i, ], sd = phi[i])
    }
    
  } else if (family == "negbinomial") {
    for (i in 1:n_iter) {
      mu_i <- exp(pred_matrix[i, ])
      size_i <- phi[i]
      y_rep[i, ] <- rnbinom(N, size = size_i, mu = mu_i)
    }
    
  } else if (family == "zinb") {
    for (i in 1:n_iter) {
      mu_i <- exp(pred_matrix[i, ])
      size_i <- phi[i]
      # zi is scalar: same zero-inflation probability for all observations
      is_zero <- rbinom(N, size = 1, prob = zi)
      y_nb <- rnbinom(N, size = size_i, mu = mu_i)
      y_rep[i, ] <- ifelse(is_zero == 1, 0, y_nb)
    }
  }
  
  return(y_rep)
}
bayesplot::color_scheme_set("blue")
## Survival Model full model----
fit_surv_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/tsyih2vbg04zf9odu2goe/fit_surv_abio_bio_endo.rds?rlkey=tp4no6tfv5mb6f85jwtkqfuyg&dl=1"))
post_surv <- rstan::extract(fit_surv_abio_bio_endo)
pred_surv <- post_surv$predS
y_surv <- demography_surv_ppt$y
y_rep_surv <- simulate_ppc(pred_surv, family = "bernoulli")
p_surv <- ppc_dens_overlay(y_surv, y_rep_surv) + ggtitle("Survival")

## Growth Model full model----
fit_grow_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/4r5062xbfc66gqh5l6xbz/fit_grow_abio_bio_endo.rds?rlkey=spnn4nj0zvzfsss1kn3qsnocj&dl=1"))
post_grow <- rstan::extract(fit_grow_abio_bio_endo)
pred_grow <- post_grow$predG
sigma_grow <- post_grow$sigma
y_grow <- demography_grow_ppt$y
y_rep_grow <- simulate_ppc(pred_grow, phi = sigma_grow, family = "normal")
p_grow <- ppc_dens_overlay(y_grow, y_rep_grow) + ggtitle("Growth")


## Flowering Model full model ----
fit_flow_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/5717xz8nt6sph3neq6jj9/fit_flow_abio_bio_endo.rds?rlkey=p4s7391sdqgepd89x53tbgw82&dl=1"))
post_flow <- rstan::extract(fit_flow_abio_bio_endo)
pred_flow <- post_flow$predS      # linear predictor
y_flow <- demography_flow_ppt$y
y_rep_flow <- simulate_ppc(pred_flow, family = "bernoulli")
p_flow <- ppc_dens_overlay(y_flow, y_rep_flow) + ggtitle("Flowering") 

## Inflorescence Model full model ----
# fit_inf_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/rnkijsri04jtrczshfmp6/fit_inf_abio_bio_endo.rds?rlkey=kyxj4f6mpwdp5m78wi0p6v64r&dl=1"))
# post_inf <- rstan::extract(fit_inf_abio_bio_endo)
# pred_inf <- post_inf$predF      # linear predictor
# phi_inf <- post_inf$phi         # dispersion
# zi_inf <- post_inf$zi           # scalar zero-inflation
# y_inf <- demography_inf_ppt$y
# y_rep_inf <- simulate_ppc(pred_inf, phi = phi_inf, zi = zi_inf, family = "zinb")
# p_inf <- ppc_dens_overlay(y_inf, y_rep_inf) + ggtitle("Inflorescence")+xlim(0, 75) 

## Spikelet Model full model ----
fit_spik_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/pebuc3ysvv9rrr2dvihl2/fit_spik_abio_bio_endo.rds?rlkey=f1l4q9ucvk4h4236600zoyjci&dl=1"))
post_spik <- rstan::extract(fit_spik_abio_bio_endo)
pred_spik <- post_spik$predF
#zi_spik <- post_spik$zi           # scalar zero-inflation
phi_spik <- post_spik$phi
y_spik <- demography_spik_ppt$y
y_rep_spik <- simulate_ppc(pred_spik, phi = phi_spik,family = "negbinomial")
p_spik <- ppc_dens_overlay(y_spik, y_rep_spik) + ggtitle("Spikelet")

# Combine all PPC plots full model ----
combined_plot <- (p_surv | p_grow) / (p_flow | p_spik) +
  plot_annotation(title = "") +
  plot_layout(guides = "collect") &
  theme_light() +
  theme(legend.position = "bottom")

print(combined_plot)

# Save to PDF
ggsave(filename = "/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/combined_ppc_plots.pdf", plot = combined_plot, width = 7, height = 6)

# Monotonic models
## Survival Model linear ----
fit_surv_linear <- readRDS(url("https://www.dropbox.com/scl/fi/khyez2xegn8j5elkchaf6/fit_surv_abio_bio_endo_linear.rds?rlkey=zh4hx9czjov9aivlmaycfcuq7&dl=1"))
post_surv_linear <- rstan::extract(fit_surv_linear)
pred_surv_linear <- post_surv_linear$predS
y_surv_linear <- demography_surv_ppt$y
y_rep_surv_linear <- simulate_ppc(pred_surv_linear, family = "bernoulli")
p_surv_linear <- ppc_dens_overlay(y_surv_linear, y_rep_surv_linear) + ggtitle("Survival")

## Growth Model linear ----
fit_grow_linear <- readRDS(url("https://www.dropbox.com/scl/fi/o62tvjf8aqqz15gjxnrjn/fit_grow_abio_bio_endo_linear.rds?rlkey=xg1s6u5ctsluampm1l2zy1wqn&dl=1"))
post_grow_linear <- rstan::extract(fit_grow_linear)
pred_grow_linear <- post_grow_linear$predG
sigma_grow_linear <- post_grow_linear$sigma
y_grow_linear <- demography_grow_ppt$y
y_rep_grow_linear <- simulate_ppc(pred_grow_linear, phi = sigma_grow_linear, family = "normal")
p_grow_linear <- ppc_dens_overlay(y_grow_linear, y_rep_grow_linear) + ggtitle("Growth")

## Flowering Model linear ----
fit_flow_linear <- readRDS(url("https://www.dropbox.com/scl/fi/1v4f4thyh826qcuiiyhub/fit_flow_abio_bio_endo_linear.rds?rlkey=raj4ls5dcqkeeexvcqj8b495m&dl=1"))
post_flow_linear <- rstan::extract(fit_flow_linear)
pred_flow_linear <- post_flow_linear$predS
y_flow_linear <- demography_flow_ppt$y
y_rep_flow_linear <- simulate_ppc(pred_flow_linear, family = "bernoulli")
p_flow_linear <- ppc_dens_overlay(y_flow_linear, y_rep_flow_linear) + ggtitle("Inflorescence")+xlim(0,100)

# ## Inflorescence Model linear ----
# fit_inf_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6ngnypika10yc0jrvr531/fit_inf_abio_bio_endo_linear.rds?rlkey=822f09t91dd2jw4r8svwmdh29&dl=1"))
# post_inf_linear <- rstan::extract(fit_inf_linear)
# pred_inf_linear <- post_inf_linear$predF
# zi_inf_linear <- post_inf_linear$zi           # scalar zero-inflation
# phi_inf_linear <- post_inf_linear$phi
# y_inf_linear <- demography_inf_ppt$y
# y_rep_inf_linear <- simulate_ppc(pred_inf_linear, zi =zi_inf_linear,phi = phi_inf_linear, family = "zinb")
# p_inf_linear <- ppc_dens_overlay(y_inf_linear, y_rep_inf_linear[1:500, ]) + ggtitle("Inflorescence")+xlim(0,100)

## Spikelet Model linear ----
fit_spik_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6fcebl4lw8mu94fz62hnh/fit_spik_abio_bio_endo_linear.rds?rlkey=zy25y44zocugs6shh68lwpy1q&dl=1"))
post_spik_linear <- rstan::extract(fit_spik_linear)
pred_spik_linear <- post_spik_linear$predF
zi_spik_linear <- post_spik_linear$zi           # scalar zero-inflation
phi_spik_linear <- post_spik_linear$phi
y_spik_linear <- demography_spik_ppt$y
y_rep_spik_linear <- simulate_ppc(pred_spik_linear,zi =y_spik_linear, phi = phi_spik_linear, family = "negbinomial")
p_spik_linear <- ppc_dens_overlay(y_spik_linear, y_rep_spik_linear) + ggtitle("Spikelet")

## Combine all PPC plots linear ----
combined_plot_linear <- (p_surv_linear | p_grow_linear) / (p_flow_linear | p_spik_linear) +
  plot_annotation(title = "") +
  plot_layout(guides = "collect") &
  theme_light() +
  theme(legend.position = "bottom")

print(combined_plot_linear)

## Save to PDF linear ----
ggsave(
  filename = "/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/combined_ppc_plots_linear.pdf",
  plot = combined_plot_linear,
  width = 7,
  height = 6
)

