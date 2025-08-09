# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis on  vital rate models (survival, growth, flowering and spikelet).
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
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
# unique(datini$Site)

# calculate the average spikelet and inflorescence number for each census
dat23 %>%
  mutate(spikelet_23 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T)), digit = 0) -> dat23_spike
dat24 %>%
  mutate(spikelet_24 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_24 = round(rowMeans(across(attachedInf_24:brokenInf_24), na.rm = T), digit = 0)) -> dat24_spike
dat25 %>%
  mutate(spikelet_25 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_25 = round(rowMeans(across(attachedInf_25:brokenInf_25), na.rm = T), digit = 0)) -> dat25_spike

# Check for duplicate Tag_IDs
datini$Tag_ID <- as.character(datini$Tag_ID)
dat23_spike$Tag_ID <- as.character(dat23_spike$Tag_ID)
dat23_spike %>% count(Tag_ID) %>% filter(n > 1)
dat24_spike$Tag_ID <- as.character(dat24_spike$Tag_ID)
dat25_spike$Tag_ID <- as.character(dat25_spike$Tag_ID)
#view(dat24_spike)
# Here 731 was two row to I change one with 7311
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

# Combine datini and dat23
combined_data <- bind_rows(datini[,c("Site","Species","Plot","Tag_ID","Population","Endo")], datini23_spike[,c("Site","Species","Plot","Tag_ID","Population","Endo" )])%>% 
  distinct(Tag_ID, .keep_all = TRUE)

dat24_spike_sp_site_tag<-left_join(x = dat24_spike, y = combined_data, by = c("Tag_ID")) %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_24","brokenInf_24"))) %>% 
  filter(!is.na(Species))
# view(dat24_spike_sp_site_tag)

dat25_spike_sp_site_tag<- dat25_spike %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_25","brokenInf_25"))) 

dat2324 <- datini23_spike %>%
  left_join(
    dat24_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Inf_24, Tiller_24, tiller_herb_24, date_24, stroma_24, spikelet_24),
    by = "Tag_ID"
  ) %>%
  mutate(
    census_year = if_else(!is.na(date_24), 2024L, 2023L)
  )

dat2425 <- dat24_spike_sp_site_tag %>%
  left_join(
    dat25_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Inf_25, Tiller_25, tiller_herb_25, date_25, stroma_25, spikelet_25),
    by = "Tag_ID"
  ) %>%
  mutate(
    census_year = if_else(!is.na(date_25), 2025L, 2024L)
  )

## Merge the demographic data with the climatic data -----
dat2324_clim <- dat2324 %>%
  left_join(
    climate_site,
    by = c("Site" = "Site", "Species" = "Species", "census_year" = "year")
  )

dat2425_clim <- dat2425 %>%
  left_join(
    climate_site,
    by = c("Site" = "Site", "Species" = "Species", "census_year" = "year")
  )

# Change variable names
dat2324_clim %>%
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
    date_t1=date_24
  ) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Tag_ID,
    Population,
    Clone,
    GreenhouseID,
    Endo,
    tiller_t,
    tiller_t1,
    inf_t,
    inf_t1,
    spikelet_t,
    spikelet_t1,
    tiller_Herb_t,
    tiller_Herb_t1,
    cumulative_pptmean,
    census_year,
    date_t,
    date_t1
  ) -> dat2324_t_t1


dat2425_clim %>%
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
    date_t1 = date_25
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
    cumulative_pptmean,
    census_year,
    date_t,
    date_t1
  ) -> dat2425_t_t1

dat_t_t1 <- bind_rows(dat2324_t_t1, dat2425_t_t1) %>%
  mutate(
    date_t = ymd(date_t),
    date_t1 = ymd(date_t1)
  )
## Merge the demographic data with the herbivory data -----
dat_t_t1_herb <- left_join(x = dat_t_t1, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
# head(dat_t_t1_herb)
# unique(dat_t_t1_herb$Species)
# view(dat2324_t_t1_herb)

# Check the number of data for each data
dat_t_t1_herb %>%
  #filter(tiller_t1 > 0) %>%
  dplyr::select(Species, census_year, tiller_t1) %>%
  group_by(Species, census_year) %>%
  summarise(n = n(), .groups = "drop")
##unclear why 2023 is kep here as a transition year, and why it has such small sample size

## Create new variables
dat_t_t1_herb %>%
  mutate(
    site_year=interaction(Site,census_year),
    surv1 = 1 * (!is.na(tiller_t) & !is.na(tiller_t1)),
    site_species_plot = interaction(Site, Species, Plot),
    grow = (log(tiller_t1 + 1) - log(tiller_t + 1))
  ) -> demography_climate

{ ##TOM
  ##look at "mortality"
demography_climate %>% filter(surv1==0) %>% View
demography_climate  %>% 
  filter(surv1==0) %>% 
  summarise(n())
## most of these have NA for tiller_t1 because they died the *preceding year*
demography_climate  %>% 
  filter(surv1==0) %>% 
  filter(tiller_t==0) %>% 
  summarise(n())
## the others have NA because there are in plots that were destroyed, mainly KER but also the flooded plot in HUN
demography_climate %>% 
  filter(surv1==0) %>% 
  filter(tiller_t>0) %>% 
  group_by(Species,Site,census_year) %>% 
  summarise(n())
## I also noticed that TagIDs are often duplicated within years, which should not happen
demography_climate %>% 
  group_by(Tag_ID,census_year) %>% 
  summarise(tag_rep = n()) %>% 
  filter(tag_rep>1)
## look at one of these as an example -- here, the 2025 data are duplicated. In other cases it's 2024.
demography_climate %>% filter(Tag_ID==119) %>% View

##once the duplicate tag issues are fixed, I think the right way to estimate mortality and growth would look something like this:
dat_t_t1_herb %>%
  mutate(
    site_year=interaction(Site,census_year),
    ##if the plant was dead or size was NA at the start of the transition year, survival is NA
    ##if the plant was alive at the start of the transition year, it survived if tillers_t1>0
    surv1 = ifelse(tiller_t>0,tiller_t1>0,NA),
    site_species_plot = interaction(Site, Species, Plot),
    ##if the plant was alive at the start of the transition year and it survived, growth is the log ratio of tiller counts, else NA
    ##note there that growth is conditional on survival, which I think is how it should be
    grow = ifelse(tiller_t>0 & tiller_t1>0,log(tiller_t1/tiller_t),NA)
  ) -> demography_climate
}

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
demography_climate %>%
  subset(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site,site_species_plot, site_year, Endo, Herbivory,
    tiller_t, surv1, cumulative_pptmean
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_year = as.integer(factor(site_year)),
    Endo = as.integer(factor(Endo)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Herbivory = as.integer(factor(Herbivory))
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    surv_t1 = surv1,
    ppt = log(cumulative_pptmean),
  ) -> demography_climate_surv

### Convert into list to run in stan
demography_surv_ppt <- list(
  nSpp = demography_climate_surv$Species %>% n_distinct(),
  nSite = demography_climate_surv$Site %>% n_distinct(),
  nsite_year = demography_climate_surv$site_year %>% n_distinct(),
  nPop = demography_climate_surv$Population %>% n_distinct(),
  nPlot = demography_climate_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_surv$Species,
  site = demography_climate_surv$Site,
  site_year = demography_climate_surv$site_year,
  pop= demography_climate_surv$Population,
  plot = demography_climate_surv$site_species_plot,
  clim = as.vector(demography_climate_surv$ppt),
  endo = demography_climate_surv$Endo -1 ,
  herb = demography_climate_surv$Herbivory -1,
  size = demography_climate_surv$log_size_t0,
  y = demography_climate_surv$surv_t1,
  N = nrow(demography_climate_surv)
)

# fit_surv_ppt <- stan(
#   file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/survival.stan",
#   data = demography_surv_ppt,
#   warmup = sim_pars$warmup,
#   control = sim_pars$control,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains
# )

# rstan::check_hmc_diagnostics(fit_surv_ppt)
summary(fit_surv_ppt)$summary[, c("Rhat", "n_eff")]
posterior_surv_ppt <- as.array(fit_surv_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_ppt,
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
# saveRDS(fit_surv_ppt, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_surv_ppt.rds')


# Growth----
## Read and format survival data to build the model
demography_climate%>%
  subset(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot,site_year, site_species_plot, Endo, Herbivory,
    tiller_t, grow,cumulative_pptmean
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    site_year = as.integer(factor(site_year)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    grow = grow,
    ppt = log(cumulative_pptmean),
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

# fit_grow_ppt <- stan(
#   file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/growth.stan",
#   data = demography_grow_ppt,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

summary(fit_grow_ppt)$summary[, c("Rhat", "n_eff")]
posterior_grow_ppt <- as.array(fit_grow_ppt) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_ppt,
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
# saveRDS(fit_grow_ppt, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_grow_ppt.rds')

# Flowering----
demography_climate %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cumulative_pptmean
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    site_year=as.integer(factor(site_year)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    flow_t1 = inf_t1,
    ppt = log(cumulative_pptmean)
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

# fit_flow_ppt <- stan(
#   file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering.stan",
#   data = demography_flow_ppt,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control = sim_pars$control)

summary(fit_flow_ppt)$summary[, c("Rhat", "n_eff")]
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
saveRDS(fit_flow_ppt, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_flow_ppt.rds')

# Spikelet----
demography_climate %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cumulative_pptmean
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    site_year=as.integer(factor(site_year)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    spi_t1 = spikelet_t1,
    ppt = log(cumulative_pptmean)
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

# fit_spik_ppt <- stan(
#   file = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet.stan",
#   data = demography_spik_ppt,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains,
#   control =sim_pars$control)

summary(fit_spik_ppt)$summary[, c("Rhat", "n_eff")]
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
# saveRDS(fit_spik_ppt, '/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/fit_spik_ppt.rds')

# Posterior predictive check----
library(rstan)
library(bayesplot)
library(patchwork)  # for combined plots

set.seed(123)

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

# Load and process each model ---------------------------------------------

library(rstan)
library(bayesplot)
library(patchwork)  # for combined plots

set.seed(123)

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

## Survival Model ----
fit_survival <- readRDS(url("https://www.dropbox.com/scl/fi/9xa46n5v7u1lddxaj69cs/fit_surv_ppt.rds?rlkey=1mbkby4394s04j7qej4kvzeo2&dl=1"))
post_surv <- rstan::extract(fit_survival)
pred_surv <- post_surv$predS
y_surv <- demography_surv_ppt$y
y_rep_surv <- simulate_ppc(pred_surv, family = "bernoulli")
p_surv <- ppc_dens_overlay(y_surv, y_rep_surv[1:500, ]) + ggtitle("Survival")

## Growth Model ----
fit_growth <- readRDS(url("https://www.dropbox.com/scl/fi/d0x30lqqcxnatupsm2hej/fit_grow_ppt.rds?rlkey=er2is1le25trin73an23ztfgm&dl=1"))
post_grow <- rstan::extract(fit_growth)
pred_grow <- post_grow$predG
sigma_grow <- post_grow$sigma
y_grow <- demography_grow_ppt$y
y_rep_grow <- simulate_ppc(pred_grow, phi = sigma_grow, family = "normal")
p_grow <- ppc_dens_overlay(y_grow, y_rep_grow[1:500, ]) + ggtitle("Growth")

## Flowering Model ----
fit_flowering <- readRDS(url("https://www.dropbox.com/scl/fi/z4kh899krlz2ssz1oimxk/fit_flow_ppt.rds?rlkey=uqmm05uas6czzb9siyg7unp1r&dl=1"))
post_flow <- rstan::extract(fit_flowering)
pred_flow <- post_flow$predF
phi_flow <- post_flow$phi
y_flow <- demography_flow_ppt$y
y_rep_flow <- simulate_ppc(pred_flow, phi = phi_flow, family = "negbinomial")
p_flow <- ppc_dens_overlay(y_flow, y_rep_flow[1:500, ]) + ggtitle("Flowering")

## Spikelet Model ----
fit_spikelet <- readRDS(url("https://www.dropbox.com/scl/fi/1ar22b0jf1urcnypduybm/fit_spik_ppt.rds?rlkey=5hiwrlm6wc2iwvqxma9mx5r15&dl=1"))
post_spik <- rstan::extract(fit_spikelet)
pred_spik <- post_spik$predF
phi_spik <- post_spik$phi
y_spik <- demography_spik_ppt$y
y_rep_spik <- simulate_ppc(pred_spik, phi = phi_spik, family = "negbinomial")
p_spik <- ppc_dens_overlay(y_spik, y_rep_spik[1:500, ]) + ggtitle("Spikelet")

# Combine all PPC plots ----
combined_plot <- (p_surv | p_grow) / (p_flow | p_spik) +
  plot_annotation(title = "")
print(combined_plot)
# Save to PDF
ggsave(filename = "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/combined_ppc_plots.pdf", plot = combined_plot, width = 8, height = 6)
