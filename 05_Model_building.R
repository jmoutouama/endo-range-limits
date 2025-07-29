# Project:
# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on  vital rate models (survival, growth, flowering and spikelet).
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
library(bayesplot)
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
# Define some basic functions that we'll use later
quote_bare <- function(...) {
  substitute(alist(...)) %>%
    eval() %>%
    sapply(deparse)
}
set.seed(13)
# Demographic data -----
# Merge the demographic census
datini <- read.csv("https://www.dropbox.com/scl/fi/9lemc8zv5fdi07r7j1jqq/Initialdata.csv?rlkey=gmryfg74motuf7h49s1aul0tm&dl=1", stringsAsFactors = F)
dat23 <- read.csv("https://www.dropbox.com/scl/fi/o18136j5t9sd94irtd46i/census2023.csv?rlkey=7byyglfehzqo5e197fp1jsasw&dl=1", stringsAsFactors = F)
dat24 <- read.csv("https://www.dropbox.com/scl/fi/w6qdatnkosrjsvvub1vsf/census2024.csv?rlkey=tt1qbdaratj2zvt4igkxz89n9&dl=1", stringsAsFactors = F)
dat25<-read.csv("https://www.dropbox.com/scl/fi/x6gttyivrogzug0w3iveo/census_2025.csv?rlkey=jidxofw54907a73gdnn7t9n3i&dl=1", stringsAsFactors = F)
datherbivory <- read.csv("https://www.dropbox.com/scl/fi/mcf4v8voqi8xpmum6l9qv/herbivory.csv?rlkey=my3araheil7cc7iajp2k7hcbt&dl=1", stringsAsFactors = F)
# unique(datini$Site)
# unique(datini$dat23)
# unique(datini$dat24)
# unique(dat25$Site)
# names(dat23)

# calculate the average spikelet and inflorescence number for each census
dat23 %>%
  mutate(spikelet_23 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T)), digit = 0) -> dat23_spike

dat24 %>%
  mutate(spikelet_24 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_24 = round(rowMeans(across(attachedInf_24:brokenInf_24), na.rm = T), digit = 0)) -> dat24_spike
dat25 %>%
  mutate(spikelet_25 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_25 = round(rowMeans(across(attachedInf_25:brokenInf_25), na.rm = T), digit = 0)) -> dat25_spike

# Check for duplicate Tag_IDs
datini$Tag_ID <- as.integer(datini$Tag_ID)
dat23_spike$Tag_ID <- as.integer(dat23_spike$Tag_ID)
dat24_spike$Tag_ID <- as.integer(dat24_spike$Tag_ID)
dat25_spike$Tag_ID <- as.integer(dat25_spike$Tag_ID)

dat23_spike %>% count(Tag_ID) %>% filter(n > 1)
# Here 731 was two row to I change one with 7311
dat24_spike %>% count(Tag_ID) %>% filter(n > 1)
dat25_spike %>% count(Tag_ID) %>% filter(n > 1)

## Merge the initial data to the 23, the 23 data with the 24 data and the 24 data with the 25 -----
datini23 <- right_join(x = datini, y = dat23_spike, by = c("Tag_ID"))
#anti_join(dat23_spike, datini, by = "Tag_ID")

# Remove rows with NA in Tag_ID
datini23_spike <- datini23 %>% filter(!is.na(Tag_ID))
dat2324 <- left_join(x = datini23_spike, y = dat24_spike, by = c("Tag_ID"))
#names(dat2324)

dat2425 <- left_join(dat2324, dat25_spike, by = "Tag_ID")

# Check mismatches
mismatch_sites <- dat2425 %>%
  filter(Site.x != Site.y)

if(nrow(mismatch_sites) == 0) {
  # Safe to keep one Site column only
  dat2425 <- dat2425 %>%
    dplyr::select(-Site.y) %>%
    rename(Site = Site.x)
} else {
  warning("Sites do not match for all Tag_IDs!")
  # Decide how to handle mismatches
}

#names(dat2425)

# Change variable names
dat2324 %>%
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
    Position,
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
    date_t,
    date_t1
  ) -> dat2324_t_t1


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
    Species = Species.x,  # fixed here
    date_t = date_24,
    date_t1 = date_25
  ) %>%
  dplyr::select(-Species.y) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Position,
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
    date_t,
    date_t1
  ) -> dat2425_t_t1

dat_t_t1<-rbind(dat2324_t_t1,dat2425_t_t1)

## Merge the demographic data with the herbivory data -----
dat_t_t1_herb <- left_join(x = dat_t_t1, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
head(dat_t_t1_herb)
unique(dat_t_t1_herb$Species)
# view(dat2324_t_t1_herb)

dat_t_t1_herb %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(Species, tiller_t1) %>%
  group_by(Species) %>%
  summarise(n = sum(tiller_t1, na.rm = T))

# Climatic data ----
climate_summary <- readRDS(url("https://www.dropbox.com/scl/fi/kfh1pnchc9jt8cw77co27/prism_means.rds?rlkey=s514h0gr61knu9mdmfhjo37nu&dl=1"))
climate_summary %>%
  rename(Site = site) -> climate_site
distance_species <- readRDS(url("https://www.dropbox.com/scl/fi/x9lb4outzenaugw6e00w7/distance_species.rds?rlkey=yrz9djcncfz0omw5dxwsx81ds&dl=1"))
distance_species %>%
  rename(Site = site_code) -> distance_species_clean

## Merge the demographic data with the climatic data -----
demography_climate <- left_join(x = dat_t_t1_herb, y = climate_site, by = c("Site"))
demography_climate_distance <- left_join(x = demography_climate, y = distance_species_clean, by = c("Site", "Species"))

## Create new variables
demography_climate_distance %>%
  mutate(
    surv1 = 1 * (!is.na(demography_climate_distance$tiller_t) & !is.na(demography_climate_distance$tiller_t1)),
    site_species_plot = interaction(demography_climate_distance$Site, demography_climate_distance$Species, demography_climate_distance$Plot),
    grow = (log(demography_climate_distance$tiller_t1 + 1) - log(demography_climate_distance$tiller_t + 1))
  ) -> demography_climate_distance

# names(demography_climate)
# view(demography_climate)
# summary(demography_climate)
# Diagnostic of the response variable
ggplot(demography_climate_distance, aes(x = spikelet_t1, fill = Species, color = Species)) +
  geom_density(alpha = 0.4) + # Density plot with transparency
  labs(
    x = "Spikelet average", y = "Density",
    title = ""
  ) +
  theme_bw() +
  # scale_color_manual(values = c("red", "#00AFBB")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "#00AFBB")) +
  theme(
    legend.position = c(0.8, 0.8), # Legend at the top
    legend.background = element_rect(fill = "white", color = "white"), # Optional: outline the legend
    legend.title = element_text(size = 10), # Optional: adjust legend title size
    legend.text = element_text(size = 8) # Optional: adjust legend text size
  )

ggplot(demography_climate_distance, aes(x = grow, fill = Species, color = Species)) +
  geom_density(alpha = 0.4) + # Density plot with transparency
  labs(
    x = "Relative growth", y = "Density",
    title = ""
  ) +
  theme_bw() +
  # scale_color_manual(values = c("red", "#00AFBB")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "#00AFBB")) +
  theme(
    legend.position = c(0.8, 0.8), # Legend at the top
    legend.background = element_rect(fill = "white", color = "white"), # Optional: outline the legend
    legend.title = element_text(size = 10), # Optional: adjust legend title size
    legend.text = element_text(size = 8) # Optional: adjust legend text size
  )
# Plot the data with regression lines
filtered_data <- demography_climate_distance %>% filter(Species %in% c("ELVI", "POAU"))
ggplot(filtered_data, aes(x = spikelet_t1, y = inf_t1, color = Species)) +
  geom_point() + # Plot points
  geom_smooth(method = "lm", se = TRUE, aes(color = Species)) + # Add regression line
  labs(x = "Spikelet average", y = "Inflorescence average", title = "") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))

# Probability density plot for Inflorescence grouped by species
ggplot(filtered_data, aes(x = inf_t1, fill = Species, color = Species)) +
  geom_density(alpha = 0.4) + # Density plot with transparency
  labs(
    x = "Inflorescence average", y = "Density",
    title = ""
  ) +
  theme_bw() +
  # scale_color_manual(values = c("red", "#00AFBB")) +  # Custom colors for the species
  # scale_fill_manual(values = c("red", "#00AFBB")) +
  theme(
    legend.position = c(0.8, 0.8), # Legend at the top
    legend.background = element_rect(fill = "white", color = "black"), # Optional: outline the legend
    legend.title = element_text(size = 10), # Optional: adjust legend title size
    legend.text = element_text(size = 8) # Optional: adjust legend text size
  )

# I am not sure what prior to use.So I am going to simulate some data and found that the difference
# Define the range of tau values
x <- seq(0.001, 5, length.out = 1000)
# Compute densities for different priors
df <- data.frame(
  x = rep(x, 4),
  density = c(
    dinvgamma(x, 0.1, 0.1),
    dinvgamma(x, 2, 1),
    dhalfcauchy(x, scale = 1),
    dnorm(x, mean = 0, sd = 1) * 2
  ), # Half-Normal (mirrored normal)
  Prior = rep(c(
    "Inv-Gamma(0.1, 0.1)",
    "Inv-Gamma(2, 1)",
    "Half-Cauchy(1)",
    "Half-Normal(1)"
  ), each = length(x))
)

# Plot the priors
ggplot(df, aes(x, density, color = Prior)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(
    title = "",
    x = "tau",
    y = "Density",
    color = "Prior"
  ) +
  theme(legend.position = c(0.8, 0.8))
# The inverse gamma(0.1, 0.1) prior has a very heavy tail, meaning it gives high probability to very large values, which can cause instability.
# The inverse gamma(2, 1) is more reasonable, providing some regularization while still allowing moderate values.
# The half-Cauchy(1) prior has a fat tail but is more controlled compared to inv_gamma(0.1, 0.1).
# The half-Normal(1) prior is much more concentrated near small values, leading to stronger regularization.

## Running the stan model
sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  thin = 2,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# Survival----
## Read and format survival data to build the model
demography_climate_distance %>%
  subset(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, surv1, sum_ppt, mean_temp, mean_vpd, distance,geo_distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)),
    Herbivory = as.integer(factor(Herbivory))
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    surv_t1 = surv1,
    ppt = log(sum_ppt),
    temp = log(mean_temp),
    vpd = mean_vpd,
    distance = log(distance),
    geo_distance=log(geo_distance)
  ) -> demography_climate_distance_surv

## Separate each variable to use the same model stan
### Cumulative precipitation
demography_surv_ppt <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop= demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$site_species_plot,
  clim = as.vector(demography_climate_distance_surv$ppt),
  endo = demography_climate_distance_surv$Endo -1 ,
  herb = demography_climate_distance_surv$Herbivory -1,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_ppt <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/survival.stan",
  data = demography_surv_ppt,
  warmup = sim_pars$warmup,
  seed = 13,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)


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


### Distance from niche centroid
demography_surv_distance <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop = demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$site_species_plot,
  clim = as.vector(demography_climate_distance_surv$distance),
  endo = demography_climate_distance_surv$Endo-1,
  herb = demography_climate_distance_surv$Herbivory-1,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/survival_distance.stan",
  data = demography_surv_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

summary(fit_surv_distance)$summary[, c("Rhat", "n_eff")]
posterior_surv_distance <- as.array(fit_surv_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_surv_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()


### Distance from geographic centroid
demography_surv_geo_distance <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop = demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$site_species_plot,
  clim = as.vector(demography_climate_distance_surv$geo_distance),
  endo = demography_climate_distance_surv$Endo-1,
  herb = demography_climate_distance_surv$Herbivory-1,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_geo_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/survival_distance.stan",
  data = demography_surv_geo_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)


summary(fit_surv_geo_distance)$summary[, c("Rhat", "n_eff")]
posterior_surv_geo_distance <- as.array(fit_surv_geo_distance) 
bayesplot::mcmc_trace(posterior_surv_geo_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

## Save RDS file for further use
saveRDS(fit_surv_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_surv_ppt.rds')
saveRDS(fit_surv_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_surv_distance.rds')
saveRDS(fit_surv_geo_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_surv_geo_distance.rds')


# Growth----
## Read and format survival data to build the model
demography_climate_distance %>%
  subset(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, grow, sum_ppt, mean_temp, mean_vpd, distance,geo_distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    grow = grow,
    ppt = log(sum_ppt),
    temp = log(mean_temp),
    vpd = log(mean_vpd),
    distance = log(distance),
    geo_distance=log(geo_distance)
  ) -> demography_climate_distance_grow

## Separate each variable to use the same model stan
### Precipitation
demography_grow_ppt <- list(
  nSpp = demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$site_species_plot,
  clim = as.vector(demography_climate_distance_grow$ppt),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_ppt <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/growth.stan",
  data = demography_grow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

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


### Distance fro niche centroid
demography_grow_distance <- list(
  nSpp= demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$site_species_plot,
  clim = as.vector(demography_climate_distance_grow$distance),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/growth_distance.stan",
  data = demography_grow_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)


summary(fit_grow_distance)$summary[, c("Rhat", "n_eff")]
posterior_grow_distance <- as.array(fit_grow_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

### Distance from geographic center
demography_grow_geo_distance <- list(
  nSpp= demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$site_species_plot,
  clim = as.vector(demography_climate_distance_grow$geo_distance),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_geo_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/growth_distance.stan",
  data = demography_grow_geo_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)


summary(fit_grow_geo_distance)$summary[, c("Rhat", "n_eff")]
posterior_grow_geo_distance <- as.array(fit_grow_geo_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_grow_geo_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

## Save RDS file for further use
# saveRDS(fit_grow_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_grow_ppt.rds')
# saveRDS(fit_grow_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_grow_distance.rds')
# saveRDS(fit_grow_geo_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_grow_geo_distance.rds')

# Flowering----
demography_climate_distance %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, sum_ppt, mean_temp, mean_vpd, distance,geo_distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    flow_t1 = inf_t1,
    ppt = log(sum_ppt),
    temp = log(mean_temp),
    vpd = log(mean_vpd),
    distance = log(distance),
    geo_distance=log(geo_distance)
  ) -> demography_climate_distance_flow

## Separate each variable to use the same model stan
### Precipitation 
demography_flow_ppt <- list(
  nSpp = demography_climate_distance_flow$Species %>% n_distinct(),
  nSite = demography_climate_distance_flow$Site %>% n_distinct(),
  nPop = demography_climate_distance_flow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_flow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_flow$Species,
  site = demography_climate_distance_flow$Site,
  pop = demography_climate_distance_flow$Population,
  plot = demography_climate_distance_flow$site_species_plot,
  clim = as.vector(demography_climate_distance_flow$ppt),
  endo = demography_climate_distance_flow$Endo,
  herb = demography_climate_distance_flow$Herbivory,
  size = demography_climate_distance_flow$log_size_t0,
  y = demography_climate_distance_flow$flow_t1,
  N = nrow(demography_climate_distance_flow)
)

fit_flow_ppt <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering.stan",
  data = demography_flow_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

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


demography_flow_distance <- list(
  nSpp = demography_climate_distance_flow$Species %>% n_distinct(),
  nSite = demography_climate_distance_flow$Site %>% n_distinct(),
  nPop = demography_climate_distance_flow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_flow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_flow$Species,
  site = demography_climate_distance_flow$Site,
  pop = demography_climate_distance_flow$Population,
  plot = demography_climate_distance_flow$site_species_plot,
  clim = as.vector(demography_climate_distance_flow$distance),
  endo = demography_climate_distance_flow$Endo,
  herb = demography_climate_distance_flow$Herbivory,
  size = demography_climate_distance_flow$log_size_t0,
  y = demography_climate_distance_flow$flow_t1,
  N = nrow(demography_climate_distance_flow)
)

fit_flow_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering_distance.stan",
  data = demography_flow_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

summary(fit_flow_distance)$summary[, c("Rhat", "n_eff")]
posterior_flow_distance <- as.array(fit_flow_distance) # Converts to an array
bayesplot::mcmc_trace(fit_flow_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

demography_flow_geo_distance <- list(
  nSpp = demography_climate_distance_flow$Species %>% n_distinct(),
  nSite = demography_climate_distance_flow$Site %>% n_distinct(),
  nPop = demography_climate_distance_flow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_flow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_flow$Species,
  site = demography_climate_distance_flow$Site,
  pop = demography_climate_distance_flow$Population,
  plot = demography_climate_distance_flow$site_species_plot,
  clim = as.vector(demography_climate_distance_flow$geo_distance),
  endo = demography_climate_distance_flow$Endo,
  herb = demography_climate_distance_flow$Herbivory,
  size = demography_climate_distance_flow$log_size_t0,
  y = demography_climate_distance_flow$flow_t1,
  N = nrow(demography_climate_distance_flow)
)

fit_flow_geo_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/flowering_distance.stan",
  data = demography_flow_geo_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control)

summary(fit_flow_geo_distance)$summary[, c("Rhat", "n_eff")]
posterior_flow_geo_distance <- as.array(fit_flow_geo_distance) # Converts to an array
bayesplot::mcmc_trace(fit_flow_geo_distance,
                      pars = quote_bare(
                        b0[1], b0[2], b0[3],
                        bendo[1], bendo[2], bendo[3],
                        bherb[1], bherb[2], bherb[3],
                        bclim[1], bclim[2], bclim[3],
                        bendoclim[1], bendoclim[2], bendoclim[3],
                        bendoherb[1], bendoherb[2], bendoherb[3]
                      )
) + theme_bw()

## Save RDS file for further use
# saveRDS(fit_flow_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_flow_ppt.rds')
# saveRDS(fit_flow_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_flow_distance.rds')
# saveRDS(fit_flow_geo_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_flow_geo_distance.rds')

# Spikelet----
demography_climate_distance %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, sum_ppt, mean_temp, mean_vpd, distance,geo_distance
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(factor(Endo)) - 1,
    Herbivory = as.integer(factor(Herbivory)) - 1
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    spi_t1 = spikelet_t1,
    ppt = log(sum_ppt),
    temp = log(mean_temp),
    vpd = log(mean_vpd),
    distance = log(distance),
    geo_distance=log(geo_distance)
  ) -> demography_climate_distance_spik

### Precipitation
demography_spik_ppt <- list(
  nSpp = demography_climate_distance_spik$Species %>% n_distinct(),
  nSite = demography_climate_distance_spik$Site %>% n_distinct(),
  nPop = demography_climate_distance_spik$Population %>% n_distinct(),
  nPlot = demography_climate_distance_spik$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_spik$Species,
  site = demography_climate_distance_spik$Site,
  pop = demography_climate_distance_spik$Population,
  plot = demography_climate_distance_spik$site_species_plot,
  clim = as.vector(demography_climate_distance_spik$ppt),
  endo = demography_climate_distance_spik$Endo,
  herb = demography_climate_distance_spik$Herbivory,
  size = demography_climate_distance_spik$log_size_t0,
  y = demography_climate_distance_spik$spikelet_t1,
  N = nrow(demography_climate_distance_spik)
)

fit_spik_ppt <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/spikelet.stan",
  data = demography_spik_ppt,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control =sim_pars$control)

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

### Distance from niche centroid
demography_spik_distance <- list(
  nSpp = demography_climate_distance_spik$Species %>% n_distinct(),
  nSite = demography_climate_distance_spik$Site %>% n_distinct(),
  nPop = demography_climate_distance_spik$Population %>% n_distinct(),
  nPlot = demography_climate_distance_spik$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_spik$Species,
  site = demography_climate_distance_spik$Site,
  pop = demography_climate_distance_spik$Population,
  plot = demography_climate_distance_spik$site_species_plot,
  clim = as.vector(demography_climate_distance_spik$distance),
  endo = demography_climate_distance_spik$Endo,
  herb = demography_climate_distance_spik$Herbivory,
  size = demography_climate_distance_spik$log_size_t0,
  y = demography_climate_distance_spik$spikelet_t1,
  N = nrow(demography_climate_distance_spik)
)

fit_spik_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/Spikelet_distance.stan",
  data = demography_spik_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control
)

summary(fit_spik_distance)$summary[, c("Rhat", "n_eff")]
posterior_spik_distance <- as.array(fit_spik_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_distance,
                      pars = quote_bare(
                        b0[1], b0[2],
                        bendo[1], bendo[2],
                        bherb[1], bherb[2],
                        bclim[1], bclim[2],
                        bendoclim[1], bendoclim[2],
                        bendoherb[1], bendoherb[2]
                      )
) + theme_bw()

### Distance from geographic center
demography_spik_geo_distance <- list(
  nSpp = demography_climate_distance_spik$Species %>% n_distinct(),
  nSite = demography_climate_distance_spik$Site %>% n_distinct(),
  nPop = demography_climate_distance_spik$Population %>% n_distinct(),
  nPlot = demography_climate_distance_spik$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_spik$Species,
  site = demography_climate_distance_spik$Site,
  pop = demography_climate_distance_spik$Population,
  plot = demography_climate_distance_spik$site_species_plot,
  clim = as.vector(demography_climate_distance_spik$geo_distance),
  endo = demography_climate_distance_spik$Endo,
  herb = demography_climate_distance_spik$Herbivory,
  size = demography_climate_distance_spik$log_size_t0,
  y = demography_climate_distance_spik$spikelet_t1,
  N = nrow(demography_climate_distance_spik)
)

fit_spik_geo_distance <- stan(
  file = "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/stan/Spikelet_distance.stan",
  data = demography_spik_geo_distance,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = sim_pars$control
)

summary(fit_spik_geo_distance)$summary[, c("Rhat", "n_eff")]
posterior_spik_geo_distance <- as.array(fit_spik_geo_distance) # Converts to an array
bayesplot::mcmc_trace(posterior_spik_geo_distance,
                      pars = quote_bare(
                        b0[1], b0[2],
                        bendo[1], bendo[2],
                        bherb[1], bherb[2],
                        bclim[1], bclim[2],
                        bendoclim[1], bendoclim[2],
                        bendoherb[1], bendoherb[2]
                      )
) + theme_bw()


## Save RDS file for further use
saveRDS(fit_spik_ppt, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_spik_ppt.rds')
saveRDS(fit_spik_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_spik_distance.rds')
saveRDS(fit_spik_geo_distance, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/range limits model output/fit_spik_geo_distance.rds')
