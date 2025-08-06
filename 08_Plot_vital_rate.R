# Project:
# Purpose: Plot  vital rate models (survival, growth, flowering and spikelet) as function of climate or distance from niche center.
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


## Create new variables
dat_t_t1_herb %>%
  mutate(
    site_year=interaction(Site,census_year),
    surv1 = 1 * (!is.na(tiller_t) & !is.na(tiller_t1)),
    site_species_plot = interaction(Site, Species, Plot),
    grow = (log(tiller_t1 + 1) - log(tiller_t + 1))
  ) -> demography_climate

# names(demography_climate)
# view(demography_climate)
# summary(demography_climate)

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

fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/9xa46n5v7u1lddxaj69cs/fit_surv_ppt.rds?rlkey=1mbkby4394s04j7qej4kvzeo2&dl=1"))
# Create a new data frame for generating predictions
climate_range <- seq(min(demography_surv_ppt$clim),
                     max(demography_surv_ppt$clim),
                     length.out = 30)
endo_status <- c(0, 1) # Endophyte negative and positive
herb_status <- c(0, 1) # Herbivory no and yes
species <- 1:3 # Species 1, 2, 3

# Create a data frame with all combinations
predictions <- expand.grid(
  clim = climate_range,
  endo = endo_status,
  herb = herb_status,
  species = species
)

# Extract posterior samples
posterior_samples <- rstan::extract(fit_surv_ppt)

# Function to calculate predictions based on the posterior samples
get_predictions <- function(clim,
                            endo,
                            herb,
                            species_index,
                            posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted survival (logit scale)
  logit_preds <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  # Convert logit to probability using logistic function
  pred_probs <- 1 / (1 + exp(-logit_preds))
  return(pred_probs)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probs_matrix[i, ] <- get_predictions(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_preds <- get_predictions(0.5, 1, 0, 1, posterior_samples)
species_2_preds <- get_predictions(0.2, 0, 1, 2, posterior_samples)
species_3_preds <- get_predictions(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probs_df <- as.data.frame(pred_probs_matrix)
colnames(pred_probs_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_df <- cbind(predictions, pred_probs_df)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_df <- gather(pred_probs_df,
                             key = "Posterior_Sample",
                             value = "Pred_Survival",
                             -clim,
                             -endo,
                             -herb,
                             -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals <- pred_probs_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    lower_95 = quantile(Pred_Survival, 0.025),
    upper_95 = quantile(Pred_Survival, 0.975),
    median = quantile(Pred_Survival, 0.5),
    mean = mean(Pred_Survival) # Calculate the mean survival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_data <- data.frame(
  clim = demography_surv_ppt$clim,
  # Your climate data
  endo = demography_surv_ppt$endo,
  # Your endophyte status data
  herb = demography_surv_ppt$herb,
  # Your herbivory status data
  species = demography_surv_ppt$Spp,
  # Your species data
  y = demography_surv_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/PrSurival_ppt.pdf",
  useDingbats = F,
  height = 9,
  width = 7
)
ggplot(cred_intervals, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(
    data = observed_data,
    aes(
      x = exp(clim),
      y = y,
      color = factor(endo)
    ),
    size = 3,
    position = position_jitter(width = 0, height = 0.02)  # jitter only y-direction
  )+
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.88, 0.085),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

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

fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/d0x30lqqcxnatupsm2hej/fit_grow_ppt.rds?rlkey=er2is1le25trin73an23ztfgm&dl=1"))
posterior_samples <- rstan::extract(fit_grow_ppt)
predictions <- expand.grid(
  clim = seq(
    min(demography_grow_ppt$clim),
    max(demography_grow_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:3
)
# Function to calculate predictions based on the posterior samples
get_predictions_grow <- function(clim,
                                 endo,
                                 herb,
                                 species_index,
                                 posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predg <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  # Keep predg
  pred_probg <- predg
  return(pred_probg)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probg_matrix[i, ] <- get_predictions_grow(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_predg <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predg <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)
species_3_predg <- get_predictions_grow(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probg_df <- as.data.frame(pred_probg_matrix)
colnames(pred_probg_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_df <- cbind(predictions, pred_probg_df)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_df <- gather(pred_probg_df,
                             key = "Posterior_Sample",
                             value = "Pred_Growth",
                             -clim,
                             -endo,
                             -herb,
                             -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalg <- pred_probg_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth) # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_grow <- data.frame(
  clim = demography_grow_ppt$clim,
  # Your climate data
  endo = demography_grow_ppt$endo,
  # Your endophyte status data
  herb = demography_grow_ppt$herb,
  # Your herbivory status data
  species = demography_grow_ppt$Spp,
  # Your species data
  y = demography_grow_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_ppt.pdf",
  useDingbats = F,
  height = 9,
  width = 7
)
ggplot(cred_intervalg, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(data = observed_grow,
             aes(
               x = exp(clim),
               y = y,
               color = factor(endo)
             ),
             size = 3,
             position = position_jitter(width = 0, height = 0.02)  # jitter only y-direction
  ) + # Observed data points
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "Predicted realtive growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.12, 0.23),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

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

fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/z4kh899krlz2ssz1oimxk/fit_flow_ppt.rds?rlkey=uqmm05uas6czzb9siyg7unp1r&dl=1"))
posterior_samples <- rstan::extract(fit_flow_ppt)

predictions <- expand.grid(
  clim = seq(
    min(demography_flow_ppt$clim),
    max(demography_flow_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:3
)
# Function to calculate predictions based on the posterior samples
get_predictions_flow <- function(clim,
                                 endo,
                                 herb,
                                 species_index,
                                 posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predf <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  #  predf
  pred_probf <- exp(predf)
  return(pred_probf)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probf_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probf_matrix[i, ] <- get_predictions_flow(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_predf <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predf <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)
species_3_predf <- get_predictions_grow(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probf_df <- as.data.frame(pred_probf_matrix)
colnames(pred_probf_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probf_df <- cbind(predictions, pred_probf_df)

# Reshape the data frame so we have long format for ggplot
pred_probf_long_df <- gather(pred_probf_df,
                             key = "Posterior_Sample",
                             value = "Pred_Flow",
                             -clim,
                             -endo,
                             -herb,
                             -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalf <- pred_probf_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    lower_95 = quantile(Pred_Flow, 0.025),
    upper_95 = quantile(Pred_Flow, 0.975),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow) # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_flow <- data.frame(
  clim = demography_flow_ppt$clim,
  # Your climate data
  endo = demography_flow_ppt$endo,
  # Your endophyte status data
  herb = demography_flow_ppt$herb,
  # Your herbivory status data
  species = demography_flow_ppt$Spp,
  # Your species data
  y = demography_flow_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flow_ppt.pdf",
  useDingbats = F,
  height = 9,
  width = 7
)
ggplot(cred_intervalf, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(data = observed_flow,
             aes(
               x = exp(clim),
               y = y,
               color = factor(endo)
             ),
             size = 3,
             position = position_jitter(width = 0, height = 0.02)) + # Observed data points
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "# Inflorescences",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.1, 0.2),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

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

fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/1ar22b0jf1urcnypduybm/fit_spik_ppt.rds?rlkey=5hiwrlm6wc2iwvqxma9mx5r15&dl=1"))
posterior_samples <- rstan::extract(fit_spik_ppt)
predictions <- expand.grid(
  clim = seq(
    min(demography_spik_ppt$clim),
    max(demography_spik_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:2
)
# Function to calculate predictions based on the posterior samples
get_predictions_spk <- function(clim,
                                endo,
                                herb,
                                species_index,
                                posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predspk <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  #  predf
  pred_probspk <- exp(predspk)
  return(pred_probspk)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probspk_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probspk_matrix[i, ] <- get_predictions_spk(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_predspk <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predspk <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probspk_df <- as.data.frame(pred_probspk_matrix)
colnames(pred_probspk_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probspk_df <- cbind(predictions, pred_probspk_df)

# Reshape the data frame so we have long format for ggplot
pred_probspk_long_df <- gather(pred_probspk_df,
                               key = "Posterior_Sample",
                               value = "Pred_Spik",
                               -clim,
                               -endo,
                               -herb,
                               -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalspk <- pred_probspk_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Spik, 0.05),
    upper_90 = quantile(Pred_Spik, 0.95),
    lower_95 = quantile(Pred_Spik, 0.025),
    upper_95 = quantile(Pred_Spik, 0.975),
    median = quantile(Pred_Spik, 0.5),
    mean = mean(Pred_Spik) # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_spk <- data.frame(
  clim = demography_spik_ppt$clim,
  # Your climate data
  endo = demography_spik_ppt$endo,
  # Your endophyte status data
  herb = demography_spik_ppt$herb,
  # Your herbivory status data
  species = demography_spik_ppt$Spp,
  # Your species data
  y = demography_spik_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spik_ppt.pdf",
  useDingbats = F,
  height = 7,
  width = 7
)
ggplot(cred_intervalspk, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(data = observed_spk,
             aes(
               x = exp(clim),
               y = y,
               color = factor(endo)
             ),
             size = 3,
             position = position_jitter(width = 0, height = 0.02)) + # Observed data points
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "ELVI", "2" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "# Spikelets",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.1, 0.35),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

