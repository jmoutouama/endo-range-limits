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

fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/hi11gxhpqlrdfg389ir0w/fit_surv_ppt.rds?rlkey=22ujyjnm74c6pw9biw50uasvu&dl=1"))

# Create a new data frame for generating predictions
climate_range <- seq(min(demography_surv_ppt$clim), max(demography_surv_ppt$clim), length.out = 30)
endo_status <- c(0, 1) # Endophyte negative and positive
herb_status <- c(0, 1) # Herbivory no and yes
species <- 1:3 # Species 1, 2, 3

# Create a data frame with all combinations
predictions <- expand.grid(clim = climate_range, endo = endo_status, herb = herb_status, species = species)

# Extract posterior samples
posterior_samples <- rstan::extract(fit_surv_ppt)

# Function to calculate predictions based on the posterior samples
get_predictions <- function(clim, endo, herb, species_index, posterior_samples) {
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
    bclim2 * clim^2 +
    bendoclim2 * endo * clim^2
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
pred_probs_long_df <- gather(pred_probs_df, key = "Posterior_Sample", value = "Pred_Survival", -clim, -endo, -herb, -species)

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
  clim = demography_surv_ppt$clim, # Your climate data
  endo = demography_surv_ppt$endo, # Your endophyte status data
  herb = demography_surv_ppt$herb, # Your herbivory status data
  species = demography_surv_ppt$Spp, # Your species data
  y = demography_surv_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurival_ppt.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals, aes(x = exp(clim), y = mean, color = factor(endo))) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) + # Credible interval
  geom_point(data = observed_data, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) + # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"), herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Precipitation (mm)",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()


### Distance from niche centroid
demography_surv_distance <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop = demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$Plot,
  clim = as.vector(demography_climate_distance_surv$distance),
  endo = demography_climate_distance_surv$Endo,
  herb = demography_climate_distance_surv$Herbivory,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_distance<-readRDS(url("https://www.dropbox.com/scl/fi/gxc8edjzdjvsrtlb8zm5o/fit_surv_distance.rds?rlkey=bmtq9q0bxf6ooafq8pz3tyavp&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_surv_distance$clim), max(demography_surv_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_surv_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival (logit scale)
  logit_preds <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probs <- 1 / (1 + exp(-logit_preds))
  return(pred_probs)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probs_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probs_distance <- as.data.frame(pred_probs_matrix_distance)
colnames(pred_probs_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_distance <- cbind(predictions_distance, pred_probs_distance)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_distance <- gather(pred_probs_distance, key = "Posterior_Sample", value = "Pred_Survival", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals_distance <- pred_probs_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    lower_95 = quantile(Pred_Survival, 0.025),
    upper_95 = quantile(Pred_Survival, 0.975),
    median = quantile(Pred_Survival, 0.5),
    mean = mean(Pred_Survival)  # Calculate the mean survival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_distance <- data.frame(
  clim = demography_surv_distance$clim,  #  climate data
  endo = demography_surv_distance$endo,  #  endophyte status data
  herb = demography_surv_distance$herb,  #  herbivory status data
  species = demography_surv_distance$Spp,  #  species data
  y = demography_surv_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurival_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Mahalanobis distance",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

### Distance from geographic center
demography_surv_geo_distance <- list(
  nSpp = demography_climate_distance_surv$Species %>% n_distinct(),
  nSite = demography_climate_distance_surv$Site %>% n_distinct(),
  nPop = demography_climate_distance_surv$Population %>% n_distinct(),
  nPlot = demography_climate_distance_surv$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_surv$Species,
  site = demography_climate_distance_surv$Site,
  pop = demography_climate_distance_surv$Population,
  plot = demography_climate_distance_surv$Plot,
  clim = as.vector(demography_climate_distance_surv$geo_distance),
  endo = demography_climate_distance_surv$Endo,
  herb = demography_climate_distance_surv$Herbivory,
  size = demography_climate_distance_surv$log_size_t0,
  y = demography_climate_distance_surv$surv_t1,
  N = nrow(demography_climate_distance_surv)
)

fit_surv_geo_distance<-readRDS(url("https://www.dropbox.com/scl/fi/ah7h4cn9l9p6gwl3cl0aw/fit_surv_geo_distance.rds?rlkey=dypxz2231pwle7p8mb534uqwr&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_surv_geo_distance$clim), max(demography_surv_geo_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_surv_geo_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival (logit scale)
  logit_preds <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probs <- 1 / (1 + exp(-logit_preds))
  return(pred_probs)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probs_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probs_distance <- as.data.frame(pred_probs_matrix_distance)
colnames(pred_probs_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_distance <- cbind(predictions_distance, pred_probs_distance)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_distance <- gather(pred_probs_distance, key = "Posterior_Sample", value = "Pred_Survival", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals_distance <- pred_probs_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    lower_95 = quantile(Pred_Survival, 0.025),
    upper_95 = quantile(Pred_Survival, 0.975),
    median = quantile(Pred_Survival, 0.5),
    mean = mean(Pred_Survival)  # Calculate the mean survival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_geo_distance <- data.frame(
  clim = demography_surv_geo_distance$clim,  #  climate data
  endo = demography_surv_geo_distance$endo,  #  endophyte status data
  herb = demography_surv_geo_distance$herb,  #  herbivory status data
  species = demography_surv_geo_distance$Spp,  #  species data
  y = demography_surv_geo_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/PrSurival_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_geo_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Distance from geographc center",
    y = "Predicted survival probability",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.8),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()


# Growth----
## Read and format survival data to build the model
demography_climate_distance %>%
  subset(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, grow, sum_ppt, mean_pet, mean_spei, distance,geo_distance
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
    pet = log(mean_pet),
    spei = mean_spei,
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

fit_grow_ppt<-readRDS(url("https://www.dropbox.com/scl/fi/85pzzrgvkxwogoybr004t/fit_grow_ppt.rds?rlkey=rz2xlg00u1aqhkxeaix7wiu9e&dl=1"))

posterior_samples <- rstan::extract(fit_grow_ppt)
predictions <- expand.grid(clim = seq(min(demography_grow_ppt$clim), max(demography_grow_ppt$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Function to calculate predictions based on the posterior samples
get_predictions_grow <- function(clim, endo, herb, species_index, posterior_samples) {
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
    bclim2 * clim^2 +
    bendoclim2 * endo * clim^2
  # Keep predg
  pred_probg <- predg
  return(pred_probg)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0)  # Number of posterior samples 
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
pred_probg_long_df <- gather(pred_probg_df, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalg <- pred_probg_long_df %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_grow <- data.frame(
  clim = demography_grow_ppt$clim,  # Your climate data
  endo = demography_grow_ppt$endo,  # Your endophyte status data
  herb = demography_grow_ppt$herb,  # Your herbivory status data
  species = demography_grow_ppt$Spp,  # Your species data
  y = demography_grow_ppt$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Growth_ppt.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervalg, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_grow, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Precipitation (mm)",
    y = "Predicted realtive growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.92, 0.085),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

### SPEI
demography_grow_spei <- list(
  nSpp = demography_climate_distance_grow$Species %>% n_distinct(),
  nSite = demography_climate_distance_grow$Site %>% n_distinct(),
  nPop = demography_climate_distance_grow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_grow$Species,
  site = demography_climate_distance_grow$Site,
  pop = demography_climate_distance_grow$Population,
  plot = demography_climate_distance_grow$site_species_plot,
  clim = as.vector(demography_climate_distance_grow$spei),
  endo = demography_climate_distance_grow$Endo,
  herb = demography_climate_distance_grow$Herbivory,
  size = demography_climate_distance_grow$log_size_t0,
  y = demography_climate_distance_grow$grow,
  N = nrow(demography_climate_distance_grow)
)

fit_grow_spei <- readRDS(url("https://www.dropbox.com/scl/fi/4iuz5ay461qkjv732b5yb/fit_grow_spei.rds?rlkey=gq68ixyetb3v4o3ds7uo24cwk&dl=1"))

# Create a data frame with all combinations
predictions_spei <- expand.grid(clim = seq(min(demography_grow_spei$clim), max(demography_grow_spei$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_spei <- rstan::extract(fit_grow_spei)
# Apply the function to generate predictions for all combinations
n_posterior_samples_spei <- length(posterior_samples_spei$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix_spei <- matrix(NA, nrow = nrow(predictions_spei), ncol = n_posterior_samples_spei)
# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_spei)) {
  pred_probg_matrix_spei[i, ] <- get_predictions_grow(
    predictions_spei$clim[i], 
    predictions_spei$endo[i], 
    predictions_spei$herb[i], 
    predictions_spei$species[i], 
    posterior_samples_spei
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probg_spei <- as.data.frame(pred_probg_matrix_spei)
colnames(pred_probg_spei) <- paste("Posterior_Sample", 1:n_posterior_samples_spei)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_spei <- cbind(predictions_spei, pred_probg_spei)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_spei <- gather(pred_probg_spei, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growival probability
cred_intervalg_spei <- pred_probg_long_spei %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean growival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed growival)
observed_spei <- data.frame(
  clim = demography_grow_spei$clim,  #  climate data
  endo = demography_grow_spei$endo,  #  endophyte status data
  herb = demography_grow_spei$herb,  #  herbivory status data
  species = demography_grow_spei$Spp,  #  species data
  y = demography_grow_spei$y  # Observed growival
)

# Plot the results with credible intervals, mean growival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Growth_spei.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervalg_spei, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median growival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean growival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_spei, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "SPEI",
    y = "Predicted relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.78, 0.945),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

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

fit_grow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/fhwjeizspvvd2dbz11195/fit_grow_distance.rds?rlkey=yt863sjawv1zliwem6a9ved1h&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_grow_distance$clim), max(demography_grow_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_grow_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival
  predg <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probg <-predg
  return(pred_probg)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probg_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probg_distance <- as.data.frame(pred_probg_matrix_distance)
colnames(pred_probg_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_distance <- cbind(predictions_distance, pred_probg_distance)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_distance <- gather(pred_probg_distance, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growth
cred_intervals_distance <- pred_probg_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean survival probability
  ) %>%
  ungroup()

sum(is.na(cred_intervals_distance))
sum(is.na(observed_distance))
# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_grow_distance$clim,  #  climate data
  endo = demography_grow_distance$endo,  #  endophyte status data
  herb = demography_grow_distance$herb,  #  herbivory status data
  species = demography_grow_distance$Spp,  #  species data
  y = demography_grow_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Growth_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Mahalanobis distance",
    y = "Predicted relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.23),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()


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

fit_grow_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/m3r1hjc7sv3xvuvue79tz/fit_grow_geo_distance.rds?rlkey=5rj6htaf8fovv2ve8qz0f6hvg&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_grow_geo_distance$clim), max(demography_grow_geo_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_grow_geo_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival
  predg <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probg <-predg
  return(pred_probg)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probg_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probg_distance <- as.data.frame(pred_probg_matrix_distance)
colnames(pred_probg_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_distance <- cbind(predictions_distance, pred_probg_distance)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_distance <- gather(pred_probg_distance, key = "Posterior_Sample", value = "Pred_Growth", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growth
cred_intervals_distance <- pred_probg_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth)  # Calculate the mean survival probability
  ) %>%
  ungroup()

sum(is.na(cred_intervals_distance))
sum(is.na(observed_distance))
# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_grow_geo_distance$clim,  #  climate data
  endo = demography_grow_geo_distance$endo,  #  endophyte status data
  herb = demography_grow_geo_distance$herb,  #  herbivory status data
  species = demography_grow_geo_distance$Spp,  #  species data
  y = demography_grow_geo_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Growth_geo_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Distance from geographic center(log scale) ",
    y = "Predicted relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.23),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()



# Flowering----
demography_climate_distance %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, sum_ppt, mean_pet, mean_spei, distance,geo_distance
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
    pet = log(mean_pet),
    spei = mean_spei,
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

table(demography_climate_distance_flow$Species)

fit_flow_ppt<-readRDS(url("https://www.dropbox.com/scl/fi/ajo9euxdtxfmle2g006l3/fit_flow_ppt.rds?rlkey=78cknj1yaqs2m5cyicdzg8up4&dl=1"))


posterior_samples <- rstan::extract(fit_flow_ppt)

predictions <- expand.grid(clim = seq(min(demography_flow_ppt$clim), max(demography_flow_ppt$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Function to calculate predictions based on the posterior samples
get_predictions_flow <- function(clim, endo, herb, species_index, posterior_samples) {
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
    bclim2 * clim^2 +
    bendoclim2 * endo * clim^2
  #  predf
  pred_probf <- exp(predf)
  return(pred_probf)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0)  # Number of posterior samples 
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
pred_probf_long_df <- gather(pred_probf_df, key = "Posterior_Sample", value = "Pred_Flow", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalf <- pred_probf_long_df %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    lower_95 = quantile(Pred_Flow, 0.025),
    upper_95 = quantile(Pred_Flow, 0.975),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow)  # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_flow <- data.frame(
  clim = demography_flow_ppt$clim,  # Your climate data
  endo = demography_flow_ppt$endo,  # Your endophyte status data
  herb = demography_flow_ppt$herb,  # Your herbivory status data
  species = demography_flow_ppt$Spp,  # Your species data
  y = demography_flow_ppt$y  # Observed survival
)

unique(observed_flow$species)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Flowth_ppt.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervalf, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_flow, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Precipitation (mm)",
    y = "Inflorescences",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.88, 0.2),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

### SPEI
demography_flow_spei <- list(
  nSpp = demography_climate_distance_flow$Species %>% n_distinct(),
  nSite = demography_climate_distance_flow$Site %>% n_distinct(),
  nPop = demography_climate_distance_flow$Population %>% n_distinct(),
  nPlot = demography_climate_distance_flow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_flow$Species,
  site = demography_climate_distance_flow$Site,
  pop = demography_climate_distance_flow$Population,
  plot = demography_climate_distance_flow$site_species_plot,
  clim = as.vector(demography_climate_distance_flow$spei),
  endo = demography_climate_distance_flow$Endo,
  herb = demography_climate_distance_flow$Herbivory,
  size = demography_climate_distance_flow$log_size_t0,
  y = demography_climate_distance_flow$flow_t1,
  N = nrow(demography_climate_distance_flow)
)

fit_flow_spei <- readRDS(url("https://www.dropbox.com/scl/fi/p4t60pr13u5qtt6hnz8kv/fit_flow_spei.rds?rlkey=5on6zasr6c1qt9wzaqux7lk0d&dl=1"))

# Create a data frame with all combinations
predictions_spei <- expand.grid(clim = seq(min(demography_flow_spei$clim), max(demography_flow_spei$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_spei <- rstan::extract(fit_grow_spei)
# Apply the function to generate predictions for all combinations
n_posterior_samples_spei <- length(posterior_samples_spei$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probf_matrix_spei <- matrix(NA, nrow = nrow(predictions_spei), ncol = n_posterior_samples_spei)
# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_spei)) {
  pred_probf_matrix_spei[i, ] <- get_predictions_flow(
    predictions_spei$clim[i], 
    predictions_spei$endo[i], 
    predictions_spei$herb[i], 
    predictions_spei$species[i], 
    posterior_samples_spei
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probf_spei <- as.data.frame(pred_probf_matrix_spei)
colnames(pred_probf_spei) <- paste("Posterior_Sample", 1:n_posterior_samples_spei)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probf_spei <- cbind(predictions_spei, pred_probf_spei)

# Reshape the data frame so we have long format for ggplot
pred_probf_long_spei <- gather(pred_probf_spei, key = "Posterior_Sample", value = "Pred_Flow", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growival probability
cred_intervalf_spei <- pred_probf_long_spei %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    lower_95 = quantile(Pred_Flow, 0.025),
    upper_95 = quantile(Pred_Flow, 0.975),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow)  # Calculate the mean 
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed growival)
observed_spei <- data.frame(
  clim = demography_flow_spei$clim,  #  climate data
  endo = demography_flow_spei$endo,  #  endophyte status data
  herb = demography_flow_spei$herb,  #  herbivory status data
  species = demography_flow_spei$Spp,  #  species data
  y = demography_flow_spei$y  # Observed growival
)

# Plot the results with credible intervals, mean growival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Flow_spei.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervalf_spei, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median growival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean growival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_spei, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "SPEI",
    y = "Predicted relative growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.90, 0.945),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

### Distance from niche centroid

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

fit_flow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/y7qvz4pmy2t2j00gqvwqd/fit_flow_distance.rds?rlkey=qe47fg37dpi4z2lu1ur1lsbkh&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_flow_distance$clim), max(demography_flow_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_flow_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probf_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival
  predf <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probf <-exp(predf)
  return(pred_probf)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probf_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probf_distance <- as.data.frame(pred_probf_matrix_distance)
colnames(pred_probf_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probf_distance <- cbind(predictions_distance, pred_probf_distance)

# Reshape the data frame so we have long format for ggplot
pred_probf_long_distance <- gather(pred_probf_distance, key = "Posterior_Sample", value = "Pred_Flow", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growth
cred_intervals_distance <- pred_probf_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    lower_95 = quantile(Pred_Flow, 0.025),
    upper_95 = quantile(Pred_Flow, 0.975),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow)  
  ) %>%
  ungroup()

sum(is.na(cred_intervals_distance))
sum(is.na(observed_distance))
# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_flow_distance$clim,  #  climate data
  endo = demography_flow_distance$endo,  #  endophyte status data
  herb = demography_flow_distance$herb,  #  herbivory status data
  species = demography_flow_distance$Spp,  #  species data
  y = demography_flow_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Flow_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Mahalanobis distance",
    y = "Inflorescences",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.23),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()


### Distance from geographic center
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

fit_flow_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/b2dez0gz22p3m2ke1zs7f/fit_flow_geo_distance.rds?rlkey=9fuwynap20617l1rqffjz62wd&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_flow_geo_distance$clim), max(demography_flow_geo_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:3 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_flow_geo_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probf_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival
  predf <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probf <-exp(predf)
  return(pred_probf)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probf_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probf_distance <- as.data.frame(pred_probf_matrix_distance)
colnames(pred_probf_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probf_distance <- cbind(predictions_distance, pred_probf_distance)

# Reshape the data frame so we have long format for ggplot
pred_probf_long_distance <- gather(pred_probf_distance, key = "Posterior_Sample", value = "Pred_Flow", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growth
cred_intervals_distance <- pred_probf_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    lower_95 = quantile(Pred_Flow, 0.025),
    upper_95 = quantile(Pred_Flow, 0.975),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow)  
  ) %>%
  ungroup()

sum(is.na(cred_intervals_distance))
sum(is.na(observed_distance))
# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_flow_geo_distance$clim,  #  climate data
  endo = demography_flow_geo_distance$endo,  #  endophyte status data
  herb = demography_flow_geo_distance$herb,  #  herbivory status data
  species = demography_flow_geo_distance$Spp,  #  species data
  y = demography_flow_geo_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Flow_geo_distance.pdf", useDingbats = F, height = 9, width = 7)
ggplot(cred_intervals_distance, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Distance from geographic center(log scale)",
    y = "Inflorescences",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.23),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()



# Spikelet----
demography_climate_distance %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  subset(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, sum_ppt, mean_pet, mean_spei, distance
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
    pet = log(mean_pet),
    spei = mean_spei,
    distance = log(distance)
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

fit_spik_ppt<-readRDS(url("https://www.dropbox.com/scl/fi/9klrwl866stup2x3m2p9x/fit_spik_ppt.rds?rlkey=okno4gvn9dj5q8uhpnu36f59i&dl=1"))

posterior_samples <- rstan::extract(fit_spik_ppt)

predictions <- expand.grid(clim = seq(min(demography_spik_ppt$clim), max(demography_spik_ppt$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:2 )
# Function to calculate predictions based on the posterior samples
get_predictions_spk <- function(clim, endo, herb, species_index, posterior_samples) {
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
    bclim2 * clim^2 +
    bendoclim2 * endo * clim^2
  #  predf
  pred_probspk <- exp(predspk)
  return(pred_probspk)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0)  # Number of posterior samples 
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
pred_probspk_long_df <- gather(pred_probspk_df, key = "Posterior_Sample", value = "Pred_Spik", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalspk <- pred_probspk_long_df %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Spik, 0.05),
    upper_90 = quantile(Pred_Spik, 0.95),
    lower_95 = quantile(Pred_Spik, 0.025),
    upper_95 = quantile(Pred_Spik, 0.975),
    median = quantile(Pred_Spik, 0.5),
    mean = mean(Pred_Spik)  # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_spk<- data.frame(
  clim = demography_spik_ppt$clim,  # Your climate data
  endo = demography_spik_ppt$endo,  # Your endophyte status data
  herb = demography_spik_ppt$herb,  # Your herbivory status data
  species = demography_spik_ppt$Spp,  # Your species data
  y = demography_spik_ppt$y  # Observed survival
)

unique(observed_spk$species)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Spik_ppt.pdf", useDingbats = F, height = 7, width = 7)
ggplot(cred_intervalspk, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_spk, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "ELVI", "2" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Precipitation (mm)",
    y = "Spikelets",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.88, 0.2),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

### SPEI
demography_spik_spei <- list(
  nSpp = demography_climate_distance_spik$Species %>% n_distinct(),
  nSite = demography_climate_distance_spik$Site %>% n_distinct(),
  nPop = demography_climate_distance_spik$Population %>% n_distinct(),
  nPlot = demography_climate_distance_spik$site_species_plot %>% n_distinct(),
  Spp = demography_climate_distance_spik$Species,
  site = demography_climate_distance_spik$Site,
  pop = demography_climate_distance_spik$Population,
  plot = demography_climate_distance_spik$site_species_plot,
  clim = as.vector(demography_climate_distance_spik$spei),
  endo = demography_climate_distance_spik$Endo,
  herb = demography_climate_distance_spik$Herbivory,
  size = demography_climate_distance_spik$log_size_t0,
  y = demography_climate_distance_spik$spikelet_t1,
  N = nrow(demography_climate_distance_spik)
)

fit_spik_spei <- readRDS(url("https://www.dropbox.com/scl/fi/uj8ok9u0hp2uw6gnpqhv3/fit_spik_spei.rds?rlkey=wfbjdxnjjx2qetc7cg8rinbg2&dl=1"))

# Create a data frame with all combinations
predictions_spei <- expand.grid(clim = seq(min(demography_spik_spei$clim), max(demography_spik_spei$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:2 )
# Extract posterior samples
posterior_samples_spei <- rstan::extract(fit_spik_spei)
# Apply the function to generate predictions for all combinations
n_posterior_samples_spei <- length(posterior_samples_spei$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probspk_matrix_spei <- matrix(NA, nrow = nrow(predictions_spei), ncol = n_posterior_samples_spei)
# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_spei)) {
  pred_probspk_matrix_spei[i, ] <- get_predictions_spk(
    predictions_spei$clim[i], 
    predictions_spei$endo[i], 
    predictions_spei$herb[i], 
    predictions_spei$species[i], 
    posterior_samples_spei
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probspk_spei <- as.data.frame(pred_probspk_matrix_spei)
colnames(pred_probspk_spei) <- paste("Posterior_Sample", 1:n_posterior_samples_spei)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probspk_spei <- cbind(predictions_spei, pred_probspk_spei)

# Reshape the data frame so we have long format for ggplot
pred_probspk_long_spei <- gather(pred_probspk_spei, key = "Posterior_Sample", value = "Pred_SPK", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growival probability
cred_intervalspk_spei <- pred_probspk_long_spei %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_SPK, 0.05),
    upper_90 = quantile(Pred_SPK, 0.95),
    lower_95 = quantile(Pred_SPK, 0.025),
    upper_95 = quantile(Pred_SPK, 0.975),
    median = quantile(Pred_SPK, 0.5),
    mean = mean(Pred_SPK)  # Calculate the mean 
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed growival)
observed_spei <- data.frame(
  clim = demography_spik_spei$clim,  #  climate data
  endo = demography_spik_spei$endo,  #  endophyte status data
  herb = demography_spik_spei$herb,  #  herbivory status data
  species = demography_spik_spei$Spp,  #  species data
  y = demography_spik_spei$y  # Observed growival
)

# Plot the results with credible intervals, mean growival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Spik_spei.pdf", useDingbats = F, height = 7, width = 7)
ggplot(cred_intervalspk_spei, aes(x = clim, y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median growival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean growival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_spei, aes(x = clim, y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "ELVI", "2" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "SPEI",
    y = "Spikelets",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.90, 0.945),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

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

fit_spik_distance <- readRDS(url("https://www.dropbox.com/scl/fi/1m01hqgokktxz77trnuu3/fit_spik_distance.rds?rlkey=jxx4nt8qgwgiospq870anu9mn&dl=1"))

predictions_distance <- expand.grid(clim = seq(min(demography_spik_distance$clim), max(demography_spik_distance$clim), length.out = 30), endo = c(0, 1), herb = c(0, 1) , species =  1:2 )
# Extract posterior samples
posterior_samples_distance <- rstan::extract(fit_spik_distance)
# Apply the function to generate predictions for all combinations
n_posterior_samples_distance <- length(posterior_samples_distance$b0)  # Number of posterior samples 
# Initialize a matrix to hold predictions for each posterior sample
pred_probspk_matrix_distance <- matrix(NA, nrow = nrow(predictions_distance), ncol = n_posterior_samples_distance)
# Function to calculate predictions based on the posterior samples
get_predictions_distance <- function(clim, endo, herb, species_index, posterior_samples_distance) {
  b0 <- posterior_samples_distance$b0[, species_index]
  bendo <- posterior_samples_distance$bendo[, species_index]
  bherb <- posterior_samples_distance$bherb[, species_index]
  bclim <- posterior_samples_distance$bclim[, species_index]
  bendoclim <- posterior_samples_distance$bendoclim[, species_index]
  bendoherb <- posterior_samples_distance$bendoherb[, species_index]
  # Predicted survival
  predspk <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb 
  # Convert logit to probability using logistic function
  pred_probspk <-exp(predspk)
  return(pred_probspk)
}

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions_distance)) {
  pred_probspk_matrix_distance[i, ] <- get_predictions_distance(
    predictions_distance$clim[i], 
    predictions_distance$endo[i], 
    predictions_distance$herb[i], 
    predictions_distance$species[i], 
    posterior_samples_distance
  )
}

# Convert the matrix into a data frame with the correct structure
pred_probspk_distance <- as.data.frame(pred_probspk_matrix_distance)
colnames(pred_probspk_distance) <- paste("Posterior_Sample", 1:n_posterior_samples_distance)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probspk_distance <- cbind(predictions_distance, pred_probspk_distance)

# Reshape the data frame so we have long format for ggplot
pred_probspk_long_distance <- gather(pred_probspk_distance, key = "Posterior_Sample", value = "Pred_Spik", -clim, -endo, -herb, -species)

# Calculate credible intervals (90% and 95%) and mean growth
cred_intervals_distance <- pred_probspk_long_distance %>%
  group_by(species, endo, herb,clim) %>%
  summarise(
    lower_90 = quantile(Pred_Spik, 0.05),
    upper_90 = quantile(Pred_Spik, 0.95),
    lower_95 = quantile(Pred_Spik, 0.025),
    upper_95 = quantile(Pred_Spik, 0.975),
    median = quantile(Pred_Spik, 0.5),
    mean = mean(Pred_Spik)  
  ) %>%
  ungroup()

sum(is.na(cred_intervals_distance))
sum(is.na(observed_distance))
# observed_data should have columns: clim, endo, herb, species, y (observed survival)
observed_distance <- data.frame(
  clim = demography_spik_distance$clim,  #  climate data
  endo = demography_spik_distance$endo,  #  endophyte status data
  herb = demography_spik_distance$herb,  #  herbivory status data
  species = demography_spik_distance$Spp,  #  species data
  y = demography_spik_distance$y  # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/ELVI-endophyte-density/Figure/Spike_distance.pdf", useDingbats = F, height = 7, width = 7)
ggplot(cred_intervals_distance, aes(x = exp(clim), y = mean, color = factor(endo))) +
  #geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) +  # Plot the mean survival probability (dashed line)
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(endo)), alpha = 0.3, color = NA) +  # Credible interval
  geom_point(data = observed_distance, aes(x = exp(clim), y = y, color = factor(endo)), size = 3) +  # Observed data points
  facet_grid(species ~ herb, scales = "free_y", labeller = labeller(species = c("1" = "ELVI", "2" = "POAU"),herb = c("0" = "Unfenced", "1" = "Fenced"))) +
  labs(
    x = "Mahalanobis distance",
    y = "Spikelets",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"), labels = c("E-", "E+")) + # Change fill labels
  theme_bw() +
  theme(
    legend.position = c(0.35, 0.23),
    legend.title = element_text(size = 10), # Reduce legend title size
    legend.text = element_text(size = 12), # Adjust legend text size
    axis.title = element_text(size = 13), # Increase axis title size
    axis.text = element_text(size = 10), # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()
