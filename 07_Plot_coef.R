# Project:
# Purpose: Fit vital rate models to test the effect of grass-endophyte symbiosis and endophyte hyphal density on  vital rate models (survival, growth, flowering,fertility).
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

# Surivival models
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/x5q2v6fxnojb9msl0yhnj/fit_surv_abio_bio_endo.rds?rlkey=ojivjcdq1s30jf15is7pbezx8&dl=1"))
## Plot the coefficients
### Precipitation
posterior_samples_surv_ppt <- rstan::extract(fit_surv_ppt)
# Convert to data frame
posterior_samples_surv_ppt_df <- as.data.frame(posterior_samples_surv_ppt)
# Get the number of species
n_species <- length(posterior_samples_surv_ppt$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_ppt_coef_list <- c(
                        "bendoclim",
                        "bendoherb",
                        "bendoherbclim")
surv_ppt_long_data <- list()
for (coef in surv_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_ppt_coef_matrix <- posterior_samples_surv_ppt[[coef]]
  # Convert to long format
  surv_ppt_long_data[[coef]] <- as.data.frame(surv_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_ppt_coef_list'
}
# Combine all into one dataframe
plot_data_surv_ppt <- bind_rows(surv_ppt_long_data)
# Convert species index to numeric
plot_data_surv_ppt$species <- as.numeric(gsub("V", "", plot_data_surv_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_ppt <- plot_data_surv_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
# Convert species numeric IDs to abbreviations
summary_stats_surv_ppt$species <- factor(
  summary_stats_surv_ppt$species,
  levels = c(1, 2, 3),
  labels = c("AGHY", "ELVI", "POAU")
)

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

# Conversion factor
pi_sqrt3 <- pi / sqrt(3)

summary_stats_surv_ppt <- summary_stats_surv_ppt %>%
  mutate(
    d_estimate = mean_estimate * pi_sqrt3,
    d_lower = lower_CI * pi_sqrt3,
    d_upper = upper_CI * pi_sqrt3
  )
# Rename parameters for readability
param_labels <- c(
  "bendoclim" = "Endophyte × Climate",
  "bendoherb" = "Endophyte × Herbivory",
  "bendoherbclim" = "Endophyte × Herbivory × Climate"
)

summary_stats_surv_ppt$parameter_label <- recode(summary_stats_surv_ppt$parameter, !!!param_labels)


# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_surv_coeff.pdf",
  width = 7,
  height = 9
)
ggplot(summary_stats_surv_ppt, 
       aes(x = species, 
           y = d_estimate, 
           fill = parameter_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = d_lower, ymax = d_upper),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = c(
    "Endophyte × Climate" = "#1b9e77",
    "Endophyte × Herbivory" = "#d95f02",
    "Endophyte × Herbivory × Climate" = "#7570b3"
  )) +
  labs(
    x = "Species",
    y = "Standardized effect size",
    fill = "Interaction Term",
    #title = "Relative Strength of Factors Modifying the Endophyte Effect on Survival",
    #subtitle = "Positive = more mutualistic; Negative = more parasitic"
  ) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

dev.off()

## Growth----
fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/6zy5h2kdj56jkzm5ca8ok/fit_grow_abio_bio_endo.rds?rlkey=92hul5ki05hpbgguo9nbg6dd3&dl=1"))
posterior_samples_grow_ppt <- rstan::extract(fit_grow_ppt)

# List of interaction coefficients
grow_ppt_coef_list <- c("bendoclim", "bendoherb", "bendoherbclim")

# Convert coefficients to long format
grow_ppt_long_data <- list()
for (coef in grow_ppt_coef_list) {
  coef_matrix <- posterior_samples_grow_ppt[[coef]]
  grow_ppt_long_data[[coef]] <- as.data.frame(coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef)
}

# Combine into one dataframe
plot_data_grow_ppt <- bind_rows(grow_ppt_long_data)

# Convert species index to numeric
plot_data_grow_ppt$species <- as.numeric(gsub("V", "", plot_data_grow_ppt$species))

# Summarize posterior: mean, median, 95% credible intervals
summary_stats_grow_ppt <- plot_data_grow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975),
    .groups = "drop"
  )

# Map species numeric IDs to abbreviations
summary_stats_grow_ppt$species <- factor(
  summary_stats_grow_ppt$species,
  levels = c(1,2,3),
  labels = c("AGHY","ELVI","PAOU")
)

# Convert to standardized effect size (Cohen's d-like)
demography_climate<-readRDS(url("https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1"))
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

# Continuous predictor SD
sd_clim <- sd(demography_grow_ppt$clim)

# Response SD
sd_y <- sd(demography_grow_ppt$y)

# Binary predictors SD
sd_endo <- sqrt(mean(demography_grow_ppt$endo) * (1 - mean(demography_grow_ppt$endo)))
sd_herb <- sqrt(mean(demography_grow_ppt$herb) * (1 - mean(demography_grow_ppt$herb)))


summary_stats_grow_ppt <- summary_stats_grow_ppt %>%
  mutate(
    d_estimate = case_when(
      parameter == "bendoclim" ~ mean_estimate * sd_clim * sd_endo / sd_y,
      parameter == "bendoherb" ~ mean_estimate * sd_endo * sd_herb / sd_y,
      parameter == "bendoherbclim" ~ mean_estimate * sd_clim * sd_endo * sd_herb / sd_y
    ),
    d_lower = case_when(
      parameter == "bendoclim" ~ lower_CI * sd_clim * sd_endo / sd_y,
      parameter == "bendoherb" ~ lower_CI * sd_endo * sd_herb / sd_y,
      parameter == "bendoherbclim" ~ lower_CI * sd_clim * sd_endo * sd_herb / sd_y
    ),
    d_upper = case_when(
      parameter == "bendoclim" ~ upper_CI * sd_clim * sd_endo / sd_y,
      parameter == "bendoherb" ~ upper_CI * sd_endo * sd_herb / sd_y,
      parameter == "bendoherbclim" ~ upper_CI * sd_clim * sd_endo * sd_herb / sd_y
    )
  )


# Rename parameters for readability
param_labels <- c(
  "bendoclim" = "Endophyte × Climate",
  "bendoherb" = "Endophyte × Herbivory",
  "bendoherbclim" = "Endophyte × Herbivory × Climate"
)
summary_stats_grow_ppt$parameter_label <- recode(summary_stats_grow_ppt$parameter, !!!param_labels)

# Plot
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_grow_coeff.pdf",
  width = 7,
  height = 9
)

ggplot(summary_stats_grow_ppt, aes(x = species, y = d_estimate, fill = parameter_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = d_lower, ymax = d_upper),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = c(
    "Endophyte × Climate" = "#1b9e77",
    "Endophyte × Herbivory" = "#d95f02",
    "Endophyte × Herbivory × Climate" = "#7570b3"
  )) +
  labs(
    x = "Species",
    y = "Standardized effect size",
    fill = "Interaction Term"
  ) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

dev.off()

## Flowering----
fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/pl6444lqmvl10s8ccxbsf/fit_flow_abio_bio_endo.rds?rlkey=x34lgk5q3n1hfryd6y2se2rm2&dl=1"))
posterior_samples_flow_ppt <- rstan::extract(fit_flow_ppt)
# Convert to data frame
posterior_samples_flow_ppt_df <- as.data.frame(posterior_samples_flow_ppt)
# Get the number of species
n_species <- length(posterior_samples_flow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_ppt_coef_list <- c(
                        "bendoclim",
                        "bendoherb",
                        "bendoherbclim")
flow_ppt_long_data <- list()
for (coef in flow_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  flow_ppt_coef_matrix <- posterior_samples_flow_ppt[[coef]]
  # Convert to long format
  flow_ppt_long_data[[coef]] <- as.data.frame(flow_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'flow_ppt_coef_list'
}

# Combine all into one dataframe
plot_data_flow_ppt <- bind_rows(flow_ppt_long_data)
# Convert species index to numeric
plot_data_flow_ppt$species <- as.numeric(gsub("V", "", plot_data_flow_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_flow_ppt <- plot_data_flow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()

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

# Change species names
summary_stats_flow_ppt$species <- recode(
  summary_stats_flow_ppt$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

# Conversion constant
conversion_factor <- sqrt(3) / pi  # ≈ 0.5513

# Add Cohen's d to each posterior sample in your long-format data:
plot_data_flow_ppt <- plot_data_flow_ppt %>%
  mutate(
    cohen_d = estimate * conversion_factor,  # estimate = posterior β
    rate_ratio = exp(estimate)               # optional: keep RR alongside
  )

# Summarize for each species and parameter:
summary_stats_flow_ppt <- plot_data_flow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975),
    mean_d = mean(cohen_d),
    lower_d = quantile(cohen_d, 0.025),
    upper_d = quantile(cohen_d, 0.975)
  ) %>% ungroup()



# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_flow_coeff.pdf",
  width = 5,
  height = 12
)
# Plot standardized effect sizes (Cohen's d) for flowering
ggplot(summary_stats_flow_ppt, aes(x = species, y = mean_d, fill = parameter)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = lower_d, ymax = upper_d),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3"
  ),
  labels = c(
    "bendoclim" = "Endophyte × Climate",
    "bendoherb" = "Endophyte × Herbivory",
    "bendoherbclim" = "Endophyte × Herbivory × Climate"
  )) +
  labs(
    x = "Species",
    y = "Standardized effect size (Cohen's d)",
    fill = "Interaction Term"
  ) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

dev.off()

## Spikelet----
fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/dyr574ub0zv4rfcbla7ov/fit_spik_abio_bio_endo.rds?rlkey=agw1q5xilj21z8wmzyziqsg6u&dl=1"))

posterior_samples_spik_ppt <- rstan::extract(fit_spik_ppt)
# Convert to data frame
posterior_samples_spik_ppt_df <- as.data.frame(posterior_samples_spik_ppt)
# Get the number of species
n_species <- length(posterior_samples_spik_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
spik_ppt_coef_list <- c("b0",
                        "bendo",
                        "bherb",
                        "bclim",
                        "bendoclim",
                        "bendoherb",
                        "bclim2",
                        "bendoclim2")
spik_ppt_long_data <- list()
for (coef in spik_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  spik_ppt_coef_matrix <- posterior_samples_spik_ppt[[coef]]
  # Convert to long format
  spik_ppt_long_data[[coef]] <- as.data.frame(spik_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'spik_ppt_coef_list'
}

# Combine all into one dataframe
plot_data_spik_ppt <- bind_rows(spik_ppt_long_data)
# Convert species index to numeric
plot_data_spik_ppt$species <- as.numeric(gsub("V", "", plot_data_spik_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_spik_ppt <- plot_data_spik_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()

# Change species names
summary_stats_spik_ppt$species <- recode(
  summary_stats_spik_ppt$species,
  "1" = "ELVI",
  "2" = "PAOU"
)

# Conversion factor
conversion_factor <- sqrt(3) / pi  # ~0.5513

# Add Cohen's d to each posterior sample
plot_data_spik_ppt <- plot_data_spik_ppt %>%
  mutate(
    cohen_d = estimate * conversion_factor,
    rate_ratio = exp(estimate)  # optional, keeps IRR for reference
  )

# Summarize by species and parameter
summary_stats_spik_ppt <- plot_data_spik_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975),
    mean_d = mean(cohen_d),
    lower_d = quantile(cohen_d, 0.025),
    upper_d = quantile(cohen_d, 0.975)
  ) %>% 
  ungroup()


library(ggplot2)

# Select only the interaction terms you want to highlight
interaction_terms <- c("bendoclim", "bendoherb", "bendoclim2")  # example

ggplot(summary_stats_spik_ppt %>% filter(parameter %in% interaction_terms),
       aes(x = species, y = mean_d, fill = parameter)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = lower_d, ymax = upper_d),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoclim2" = "#7570b3"
  ),
  labels = c(
    "bendoclim" = "Endophyte × Climate",
    "bendoherb" = "Endophyte × Herbivory",
    "bendoclim2" = "Endophyte × Climate²"
  )) +
  labs(
    x = "Species",
    y = "Standardized effect size (Cohen's d)",
    fill = "Interaction Term"
  ) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )


# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_spik_coeff.pdf",
  width = 5,
  height = 12
)
ggplot(summary_stats_spik_ppt,
       aes(
         x = factor(species),
         y = mean_estimate,
         color = species
       )) +
  geom_pointrange(
    aes(ymin = lower_CI, ymax = upper_CI),
    position = position_dodge(width = 0.6),
    size = 1
  ) + # Adds the error bars with credible intervals
  facet_grid(parameter ~ .,
             scales = "free_y",
             labeller = labeller(parameter = as_labeller(
               c(
                 "b0" = "Intercept",
                 "bendo" = "Endophyte",
                 "bherb" = "Herbivory",
                 "bclim" = "Climate",
                 "bclim2" = "Climate^2",
                 "bendoclim2" = "Endophyte:Climate^2",
                 "bendoclim" = "Endophyte:Climate",
                 "bendoherb" = "Endophyte:Herbivory"
               ),
               default = label_parsed # This tells ggplot to interpret as expressions
             ))) + # Facet by parameter
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") + # Add horizontal dashed line at y = 0
  theme_classic() +
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black")) +
  scale_color_manual(values = c(
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()
