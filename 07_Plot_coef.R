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
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_surv_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))
## Plot the coefficients
### Precipitation
posterior_samples_surv_ppt <- rstan::extract(fit_surv_ppt)
# Convert to data frame
posterior_samples_surv_ppt_df <- as.data.frame(posterior_samples_surv_ppt)
# Get the number of species
n_species <- length(posterior_samples_surv_ppt$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_ppt_coef_list <- c("b0",
                        "bendo",
                        "bherb",
                        "bclim",
                        "bendoclim",
                        "bendoherb",
                        "bclim2",
                        "bendoclim2")
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
summary_stats_surv_ppt$species <- recode(
  summary_stats_surv_ppt$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_surv_coeff.pdf",
  width = 6,
  height = 12
)
ggplot(summary_stats_surv_ppt,
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
    "AGHY" = "#009E73",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

## Growth----
fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/5oduhrkn3l0cu5b9soju5/fit_grow_ppt.rds?rlkey=mpxxl4aejowhdm29pij9jv8kk&dl=1"))
posterior_samples_grow_ppt <- rstan::extract(fit_grow_ppt)
# Convert to data frame
posterior_samples_grow_ppt_df <- as.data.frame(posterior_samples_grow_ppt)
# Get the number of species
n_species <- length(posterior_samples_grow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
grow_ppt_coef_list <- c("b0",
                        "bendo",
                        "bherb",
                        "bclim",
                        "bendoclim",
                        "bendoherb",
                        "bclim2",
                        "bendoclim2")
grow_ppt_long_data <- list()
for (coef in grow_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  grow_ppt_coef_matrix <- posterior_samples_grow_ppt[[coef]]
  # Convert to long format
  grow_ppt_long_data[[coef]] <- as.data.frame(grow_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'grow_ppt_coef_list'
}

# Combine all into one dataframe
plot_data_grow_ppt <- bind_rows(grow_ppt_long_data)
# Convert species index to numeric
plot_data_grow_ppt$species <- as.numeric(gsub("V", "", plot_data_grow_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_grow_ppt <- plot_data_grow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()

# Change species names
summary_stats_grow_ppt$species <- recode(
  summary_stats_grow_ppt$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_grow_coeff.pdf",
  width = 5,
  height = 12
)
ggplot(summary_stats_grow_ppt,
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
    "AGHY" = "#009E73",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

## Flowering----
fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/1j3ln3jxk94s56c9j193q/fit_flow_ppt.rds?rlkey=ag5bdlhngtg2gsfbfbx15purr&dl=1"))
posterior_samples_flow_ppt <- rstan::extract(fit_flow_ppt)
# Convert to data frame
posterior_samples_flow_ppt_df <- as.data.frame(posterior_samples_flow_ppt)
# Get the number of species
n_species <- length(posterior_samples_flow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_ppt_coef_list <- c("b0",
                        "bendo",
                        "bherb",
                        "bclim",
                        "bendoclim",
                        "bendoherb",
                        "bclim2",
                        "bendoclim2")
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

# Change species names
summary_stats_flow_ppt$species <- recode(
  summary_stats_flow_ppt$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_flow_coeff.pdf",
  width = 5,
  height = 12
)
ggplot(summary_stats_flow_ppt,
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
    "AGHY" = "#009E73",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

## Spikelet----
fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/fg562lkkl077j8qoegb8q/fit_spik_ppt.rds?rlkey=cbfyy7tq7tja3e9mxujv0scm6&dl=1"))

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
