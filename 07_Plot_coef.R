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
fit_surv_ppt <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/hi11gxhpqlrdfg389ir0w/fit_surv_ppt.rds?rlkey=22ujyjnm74c6pw9biw50uasvu&dl=1"
  )
)
fit_surv_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/gxc8edjzdjvsrtlb8zm5o/fit_surv_distance.rds?rlkey=bmtq9q0bxf6ooafq8pz3tyavp&dl=1"
  )
)
fit_surv_geo_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/ah7h4cn9l9p6gwl3cl0aw/fit_surv_geo_distance.rds?rlkey=dypxz2231pwle7p8mb534uqwr&dl=1"
  )
)


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

### Distance from niche center
posterior_samples_surv_distance <- rstan::extract(fit_surv_distance)
# Convert to data frame
posterior_samples_surv_distance_df <- as.data.frame(posterior_samples_surv_distance)
# Get the number of species
n_species <- length(posterior_samples_surv_distance$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
surv_distance_long_data <- list()
for (coef in surv_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_distance_coef_matrix <- posterior_samples_surv_distance[[coef]]
  # Convert to long format
  surv_distance_long_data[[coef]] <- as.data.frame(surv_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_distance_coef_list'
}
# Combine all into one dataframe
plot_data_surv_distance <- bind_rows(surv_distance_long_data)
# Convert species index to numeric
plot_data_surv_distance$species <- as.numeric(gsub("V", "", plot_data_surv_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_distance <- plot_data_surv_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_surv_distance$species <- recode(
  summary_stats_surv_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)
# unique(summary_stats_surv_distance$parameter)

### Distance from geographic center
posterior_samples_surv_geo_distance <- rstan::extract(fit_surv_geo_distance)
# Convert to data frame
posterior_samples_surv_geo_distance_df <- as.data.frame(posterior_samples_surv_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_surv_geo_distance$bendo_s[1, ])
# Convert each coefficient into a long-format data frame
surv_geo_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
surv_geo_distance_long_data <- list()
for (coef in surv_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  surv_geo_distance_coef_matrix <- posterior_samples_surv_geo_distance[[coef]]
  # Convert to long format
  surv_geo_distance_long_data[[coef]] <- as.data.frame(surv_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'surv_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_surv_geo_distance <- bind_rows(surv_geo_distance_long_data)
# Convert species index to numeric
plot_data_surv_geo_distance$species <- as.numeric(gsub("V", "", plot_data_surv_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_surv_geo_distance <- plot_data_surv_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_surv_geo_distance$species <- recode(
  summary_stats_surv_geo_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)
# unique(summary_stats_surv_geo_distance$parameter)
# Create the coefficient plot with error bars (credible intervals)
pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_surv_coeff.pdf",
  width = 5,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/distance_surv_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_surv_distance,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/geodistance_surv_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_surv_geo_distance,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()



## Growth----
fit_grow_ppt <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/85pzzrgvkxwogoybr004t/fit_grow_ppt.rds?rlkey=rz2xlg00u1aqhkxeaix7wiu9e&dl=1"
  )
)
fit_grow_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/1m01hqgokktxz77trnuu3/fit_spik_distance.rds?rlkey=jxx4nt8qgwgiospq870anu9mn&dl=1"
  )
)
fit_grow_geo_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/q7fvaf40gwkww9kbxgdj8/fit_spik_geo_distance.rds?rlkey=jw5f1k57l3uimleyutuek63q5&dl=1"
  )
)

posterior_samples_grow_ppt <- rstan::extract(fit_grow_ppt)
# Convert to data frame
posterior_samples_grow_ppt_df <- as.data.frame(posterior_samples_grow_ppt)
# Get the number of species
n_species <- length(posterior_samples_grow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
grow_ppt_coef_list <- c("bendo",
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

unique(summary_stats_grow_ppt$parameter)

### Distance from climatic center
posterior_samples_grow_distance <- rstan::extract(fit_grow_distance)
# Convert to data frame
posterior_samples_grow_distance_df <- as.data.frame(posterior_samples_grow_distance)
# Get the number of species
n_species <- length(posterior_samples_grow_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
grow_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
grow_distance_long_data <- list()
for (coef in grow_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  grow_distance_coef_matrix <- posterior_samples_grow_distance[[coef]]
  # Convert to long format
  grow_distance_long_data[[coef]] <- as.data.frame(grow_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'grow_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_grow_distance <- bind_rows(grow_distance_long_data)
# Convert species index to numeric
plot_data_grow_distance$species <- as.numeric(gsub("V", "", plot_data_grow_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_grow_distance <- plot_data_grow_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_grow_distance$species <- recode(
  summary_stats_grow_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)


### Distance from geographic center
posterior_samples_grow_geo_distance <- rstan::extract(fit_grow_geo_distance)
# Convert to data frame
posterior_samples_grow_geo_distance_df <- as.data.frame(posterior_samples_grow_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_grow_geo_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
grow_geo_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
grow_geo_distance_long_data <- list()
for (coef in grow_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  grow_geo_distance_coef_matrix <- posterior_samples_grow_geo_distance[[coef]]
  # Convert to long format
  grow_geo_distance_long_data[[coef]] <- as.data.frame(grow_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'grow_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_grow_geo_distance <- bind_rows(grow_geo_distance_long_data)
# Convert species index to numeric
plot_data_grow_geo_distance$species <- as.numeric(gsub("V", "", plot_data_grow_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_grow_geo_distance <- plot_data_grow_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_grow_geo_distance$species <- recode(
  summary_stats_grow_geo_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

# Create the coefficient plot with error bars (credible intervals)

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_grow_coeff.pdf",
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/distance_grow_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_grow_distance,
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
                 "bendoclim" = "Endophyte:Climate",
                 "bendoherb" = "Endophyte:Herbivory"
               ),
               default = label_parsed # This tells ggplot to interpret as expressions
             ))) + # Facet by parameter
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme_classic() +
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black")) +
  scale_color_manual(values = c(
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/geodistance_grow_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_grow_geo_distance,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

## Flowering----
fit_flow_ppt <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/ajo9euxdtxfmle2g006l3/fit_flow_ppt.rds?rlkey=78cknj1yaqs2m5cyicdzg8up4&dl=1"
  )
)
fit_flow_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/y7qvz4pmy2t2j00gqvwqd/fit_flow_distance.rds?rlkey=qe47fg37dpi4z2lu1ur1lsbkh&dl=1"
  )
)
fit_flow_geo_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/b2dez0gz22p3m2ke1zs7f/fit_flow_geo_distance.rds?rlkey=9fuwynap20617l1rqffjz62wd&dl=1"
  )
)

posterior_samples_flow_ppt <- rstan::extract(fit_flow_ppt)
# Convert to data frame
posterior_samples_flow_ppt_df <- as.data.frame(posterior_samples_flow_ppt)
# Get the number of species
n_species <- length(posterior_samples_flow_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_ppt_coef_list <- c("bendo",
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

### Distance from climatic center
posterior_samples_flow_distance <- rstan::extract(fit_flow_distance)
# Convert to data frame
posterior_samples_flow_distance_df <- as.data.frame(posterior_samples_flow_distance)
# Get the number of species
n_species <- length(posterior_samples_flow_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
flow_distance_long_data <- list()
for (coef in flow_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  flow_distance_coef_matrix <- posterior_samples_flow_distance[[coef]]
  # Convert to long format
  flow_distance_long_data[[coef]] <- as.data.frame(flow_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'flow_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_flow_distance <- bind_rows(flow_distance_long_data)
# Convert species index to numeric
plot_data_flow_distance$species <- as.numeric(gsub("V", "", plot_data_flow_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_flow_distance <- plot_data_flow_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_flow_distance$species <- recode(
  summary_stats_flow_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

### Distance from geographic center
posterior_samples_flow_geo_distance <- rstan::extract(fit_flow_geo_distance)
# Convert to data frame
posterior_samples_flow_geo_distance_df <- as.data.frame(posterior_samples_flow_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_flow_geo_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
flow_geo_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
flow_geo_distance_long_data <- list()
for (coef in flow_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  flow_geo_distance_coef_matrix <- posterior_samples_flow_geo_distance[[coef]]
  # Convert to long format
  flow_geo_distance_long_data[[coef]] <- as.data.frame(flow_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'flow_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_flow_geo_distance <- bind_rows(flow_geo_distance_long_data)
# Convert species index to numeric
plot_data_flow_geo_distance$species <- as.numeric(gsub("V", "", plot_data_flow_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_flow_geo_distance <- plot_data_flow_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_flow_geo_distance$species <- recode(
  summary_stats_flow_geo_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_flow_coeff.pdf",
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/distance_flow_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_flow_distance,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/geodistance_flow_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_flow_geo_distance,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()


## Spikelet----
fit_spik_ppt <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/9klrwl866stup2x3m2p9x/fit_spik_ppt.rds?rlkey=okno4gvn9dj5q8uhpnu36f59i&dl=1"
  )
)
fit_spik_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/1m01hqgokktxz77trnuu3/fit_spik_distance.rds?rlkey=jxx4nt8qgwgiospq870anu9mn&dl=1"
  )
)
fit_spik_geo_distance <- readRDS(
  url(
    "https://www.dropbox.com/scl/fi/q7fvaf40gwkww9kbxgdj8/fit_spik_geo_distance.rds?rlkey=jw5f1k57l3uimleyutuek63q5&dl=1"
  )
)

posterior_samples_spik_ppt <- rstan::extract(fit_spik_ppt)
# Convert to data frame
posterior_samples_spik_ppt_df <- as.data.frame(posterior_samples_spik_ppt)
# Get the number of species
n_species <- length(posterior_samples_spik_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
spik_ppt_coef_list <- c("bendo",
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
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

### Distance from climatic center
posterior_samples_spik_distance <- rstan::extract(fit_spik_distance)
# Convert to data frame
posterior_samples_spik_distance_df <- as.data.frame(posterior_samples_spik_distance)
# Get the number of species
n_species <- length(posterior_samples_spik_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
spik_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
spik_distance_long_data <- list()
for (coef in spik_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  spik_distance_coef_matrix <- posterior_samples_spik_distance[[coef]]
  # Convert to long format
  spik_distance_long_data[[coef]] <- as.data.frame(spik_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'spik_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_spik_distance <- bind_rows(spik_distance_long_data)
# Convert species index to numeric
plot_data_spik_distance$species <- as.numeric(gsub("V", "", plot_data_spik_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_spik_distance <- plot_data_spik_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_spik_distance$species <- recode(
  summary_stats_spik_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

### Distance from geographic center
posterior_samples_spik_geo_distance <- rstan::extract(fit_spik_geo_distance)
# Convert to data frame
posterior_samples_spik_geo_distance_df <- as.data.frame(posterior_samples_spik_geo_distance)
# Get the number of species
n_species <- length(posterior_samples_spik_geo_distance$bendo[1, ])
# Convert each coefficient into a long-format data frame
spik_geo_distance_coef_list <- c("b0", "bendo", "bherb", "bclim", "bendoclim", "bendoherb")
spik_geo_distance_long_data <- list()
for (coef in spik_geo_distance_coef_list) {
  # Extract the coefficient matrix for the current parameter
  spik_geo_distance_coef_matrix <- posterior_samples_spik_geo_distance[[coef]]
  # Convert to long format
  spik_geo_distance_long_data[[coef]] <- as.data.frame(spik_geo_distance_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'spik_geo_distance_coef_list'
}
# Combine all into one dataframe
plot_data_spik_geo_distance <- bind_rows(spik_geo_distance_long_data)
# Convert species index to numeric
plot_data_spik_geo_distance$species <- as.numeric(gsub("V", "", plot_data_spik_geo_distance$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_spik_geo_distance <- plot_data_spik_geo_distance %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975)
  ) %>%
  ungroup()
# Change species names
summary_stats_spik_geo_distance$species <- recode(
  summary_stats_spik_geo_distance$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/ppt_spik_coeff.pdf",
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()

pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/distance_spik_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_spik_distance,
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
                 "bendoclim" = "Endophyte:Climate",
                 "bendoherb" = "Endophyte:Herbivory"
               ),
               default = label_parsed # This tells ggplot to interpret as expressions
             ))) + # Facet by parameter
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") + # Add horizontal dashed line at y = 0
  theme_bw() +
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme_classic() +
  labs(x = "Species", y = "Coefficient estimate", title = "") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black")) +
  scale_color_manual(values = c(
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()


pdf(
  "/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/endo-range-limits/Figure/geodistance_spik_coeff.pdf",
  width = 5,
  height = 9
)
ggplot(summary_stats_spik_geo_distance,
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
    "AGHY" = "#000000",
    "ELVI" = "#D55E00",
    "PAOU" = "#0072B2"
  )) # Assign unique colors for species
dev.off()
