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
library(dplyr)
library(readr)
library(ggpubr)

# Surivival models
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/khyez2xegn8j5elkchaf6/fit_surv_abio_bio_endo_linear.rds?rlkey=zh4hx9czjov9aivlmaycfcuq7&dl=1"))
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



# Rename parameters for readability
param_labels <- c(
  "bendoclim" = "Endophyte × Climate",
  "bendoherb" = "Endophyte × Herbivory",
  "bendoherbclim" = "Endophyte × Herbivory × Climate"
)

summary_stats_surv_ppt$parameter_label <- recode(summary_stats_surv_ppt$parameter, !!!param_labels)


panel_labels <- data.frame(
  label = c("(a)", "(b)", "(c)", "(d)"),
  x = c(0.5, 0.5, 0.5, 0.5),  # x positions for each label
  y = c(1.5, 0.7, 1.2, 0.48)       # y positions for each label
)

# Create the coefficient plot with error bars (credible intervals)

Fig6a<-ggplot(summary_stats_surv_ppt, 
       aes(x = species, 
           y = median_estimate, 
           fill = parameter_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
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
    fill = "Drivers",
    #title = "Relative Strength of Factors Modifying the Endophyte Effect on Survival",
    #subtitle = "Positive = more mutualistic; Negative = more parasitic"
  ) +
  theme_light(base_size = 12) +
  geom_text(
    data = panel_labels[1, ],  # first label "a"
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  )+
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

## Growth----
fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/o62tvjf8aqqz15gjxnrjn/fit_grow_abio_bio_endo_linear.rds?rlkey=xg1s6u5ctsluampm1l2zy1wqn&dl=1"))
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

# Rename parameters for readability
param_labels <- c(
  "bendoclim" = "Endophyte × Climate",
  "bendoherb" = "Endophyte × Herbivory",
  "bendoherbclim" = "Endophyte × Herbivory × Climate"
)
summary_stats_grow_ppt$parameter_label <- recode(summary_stats_grow_ppt$parameter, !!!param_labels)

# Plot

Fig6b<-ggplot(summary_stats_grow_ppt, aes(x = species, y = median_estimate, fill = parameter_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
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
    fill = "Drivers"
  ) +
  geom_text(
    data = panel_labels[2, ],  # first label "a"
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  )+
  theme_light(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

## Flowering----
fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/1v4f4thyh826qcuiiyhub/fit_flow_abio_bio_endo_linear.rds?rlkey=raj4ls5dcqkeeexvcqj8b495m&dl=1"))
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


# Change species names
summary_stats_flow_ppt$species <- recode(
  summary_stats_flow_ppt$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)

# Plot standardized effect sizes (Cohen's d) for flowering
Fig6c <- ggplot(summary_stats_flow_ppt, aes(x = species, y = mean_estimate, fill = parameter)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
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
    y = "Standardized effect size",
    fill = "Drivers"   # <-- changed here
  ) +
  geom_text(
    data = panel_labels[3, ],  # first label "a"
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  )+
  theme_light(base_size = 12) +
  theme(
    legend.position = c(0.5, 0.9),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )


## Spikelet----
fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/6fcebl4lw8mu94fz62hnh/fit_spik_abio_bio_endo_linear.rds?rlkey=zy25y44zocugs6shh68lwpy1q&dl=1"))

posterior_samples_spik_ppt <- rstan::extract(fit_spik_ppt)
# Convert to data frame
posterior_samples_spik_ppt_df <- as.data.frame(posterior_samples_spik_ppt)
# Get the number of species
n_species <- length(posterior_samples_spik_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
spik_ppt_coef_list <- c(
                        "bendoclim",
                        "bendoherb",
                        "bendoherbclim")
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



# Select only the interaction terms you want to highlight
Drivers <- c("bendoclim", "bendoherb", "bendoherbclim")  # example

Fig6d<-ggplot(summary_stats_spik_ppt %>% filter(parameter %in% Drivers),
       aes(x = species, y = median_estimate, fill = parameter)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
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
    y = "Standardized effect size",
    fill = "Drivers"
  ) +
  geom_text(
    data = panel_labels[4, ],  # first label "a"
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  )+
  theme_light(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

#Combine with a shared legend

# Assuming your plots are named Fig6a, Fig6b, Fig6c, Fig6d
Figure6 <- ggarrange(
  Fig6a, Fig6b, Fig6c, Fig6d,
  ncol = 2, nrow = 2,       # 2x2 layout
  #labels = c("a", "b", "c", "d"),
  label.x = 0.02,            # optional: adjust label position
  label.y = 0.95,
  common.legend = TRUE,      # share a single legend
  legend = "bottom"          # legend position
)
ggsave("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Figure6.pdf", Figure6, width = 10, height = 8)


# Add a new column for each data frame
surv_df <- summary_stats_surv_ppt %>%
  mutate(Fitness_component = "Survival") %>%
  rename(Parameter = parameter_label, Estimate = median_estimate, Lower_CI = lower_CI, Upper_CI = upper_CI)

grow_df <- summary_stats_grow_ppt %>%
  mutate(Fitness_component = "Growth") %>%
  rename(Parameter = parameter_label, Estimate = median_estimate, Lower_CI = lower_CI, Upper_CI = upper_CI)

flow_df <- summary_stats_flow_ppt %>%
  mutate(Fitness_component = "Flowering") %>%
  rename(Parameter = parameter, Estimate = mean_estimate, Lower_CI = lower_CI, Upper_CI = upper_CI)

spike_df <- summary_stats_spik_ppt %>%
  mutate(Fitness_component = "Spikelet") %>%
  rename(Parameter = parameter, Estimate = median_estimate, Lower_CI = lower_CI, Upper_CI = upper_CI)

# Combine all data frames
combined_summary <- bind_rows(surv_df, grow_df, flow_df, spike_df) %>%
  dplyr::select(Fitness_component, species, Parameter, Estimate, Lower_CI, Upper_CI)

# Export to CSV
write_csv(combined_summary, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/combined_summary_table_drivers.csv")
