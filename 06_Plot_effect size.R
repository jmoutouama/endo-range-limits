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
    prob_gt0 = mean(estimate > 0),   # probability effect is positive
    prob_lt0 = mean(estimate < 0),   # probability effect is negative
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    OR_mean   = exp(mean_estimate),
    OR_median = exp(median_estimate),
    OR_lower  = exp(lower_CI),
    OR_upper  = exp(upper_CI),
    strong_effect = (prob_gt0 > 0.8 | prob_lt0 > 0.8),
    # Fill color: use driver color if strong, else white
    fill_color = ifelse(strong_effect, parameter, "empty")
  )

# Change species names
# Convert species numeric IDs to abbreviations
summary_stats_surv_ppt$species <- factor(
  summary_stats_surv_ppt$species,
  levels = c(1, 2, 3),
  labels = c("AGHY", "ELVI", "POAU")
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
  x = c(-1, -0.9, 0, 0.55),  # x positions for each label
  y = c(3.35, 3.38, 3.34, 2.38)       # y positions for each label
)

# Create the forest + density plot
Fig6a <- ggplot() +
  geom_tile(
    data = summary_stats_surv_ppt %>%
      distinct(species) %>%
      mutate(
        y_id = as.numeric(factor(species)),
        fill_col = rep(c("gray90", "gray90","gray90"), length.out = n())
      ),
    aes(y = species, x = 0, fill = fill_col),
    width = Inf,
    height = 0.98,
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  # Horizontal error bars (95% CI)
  geom_errorbar(
    data = summary_stats_surv_ppt,
    aes(y = species, xmin = OR_lower, xmax = OR_upper, color = parameter_label),
    position = position_dodge(width = 0.4),
    height = 0.4, linewidth = 0.8
  ) +
  # Median points with fill reflecting strong effect
  geom_point(
    data = summary_stats_surv_ppt,
    aes(x = OR_mean, y = species, color = parameter_label, fill = fill_color),
    shape = 21,  # circle with border
    position = position_dodge(width = 0.4),
    size = 5,
    stroke = 1.2
  ) +
  xlim(-1,7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
  # Border color
  scale_color_manual(values = c(
    "Endophyte × Climate" = "#1b9e77",
    "Endophyte × Herbivory" = "#d95f02",
    "Endophyte × Herbivory × Climate" = "#7570b3"
  ), 
  guide = "none") +
  # Fill color
  scale_fill_manual(
    values = c(
      "bendoclim" = "#1b9e77",
      "bendoherb" = "#d95f02",
      "bendoherbclim" = "#7570b3",
      "empty" = "white"
    )
  ) +
  labs(
    x = "Effect size (Odds Ratio)",
    y = NULL
  ) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw(base_size = 12) +
  geom_text(
    data = panel_labels[1, ],
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0, "lines")
  )

# selected_params <- c("bendoclim", "bendoherb", "bendoherbclim")
# 
# # Convert posterior samples to long format
# surv_ridge_data <- bind_rows(
#   lapply(selected_params, function(param) {
#     mat <- posterior_samples_surv_ppt[[param]]
#     n_spp <- ncol(mat)
#     df <- as.data.frame(mat)
#     df_long <- df %>%
#       setNames(as.character(1:n_spp)) %>%
#       pivot_longer(cols = everything(),
#                    names_to = "species_idx",
#                    values_to = "logit_estimate") %>%
#       mutate(
#         species = factor(as.numeric(species_idx), levels = 1:3, labels = c("AGHY", "ELVI", "POAU")),
#         parameter = param,
#         OR = exp(logit_estimate)
#       )
#     df_long
#   })
# )
# 
# # Keep only one parameter for example (or can combine later)
# surv_ridge_data <- surv_ridge_data %>% filter(parameter %in% selected_params)
# 
# # Compute density for each species
# library(tidyr)
# library(dplyr)
# library(tidyr)
# library(ggridges)
# 
# # Compute density manually and extract x, y
# dens_data <- surv_ridge_data %>%
#   filter(parameter %in% selected_params) %>%
#   group_by(species, parameter) %>%
#   group_modify(~ {
#     d <- density(.x$OR, from = min(.x$OR), to = max(.x$OR), n = 1000)
#     tibble(
#       x = d$x,
#       y = d$y,
#       OR_category = ifelse(d$x > 1, "OR > 1", "OR < 1")
#     )
#   }) %>%
#   ungroup()
# 
# # Plot
# ggplot(dens_data, aes(x = x, y = species, height = y, fill = OR_category)) +
#   geom_ridgeline(scale = 1, alpha = 0.8, color = NA) +
#   facet_wrap(~parameter, scales = "free_x") +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
#   scale_fill_manual(values = c("OR > 1" = "#E69F00", "OR < 1" = "#56B4E9")) +
#   labs(
#     x = "Odds ratio (effect size)",
#     y = "Species",
#     fill = "Effect direction"
#   ) +
#   theme_light(base_size = 12) +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank()
#   )

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
# Update summary stats with strong/weak effect
summary_stats_grow_ppt <- plot_data_grow_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    prob_gt0 = mean(estimate > 0),
    prob_lt0 = mean(estimate < 0),
    lower_CI = quantile(estimate, 0.025),
    upper_CI = quantile(estimate, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    strong_effect = (prob_gt0 > 0.8 | prob_lt0 > 0.8),
    # Fill color: if strong, use driver color; else white
    fill_color = ifelse(strong_effect, parameter, "empty")
  )

# Update summary stats with strong/weak effect

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
#view(summary_stats_grow_ppt)
# Plot

# # Define fill colors (driver colors + white for weak effect)
# fill_values <- c(
#   "bendoclim" = "#1b9e77",
#   "bendoherb" = "#d95f02",
#   "bendoherbclim" = "#7570b3",
#   "empty" = "white"
# )

# Plot
Fig6b <- ggplot() +
  geom_tile(
    data = summary_stats_grow_ppt %>%
      distinct(species) %>%
      mutate(
        y_id = as.numeric(factor(species)),
        fill_col = rep(c("gray90", "gray90","gray90"), length.out = n())
      ),
    aes(y = species, x = 0, fill = fill_col),
    width = Inf,
    height = 0.98,  # slightly less than 1 to leave gaps
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_errorbar(
    data = summary_stats_grow_ppt,
    aes(y = species, xmin = lower_CI, xmax = upper_CI, color = parameter),
    position = position_dodge(width = 0.4), height = 0.4, linewidth = 1
  ) +
  geom_point(
    data = summary_stats_grow_ppt,
    aes(x = mean_estimate, y = species, color = parameter, fill = fill_color),
    shape = 21, size = 5, stroke = 1.2, position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3"
  ), guide = "none") +
  scale_fill_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3",
    "empty" = "white")) +  # hide fill legend
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  xlim(-1,1) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Effect size (regression coefficient)", y = NULL) +
  geom_text(
    data = panel_labels[2, ],
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )



## Inflorescence----
fit_inf_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/6ngnypika10yc0jrvr531/fit_inf_abio_bio_endo_linear.rds?rlkey=822f09t91dd2jw4r8svwmdh29&dl=1"))
posterior_samples_inf_ppt <- rstan::extract(fit_inf_ppt)
# Convert to data frame
posterior_samples_inf_ppt_df <- as.data.frame(posterior_samples_inf_ppt)
# Get the number of species
n_species <- length(posterior_samples_inf_ppt$bendo[1, ])
# Convert each coefficient into a long-format data frame
inf_ppt_coef_list <- c(
                        "bendoclim",
                        "bendoherb",
                        "bendoherbclim")
inf_ppt_long_data <- list()
for (coef in inf_ppt_coef_list) {
  # Extract the coefficient matrix for the current parameter
  inf_ppt_coef_matrix <- posterior_samples_inf_ppt[[coef]]
  # Convert to long format
  inf_ppt_long_data[[coef]] <- as.data.frame(inf_ppt_coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef) # Use 'coef' instead of 'inf_ppt_coef_list'
}

# Combine all into one dataframe
plot_data_inf_ppt <- bind_rows(inf_ppt_long_data)
# Convert species index to numeric
plot_data_inf_ppt$species <- as.numeric(gsub("V", "", plot_data_inf_ppt$species))
# Calculate the mean, median, and 95% credible intervals for each species and coefficient
summary_stats_inf_ppt <- plot_data_inf_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    median_estimate = median(estimate),
    prob_gt0 = mean(estimate > 0),
    prob_lt0 = mean(estimate < 0),
    # Compute IRR percentiles directly from posterior draws
    IRR_lower = quantile(exp(estimate), 0.025),
    IRR_median = quantile(exp(estimate), 0.5),
    IRR_upper = quantile(exp(estimate), 0.975),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    # Flag strong effects
    strong_effect = (prob_gt0 > 0.8 | prob_lt0 > 0.8),
    # Fill color: use parameter color if strong, else white
    fill_color = ifelse(strong_effect, parameter, "empty")
  )

# Change species names
summary_stats_inf_ppt$species <- recode(
  summary_stats_inf_ppt$species,
  "1" = "AGHY",
  "2" = "ELVI",
  "3" = "PAOU"
)
param_labels <- c(
  "bendoclim" = "Endophyte × Climate",
  "bendoherb" = "Endophyte × Herbivory",
  "bendoherbclim" = "Endophyte × Herbivory × Climate"
)
summary_stats_inf_ppt$parameter_label <- recode(summary_stats_inf_ppt$parameter, !!!param_labels)

# Plot standardized effect sizes (Cohen's d) for infering
Fig6c <- ggplot() +
  geom_tile(
    data = summary_stats_inf_ppt %>%
      distinct(species) %>%
      mutate(
        y_id = as.numeric(factor(species)),
        fill_col = rep(c("gray90", "gray90","gray90"), length.out = n())
      ),
    aes(y = species, x = 0, fill = fill_col),
    width = Inf,
    height = 0.98,  # slightly less than 1 to leave gaps
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_errorbarh(
    data = summary_stats_inf_ppt,
    aes(y = species, xmin = IRR_lower, xmax = IRR_upper, color = parameter),
    position = position_dodge(width = 0.4), height = 0.4, linewidth = 1
  ) +
  geom_point(
    data = summary_stats_inf_ppt,
    aes(x = IRR_median, y = species, color = parameter, fill = fill_color),
    shape = 21, size = 5, stroke = 1.2, position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3"
  ), guide = "none") +
  scale_fill_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3",
    "empty" = "white")) +  # hide fill legend
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
  xlim(-0.15,4) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Effect size (Incidence Rate Ratios)", y = NULL) +
  geom_text(
    data = panel_labels[3, ],
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )



## Spikelet----
fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/6fcebl4lw8mu94fz62hnh/fit_spik_abio_bio_endo_linear.rds?rlkey=zy25y44zocugs6shh68lwpy1q&dl=1"))

posterior_samples_spik_ppt <- rstan::extract(fit_spik_ppt)

spik_ppt_coef_list <- c("bendoclim", "bendoherb", "bendoherbclim")
spik_ppt_long_data <- list()

for (coef in spik_ppt_coef_list) {
  coef_matrix <- posterior_samples_spik_ppt[[coef]]
  spik_ppt_long_data[[coef]] <- as.data.frame(coef_matrix) %>%
    pivot_longer(cols = everything(),
                 names_to = "species",
                 values_to = "estimate") %>%
    mutate(parameter = coef)
}

plot_data_spik_ppt <- bind_rows(spik_ppt_long_data)
plot_data_spik_ppt$species <- as.numeric(gsub("V", "", plot_data_spik_ppt$species))

# Compute summary statistics and IRRs
summary_stats_spik_ppt <- plot_data_spik_ppt %>%
  group_by(parameter, species) %>%
  summarize(
    median_estimate = median(estimate),
    prob_gt0 = mean(estimate > 0),
    prob_lt0 = mean(estimate < 0),
    IRR_lower = quantile(exp(estimate), 0.025),
    IRR_median = quantile(exp(estimate), 0.5),
    IRR_upper = quantile(exp(estimate), 0.975),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    strong_effect = (prob_gt0 > 0.8 | prob_lt0 > 0.8),
    fill_color = ifelse(strong_effect, parameter, "empty")
  )

# Change species names
summary_stats_spik_ppt$species <- recode(
  summary_stats_spik_ppt$species,
  "1" = "ELVI",
  "2" = "PAOU"
)

# Parameter labels
param_labels <- c(
  "bendoclim" = "Endophyte × Climate",
  "bendoherb" = "Endophyte × Herbivory",
  "bendoherbclim" = "Endophyte × Herbivory × Climate"
)
summary_stats_spik_ppt$parameter_label <- recode(summary_stats_spik_ppt$parameter, !!!param_labels)

# Plot
Fig6d <- ggplot() +
  geom_tile(
    data = summary_stats_spik_ppt %>%
      distinct(species) %>%
      mutate(
        y_id = as.numeric(factor(species)),
        fill_col = rep(c("gray90", "gray90"), length.out = n())
      ),
    aes(y = species, x = 0.5, fill = fill_col),
    width = Inf,
    height = 0.98,  # slightly less than 1 to leave gaps
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_errorbarh(
    data = summary_stats_spik_ppt,
    aes(y = species, xmin = IRR_lower, xmax = IRR_upper, color = parameter),
    position = position_dodge(width = 0.4), height = 0.4, linewidth = 1
  ) +
  geom_point(
    data = summary_stats_spik_ppt,
    aes(x = IRR_median, y = species, color = parameter, fill = fill_color),
    shape = 21, size = 5, stroke = 1.2, position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3"
  ),guide = "none") +
  scale_fill_manual(values = c(
    "bendoclim" = "#1b9e77",
    "bendoherb" = "#d95f02",
    "bendoherbclim" = "#7570b3",
    "empty" = "white")) +  # hide fill legend
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
  xlim(0.5, 2) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Effect size (Incidence Rate Ratios)", y = NULL) +
  geom_text(
    data = panel_labels[4, ],
    aes(x = x, y = y, label = label),
    fontface = "plain", size = 4.5, hjust = 0,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

#Combine with a shared legend

# Assuming your plots are named Fig6a, Fig6b, Fig6c, Fig6d
Figure6 <- ggarrange(
  Fig6a, Fig6b, Fig6c, Fig6d,
  ncol = 2, nrow = 2,       # 2x2 layout
  #labels = c("a", "b", "c", "d"),
  # label.x = 0.02,            # optional: adjust label position
  # label.y = 0.95,
  common.legend = TRUE,      # share a single legend
  legend = "top"          # legend position
)
ggsave("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Figure6.pdf", Figure6, width = 10, height = 8)


# Add a new column for each data frame
surv_df <- summary_stats_surv_ppt %>%
  mutate(Fitness_component = "Survival") %>%
  rename(Parameter = parameter_label, Estimate =OR_median , Lower_CI = OR_lower, Upper_CI = OR_upper, Pr=prob_1)

grow_df <- summary_stats_grow_ppt %>%
  mutate(Fitness_component = "Growth") %>%
  rename(Parameter = parameter_label, Estimate = median_estimate, Lower_CI = lower_CI, Upper_CI = upper_CI,Pr=prob_gt0)

flow_df <- summary_stats_flow_ppt %>%
  mutate(Fitness_component = "Flowering") %>%
  rename(Parameter = parameter, Estimate = IRR_median, Lower_CI = IRR_lower, Upper_CI = IRR_upper,Pr=prob_gt1)

spike_df <- summary_stats_spik_ppt %>%
  mutate(Fitness_component = "Spikelet") %>%
  rename(Parameter = parameter, Estimate = IRR_median, Lower_CI = IRR_lower, Upper_CI = IRR_upper, Pr=prob_gt1)

# Combine all data frames
combined_summary <- bind_rows(surv_df, grow_df, flow_df, spike_df) %>%
  dplyr::select(Fitness_component, species, Parameter, Estimate, Lower_CI, Upper_CI,Pr)

# Export to CSV
write_csv(combined_summary, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/combined_summary_table_drivers.csv")
