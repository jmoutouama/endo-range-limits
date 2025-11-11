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
library(ggh4x)
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
demography_climate<-readRDS(url("https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1"))
climate_max<-readRDS(url("https://www.dropbox.com/scl/fi/v9lp1wvc2eeatgf0uf7ol/prism_max_2023_2025.rds?rlkey=debj1y03k9buxank3s5wadbcc&dl=1"))
climate_scaled<-readRDS(url("https://www.dropbox.com/scl/fi/irecsnoh3xrq6g8d5cysa/climate_site_scaled.rds?rlkey=63r7ugrtkuo5ncmqfoywbidps&dl=1"))

yearly_prism_max <- climate_max %>%
  group_by(species, year) %>%
  summarise(
    ppt_sum = sum(ppt, na.rm = TRUE),
    tmean_mean = mean(tmean, na.rm = TRUE)
  ) %>%
  ungroup()
yearly_prism_max_avg <- yearly_prism_max %>%
  group_by(species) %>%
  summarise(
    ppt_avg = mean(ppt_sum, na.rm = TRUE),
    tmean_avg = mean(tmean_mean, na.rm = TRUE)
  ) %>%
  ungroup()

# Optional: also get corresponding real precipitation values for reference
ppt_mean <- mean(climate_scaled$ppt_log)
ppt_sd   <- sd(climate_scaled$ppt_log)

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
    log_size_t0 = log(tiller_t),
    surv_t1     = as.integer(surv1),
    ppt         = ppt_scaled
  )

# summary(demography_climate_surv$cum_ppt)
# any(demography_climate_surv$cum_ppt <= 0)
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

# Load model and data
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/khyez2xegn8j5elkchaf6/fit_surv_abio_bio_endo_linear.rds?rlkey=zh4hx9czjov9aivlmaycfcuq7&dl=1"))

# Prediction grid
climate_range <- seq(min(demography_surv_ppt$clim),
                     max(demography_surv_ppt$clim),
                     length.out = 30)
predictions <- expand.grid(
  clim = climate_range,
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:3
)
climate_range_mm <- exp(climate_range * ppt_sd + ppt_mean)

# Extract posterior samples
posterior_samples_survival <- rstan::extract(fit_surv_ppt)
# Prediction function
get_predictions_survival <- function(clim, endo, herb, species_index, posterior_samples_survival) {
  with(posterior_samples_survival, {
    logit_preds <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bendoherbclim[, species_index] * endo * herb * clim 
    1 / (1 + exp(-logit_preds))
  })
}

# Generate predictions
n_post_survival <- nrow(posterior_samples_survival$b0)
pred_probs_matrix_survival <- matrix(NA, nrow = nrow(predictions), ncol = n_post_survival)

for (i in seq_len(nrow(predictions))) {
  pred_probs_matrix_survival[i, ] <- get_predictions_survival(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples_survival
  )
}

# Combine with predictors
pred_probs_df_survival <- cbind(predictions, as.data.frame(pred_probs_matrix_survival))
pred_probs_long_df_survival <- pred_probs_df_survival %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Posterior_Sample",
    values_to = "Pred_Survival"
  )

# Credible intervals
cred_intervals_survival <- pred_probs_long_df_survival %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    median = quantile(Pred_Survival, 0.5),
    mean = mean(Pred_Survival),
    .groups = "drop"
  ) %>%
  mutate(panel = "Pr (survival)")

# Observed data
observed_data_survival <- demography_surv_ppt %>% 
  data.frame(
    clim = .$clim,
    endo = .$endo,
    herb = .$herb,
    species = .$Spp,
    plot = .$plot,
    y = .$y
  ) %>% 
  group_by(plot, species, herb, clim, endo) %>% 
  summarise(y_plot_mean = mean(y, na.rm = TRUE), .groups = "drop") %>%
  mutate(panel = "Pr (survival)")

# Differences (E+ - E-)
diff_df_survival <- pred_probs_long_df_survival %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Survival[endo == 1]) - mean(Pred_Survival[endo == 0]),
    .groups = "drop"
  )

diff_ci_survival <- diff_df_survival %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    mean = mean(diff),
    .groups = "drop"
  ) %>%
  mutate(panel = "Δ (E+ - E-)")

# Merge panels
plot_data_survival <- bind_rows(cred_intervals_survival, diff_ci_survival)
# Ensure species column is factor with parseable labels for italics
plot_data_survival$species <- factor(
  plot_data_survival$species,
  levels = c("1","2","3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

observed_data_survival$species <- factor(
  observed_data_survival$species,
  levels = c("1","2","3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# Add back-transformed climate to plot_data_survival
plot_data_survival <- plot_data_survival %>%
  mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

observed_data_survival <- observed_data_survival %>%
  mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))


# Create label table (a–f for 3 species × 2 herbivory)
panel_labels <- data.frame(
  species = rep(c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  ), each = 2),
  herb = rep(c(0, 1), times = 3),
  label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
  panel = "Pr (survival)"   # only place labels on upper panels
)
# Trim panel names
plot_data_survival$panel <- trimws(as.character(plot_data_survival$panel))
observed_data_survival$panel <- trimws(as.character(observed_data_survival$panel))

# Plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/PrSurvival_diff.pdf",
  width = 6, height = 7
)
ggplot(plot_data_survival) +
  # Survival panel
  geom_line(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_survival, panel == "Pr (survival)"),
    aes(x = climate_mm, y = y_plot_mean, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +
  
  # Δ panel
  geom_line(
    data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  geom_hline(
    data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0), linetype = "dashed", color = "black"
  ) +
  
  # Facets
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Unfenced", "1" = "Fenced")
    )
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      panel == "Δ (E+ - E-)" ~ scale_y_continuous(limits = c(-0.3, 0.3), expand = c(0,0)),
      panel == "Pr (survival)" ~ scale_y_continuous(expand = c(0,0))
    )
  ) +
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.058, 0.23),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1)
  ) +
  # Add facet labels (a–f) only on top panels
  geom_text(
    data = panel_labels,
    aes(x = 490, y = 0.9, label = label),
    fontface = "plain", size = 3.5, hjust = 0,
    inherit.aes = FALSE
  )

dev.off()

# Quantiles of climate for survival
clim_quantiles_surv <- quantile(demography_surv_ppt$clim, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

# Function to compute Δ(E+ − E−) for survival 
compute_delta_surv <- function(clim_val, posterior_samples_surv, herb_values = c(0,1)) {
  n_species <- dim(posterior_samples_surv$b0)[2]
  n_post <- dim(posterior_samples_surv$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      # Posterior predictions for E+ (endo present) and E− (endo absent)
      pred_Eplus <- 1 / (1 + exp(-(
        posterior_samples_surv$b0[, sp] +
          posterior_samples_surv$bendo[, sp] * 1 +
          posterior_samples_surv$bherb[, sp] * h +
          posterior_samples_surv$bclim[, sp] * clim_val +
          posterior_samples_surv$bendoclim[, sp] * clim_val * 1 +
          posterior_samples_surv$bendoherb[, sp] * 1 * h +
          posterior_samples_surv$bendoherbclim[, sp] * 1 * h * clim_val 
      )))
      
      pred_Eminus <- 1 / (1 + exp(-(
        posterior_samples_surv$b0[, sp] +
          posterior_samples_surv$bendo[, sp] * 0 +
          posterior_samples_surv$bherb[, sp] * h +
          posterior_samples_surv$bclim[, sp] * clim_val +
          posterior_samples_surv$bendoclim[, sp] * clim_val * 0 +
          posterior_samples_surv$bendoherb[, sp] * 0 * h +
          posterior_samples_surv$bendoherbclim[, sp] * 0 * h * clim_val 
      )))
      
      delta <- pred_Eplus - pred_Eminus
      
      data.frame(
        species = sp,
        herb = h,
        clim = clim_val,
        Posterior_Sample = 1:n_post,
        delta = delta
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  return(result)
}

# Compute Δ for all climate quantiles ---
delta_surv_quantiles <- lapply(clim_quantiles_surv, function(cl) compute_delta_surv(cl, posterior_samples_survival)) %>%
  bind_rows()

# Summarize Δ
delta_surv_summary <- delta_surv_quantiles %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  )

# Relabel species and herbivore treatments
delta_surv_summary$species <- factor(delta_surv_summary$species, levels = 1:3,
                                     labels = c(
                                       "Agrostis hyemalis",
                                       "Elymus virginicus",
                                       "Poa autumnalis"
                                     ))
delta_surv_summary$herb <- factor(delta_surv_summary$herb, levels = c(0,1),
                                  labels = c("Unfenced", "Fenced"))

# Update delta_surv_summary with back-transformed climate in mm
delta_surv_summary <- delta_surv_summary %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))

# Prepare data for plotting
delta_long_surv <- delta_surv_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  tidyr::pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           "median_delta" = "Median Δ (E+ − E−)",
                           "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = dplyr::case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Plot for survival
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/PrSurvival_diff_stat.pdf",
  width = 7, height = 6
)

ggplot(delta_long_surv, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  # Horizontal dashed lines
  geom_hline(
    data = delta_long_surv %>% filter(metric == "Median Δ (E+ − E−)"),
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    data = delta_long_surv %>% filter(metric == "Pr (Δ > 0)"),
    aes(yintercept = 0.5),
    linetype = "dashed",
    color = "black"
  ) +
  # Vertical lines for species-specific average precipitation
  geom_vline(
    data = yearly_prism_max_avg %>% 
      mutate(species_label = case_when(
        species == "aghy" ~ "italic('Agrostis hyemalis')",
        species == "elvi" ~ "italic('Elymus virginicus')",
        species == "poau" ~ "italic('Poa autumnalis')"
      )),
    aes(xintercept = ppt_avg),
    linetype = "dotted",
    color = "#0072B2",
    size = 0.8
  ) +
  # Facets
  facet_grid(metric ~ species_label, scales = "free_y",
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Unfenced" = "#E69F00", "Fenced" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore exclusion",
    title = ""
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black", face = "plain"),
    strip.text.y = element_text(size = 10, color = "black", face = "plain"),
    strip.background = element_rect(color="black", fill="grey80", size=0.1, linetype="solid")
  )

dev.off()

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

fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/o62tvjf8aqqz15gjxnrjn/fit_grow_abio_bio_endo_linear.rds?rlkey=xg1s6u5ctsluampm1l2zy1wqn&dl=1"))
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
# Extract posterior samples
posterior_samples_grow <- rstan::extract(fit_grow_ppt)
# Prediction function for growth
get_predictions_grow <- function(clim, endo, herb, species_index, posterior_samples_grow) {
  with(posterior_samples_grow, {
    mu_preds <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bendoherbclim[, species_index] * endo * herb * clim
    mu_preds  
  })
}

# Generate predictions
n_post_grow <- nrow(posterior_samples_grow$b0)
pred_matrix_grow <- matrix(NA, nrow = nrow(predictions), ncol = n_post_grow)

for (i in seq_len(nrow(predictions))) {
  pred_matrix_grow[i, ] <- get_predictions_grow(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples_grow
  )
}

# Compute credible intervals
pred_grow_df <- cbind(predictions, as.data.frame(pred_matrix_grow))
pred_grow_long <- pred_grow_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Posterior_Sample",
    values_to = "Pred_Growth"
  )

# Credible intervals
cred_intervals_grow <- pred_grow_long %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth),
    .groups = "drop"
  ) %>%
  mutate(panel = "Growth")
# Compute E+ − E− difference panels
diff_df_grow <- pred_grow_long %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Growth[endo == 1]) - mean(Pred_Growth[endo == 0]),
    .groups = "drop"
  )

diff_ci_grow <- diff_df_grow %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    mean = mean(diff),
    .groups = "drop"
  ) %>%
  mutate(panel = "Δ (E+ - E-)")
# Combine both panels
plot_data_grow <- bind_rows(cred_intervals_grow, diff_ci_grow)

plot_data_grow$species <- factor(
  plot_data_grow$species,
  levels = c("1","2","3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

observed_data_grow <- demography_grow_ppt %>%
  data.frame(
    clim = .$clim,
    endo = .$endo,
    herb = .$herb,
    species = .$Spp,
    plot = .$plot,
    y = .$y
  ) %>%
  mutate(
    panel = "Growth",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

observed_data_grow$species <- factor(
  observed_data_grow$species,
  levels = c("1", "2", "3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# Back-transform climate for plotting
plot_data_grow <- plot_data_grow %>%
  mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

y_limits <- plot_data_grow %>%
  dplyr::filter(panel == "Δ (E+ - E-)") %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    ymin = min(lower_90, na.rm = TRUE),
    ymax = max(upper_90, na.rm = TRUE)
  )

panel_labels_grow <- data.frame(
  species = rep(c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  ), each = 2),
  herb = rep(c(0, 1), times = 3),
  label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
  panel = "Growth"   # only place labels on upper panels
)

# Plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_diff.pdf",
  width = 7, height = 8.75
)
ggplot(plot_data_grow) +
  # Upper panel: predicted growth
  geom_line(
    data = subset(plot_data_grow, panel == "Growth"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_grow, panel == "Growth"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  # Observed points
  geom_point(
    data = subset(observed_data_grow, panel == "Growth"),
    aes(x =climate_mm, y = y, color = factor(endo)),
    size = 0.75,
    position = position_jitter(width = 0, height = 0.01)
  ) +
  # Lower panel: Δ(E+ − E−) differences
  geom_line(
    data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
    aes(x =climate_mm, y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  geom_hline(
    data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0), linetype = "dashed", color = "black"
  ) +
  # Facets
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Unfenced", "1" = "Fenced")
    )
  ) +
  # Dynamic y-axis scales
  ggh4x::facetted_pos_scales(
    y = list(
      # Lower panels (Δ(E+ − E−)) – only 0 tick, no labels
      panel == "Δ (E+ - E-)" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(
          breaks = 0,
          labels = 0,
          minor_breaks = NULL,
          limits = c(
            y_limits$ymin[y_limits$species == "italic('Agrostis hyemalis')"],
            y_limits$ymax[y_limits$species == "italic('Agrostis hyemalis')"]
          ),
          expand = c(0, 0)
        ),
      panel == "Δ (E+ - E-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(
          breaks = 0,
          labels = 0,
          minor_breaks = NULL,
          limits = c(
            y_limits$ymin[y_limits$species == "italic('Elymus virginicus')"],
            y_limits$ymax[y_limits$species == "italic('Elymus virginicus')"]
          ),
          expand = c(0, 0)
        ),
      panel == "Δ (E+ - E-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(
          breaks = 0,
          labels = 0,
          minor_breaks = NULL,
          limits = c(
            y_limits$ymin[y_limits$species == "italic('Poa autumnalis')"],
            y_limits$ymax[y_limits$species == "italic('Poa autumnalis')"]
          ),
          expand = c(0, 0)
        ),
      
      # Upper panels (Growth) – manually reduced per species
      panel == "Growth" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(limits = c(-3.5, 1), expand = c(0, 0)),
      panel == "Growth" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(-2, 1.5), expand = c(0, 0)),
      panel == "Growth" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(-4, 2), expand = c(0, 0))
    )
  ) +
  # Labels and theme
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.4, 0.2),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 12, color = "black", face = "plain"),
    strip.text.y = element_text(size = 8.5, color = "black", face = "plain"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
  )+
  geom_text(
    data = panel_labels_grow,
    aes(x = 490, y = 0.8, label = label),
    fontface = "plain", size = 3.5, hjust = 0,
    inherit.aes = FALSE
  )

dev.off()


# Plot
# Cairo::CairoPDF(
#   "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_diff_h.pdf",
#   width = 11, height = 6.5
# )
# 
# ggplot(plot_data_grow) +
#   # Upper panel: predicted growth
#   geom_line(
#     data = subset(plot_data_grow, panel == "Growth"),
#     aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
#     size = 0.5
#   ) +
#   geom_ribbon(
#     data = subset(plot_data_grow, panel == "Growth"),
#     aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
#     alpha = 0.3, color = NA
#   ) +
#   # Observed points
#   geom_point(
#     data = subset(observed_data_grow, panel == "Growth"),
#     aes(x =climate_mm, y = y, color = factor(endo)),
#     size = 0.75,
#     position = position_jitter(width = 0, height = 0.01)
#   ) +
#   # Lower panel: Δ(E+ − E−) differences
#   geom_line(
#     data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
#     aes(x =climate_mm, y = mean), color = "black", size = 0.5
#   ) +
#   geom_ribbon(
#     data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
#     aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
#     fill = "#9B6B96", alpha = 0.5
#   ) +
#   geom_hline(
#     data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
#     aes(yintercept = 0), linetype = "dashed", color = "black"
#   ) +
#   # Facets
#   ggh4x::facet_nested(
#     panel ~ species + herb,
#     scales = "free_y",
#     space = "free_y",
#     labeller = labeller(
#       species = label_parsed,
#       herb = c("0" = "Unfenced", "1" = "Fenced")
#     )
#   ) +
#   ggh4x::facetted_pos_scales(
#     y = list(
#       # Lower panels – Δ(E+ − E−)
#       panel == "Δ (E+ - E-)" & species == "italic('Agrostis hyemalis')" ~
#         scale_y_continuous(
#           breaks = 0, labels = 0,
#           limits = c(
#             y_limits$ymin[y_limits$species == "italic('Agrostis hyemalis')"],
#             y_limits$ymax[y_limits$species == "italic('Agrostis hyemalis')"]
#           ),
#           expand = c(0, 0)
#         ),
#       panel == "Δ (E+ - E-)" & species == "italic('Elymus virginicus')" ~
#         scale_y_continuous(
#           breaks = 0, labels = 0,
#           limits = c(
#             y_limits$ymin[y_limits$species == "italic('Elymus virginicus')"],
#             y_limits$ymax[y_limits$species == "italic('Elymus virginicus')"]
#           ),
#           expand = c(0, 0)
#         ),
#       panel == "Δ (E+ - E-)" & species == "italic('Poa autumnalis')" ~
#         scale_y_continuous(
#           breaks = 0, labels = 0,
#           limits = c(
#             y_limits$ymin[y_limits$species == "italic('Poa autumnalis')"],
#             y_limits$ymax[y_limits$species == "italic('Poa autumnalis')"]
#           ),
#           expand = c(0, 0)
#         ),
#       
#       # Upper panels – Growth
#       panel == "Growth" & species == "italic('Agrostis hyemalis')" ~
#         scale_y_continuous(limits = c(-4, 3), expand = c(0, 0)),
#       panel == "Growth" & species == "italic('Elymus virginicus')" ~
#         scale_y_continuous(limits = c(-4, 3), expand = c(0, 0)),
#       panel == "Growth" & species == "italic('Poa autumnalis')" ~
#         scale_y_continuous(limits = c(-13, 7), expand = c(0, 0))
#     )
#   )+
#   # Labels and theme
#   labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
#   scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
#   scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
#   theme_light() +
#   theme(
#     legend.position = c(0.05, 0.85),
#     legend.title = element_text(size = 8),
#     legend.text = element_text(size = 8),
#     panel.spacing.y = unit(0.0, "cm"),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 6),
#     axis.line.y = element_line(color = "black", size = 0.1),
#     axis.line.x = element_line(color = "black", size = 0.1),
#     text = element_text(family = "Arial"),
#     strip.text.x = element_text(size = 14, color = "black", face = "plain"),
#     strip.text.y = element_text(size = 14, color = "black", face = "plain"),
#     strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
#   )+
#   geom_text(
#     data = panel_labels_grow,
#     aes(x = 490, y = 2.75, label = label),
#     fontface = "plain", size = 3.5, hjust = 0,
#     inherit.aes = FALSE
#   )
# 
# 
# dev.off()

# Climate quantiles for grow 
clim_quantiles_grow <- quantile(demography_grow_ppt$clim, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

#Function to compute Δ(E+ − E−) for grow 
compute_delta_grow <- function(clim_val, posterior_samples_grow, herb_values = c(0, 1)) {
  n_species <- dim(posterior_samples_grow$b0)[2]
  n_post <- dim(posterior_samples_grow$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      # Predicted mean growth for E+ (endophyte present)
      pred_Eplus <- posterior_samples_grow$b0[, sp] +
        posterior_samples_grow$bendo[, sp] * 1 +
        posterior_samples_grow$bherb[, sp] * h +
        posterior_samples_grow$bclim[, sp] * clim_val +
        posterior_samples_grow$bendoclim[, sp] * clim_val * 1 +
        posterior_samples_grow$bendoherb[, sp] * 1 * h +
        posterior_samples_grow$bendoherbclim[, sp] * 1 * h * clim_val 
      
      
      # Predicted mean growth for E− (endophyte absent)
      pred_Eminus <- posterior_samples_grow$b0[, sp] +
        posterior_samples_grow$bendo[, sp] * 0 +
        posterior_samples_grow$bherb[, sp] * h +
        posterior_samples_grow$bclim[, sp] * clim_val +
        posterior_samples_grow$bendoclim[, sp] * clim_val * 0 +
        posterior_samples_grow$bendoherb[, sp] * 0 * h +
        posterior_samples_grow$bendoherbclim[, sp] * 0 * h * clim_val 

      
      # Δ(E+ − E−)
      delta <- pred_Eplus - pred_Eminus
      
      data.frame(
        species = sp,
        herb = h,
        clim = clim_val,
        Posterior_Sample = 1:n_post,
        delta = delta
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  return(result)
}

#  Compute Δ for all climate quantiles
delta_grow_quantiles <- lapply(clim_quantiles_grow, function(cl) compute_delta_grow(cl, posterior_samples_grow)) %>%
  bind_rows()

#  Summarize Δ 
delta_grow_summary <- delta_grow_quantiles %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  )

# Relabel species and herbivore treatments 
delta_grow_summary$species <- factor(delta_grow_summary$species, levels = 1:3,
                                     labels = c(
                                       "Agrostis hyemalis",
                                       "Elymus virginicus",
                                       "Poa autumnalis"
                                     ))
delta_grow_summary$herb <- factor(delta_grow_summary$herb, levels = c(0,1),
                                  labels = c("Unfenced", "Fenced"))

# Add exponentiated climate for plotting 
delta_grow_summary <- delta_grow_summary %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))


# Prepare long-format data for growth
delta_long_grow <- delta_grow_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  tidyr::pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           "median_delta" = "Median Δ (E+ − E−)",
                           "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = dplyr::case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Create the growth plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_diff_stat.pdf",
  width = 7, height = 6
)
ggplot(delta_long_grow, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  # Horizontal dashed line at y = 0 for Median Δ
  geom_hline(
    data = delta_long_grow %>% filter(metric == "Median Δ (E+ − E−)"),
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  # Horizontal dashed line at y = 0.5 for Pr(Δ > 0)
  geom_hline(
    data = delta_long_grow %>% filter(metric == "Pr (Δ > 0)"),
    aes(yintercept = 0.5),
    linetype = "dashed",
    color = "black"
  ) +
  # Vertical lines for species-specific average precipitation
  geom_vline(
    data = yearly_prism_max_avg %>% 
      mutate(species_label = case_when(
        species == "aghy" ~ "italic('Agrostis hyemalis')",
        species == "elvi" ~ "italic('Elymus virginicus')",
        species == "poau" ~ "italic('Poa autumnalis')"
      )),
    aes(xintercept = ppt_avg),
    linetype = "dotted",
    color = "#0072B2",
    size = 0.8
  ) +
  facet_grid(metric ~ species_label, scales = "free_y",
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Unfenced" = "#E69F00", "Fenced" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore exclusion",
    title = ""
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black", face = "plain"),
    strip.text.y = element_text(size = 10, color = "black", face = "plain"),
    strip.background = element_rect(color="black", fill="grey80", size=0.1, linetype="solid")
  )

dev.off()

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
    flow_t1 = inf_t1,
    ppt = ppt_scaled
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

fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/1v4f4thyh826qcuiiyhub/fit_flow_abio_bio_endo_linear.rds?rlkey=raj4ls5dcqkeeexvcqj8b495m&dl=1"))
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
# Extract posterior samples
posterior_samples_flow <- rstan::extract(fit_flow_ppt)
# Prediction function for flowering
get_predictions_flow <- function(clim, endo, herb, species_index, posterior_samples_flow) {
  with(posterior_samples_flow, {
    # Linear predictor
    eta <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bendoherbclim[, species_index] * endo * herb * clim
    
    mu <- exp(eta)   # negative binomial mean on response scale
    
    # Apply zero-inflation (global parameter, not species-specific)
    (1 - zi) * mu
  })
}

# Generate predictions
n_post_flow <- nrow(posterior_samples_flow$b0)
pred_matrix_flow <- matrix(NA, nrow = nrow(predictions), ncol = n_post_flow)

for (i in seq_len(nrow(predictions))) {
  pred_matrix_flow[i, ] <- get_predictions_flow(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples_flow
  )
}

# Combine predictions with original predictors
pred_flow_df <- cbind(predictions, as.data.frame(pred_matrix_flow))

# Pivot longer for tidy format
pred_flow_long <- pred_flow_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Posterior_Sample",
    values_to = "Pred_Flow"
  )

# Compute credible intervals per species × endo × herb × clim
cred_intervals_flow <- pred_flow_long %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow),
    .groups = "drop"
  ) %>%
  mutate(panel = "Flowering")

# Compute Δ (E+ − E−) for flowering
diff_df_flow <- pred_flow_long %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Flow[endo == 1]) - mean(Pred_Flow[endo == 0]),
    .groups = "drop"
  )

diff_ci_flow <- diff_df_flow %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    mean = mean(diff),
    .groups = "drop"
  ) %>%
  mutate(panel = "Δ (E+ - E-)")

# Combine credible intervals and difference panels
plot_data_flow <- bind_rows(cred_intervals_flow, diff_ci_flow)

# Relabel species for plotting in italics
plot_data_flow$species <- factor(
  plot_data_flow$species,
  levels = c("1","2","3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# Change panel name for upper panel
plot_data_flow <- plot_data_flow %>%
  mutate(panel = ifelse(panel == "Flowering", "#Inflorescences", panel))

observed_data_flow <- demography_flow_ppt %>%  # or your dataset for flowering
  data.frame(
    clim = .$clim,
    endo = .$endo,
    herb = .$herb,
    species = .$Spp,
    plot = .$plot,
    y = .$y
  ) %>%
  # group_by(plot, species, herb, clim, endo) %>%
  # summarise(
  #   y_plot_mean = mean(y, na.rm = TRUE),
  #   .groups = "drop"
  # ) %>%
  mutate(panel = "#Inflorescences",
         climate_mm = exp(clim * ppt_sd + ppt_mean))

# Relabel species (matching your plot_data_flow)
observed_data_flow$species <- factor(
  observed_data_flow$species,
  levels = c("1", "2", "3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# species_levels <- c(
#   "italic('Agrostis hyemalis')",
#   "italic('Elymus virginicus')",
#   "italic('Poa autumnalis')"
# )

# plot_data_flow$species <- factor(plot_data_flow$species, levels = species_levels)
# observed_data_flow$species <- factor(observed_data_flow$species, levels = species_levels)
# panel_labels_flow$species <- factor(panel_labels_flow$species, levels = species_levels)


# observed_data_flow <- observed_data_flow %>%
#   mutate(panel = ifelse(panel == "Flowering", "#Inflorescences", panel))
# Back-transform climate for plotting
plot_data_flow <- plot_data_flow %>%
  mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

# Δ-panel limits for each species (same, only Δ panel)
y_limits_flow <- plot_data_flow %>%
  filter(panel == "Δ (E+ - E-)") %>%
  group_by(species) %>%
  summarise(
    ymin = min(lower_90, na.rm = TRUE),
    ymax = max(upper_90, na.rm = TRUE),
    .groups = "drop"
  )

panel_labels_flow <- data.frame(
  species = rep(c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  ), each = 2),
  herb = rep(c(0, 1), times = 3),
  label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
  panel = "#Inflorescences"   # only place labels on upper panels
)

#Plot with updated panel labels
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flowering_diff.pdf",
  width = 7, height = 10
)
ggplot(plot_data_flow) +
  # Upper panel: predicted flowering counts (#Inflorescences)
  geom_line(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_flow, panel == "#Inflorescences"),
    aes(x = climate_mm, y = y, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +

  # Lower panel: Δ (E+ - E-) differences
  geom_line(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +

  # Facets: species vertically, herb horizontally
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Unfenced", "1" = "Fenced")
    )
  ) +

  # Facetted scales
  ggh4x::facetted_pos_scales(
    y = list(
      # Lower panels – Δ (E+ - E-) only
      panel == "Δ (E+ - E-)" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(
          breaks = 0, labels = 0, minor_breaks = NULL,
          limits = c(y_limits_flow$ymin[y_limits_flow$species == "italic('Agrostis hyemalis')"],
                     y_limits_flow$ymax[y_limits_flow$species == "italic('Agrostis hyemalis')"]),
          expand = c(0, 0)
        ),
      panel == "Δ (E+ - E-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(
          breaks = 0, labels = 0, minor_breaks = NULL,
          limits = c(y_limits_flow$ymin[y_limits_flow$species == "italic('Elymus virginicus')"],
                     y_limits_flow$ymax[y_limits_flow$species == "italic('Elymus virginicus')"]),
          expand = c(0, 0)
        ),
      panel == "Δ (E+ - E-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(
          breaks = 0, labels = 0, minor_breaks = NULL,
          limits = c(y_limits_flow$ymin[y_limits_flow$species == "italic('Poa autumnalis')"],
                     y_limits_flow$ymax[y_limits_flow$species == "italic('Poa autumnalis')"]),
          expand = c(0, 0)
        ),

      # Upper panels – #Inflorescences custom limits per species
      panel == "#Inflorescences" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(limits = c(0, 15), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 7), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 25), expand = c(0, 0))
    )
  ) +
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.075, 0.4),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black", face = "plain"),
    strip.text.y = element_text(size = 8, color = "black", face = "plain"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
  )

dev.off()


Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flowering_diff_h.pdf",
  width = 9, height = 5
)

ggplot(plot_data_flow) +
  geom_line(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_flow, panel == "#Inflorescences"),
    aes(x = climate_mm, y = y, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +
  geom_line(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, y = mean), color = "black", size = 0.5
  ) +
  # Add a horizontal dashed line at y = 0 in the Δ-panel
  geom_hline(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0),
    color = "black",
    linetype = "dashed",
    size = 0.3
  )+
  geom_ribbon(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  ggh4x::facet_nested(
    panel ~ species + herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Unfenced", "1" = "Fenced")
    )
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      panel == "Δ (E+ - E-)" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_flow$ymin[y_limits_flow$species == "italic('Agrostis hyemalis')"],
          y_limits_flow$ymax[y_limits_flow$species == "italic('Agrostis hyemalis')"]
        ), expand = c(0, 0)),
      panel == "Δ (E+ - E-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_flow$ymin[y_limits_flow$species == "italic('Elymus virginicus')"],
          y_limits_flow$ymax[y_limits_flow$species == "italic('Elymus virginicus')"]
        ), expand = c(0, 0)),
      panel == "Δ (E+ - E-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_flow$ymin[y_limits_flow$species == "italic('Poa autumnalis')"],
          50
        ), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(limits = c(0, 30), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 15), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 70), expand = c(0, 0))
    )
  ) +
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.05, 0.8),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 12, color = "black", face = "plain"),
    strip.text.y = element_text(size = 12, color = "black", face = "plain"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
  )+
  geom_text(
    data = panel_labels_flow,
    aes(x = 490, y = 19, label = label),
    fontface = "plain", size = 3.5, hjust = 0,
    inherit.aes = FALSE
  )

dev.off()

# Climate quantiles for flowering
clim_quantiles_flow <- quantile(demography_flow_ppt$clim, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

# Function to compute Δ(E+ − E−) for flowering
compute_delta_flow <- function(clim_val, posterior_samples_flow, herb_values = c(0, 1)) {
  n_species <- dim(posterior_samples_flow$b0)[2]
  n_post <- dim(posterior_samples_flow$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      # Linear predictor for E+
      eta_Eplus <- posterior_samples_flow$b0[, sp] +
        posterior_samples_flow$bendo[, sp] * 1 +
        posterior_samples_flow$bherb[, sp] * h +
        posterior_samples_flow$bclim[, sp] * clim_val +
        posterior_samples_flow$bendoclim[, sp] * clim_val * 1 +
        posterior_samples_flow$bendoherb[, sp] * 1 * h +
        posterior_samples_flow$bendoherbclim[, sp] * 1 * h * clim_val 
        #posterior_samples_flow$bclim2[, sp] * clim_val^2
      
      # Linear predictor for E−
      eta_Eminus <- posterior_samples_flow$b0[, sp] +
        posterior_samples_flow$bendo[, sp] * 0 +
        posterior_samples_flow$bherb[, sp] * h +
        posterior_samples_flow$bclim[, sp] * clim_val +
        posterior_samples_flow$bendoclim[, sp] * clim_val * 0 +
        posterior_samples_flow$bendoherb[, sp] * 0 * h +
        posterior_samples_flow$bendoherbclim[, sp] * 0 * h * clim_val 
        #posterior_samples_flow$bclim2[, sp] * clim_val^2
      
      # Predicted mean flowering counts (log link → exp)
      pred_Eplus  <- exp(eta_Eplus)
      pred_Eminus <- exp(eta_Eminus)
      
      # Δ(E+ − E−): difference in expected flower count
      delta <- pred_Eplus - pred_Eminus
      
      data.frame(
        species = sp,
        herb = h,
        clim = clim_val,
        Posterior_Sample = 1:n_post,
        delta = delta
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  return(result)
}

# Compute Δ for all climate quantiles (flowering)
delta_flow_quantiles <- lapply(clim_quantiles_flow, function(cl) compute_delta_flow(cl, posterior_samples_flow)) %>%
  bind_rows()

# Summarize Δ
delta_flow_summary <- delta_flow_quantiles %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  )

# Relabel species and herbivore treatments
delta_flow_summary$species <- factor(delta_flow_summary$species, levels = 1:3,
                                     labels = c(
                                       "Agrostis hyemalis",
                                       "Elymus virginicus",
                                       "Poa autumnalis"
                                     ))
delta_flow_summary$herb <- factor(delta_flow_summary$herb, levels = c(0,1),
                                  labels = c("Unfenced", "Fenced"))

# Add exponentiated climate for plotting
# Add exponentiated climate for plotting 
delta_flow_summary <- delta_flow_summary %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))


# Prepare long-format data for flowering
delta_long_flow <- delta_flow_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  tidyr::pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           "median_delta" = "Median Δ (E+ − E−)",
                           "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = dplyr::case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Create the flowering plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flowering_diff_stat.pdf",
  width = 7, height = 5
)
ggplot(delta_long_flow, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  # Horizontal dashed line at y = 0 for Median Δ
  geom_hline(
    data = delta_long_flow %>% filter(metric == "Median Δ (E+ − E−)"),
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  # Horizontal dashed line at y = 0.5 for Pr(Δ > 0)
  geom_hline(
    data = delta_long_flow %>% filter(metric == "Pr (Δ > 0)"),
    aes(yintercept = 0.5),
    linetype = "dashed",
    color = "black"
  ) +
  # Vertical lines for species-specific average precipitation
  geom_vline(
    data = yearly_prism_max_avg %>% 
      mutate(species_label = case_when(
        species == "aghy" ~ "italic('Agrostis hyemalis')",
        species == "elvi" ~ "italic('Elymus virginicus')",
        species == "poau" ~ "italic('Poa autumnalis')"
      )),
    aes(xintercept = ppt_avg),
    linetype = "dotted",
    color = "#0072B2",
    size = 0.6
  ) +
  facet_grid(metric ~ species_label, scales = "free_y", 
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Unfenced" = "#E69F00", "Fenced" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore exclusion",
    title = ""
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black", face = "plain"),
    strip.text.y = element_text(size = 10, color = "black", face = "plain"),
    strip.background = element_rect(color="black", fill="grey80", size=0.1, linetype="solid")
  )

dev.off()

# Spikelet----
demography_climate %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  filter(tiller_t1 > 0,inf_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cum_ppt,ppt_scaled
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


fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/6fcebl4lw8mu94fz62hnh/fit_spik_abio_bio_endo_linear.rds?rlkey=zy25y44zocugs6shh68lwpy1q&dl=1"))
posterior_samples_spik <- rstan::extract(fit_spik_ppt)
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

# Prediction function
get_predictions_spik <- function(clim, endo, herb, species_index, posterior_samples_spik) {
  with(posterior_samples_spik, {
    eta <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bendoherbclim[, species_index] * endo * herb * clim 
    
    mu <- exp(eta)
    mu
  })
}

# Generate posterior predictions 
n_post_spik <- nrow(posterior_samples_spik$b0)
pred_matrix_spik <- matrix(NA, nrow = nrow(predictions), ncol = n_post_spik)

for (i in seq_len(nrow(predictions))) {
  pred_matrix_spik[i, ] <- get_predictions_spik(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples_spik
  )
}
# Summarize posterior predictions 
pred_spik_df <- cbind(predictions, as.data.frame(pred_matrix_spik))
pred_spik_long <- pred_spik_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Posterior_Sample",
    values_to = "Prediction"
  )

plot_data_spik <- pred_spik_long %>%
  group_by(species, clim, endo, herb) %>%
  summarise(
    mean = mean(Prediction),
    lower_90 = quantile(Prediction, 0.05),
    upper_90 = quantile(Prediction, 0.95),
    .groups = "drop"
  )

# Δ (E+ - E-) calculation
delta_spik <- plot_data_spik %>%
  pivot_wider(names_from = endo, values_from = c(mean, lower_90, upper_90)) %>%
  mutate(
    mean = mean_1 - mean_0,
    lower_90 = lower_90_1 - upper_90_0,
    upper_90 = upper_90_1 - lower_90_0,
    panel = "Δ (E+ - E-)"
  ) %>%
  dplyr::select(species, clim, herb, mean, lower_90, upper_90, panel)

plot_data_spik <- plot_data_spik %>%
  mutate(panel = "#Spikelets") %>%
  bind_rows(delta_spik)

# Species labels
plot_data_spik <- plot_data_spik %>%
  mutate(
    species = factor(species,
                     labels = c("italic('Elymus virginicus')",
                                "italic('Poa autumnalis')"))
  )

# Observed data for spikelets 
observed_data_spik <- demography_spik_ppt %>%
  data.frame(
    clim = .$clim,
    endo = .$endo,
    herb = .$herb,
    species = .$Spp,
    plot = .$plot,
    y = .$y
  ) %>%
  group_by(plot, species, herb, clim, endo) %>%
  # summarise(
  #   y_plot_mean = mean(y, na.rm = TRUE),
  #   .groups = "drop"
  # ) %>%
  mutate(panel = "#Spikelets",
         climate_mm = exp(clim * ppt_sd + ppt_mean))

# Relabel species
observed_data_spik$species <- factor(
  observed_data_spik$species,
  levels = c("1", "2"),
  labels = c(
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# Back-transform climate for plotting
plot_data_spik <- plot_data_spik %>%
  mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

# y-limits for Δ panels 
y_limits_spik <- plot_data_spik %>%
  filter(panel == "Δ (E+ - E-)") %>%
  group_by(species) %>%
  summarise(
    ymin = min(lower_90, na.rm = TRUE),
    ymax = max(upper_90, na.rm = TRUE),
    .groups = "drop"
  )
panel_labels_spik <- data.frame(
  species = rep(c(
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  ), each = 2),
  herb = rep(c(0, 1), times = 2),
  label = c("(a)", "(b)", "(c)", "(d)"),
  panel = "#Spikelets"   # only place labels on upper panels
)

# --- Plot (horizontal layout for paper figure) ---
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spikelet_diff.pdf",
  width = 6, height = 6
)

ggplot(plot_data_spik) +
  # Upper panel: predicted spikelets
  geom_line(
    data = subset(plot_data_spik, panel == "#Spikelets"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_spik, panel == "#Spikelets"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90,
        fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  
  # Lower panel: Δ (E+ - E-)
  geom_line(
    data = subset(plot_data_spik, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, y = mean),
    color = "black", size = 0.5
  ) +
  geom_hline(
    data = subset(plot_data_spik, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0),
    color = "black", linetype = "dashed", size = 0.3
  ) +
  geom_point(
    data = subset(observed_data_spik, panel == "#Spikelets"),
    aes(x = climate_mm, y = y, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +
  geom_ribbon(
    data = subset(plot_data_spik, panel == "Δ (E+ - E-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  
  # Facets
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Unfenced", "1" = "Fenced")
    )
  ) +
  
  # Facetted y-scales
  ggh4x::facetted_pos_scales(
    y = list(
      panel == "Δ (E+ - E-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_spik$ymin[y_limits_spik$species == "italic('Elymus virginicus')"],
          y_limits_spik$ymax[y_limits_spik$species == "italic('Elymus virginicus')"]
        ), expand = c(0, 0)),
      panel == "Δ (E+ - E-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_spik$ymin[y_limits_spik$species == "italic('Poa autumnalis')"],
          y_limits_spik$ymax[y_limits_spik$species == "italic('Poa autumnalis')"]
        ), expand = c(0, 0)),
      panel == "#Spikelets" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 50), expand = c(0, 0)),
      panel == "#Spikelets" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 60), expand = c(0, 0))
    )
  ) +
  labs(
    x = "Precipitation (mm)",
    y = "",
    color = "Endophyte",
    fill = "Endophyte"
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.12, 0.47),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 12, color = "black", face = "plain"),
    strip.text.y = element_text(size = 10, color = "black", face = "plain"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
  )+
  geom_text(
    data = panel_labels_spik,
    aes(x = 490, y = 42, label = label),
    fontface = "plain", size = 3.5, hjust = 0,
    inherit.aes = FALSE
  )

dev.off()

# Climate quantiles for spike
clim_quantiles_spik <- quantile(demography_spik_ppt$clim, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

# Function to compute Δ(E+ − E−) for spike
compute_delta_spik <- function(clim_val, posterior_samples_spik, herb_values = c(0, 1)) {
  n_species <- dim(posterior_samples_spik$b0)[2]
  n_post <- dim(posterior_samples_spik$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      # Linear predictor for E+
      eta_Eplus <- posterior_samples_spik$b0[, sp] +
        posterior_samples_spik$bendo[, sp] * 1 +
        posterior_samples_spik$bherb[, sp] * h +
        posterior_samples_spik$bclim[, sp] * clim_val +
        posterior_samples_spik$bendoclim[, sp] * clim_val * 1 +
        posterior_samples_spik$bendoherb[, sp] * 1 * h +
        posterior_samples_spik$bendoherbclim[, sp] * 1 * h * clim_val 
      
      # Linear predictor for E−
      eta_Eminus <- posterior_samples_spik$b0[, sp] +
        posterior_samples_spik$bendo[, sp] * 0 +
        posterior_samples_spik$bherb[, sp] * h +
        posterior_samples_spik$bclim[, sp] * clim_val +
        posterior_samples_spik$bendoclim[, sp] * clim_val * 0 +
        posterior_samples_spik$bendoherb[, sp] * 0 * h +
        posterior_samples_spik$bendoherbclim[, sp] * 0 * h * clim_val 
      
      # Predicted mean spike counts (log link → exp)
      pred_Eplus  <- exp(eta_Eplus)
      pred_Eminus <- exp(eta_Eminus)
      
      # Δ(E+ − E−): difference in expected spike number
      delta <- pred_Eplus - pred_Eminus
      
      data.frame(
        species = sp,
        herb = h,
        clim = clim_val,
        Posterior_Sample = 1:n_post,
        delta = delta
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  return(result)
}

# Compute Δ for all climate quantiles
delta_spik_quantiles <- lapply(clim_quantiles_spik, function(cl) compute_delta_spik(cl, posterior_samples_spik)) %>%
  bind_rows()

# Summarize Δ
delta_spik_summary <- delta_spik_quantiles %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  )

# Relabel species and herbivore treatments
delta_spik_summary$species <- factor(delta_spik_summary$species, levels = 1:2,
                                     labels = c(
                                       "Elymus virginicus",
                                       "Poa autumnalis"
                                     ))
delta_spik_summary$herb <- factor(delta_spik_summary$herb, levels = c(0,1),
                                  labels = c("Unfenced", "Fenced"))

# Add exponentiated climate for plotting
delta_spik_summary <- delta_spik_summary %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))

# Prepare long-format data for plotting
delta_long_spik <- delta_spik_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  tidyr::pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           "median_delta" = "Median Δ (E+ − E−)",
                           "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = dplyr::case_when(
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Optional: plot Δ for spike
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spike_diff_stat.pdf",
  width = 6, height = 5
)

ggplot(
  delta_long_spik %>%
    filter(!is.na(species_label)),  # ✅ remove NA panels
  aes(x = clim_mm, y = value, color = herb, group = herb)
) +
  geom_line(size = 0.5) +
  geom_hline(
    data = delta_long_spik %>%
      filter(metric == "Median Δ (E+ − E−)", !is.na(species_label)),
    aes(yintercept = 0),
    linetype = "dashed", color = "black"
  ) +
  geom_hline(
    data = delta_long_spik %>%
      filter(metric == "Pr (Δ > 0)", !is.na(species_label)),
    aes(yintercept = 0.5),
    linetype = "dashed", color = "black"
  ) +
  geom_vline(
    data = yearly_prism_max_avg %>%
      filter(species!="aghy") %>% 
      mutate(species_label = case_when(
        species == "elvi" ~ "italic('Elymus virginicus')",
        species == "poau" ~ "italic('Poa autumnalis')"
      )),
    aes(xintercept = ppt_avg),
    linetype = "dotted",
    color = "#0072B2",
    size = 0.8
  ) +
  facet_grid(
    metric ~ species_label,
    scales = "free_y",
    labeller = labeller(species_label = label_parsed, metric = label_value)
  ) +
  scale_color_manual(values = c("Unfenced" = "#E69F00", "Fenced" = "#009E73")) +
  labs(x = "Precipitation (mm)", y = NULL, color = "Herbivore exclusion") +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.spacing.y = unit(0.0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.line.y = element_line(color = "black", size = 0.1),
    axis.line.x = element_line(color = "black", size = 0.1),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black", face = "plain"),
    strip.text.y = element_text(size = 10, color = "black", face = "plain"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
  )

dev.off()


