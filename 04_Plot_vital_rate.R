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

# Survival----
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

# Load model and data
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/x5q2v6fxnojb9msl0yhnj/fit_surv_abio_bio_endo.rds?rlkey=ojivjcdq1s30jf15is7pbezx8&dl=1"))

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
      bendoherbclim[, species_index] * endo * herb * clim +
      bclim2[, species_index] * clim^2
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

# Trim panel names
plot_data_survival$panel <- trimws(as.character(plot_data_survival$panel))
observed_data_survival$panel <- trimws(as.character(observed_data_survival$panel))

# 9. Plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/PrSurvival_diff.pdf",
  width = 5, height = 7
)

ggplot(plot_data_survival) +
  # Survival panel
  geom_line(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = exp(clim), y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_survival, panel == "Pr (survival)"),
    aes(x = exp(clim), y = y_plot_mean, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +
  
  # Δ panel
  geom_line(
    data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  geom_hline(
    data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0), linetype = "dashed", color = "black"
  ) +
  # Optional: add "0" label on lower panels
  # geom_text(
  #   data = subset(plot_data_survival, panel == "Δ (E+ - E-)"),
  #   aes(x = min(exp(clim)), y = 0, label = "0"),
  #   inherit.aes = FALSE,
  #   hjust = 1.1, vjust = 0.5,
  #   size = 2.5
  # ) +
  # Facets with heights
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
      # Agrostis hyemalis 
      panel == "Δ (E+ - E-)" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(
          breaks = 0,
          labels = 0,
          minor_breaks = NULL,
          limits = c(-0.3, 0.3),
          expand = c(0,0)
        ),
      
      # Elymus virginicus 
      panel == "Δ (E+ - E-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(
          breaks = 0,
          labels = 0,
          minor_breaks = NULL,
          limits = c(-0.3, 0.3),
          expand = c(0,0)
        ),
      
      # Poa autumnalis
      panel == "Δ (E+ - E-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(
          breaks = 0,
          labels = 0,
          minor_breaks = NULL,
          limits = c(-0.3, 0.3),
          expand = c(0,0)
        ),
      
      # Upper panels (same for all)
      panel == "Pr (survival)" ~ scale_y_continuous(expand = c(0,0))
    )
  )+
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.075, 0.27),
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
    strip.background = element_rect(color="black", fill="grey80", size=0.1, linetype="solid")
  )

dev.off()

# Growth----
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

fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/6zy5h2kdj56jkzm5ca8ok/fit_grow_abio_bio_endo.rds?rlkey=92hul5ki05hpbgguo9nbg6dd3&dl=1"))
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
      bendoherbclim[, species_index] * endo * herb * clim +
      bclim2[, species_index] * clim^2
    mu_preds  # assuming Gaussian model, not logistic
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
  group_by(plot, species, herb, clim, endo) %>%
  summarise(
    y_plot_mean = mean(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(panel = "Growth")

observed_data_grow$species <- factor(
  observed_data_grow$species,
  levels = c("1", "2", "3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

y_limits <- plot_data_grow %>%
  dplyr::filter(panel == "Δ (E+ - E-)") %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    ymin = min(lower_90, na.rm = TRUE),
    ymax = max(upper_90, na.rm = TRUE)
  )


# Plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_diff.pdf",
  width = 6, height = 9.5
)

ggplot(plot_data_grow) +
  # Upper panel: predicted growth
  geom_line(
    data = subset(plot_data_grow, panel == "Growth"),
    aes(x = exp(clim), y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_grow, panel == "Growth"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  # Observed points
  geom_point(
    data = subset(observed_data_grow, panel == "Growth"),
    aes(x = exp(clim), y = y_plot_mean, color = factor(endo)),
    size = 0.75,
    position = position_jitter(width = 0, height = 0.01)
  ) +
  # Lower panel: Δ(E+ − E−) differences
  geom_line(
    data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_grow, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
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
  # Dynamic y-axis scales
  ggh4x::facetted_pos_scales(
    y = list(
      # Lower panels (Δ(Growth)) – only 0 tick, no labels
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
      # Upper panels (Growth) – regular y-axis
      panel == "Growth" ~ scale_y_continuous(expand = c(0, 0))
    )
  ) +
  # Labels and theme
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("E-", "E+")) +
  theme_light() +
  theme(
    legend.position = c(0.4, 0.2),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 8),
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

# Flowering----
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

fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/pl6444lqmvl10s8ccxbsf/fit_flow_abio_bio_endo.rds?rlkey=x34lgk5q3n1hfryd6y2se2rm2&dl=1"))
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

# Prediction function for flowering
get_predictions_flow <- function(clim, endo, herb, species_index, posterior_samples_flow) {
  with(posterior_samples_flow, {
    # Linear predictor (log link)
    eta <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bendoherbclim[, species_index] * endo * herb * clim +
      bclim2[, species_index] * clim^2
    
    mu <- exp(eta)  # inverse log link for negative binomial
    mu
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

# Change panel name for upper panel
plot_data_flow <- plot_data_flow %>%
  mutate(panel = ifelse(panel == "Flowering", "#Inflorescences", panel))

observed_data_flow <- observed_data_flow %>%
  mutate(panel = ifelse(panel == "Flowering", "#Inflorescences", panel))

# Δ-panel limits for each species (same, only Δ panel)
y_limits_flow <- plot_data_flow %>%
  filter(panel == "Δ (E+ - E-)") %>%
  group_by(species) %>%
  summarise(
    ymin = min(lower_90, na.rm = TRUE),
    ymax = max(upper_90, na.rm = TRUE),
    .groups = "drop"
  )

# Plot with updated panel labels
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flowering_diff.pdf",
  width = 6, height = 12
)

ggplot(plot_data_flow) +
  # Upper panel: predicted flowering counts (#Inflorescences)
  geom_line(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = exp(clim), y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_flow, panel == "#Inflorescences"),
    aes(x = exp(clim), y = y_plot_mean, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +
  
  # Lower panel: Δ (E+ - E-) differences
  geom_line(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
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
          limits = c(y_limits_flow$ymin[y_limits_flow$species == "italic('Poa autumnalis')"], 50),
          expand = c(0, 0)
        ),
      
      # Upper panels – #Inflorescences custom limits per species
      panel == "#Inflorescences" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(limits = c(0, 20), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 7), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
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
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flowering_diff_v.pdf",
  width = 10, height = 5
)

ggplot(plot_data_flow) +
  geom_line(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = exp(clim), y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_flow, panel == "#Inflorescences"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_flow, panel == "#Inflorescences"),
    aes(x = exp(clim), y = y_plot_mean, color = factor(endo)),
    size = 0.75, position = position_jitter(width = 0, height = 0.01)
  ) +
  geom_line(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), y = mean), color = "black", size = 0.5
  ) +
  # Add a horizontal dashed line at y = 0 in the Δ-panel
  geom_hline(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0),
    color = "black",
    linetype = "dashed",
    size = 0.3,
    inherit.aes = FALSE
  )+
  geom_ribbon(
    data = subset(plot_data_flow, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
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
        scale_y_continuous(limits = c(0, 20), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 7), expand = c(0, 0)),
      panel == "#Inflorescences" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
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
    strip.text.x = element_text(size = 10, color = "black", face = "plain"),
    strip.text.y = element_text(size = 10, color = "black", face = "plain"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.1, linetype = "solid")
  )

dev.off()

# Spikelet----
demography_climate %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cum_ppt
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
    ppt = log(cum_ppt)
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


fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/fg562lkkl077j8qoegb8q/fit_spik_ppt.rds?rlkey=cbfyy7tq7tja3e9mxujv0scm6&dl=1"))
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
  plot=demography_spik_ppt$plot,
  # Your herbivory status data
  species = demography_spik_ppt$Spp,
  # Your species data
  y = demography_spik_ppt$y # Observed survival
)

# Compute mean spikelet per plot, grouped by endo
observed_plot_spk <- observed_spk %>%
  group_by(plot, species, herb, clim, endo) %>%  # include endo now
  summarise(
    y_plot_mean = mean(y, na.rm = TRUE),
    .groups = "drop"
  )

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spik_ppt.pdf",
  useDingbats = F,
  height = 5,
  width = 6.5
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
  geom_point(data = observed_plot_spk,
             aes(
               x = exp(clim),
               y = y_plot_mean,
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
  theme_bw() +
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

# Calculate differences (E+ - E-) for spikes ---
diff_spk <- pred_probspk_long_df %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Spik[endo == 1]) - mean(Pred_Spik[endo == 0]),
    .groups = "drop"
  )

diff_spk_ci <- diff_spk %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    mean = mean(diff),
    median = median(diff),
    .groups = "drop"
  ) %>%
  mutate(panel = "Δ (E+ - E-)")

# Prepare predicted panel if needed ---
cred_intervalspk <- cred_intervalspk %>%
  mutate(panel = "Predicted spikes")

# observed data panel
observed_spk <- observed_spk %>%
  mutate(panel = "Predicted spikes")

# Combine into one data frame
plot_spk_data <- bind_rows(cred_intervalspk, diff_spk_ci) %>%
  mutate(panel = factor(panel, levels = c("Predicted spikes", "Δ (E+ - E-)")))

# Plot only Δ (E+ − E-) ---
Cairo::CairoPDF("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spike_delta_only.pdf",
                width = 5.5, height = 5)

ggplot(subset(plot_spk_data, panel == "Δ (E+ - E-)")) +
  geom_line(aes(x = exp(clim), y = mean), color = "black", size = 1) +
  geom_ribbon(aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
              fill = "#9B6B96", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(species ~ herb, scales = "free_y",
             labeller = labeller(
               species = c("1"="ELVI", "2" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(x = "Precipitation (mm)",
       y = "Δ Spikelets (E+ - E-)") +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 13)
  )

dev.off()

library(dplyr)
library(ggplot2)
library(ggh4x)

# Optional: define panel heights
heights_spk <- c(
  "1.Predicted spikes" = 2.5, "1.Δ (E+ - E-)" = 2,
  "2.Predicted spikes" = 2.5, "2.Δ (E+ - E-)" = 2
)

# Ensure no trailing spaces in panel
plot_spk_data <- plot_spk_data %>%
  mutate(panel = trimws(panel))

observed_spk <- observed_spk %>%
  mutate(panel = trimws(panel))


# Update panel names
plot_spk_data <- plot_spk_data %>%
  mutate(panel = recode(panel, "Predicted spikes" = "Spikelets"))

observed_spk <- observed_plot_spk %>%
  mutate(panel = "Spikelets")  # match panel names

# Ensure panel is character
plot_spk_data$panel <- as.character(plot_spk_data$panel)

# Save PDF
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spik_ppt_diff.pdf",
  height = 7,
  width = 6.5,
  useDingbats = FALSE
)

ggplot(plot_spk_data) +
  # Predicted spikes panel
  geom_line(
    data = subset(plot_spk_data, panel == "Spikelets"),
    aes(x = exp(clim), y = mean, color = factor(endo), group = endo),
    size = 1
  ) +
  geom_ribbon(
    data = subset(plot_spk_data, panel == "Spikelets"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = observed_spk,
    aes(x = exp(clim), y = y_plot_mean, color = factor(endo)),
    size = 2.5,
    position = position_jitter(width = 0, height = 0.02)
  ) +
  
  # Δ panel
  geom_line(
    data = subset(plot_spk_data, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), y = mean),
    color = "black", size = 1
  ) +
  geom_ribbon(
    data = subset(plot_spk_data, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  geom_hline(
    data = subset(plot_spk_data, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0),
    linetype = "dashed", color = "black"
  ) +
  
  # Facets: species × panel stacked above herb
  facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = c("1" = "ELVI", "2" = "POAU"),
      herb = c("0" = "Unfenced", "1" = "Fenced")
    )
  ) +
  facetted_pos_scales(
    y = list(
      panel == "Predicted spikes" ~ scale_y_continuous(limits = c(0, max(plot_spk_data$mean + 1)), expand = c(0, 0)),
      panel == "Δ (E+ - E-)"    ~ scale_y_continuous(limits = c(-10, 10), expand = c(0, 0))
    )
  ) +
  
  # Labels and colors
  labs(x = "Precipitation (mm)", y = "# Spikelets", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) +
  
  theme_light() +
  theme(
    legend.position = c(0.08, 0.40),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    ggh4x.facet.nest.heights = heights_spk,
    panel.spacing.y = unit(0.00, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "grey40", face = "bold"),
    strip.text.y = element_text(size = 10, color = "grey40", face = "bold"),
    strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid")
  )

dev.off()

