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
library(xtable)
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
climate_max<-readRDS(url("https://www.dropbox.com/scl/fi/h8m15t64sxkkghgsnsypp/prism_edge_yr_means.rds?rlkey=ljjf16kqvvoe15qmhrqh8qfah&dl=1"))
climate_scaled<-readRDS(url("https://www.dropbox.com/scl/fi/irecsnoh3xrq6g8d5cysa/climate_site_scaled.rds?rlkey=63r7ugrtkuo5ncmqfoywbidps&dl=1"))

# yearly_prism_max <- climate_max %>%
#   group_by(species, year) %>%
#   summarise(
#     ppt_sum = sum(ppt, na.rm = TRUE),
#     tmean_mean = mean(tmean, na.rm = TRUE)
#   ) %>%
#   ungroup()
# yearly_prism_max_avg <- yearly_prism_max %>%
#   group_by(species) %>%
#   summarise(
#     ppt_avg = mean(ppt_sum, na.rm = TRUE),
#     tmean_avg = mean(tmean_mean, na.rm = TRUE)
#   ) %>%
#   ungroup()

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
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/mh5es9xqo4t608h12zg4q/fit_surv_abio_bio_endo_linear.rds?rlkey=akzrlhtqbrx3sut9h58aidp0v&dl=1"))

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
      
      # Two-way interactions
      bendoclim[, species_index] * endo * clim +
      bendoherb[, species_index] * endo * herb +
      bherbclim[, species_index] * herb * clim +
      
      # Three-way interaction
      bendoherbclim[, species_index] * endo * herb * clim 
    
    plogis(logit_preds)
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
  summarise(
    y_plot_mean = mean(y, na.rm = TRUE),
    n_obs = sum(!is.na(y)),  # counts non-missing observations
    .groups = "drop"
  ) %>%
  mutate(panel = "Pr (survival)")

# Differences (S+ - S-)
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
  mutate(panel = "Δ (S+ - S-)")

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
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo)
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.2, color = NA
  ) +
  geom_point(
    data = observed_data_survival,
    aes(
      x = climate_mm,
      y = pmin(pmax(y_plot_mean, 0.01), 0.99),  # clamp to 0.01–0.99
      color = factor(endo),
      size = n_obs
    ),
    alpha = 0.3,
    position = position_jitter(width = 0.05, height = 0),
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(0.5, 3)) +
  
  # Δ panel
  geom_line(
    data = subset(plot_data_survival, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, y = mean),
    color = "black",
    linewidth = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.6
  ) +
  geom_hline(
    data = subset(plot_data_survival, panel == "Δ (S+ - S-)"),
    aes(yintercept = 0),
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  
  # Facets
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Herbivory access", "1" = "Herbivory exclusion")
    )
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      panel == "Δ (S+ - S-)" ~ scale_y_continuous(
        limits = c(-0.3, 0.35),
        expand = c(0,0),
        breaks = 0,
        labels = 0,
        minor_breaks = NULL
      ),
      panel == "Pr (survival)" ~ scale_y_continuous(
        limits = c(0, 1),
        expand = c(0,0)
      )
    )
  ) +
  
  # Labels and colors
  labs(x = "Precipitation (mm)", y = "", color = "Symbiont", fill = "Symbiont") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("S-", "S+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("S-", "S+")) +
  
  # Theme
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    axis.line = element_line(color = "black", linewidth = 0.1),
    legend.position = c(0.415, 0.20),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.2, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", linewidth = 0.2)
  ) +
  
  # Add facet labels (a–f)
  geom_text(
    data = panel_labels,
    aes(x = 490, y = 0.7, label = label),
    fontface = "plain",
    size = 3.5,
    hjust = 0,
    inherit.aes = FALSE
  )
dev.off()

# Function to compute Δ(S+ − S−) for survival 
compute_delta_surv <- function(clim_val, posterior_samples_surv, herb_values = c(0,1)) {
  n_species <- dim(posterior_samples_surv$b0)[2]
  n_post <- dim(posterior_samples_surv$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      
      # Posterior predictions for S+ (endo present) and S− (endo absent)
      pred_Splus <- 1 / (1 + exp(-(
        posterior_samples_surv$b0[, sp] +
          posterior_samples_surv$bendo[, sp] * 1 +
          posterior_samples_surv$bherb[, sp] * h +
          posterior_samples_surv$bclim[, sp] * clim_val +
          
          # Two-way interactions
          posterior_samples_surv$bendoclim[, sp] * 1 * clim_val +
          posterior_samples_surv$bendoherb[, sp] * 1 * h +
          posterior_samples_surv$bherbclim[, sp] * h * clim_val +
          
          # Three-way interaction
          posterior_samples_surv$bendoherbclim[, sp] * 1 * h * clim_val
      )))
      
      pred_Sminus <- 1 / (1 + exp(-(
        posterior_samples_surv$b0[, sp] +
          posterior_samples_surv$bendo[, sp] * 0 +
          posterior_samples_surv$bherb[, sp] * h +
          posterior_samples_surv$bclim[, sp] * clim_val +
          
          # Two-way interactions
          posterior_samples_surv$bendoclim[, sp] * 0 * clim_val +
          posterior_samples_surv$bendoherb[, sp] * 0 * h +
          posterior_samples_surv$bherbclim[, sp] * h * clim_val +
          
          # Three-way interaction
          posterior_samples_surv$bendoherbclim[, sp] * 0 * h * clim_val
      )))
      
      delta <- pred_Splus - pred_Sminus
      
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

# Species-specific climate ranges for predictions
climate_range_per_species <- lapply(1:3, function(sp) {
  sp_clim <- demography_surv_ppt$clim[demography_surv_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})
names(climate_range_per_species) <- 1:3
# Compute Δ survival across species-specific ranges
delta_surv_species_range <- lapply(1:3, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_surv(cl, posterior_samples_survival, herb_values = c(0,1)) %>%
      filter(species == sp)
  }) %>% bind_rows()
}) %>% bind_rows()

# Summarize Δ for plotting
delta_surv_summary <- delta_surv_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  ) %>%
  mutate(
    species = factor(species, levels = 1:3,
                     labels = c("Agrostis hyemalis",
                                "Elymus virginicus",
                                "Poa autumnalis")),
    herb = factor(herb, levels = c(0,1),
                  labels = c("Herbivory access", "Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  ) %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()

for(sp in unique(delta_surv_summary$species)) {
  for(h in unique(delta_surv_summary$herb)) {
    df_sub <- delta_surv_summary %>%
      filter(species == sp, herb == h)
    
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}

delta_surv_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm = round(clim_mm, 0),
    median_delta = round(median_delta, 3),
    lower_90 = round(lower_90, 3),
    upper_90 = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

# Prepare long format for plotting
delta_long_surv <- delta_surv_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    "median_delta" = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Plot Δ survival
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/PrSurvival_diff_stat_species.pdf",
  width = 7, height = 6
)

ggplot(delta_long_surv, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  geom_hline(
    data = delta_long_surv %>% filter(metric == "Median Δ (S+ − S−)"),
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
  facet_grid(metric ~ species_label, scales = "free",
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Herbivory access" = "#E69F00", 
                                "Herbivory exclusion" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore treatment"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = "bottom",
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  )

dev.off()


# Select only the essential columns for reporting
# Map full species names to abbreviated form
species_abbrev <- c(
  "Agrostis hyemalis"   = "A. hyemalis",
  "Elymus virginicus"   = "E. virginicus",
  "Poa autumnalis"      = "P. autumnalis"
)

# Prepare data
delta_surv_latex <- delta_surv_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),  # italic for LaTeX
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm = clim_mm,
    Median_delta = median_delta,
    Lower_90 = lower_90,
    Upper_90 = upper_90,
    Prob_delta_gt0 = prob_delta_gt0
  )

# Convert to xtable with descriptive name
delta_surv_xt <- xtable(
  delta_surv_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)

# Custom column names for LaTeX header (pass LaTeX commands as text)
colnames(delta_surv_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",  # multi-line header in LaTeX
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)

# Print LaTeX table
print(
  delta_surv_xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  floating = FALSE
)

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

fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/mu3gpry42ad7fkfbtlvne/fit_grow_abio_bio_endo_linear.rds?rlkey=pz8uiqevdm5ogy7qvm369qyvb&dl=1"))
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
      
      # Two-way interactions
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bherbclim[, species_index] * herb * clim +
      
      # Three-way interaction
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
# Compute S+ − S− difference panels
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
  mutate(panel = "Δ (S+ - S-)")
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
    n_obs = sum(!is.na(y)),  # count of non-missing observations per plot
    .groups = "drop"
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
  dplyr::filter(panel == "Δ (S+ - S-)") %>%
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
  panel = "Growth",
  # manually tweak y positions
  ymax = c(1, 1, 1.25, 1.25, 3.3, 3.3)   # tweak last one slightly higher
)

# Plot
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_diff.pdf",
  width = 7, height = 9
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
    alpha = 0.2, color = NA
  ) +
  # Observed points
  geom_point(
    data = subset(observed_data_grow, panel == "Growth"),
    aes(
      x = climate_mm,
      y = y_plot_mean,
      color = factor(endo),
      size = n_obs  # map point size to sample size
    ),
    alpha = 0.3,
    position = position_jitter(width = 5, height = 0.05),
    show.legend = FALSE  # hide size legend
  ) +
  # optional: control relative sizes
  scale_size_continuous(range = c(0.5, 3))+

  # Lower panel: Δ(S+ − S−) differences
  geom_line(
    data = subset(plot_data_grow, panel == "Δ (S+ - S-)"),
    aes(x =climate_mm, y = mean), color = "black", size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_grow, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.6
  ) +
  geom_hline(
    data = subset(plot_data_grow, panel == "Δ (S+ - S-)"),
    aes(yintercept = 0), linetype = "dashed", color = "black"
  ) +
  # Facets
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Herbivory access", "1" = "Herbivory exclusion")
    )
  ) +
  # Dynamic y-axis scales
  ggh4x::facetted_pos_scales(
    y = list(
      # Lower panels (Δ(S+ − S−)) – only 0 tick, no labels
      panel == "Δ (S+ - S-)" & species == "italic('Agrostis hyemalis')" ~
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
      panel == "Δ (S+ - S-)" & species == "italic('Elymus virginicus')" ~
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
      panel == "Δ (S+ - S-)" & species == "italic('Poa autumnalis')" ~
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
        scale_y_continuous(limits = c(-2.5, 1), expand = c(0, 0)),
      panel == "Growth" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(-1.1, 1.25), expand = c(0, 0)),
      panel == "Growth" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(-4, 3), expand = c(0, 0))
    )
  ) +
  # Labels and theme
  labs(x = "Precipitation (mm)", y = "Log ratio of size", color = "Symbiont", fill = "Symbiont") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("S-", "S+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("S-", "S+")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = c(0.37, 0.23),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.2, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.ticks.y = element_line(color = "black", size = 0.2),
    #  axis.line.y = element_line(color = "black", size = 0.01),
    # axis.line.x = element_line(color = "black", size = 0.01),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  ) +
  # Add facet labels (a–f) only on top panels
  geom_text(
    data = panel_labels_grow,
    aes(x = 490, y = ymax * 0.80, label = label),
    hjust = 0,
    size = 3.5,
    inherit.aes = FALSE
  )

dev.off()

# Climate quantiles for grow 
#clim_quantiles_grow <- quantile(demography_grow_ppt$clim, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

#Function to compute Δ(S+ − S−) for grow 
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

      
      # Δ(S+ − S−)
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

# Species-specific climate ranges for predictions
climate_range_per_species <- lapply(1:3, function(sp) {
  sp_clim <- demography_grow_ppt$clim[demography_grow_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})
names(climate_range_per_species) <- 1:3

# Compute Δ growival across species-specific ranges
delta_grow_species_range <- lapply(1:3, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_grow(cl, posterior_samples_grow, herb_values = c(0,1)) %>%
      filter(species == sp)
  }) %>% bind_rows()
}) %>% bind_rows()

# Summarize Δ for plotting
delta_grow_summary <- delta_grow_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  ) %>%
  mutate(
    species = factor(species, levels = 1:3,
                     labels = c("Agrostis hyemalis",
                                "Elymus virginicus",
                                "Poa autumnalis")),
    herb = factor(herb, levels = c(0,1),
                  labels = c("Herbivory access", "Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  ) %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()

for(sp in unique(delta_grow_summary$species)) {
  for(h in unique(delta_grow_summary$herb)) {
    df_sub <- delta_grow_summary %>%
      filter(species == sp, herb == h)
    
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}

delta_grow_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm = round(clim_mm, 0),
    median_delta = round(median_delta, 3),
    lower_90 = round(lower_90, 3),
    upper_90 = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

# Prepare long format for plotting
delta_long_grow <- delta_grow_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    "median_delta" = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Plot Δ survival
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_diff_stat_species.pdf",
  width = 7, height = 6
)

ggplot(delta_long_grow, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  geom_hline(
    data = delta_long_grow %>% filter(metric == "Median Δ (S+ − S−)"),
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    data = delta_long_grow %>% filter(metric == "Pr (Δ > 0)"),
    aes(yintercept = 0.5),
    linetype = "dashed",
    color = "black"
  ) +
  facet_grid(metric ~ species_label, scales = "free",
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Herbivory access" = "#E69F00", 
                                "Herbivory exclusion" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore treatment"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = "bottom",
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  )

dev.off()


# Select only the essential columns for reporting
# Prepare data
delta_grow_latex <- delta_grow_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),  # italic for LaTeX
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm = clim_mm,
    Median_delta = median_delta,
    Lower_90 = lower_90,
    Upper_90 = upper_90,
    Prob_delta_gt0 = prob_delta_gt0
  )

# Convert to xtable with descriptive name
delta_grow_xt <- xtable(
  delta_grow_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)

# Custom column names for LaTeX header (pass LaTeX commands as text)
colnames(delta_grow_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",  # multi-line header in LaTeX
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)

# Print LaTeX table
print(
  delta_grow_xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  floating = FALSE
)

# Inflorescence----
demography_climate %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cum_ppt,ppt_scaled,inf_t1
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
    inf_t1 = inf_t1,
    ppt = ppt_scaled
  ) -> demography_climate_inf

## Separate each variable to use the same model stan
### Precipitation 
demography_inf_ppt <- list(
  nSpp = demography_climate_inf$Species %>% n_distinct(),
  nSite = demography_climate_inf$Site %>% n_distinct(),
  nsite_year = demography_climate_inf$site_year %>% n_distinct(),
  nPop = demography_climate_inf$Population %>% n_distinct(),
  nPlot = demography_climate_inf$site_species_plot %>% n_distinct(),
  Spp = demography_climate_inf$Species,
  site = demography_climate_inf$Site,
  site_year = demography_climate_inf$site_year,
  pop = demography_climate_inf$Population,
  plot = demography_climate_inf$site_species_plot,
  clim = as.vector(demography_climate_inf$ppt),
  endo = demography_climate_inf$Endo,
  herb = demography_climate_inf$Herbivory,
  size = demography_climate_inf$log_size_t0,
  y = demography_climate_inf$inf_t1,
  N = nrow(demography_climate_inf)
)

fit_inf_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/5lfkgq6d5a2vzx5h2t9yr/fit_inf_abio_bio_endo_linear.rds?rlkey=pfqvm7sg8un5c14slypn1aqa9&dl=1"))
predictions <- expand.grid(
  clim = seq(
    min(demography_inf_ppt$clim),
    max(demography_inf_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:3
)
# Extract posterior samples
posterior_samples_inf <- rstan::extract(fit_inf_ppt)
# Prediction function for infering
get_predictions_inf <- function(clim, endo, herb, species_index, posterior_samples_inf) {
  with(posterior_samples_inf, {
    # Linear predictor
    eta <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      
      # Two-way interactions
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bherbclim[, species_index] * herb * clim +
      
      # Three-way interaction
      bendoherbclim[, species_index] * endo * herb * clim
    
    # Negative binomial mean on response scale
    mu <- exp(eta)
    
    # Apply zero-inflation (global parameter, not species-specific)
    (1 - zi) * mu
  })
}

# Generate predictions
n_post_inf <- nrow(posterior_samples_inf$b0)
pred_matrix_inf <- matrix(NA, nrow = nrow(predictions), ncol = n_post_inf)

for (i in seq_len(nrow(predictions))) {
  pred_matrix_inf[i, ] <- get_predictions_inf(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples_inf
  )
}

# Combine predictions with original predictors
pred_inf_df <- cbind(predictions, as.data.frame(pred_matrix_inf))

# Pivot longer for tidy format
pred_inf_long <- pred_inf_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Posterior_Sample",
    values_to = "Pred_Inf"
  )

# Compute credible intervals per species × endo × herb × clim
cred_intervals_inf <- pred_inf_long %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Inf, 0.05),
    upper_90 = quantile(Pred_Inf, 0.95),
    median = quantile(Pred_Inf, 0.5),
    mean = mean(Pred_Inf),
    .groups = "drop"
  ) %>%
  mutate(panel = "Inflorescences")

# Compute Δ (S+ − S−) for infering
diff_df_inf <- pred_inf_long %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Inf[endo == 1]) - mean(Pred_Inf[endo == 0]),
    .groups = "drop"
  )

diff_ci_inf <- diff_df_inf %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    mean = mean(diff),
    .groups = "drop"
  ) %>%
  mutate(panel = "Δ (S+ - S-)")

# Combine credible intervals and difference panels
plot_data_inf <- bind_rows(cred_intervals_inf, diff_ci_inf)

# Relabel species for plotting in italics
plot_data_inf$species <- factor(
  plot_data_inf$species,
  levels = c("1","2","3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# Change panel name for upper panel
plot_data_inf <- plot_data_inf %>%
  mutate(panel = ifelse(panel == "Inflorescences", "Inflorescences", panel))

observed_data_inf <- demography_inf_ppt %>%  # or your dataset for inflorescences
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
    n_obs = sum(!is.na(y)),  # count of non-missing observations per plot
    .groups = "drop"
  ) %>%
  mutate(
    panel = "Inflorescences",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

# Relabel species (matching your plot_data_inf)
observed_data_inf$species <- factor(
  observed_data_inf$species,
  levels = c("1", "2", "3"),
  labels = c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  )
)

# observed_data_inf <- observed_data_inf %>%
#   mutate(panel = ifelse(panel == "Inflorescence", "#Inflorescences", panel))
# Back-transform climate for plotting
plot_data_inf <- plot_data_inf %>%
  mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

# Δ-panel limits for each species (same, only Δ panel)
y_limits_inf <- plot_data_inf %>%
  filter(panel == "Δ (S+ - S-)") %>%
  group_by(species) %>%
  summarise(
    ymin = min(lower_90, na.rm = TRUE),
    ymax = max(upper_90, na.rm = TRUE),
    .groups = "drop"
  )

panel_labels_inf <- data.frame(
  species = rep(c(
    "italic('Agrostis hyemalis')",
    "italic('Elymus virginicus')",
    "italic('Poa autumnalis')"
  ), each = 2),
  herb = rep(c(0, 1), times = 3),
  label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
  panel = "Inflorescences",
  ymax = rep(c(40, 8, 110), each = 2)   # 👈 from your facetted scales
)

#Plot with updated panel labels
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Inflorescence_diff_v.pdf",
  width = 7, height = 9
)
ggplot(plot_data_inf) +
  # Upper panel: predicted inflorescences
  geom_line(
    data = subset(plot_data_inf, panel == "Inflorescences"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    linewidth = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_inf, panel == "Inflorescences"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90, fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = subset(observed_data_inf, panel == "Inflorescences"),
    aes(
      x = climate_mm,
      y = y_plot_mean,
      color = factor(endo),
      size = n_obs
    ),
    alpha = 0.3,
    position = position_jitter(width = 0, height = 0),
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(0.5, 3)) +
  
  # Lower panel: Δ (S+ − S−)
  geom_line(
    data = subset(plot_data_inf, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, y = mean), color = "black", linewidth = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_inf, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.6
  ) +
  geom_hline(
    data = subset(plot_data_inf, panel == "Δ (S+ - S-)"),
    aes(yintercept = 0), linetype = "dashed", color = "black"
  ) +
  
  # Facets: species vertically, herb horizontally
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Herbivory access", "1" = "Herbivory exclusion")
    )
  ) +
  
  # Facetted scales
  ggh4x::facetted_pos_scales(
    y = list(
      # Δ panels
      panel == "Δ (S+ - S-)" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(
          breaks = 0, labels = 0, minor_breaks = NULL,
          limits = c(y_limits_inf$ymin[y_limits_inf$species == "italic('Agrostis hyemalis')"],
                     y_limits_inf$ymax[y_limits_inf$species == "italic('Agrostis hyemalis')"]),
          expand = c(0, 0)
        ),
      panel == "Δ (S+ - S-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(
          breaks = 0, labels = 0, minor_breaks = NULL,
          limits = c(-1.5, 2),
          expand = c(0, 0)
        ),
      panel == "Δ (S+ - S-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(
          breaks = 0, labels = 0, minor_breaks = NULL,
          limits = c(y_limits_inf$ymin[y_limits_inf$species == "italic('Poa autumnalis')"],
                     y_limits_inf$ymax[y_limits_inf$species == "italic('Poa autumnalis')"]),
          expand = c(0, 0)
        ),
      
      # Upper panels – #Inflorescences per species
      panel == "Inflorescences" & species == "italic('Agrostis hyemalis')" ~
        scale_y_continuous(limits = c(0, 40)),
      panel == "Inflorescences" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 8)),
      panel == "Inflorescences" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 110))
    )
  ) +
  
  # Labels and colors
  labs(x = "Precipitation (mm)", y = "", color = "Symbiont", fill = "Symbiont") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("S-", "S+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"), labels = c("S-", "S+")) +
  
  # Theme
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    axis.line = element_line(color = "black", linewidth = 0.1),
    legend.position = c(0.13, 0.26),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.2, "cm"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 6),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", linewidth = 0.2)
  ) +
  
  # Add facet labels (a–f) on top panels
  geom_text(
    data = panel_labels_inf,
    aes(x = 490, y = ymax * 0.93, label = label),
    hjust = 0,
    size = 3.5,
    inherit.aes = FALSE
  )
dev.off()

# Climate quantiles for infering
# Function to compute Δ(S+ − S−) for infering
compute_delta_inf <- function(clim, posterior_samples_inf, herb_values = c(0, 1)) {
  n_species <- dim(posterior_samples_inf$b0)[2]
  n_post <- dim(posterior_samples_inf$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      # Linear predictor for S+ (endo present)
      eta_Splus <- posterior_samples_inf$b0[, sp] +
        posterior_samples_inf$bendo[, sp] * 1 +
        posterior_samples_inf$bherb[, sp] * h +
        posterior_samples_inf$bclim[, sp] * clim +
        
        # Two-way interactions
        posterior_samples_inf$bendoclim[, sp] * 1 * clim +
        posterior_samples_inf$bendoherb[, sp] * 1 * h +
        posterior_samples_inf$bherbclim[, sp] * h * clim +
        
        # Three-way interaction
        posterior_samples_inf$bendoherbclim[, sp] * 1 * h * clim
      
      # Linear predictor for S− (endo absent)
      eta_Sminus <- posterior_samples_inf$b0[, sp] +
        posterior_samples_inf$bendo[, sp] * 0 +
        posterior_samples_inf$bherb[, sp] * h +
        posterior_samples_inf$bclim[, sp] * clim +
        
        # Two-way interactions
        posterior_samples_inf$bendoclim[, sp] * 0 * clim +
        posterior_samples_inf$bendoherb[, sp] * 0 * h +
        posterior_samples_inf$bherbclim[, sp] * h * clim +
        
        # Three-way interaction
        posterior_samples_inf$bendoherbclim[, sp] * 0 * h * clim
      
      # Predicted mean inflorescence counts (log link → exp)
      pred_Splus  <- exp(eta_Splus)
      pred_Sminus <- exp(eta_Sminus)
      
      # Δ(S+ − S−): difference in expected inflorescence count
      delta <- pred_Splus - pred_Sminus
      
      data.frame(
        species = sp,
        herb = h,
        clim = clim,
        Posterior_Sample = 1:n_post,
        delta = delta
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  return(result)
}
# Species-specific climate ranges for predictions
climate_range_per_species <- lapply(1:3, function(sp) {
  sp_clim <- demography_inf_ppt$clim[demography_inf_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})
names(climate_range_per_species) <- 1:3

# Compute Δ infival across species-specific ranges
delta_inf_species_range <- lapply(1:3, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_inf(cl, posterior_samples_inf, herb_values = c(0,1)) %>%
      filter(species == sp)
  }) %>% bind_rows()
}) %>% bind_rows()

# Summarize Δ for plotting
delta_inf_summary <- delta_inf_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  ) %>%
  mutate(
    species = factor(species, levels = 1:3,
                     labels = c("Agrostis hyemalis",
                                "Elymus virginicus",
                                "Poa autumnalis")),
    herb = factor(herb, levels = c(0,1),
                  labels = c("Herbivory access", "Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  ) %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()

for(sp in unique(delta_inf_summary$species)) {
  for(h in unique(delta_inf_summary$herb)) {
    df_sub <- delta_inf_summary %>%
      filter(species == sp, herb == h)
    
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}

delta_inf_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm = round(clim_mm, 0),
    median_delta = round(median_delta, 3),
    lower_90 = round(lower_90, 3),
    upper_90 = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

# Prepare long format for plotting
delta_long_inf <- delta_inf_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    "median_delta" = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Plot Δ survival
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Inflorescence_diff_stat_species.pdf",
  width = 7, height = 6
)

ggplot(delta_long_inf, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  geom_hline(
    data = delta_long_inf %>% filter(metric == "Median Δ (S+ − S−)"),
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    data = delta_long_inf %>% filter(metric == "Pr (Δ > 0)"),
    aes(yintercept = 0.5),
    linetype = "dashed",
    color = "black"
  ) +
  facet_grid(metric ~ species_label, scales = "free",
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Herbivory access" = "#E69F00", 
                                "Herbivory exclusion" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore treatment"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = "bottom",
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  )

dev.off()


# Select only the essential columns for reporting
# Prepare data
delta_inf_latex <- delta_inf_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),  # italic for LaTeX
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm = clim_mm,
    Median_delta = median_delta,
    Lower_90 = lower_90,
    Upper_90 = upper_90,
    Prob_delta_gt0 = prob_delta_gt0
  )

# Convert to xtable with descriptive name
delta_inf_xt <- xtable(
  delta_inf_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)

# Custom column names for LaTeX header (pass LaTeX commands as text)
colnames(delta_inf_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",  # multi-line header in LaTeX
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)

# Print LaTeX table
print(
  delta_inf_xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  floating = FALSE
)

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

fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/7ivmicuigz1pahg4vxa7q/fit_spik_abio_bio_endo_linear.rds?rlkey=8h3js8dnaue95ojom8evrkmkl&dl=1"))
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
    
    # Linear predictor
    eta <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      
      # Two-way interactions
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bherbclim[, species_index] * herb * clim +
      
      # Three-way interaction
      bendoherbclim[, species_index] * endo * herb * clim
    
    # Predicted mean spikelet counts (log link → exp)
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

# Δ (S+ - S-) calculation
delta_spik <- plot_data_spik %>%
  pivot_wider(names_from = endo, values_from = c(mean, lower_90, upper_90)) %>%
  mutate(
    mean = mean_1 - mean_0,
    lower_90 = lower_90_1 - upper_90_0,
    upper_90 = upper_90_1 - lower_90_0,
    panel = "Δ (S+ - S-)"
  ) %>%
  dplyr::select(species, clim, herb, mean, lower_90, upper_90, panel)

plot_data_spik <- plot_data_spik %>%
  mutate(panel = "Spikelets") %>%
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
  summarise(
    y_plot_mean = mean(y, na.rm = TRUE),
    n_obs = sum(!is.na(y)),   # sample size per plot
    .groups = "drop"
  ) %>%
  mutate(
    panel = "Spikelets",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

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
  filter(panel == "Δ (S+ - S-)") %>%
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
  panel = "Spikelets"   # only place labels on upper panels
)

# Plot (horizontal layout for paper figure) 
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spikelet_diff.pdf",
  width = 10, height = 6
)
ggplot(plot_data_spik) +
  # Upper panel: predicted spikelets
  geom_line(
    data = subset(plot_data_spik, panel == "Spikelets"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo),
    size = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_spik, panel == "Spikelets"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90,
        fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  
  # Lower panel: Δ (S+ - S-)
  geom_line(
    data = subset(plot_data_spik, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, y = mean),
    color = "black", size = 0.5
  ) +
  geom_hline(
    data = subset(plot_data_spik, panel == "Δ (S+ - S-)"),
    aes(yintercept = 0),
    color = "black", linetype = "dashed", size = 0.3
  ) +
  geom_point(
    data = subset(observed_data_spik, panel == "Spikelets"),
    aes(
      x = climate_mm,
      y = y_plot_mean,
      color = factor(endo),
      size = n_obs   # map point size to sample size
    ),
    alpha = 0.3,
    position = position_jitter(width = 5, height = 0.05),
    show.legend = FALSE   # hide size legend
  ) +
  scale_size_continuous(range = c(0.5, 3))+
  geom_ribbon(
    data = subset(plot_data_spik, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  
  # Facets
  ggh4x::facet_nested(
    panel ~ species + herb,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      species = label_parsed,
      herb = c("0" = "Herbivory access", "1" = "Herbivory exclusion")
    )
  ) +
  # Facetted y-scales
  ggh4x::facetted_pos_scales(
    y = list(
      panel == "Δ (S+ - S-)" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_spik$ymin[y_limits_spik$species == "italic('Elymus virginicus')"],
          y_limits_spik$ymax[y_limits_spik$species == "italic('Elymus virginicus')"]
        ), expand = c(0, 0)),
      panel == "Δ (S+ - S-)" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(breaks = 0, labels = 0, limits = c(
          y_limits_spik$ymin[y_limits_spik$species == "italic('Poa autumnalis')"],
          y_limits_spik$ymax[y_limits_spik$species == "italic('Poa autumnalis')"]
        ), expand = c(0, 0)),
      panel == "Spikelets" & species == "italic('Elymus virginicus')" ~
        scale_y_continuous(limits = c(0, 50), expand = c(0, 0)),
      panel == "Spikelets" & species == "italic('Poa autumnalis')" ~
        scale_y_continuous(limits = c(0, 60), expand = c(0, 0))
    )
  ) +
  labs(
    x = "Precipitation (mm)",
    y = "",
    color = "Symbiont",
    fill = "Symbiont"
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("S-", "S+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("S-", "S+")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = c(0.1, 0.88),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.2, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.ticks.y = element_line(color = "black", size = 0.2),
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 10, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  )+
  geom_text(
    data = panel_labels_spik,
    aes(x = 490, y = 42, label = label),
    fontface = "plain", size = 3.5, hjust = 0,
    inherit.aes = FALSE
  )
dev.off()

# Function to compute Δ(S+ − S−) for spike
compute_delta_spik <- function(clim, posterior_samples_spik, herb_values = c(0, 1)) {
  n_species <- dim(posterior_samples_spik$b0)[2]
  n_post <- dim(posterior_samples_spik$b0)[1]
  
  result <- lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      # Linear predictor for S+ (symbiont present)
      eta_Splus <- posterior_samples_spik$b0[, sp] +
        posterior_samples_spik$bendo[, sp] * 1 +
        posterior_samples_spik$bherb[, sp] * h +
        posterior_samples_spik$bclim[, sp] * clim +
        
        # Two-way interactions
        posterior_samples_spik$bendoclim[, sp] * clim * 1 +
        posterior_samples_spik$bendoherb[, sp] * 1 * h +
        posterior_samples_spik$bherbclim[, sp] * h * clim +
        
        # Three-way interaction
        posterior_samples_spik$bendoherbclim[, sp] * 1 * h * clim
      
      # Linear predictor for S− (symbiont absent)
      eta_Sminus <- posterior_samples_spik$b0[, sp] +
        posterior_samples_spik$bendo[, sp] * 0 +
        posterior_samples_spik$bherb[, sp] * h +
        posterior_samples_spik$bclim[, sp] * clim +
        
        # Two-way interactions
        posterior_samples_spik$bendoclim[, sp] * clim * 0 +
        posterior_samples_spik$bendoherb[, sp] * 0 * h +
        posterior_samples_spik$bherbclim[, sp] * h * clim +
        
        # Three-way interaction
        posterior_samples_spik$bendoherbclim[, sp] * 0 * h * clim
      
      # Predicted mean spike counts (log link → exp)
      pred_Splus  <- exp(eta_Splus)
      pred_Sminus <- exp(eta_Sminus)
      
      # Δ(S+ − S−): difference in expected spike number
      delta <- pred_Splus - pred_Sminus
      
      data.frame(
        species = sp,
        herb = h,
        clim = clim,
        Posterior_Sample = 1:n_post,
        delta = delta
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  return(result)
}
# Species-specific climate ranges for predictions
climate_range_per_species <- lapply(1:2, function(sp) {
  sp_clim <- demography_spik_ppt$clim[demography_spik_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})
names(climate_range_per_species) <- 1:2

# Compute Δ spikival across species-specific ranges
delta_spik_species_range <- lapply(1:2, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_spik(cl, posterior_samples_spik, herb_values = c(0,1)) %>%
      filter(species == sp)
  }) %>% bind_rows()
}) %>% bind_rows()

# Summarize Δ for plotting
delta_spik_summary <- delta_spik_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta = median(delta),
    lower_90 = quantile(delta, 0.05),
    upper_90 = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  ) %>%
  mutate(
    species = factor(species, levels = 1:3,
                     labels = c("Agrostis hyemalis",
                                "Elymus virginicus",
                                "Poa autumnalis")),
    herb = factor(herb, levels = c(0,1),
                  labels = c("Herbivory access", "Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  ) %>%
  mutate(clim_mm = exp(clim * ppt_sd + ppt_mean))

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()

for(sp in unique(delta_spik_summary$species)) {
  for(h in unique(delta_spik_summary$herb)) {
    df_sub <- delta_spik_summary %>%
      filter(species == sp, herb == h)
    
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}

delta_spik_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm = round(clim_mm, 0),
    median_delta = round(median_delta, 3),
    lower_90 = round(lower_90, 3),
    upper_90 = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

# Prepare long format for plotting
delta_long_spik <- delta_spik_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(
    cols = c(median_delta, prob_delta_gt0),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    "median_delta" = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')",
      species == "Poa autumnalis" ~ "italic('Poa autumnalis')"
    )
  )

# Plot Δ survival
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spikelet_diff_stat.pdf",
  width = 7, height = 6
)

ggplot(delta_long_spik, aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 0.5) +
  geom_hline(
    data = delta_long_spik %>% filter(metric == "Median Δ (S+ − S−)"),
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    data = delta_long_spik %>% filter(metric == "Pr (Δ > 0)"),
    aes(yintercept = 0.5),
    linetype = "dashed",
    color = "black"
  ) +
  facet_grid(metric ~ species_label, scales = "free",
             labeller = labeller(
               species_label = label_parsed,
               metric = label_value
             )) +
  scale_color_manual(values = c("Herbivory access" = "#E69F00", 
                                "Herbivory exclusion" = "#009E73")) +
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore treatment"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = "bottom",
    text = element_text(family = "Arial"),
    strip.text.x = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  )

dev.off()


# Select only the essential columns for reporting
# Prepare data
delta_spik_latex <- delta_spik_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),  # italic for LaTeX
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm = clim_mm,
    Median_delta = median_delta,
    Lower_90 = lower_90,
    Upper_90 = upper_90,
    Prob_delta_gt0 = prob_delta_gt0
  )

# Convert to xtable with descriptive name
delta_spik_xt <- xtable(
  delta_spik_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)

# Custom column names for LaTeX header (pass LaTeX commands as text)
colnames(delta_spik_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",  # multi-line header in LaTeX
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)

# Print LaTeX table
print(
  delta_spik_xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  floating = FALSE
)

delta_long_surv$trait <- "Survival"
delta_long_grow$trait <- "Growth"
delta_long_inf$trait <- "Inflorescence"

# Combine all three datasets
delta_long_all <- dplyr::bind_rows(delta_long_surv, delta_long_grow, delta_long_inf)
delta_long_all <- delta_long_all %>%
  mutate(
    species_label = case_when(
      species == "Agrostis hyemalis"   ~ "italic('A. hyemalis')",
      species == "Elymus virginicus"   ~ "italic('E. virginicus')",
      species == "Poa autumnalis"      ~ "italic('P. autumnalis')",
      TRUE ~ as.character(species)
    )
  )

# climate_max <- climate_max %>%
#   mutate(species_label = case_when(
#     species == "aghy" ~ "italic('A. hyemalis')",
#     species == "elvi" ~ "italic('E. virginicus')",
#     species == "poau" ~ "italic('P. autumnalis')"
#   ))

# Plot
# Plot
# p_all <- ggplot(delta_long_all, aes(x = clim_mm, y = value, color = herb, group = herb)) +
#   geom_line(size = 1) +
#   geom_hline(
#     data = delta_long_all %>% filter(metric == "Median Δ (E+ − E−)"),
#     aes(yintercept = 0),
#     linetype = "dashed",
#     color = "black"
#   ) +
#   geom_hline(
#     data = delta_long_all %>% filter(metric == "Pr (Δ > 0)"), 
#     aes(yintercept = 0.5), 
#     linetype = "dashed", 
#     color = "black"
#   ) +
#   geom_vline(
#     data = climate_max,
#     aes(xintercept = mean_annual_ppt),
#     linetype = "dashed",
#     color = "#0072B2",
#     size = 0.8
#   ) +
#   facet_grid(
#     metric ~ species_label + trait,
#     scales = "free_y",
#     labeller = labeller(
#       species_label = label_parsed,
#       metric = label_value,
#       trait = label_value
#     )
#   ) +
#   scale_color_manual(values = c(
#     "Herbivory access" = "#E69F00", 
#     "Herbivory exclusion" = "#009E73"
#   )) +
#   scale_fill_manual(values = c(
#     "Herbivory access" = "#E69F00", 
#     "Herbivory exclusion" = "#009E73"
#   )) +
#   labs(
#     x = "Precipitation (mm)", 
#     y = NULL, 
#     color = "Herbivore treatment"
#   ) +
#   theme_classic(base_size = 18, base_family = "Arial") +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#     axis.line = element_line(color = "black", size = 0.1),
#     legend.position = "bottom",
#     legend.title = element_text(size = 30),
#     legend.text  = element_text(size = 30),
#     axis.title.x = element_text(size = 28),
#     axis.text    = element_text(size = 20),
#     strip.text.x = element_text(size = 32, color = "black"),
#     strip.text.y = element_text(size = 32, color = "black"),
#     #  axis.line.y = element_line(color = "black", size = 0.01),
#     # axis.line.x = element_line(color = "black", size = 0.01),
#     text = element_text(family = "Arial"),
#     strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
#   ) 
# 
# 
# 
# Cairo::CairoPDF("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/All_traits_diff_stat.pdf", width = 30, height = 15)
# print(p_all)
# dev.off()

# Filter only the lower panel (Pr (Δ > 0))
p_lower <- delta_long_all %>%
  filter(metric == "Pr (Δ > 0)") %>%
  ggplot(aes(x = clim_mm, y = value, color = herb, group = herb)) +
  geom_line(size = 1) +
  geom_hline(
    yintercept = 0.5, 
    linetype = "dashed", 
    color = "grey50"
  ) +
  facet_grid(
    . ~ species_label + trait,
    scales = "free",
    labeller = labeller(
      species_label = label_parsed,
      trait = label_value
    )
  ) +
  scale_color_manual(values = c(
    "Herbivory access" = "#E69F00", 
    "Herbivory exclusion" = "#009E73"
  )) +
  labs(
    x = "Precipitation (mm)", 
    y = "P(Δ > 0)", 
    color = "Herbivore treatment"
  ) +
  theme_classic(base_size = 18, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.line = element_line(color = "black", size = 0.1),
    legend.position = c(0.998, 0.95),       # top-right inside last panel
    legend.justification = c(1, 1),        # align top-right corner
    legend.background = element_rect(fill = alpha('white', 0.7), color = "black"),
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.text    = element_text(size = 20),
    strip.text.x = element_text(size = 32, color = "black"),
    text = element_text(family = "Arial"),
    strip.background = element_rect(color = "black", fill = "grey80", size = 0.2)
  )

Cairo::CairoPDF("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/All_traits_diff_stat_lower.pdf", width = 34, height = 12)
print(p_lower)
dev.off()
