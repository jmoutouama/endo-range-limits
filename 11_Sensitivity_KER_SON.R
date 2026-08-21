# Project: Endo range limits
# Purpose: Sensitivity analysis for reviewer — refit the full Bayesian models
#          (survival, growth, inflorescence, spikelet) EXCLUDING sites KER and
#          SON, then compare Delta (S+ - S-), the endophyte x precipitation
#          interaction, Pr(Delta>0), and 90% CrIs against the full-data models.
#
#          Primary analysis is NOT changed. This is purely a supplementary
#          check: "are the results robust to excluding KER and SON?"
#
# Mirrors data prep / model specs from 03_Model_building.R and the
# compute_delta_* link functions from 05_Plot_vital_rate.R exactly, so the
# comparison is apples-to-apples.
#
# Authors: Jacob Moutouama
# Date last modified (2026-08-18):

rm(list = ls())

# ── Packages ────────────────────────────────────────────────────────────────
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(bayesplot)
library(ggpubr)
library(Cairo)

# ── Shared conventions (match 05_Plot_vital_rate.R / 09_Plot_effect_size.R) ──
MODEL_COLORS <- c("Full data" = "grey35", "Excl. KER/SON" = "#d95f02")

SPP3_LABELS <- c("1" = "A. hyemalis", "2" = "E. virginicus", "3" = "P. autumnalis")
SPP2_LABELS <- c("1" = "E. virginicus", "2" = "P. autumnalis")

HERB_LABELS <- c("0" = "Herbivory access", "1" = "Herbivory exclusion")

sens_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 9.5, face = "bold"),
    strip.background   = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    panel.border       = element_rect(color = "black", linewidth = 0.3),
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title         = element_text(size = 9.5),
    axis.text          = element_text(size = 7.5),
    text               = element_text(family = "Arial")
  )


# ══════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA & EXCLUDE KER / SON
# ══════════════════════════════════════════════════════════════════════════
# Same object produced/saved at the end of the data-wrangling section of
# 03_Model_building.R
demography_climate_full <- readRDS(url(
  "https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1"
))

# Precipitation back-transform constants (ppt_scaled = (log(cum_ppt) - mean)/sd)
# — needed to plot Delta on the natural mm axis, as in 05_Plot_vital_rate.R
climate_scaled_ref <- readRDS(url(
  "https://www.dropbox.com/scl/fi/irecsnoh3xrq6g8d5cysa/climate_site_scaled.rds?rlkey=63r7ugrtkuo5ncmqfoywbidps&dl=1"
))
ppt_mean <- mean(climate_scaled_ref$ppt_log)
ppt_sd   <- sd(climate_scaled_ref$ppt_log)
ppt_mm   <- function(ppt_scaled_val) exp(ppt_scaled_val * ppt_sd + ppt_mean)
EXCLUDE_SITES <- c("KER", "SON")

demography_climate_sens <- demography_climate_full %>%
  filter(!Site %in% EXCLUDE_SITES)

cat("Full data:", nrow(demography_climate_full), "rows,",
    n_distinct(demography_climate_full$Site), "sites\n")
cat("Sensitivity data (excl. KER/SON):", nrow(demography_climate_sens), "rows,",
    n_distinct(demography_climate_sens$Site), "sites\n")

sim_pars <- list(
  warmup  = 1000,
  iter    = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains  = 3
)

stan_dir <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/stan"
out_dir  <- "/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output"

# ══════════════════════════════════════════════════════════════════════════
# 2. BUILD STAN DATA LISTS (sensitivity data) — same recipe as 03, applied
#    to demography_climate_sens instead of the full data
# ══════════════════════════════════════════════════════════════════════════

## Survival ------------------------------------------------------------------
demography_climate_surv_sens <- demography_climate_sens %>%
  filter(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, site_species_plot, site_year, Endo, Herbivory,
    tiller_t, surv1, cum_ppt, ppt_scaled
  ) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_year         = as.integer(factor(site_year)),
    site_species_plot = as.integer(factor(site_species_plot)),
    log_size_t0       = ppt_scaled,
    surv_t1           = as.integer(surv1),
    ppt               = ppt_scaled
  )

demography_surv_ppt_sens <- list(
  nSpp       = n_distinct(demography_climate_surv_sens$Species),
  nSite      = n_distinct(demography_climate_surv_sens$Site),
  nsite_year = n_distinct(demography_climate_surv_sens$site_year),
  nPop       = n_distinct(demography_climate_surv_sens$Population),
  nPlot      = n_distinct(demography_climate_surv_sens$site_species_plot),
  Spp        = demography_climate_surv_sens$Species,
  site       = demography_climate_surv_sens$Site,
  site_year  = demography_climate_surv_sens$site_year,
  pop        = demography_climate_surv_sens$Population,
  plot       = demography_climate_surv_sens$site_species_plot,
  clim       = as.vector(demography_climate_surv_sens$ppt),
  clim2      = as.vector(demography_climate_surv_sens$ppt^2),
  endo       = demography_climate_surv_sens$Endo,
  herb       = demography_climate_surv_sens$Herbivory,
  size       = demography_climate_surv_sens$log_size_t0,
  y          = demography_climate_surv_sens$surv_t1,
  N          = nrow(demography_climate_surv_sens)
)

## Growth ---------------------------------------------------------------------
demography_climate_grow_sens <- demography_climate_sens %>%
  filter(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot, site_year, site_species_plot, Endo, Herbivory,
    tiller_t, grow, cum_ppt, ppt_scaled
  ) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    site_year         = as.integer(factor(site_year)),
    Endo              = as.integer(Endo),
    Herbivory         = as.integer(Herbivory),
    log_size_t0       = log(tiller_t),
    ppt               = ppt_scaled
  )

demography_grow_ppt_sens <- list(
  nSpp       = n_distinct(demography_climate_grow_sens$Species),
  nSite      = n_distinct(demography_climate_grow_sens$Site),
  nsite_year = n_distinct(demography_climate_grow_sens$site_year),
  nPop       = n_distinct(demography_climate_grow_sens$Population),
  nPlot      = n_distinct(demography_climate_grow_sens$site_species_plot),
  Spp        = demography_climate_grow_sens$Species,
  site       = demography_climate_grow_sens$Site,
  site_year  = demography_climate_grow_sens$site_year,
  pop        = demography_climate_grow_sens$Population,
  plot       = demography_climate_grow_sens$site_species_plot,
  clim       = as.vector(demography_climate_grow_sens$ppt),
  clim2      = as.vector(demography_climate_grow_sens$ppt^2),
  endo       = demography_climate_grow_sens$Endo,
  herb       = demography_climate_grow_sens$Herbivory,
  size       = demography_climate_grow_sens$log_size_t0,
  y          = demography_climate_grow_sens$grow,
  N          = nrow(demography_climate_grow_sens)
)

## Inflorescence ----------------------------------------------------------------
demography_climate_inf_sens <- demography_climate_sens %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cum_ppt, ppt_scaled
  ) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    site_year         = as.integer(factor(site_year)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo              = as.integer(Endo),
    Herbivory         = as.integer(Herbivory),
    log_size_t0       = log(tiller_t),
    ppt               = ppt_scaled
  )

demography_inf_ppt_sens <- list(
  nSpp       = n_distinct(demography_climate_inf_sens$Species),
  nSite      = n_distinct(demography_climate_inf_sens$Site),
  nsite_year = n_distinct(demography_climate_inf_sens$site_year),
  nPop       = n_distinct(demography_climate_inf_sens$Population),
  nPlot      = n_distinct(demography_climate_inf_sens$site_species_plot),
  Spp        = demography_climate_inf_sens$Species,
  site       = demography_climate_inf_sens$Site,
  site_year  = demography_climate_inf_sens$site_year,
  pop        = demography_climate_inf_sens$Population,
  plot       = demography_climate_inf_sens$site_species_plot,
  clim       = as.vector(demography_climate_inf_sens$ppt),
  clim2      = as.vector(demography_climate_inf_sens$ppt^2),
  endo       = demography_climate_inf_sens$Endo,
  herb       = demography_climate_inf_sens$Herbivory,
  size       = demography_climate_inf_sens$log_size_t0,
  y          = demography_climate_inf_sens$inf_t1,
  N          = nrow(demography_climate_inf_sens)
)

## Spikelet (ELVI & POAU only) --------------------------------------------------
demography_climate_spik_sens <- demography_climate_sens %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  filter(tiller_t1 > 0, inf_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cum_ppt, ppt_scaled, inf_t1
  ) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    site_year         = as.integer(factor(site_year)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo              = as.integer(Endo),
    Herbivory         = as.integer(Herbivory),
    log_size_t0       = log(tiller_t),
    ppt               = ppt_scaled
  )

demography_spik_ppt_sens <- list(
  nSpp       = n_distinct(demography_climate_spik_sens$Species),
  nSite      = n_distinct(demography_climate_spik_sens$Site),
  nsite_year = n_distinct(demography_climate_spik_sens$site_year),
  nPop       = n_distinct(demography_climate_spik_sens$Population),
  nPlot      = n_distinct(demography_climate_spik_sens$site_species_plot),
  Spp        = demography_climate_spik_sens$Species,
  site       = demography_climate_spik_sens$Site,
  site_year  = demography_climate_spik_sens$site_year,
  pop        = demography_climate_spik_sens$Population,
  plot       = demography_climate_spik_sens$site_species_plot,
  clim       = as.vector(demography_climate_spik_sens$ppt),
  clim2      = as.vector(demography_climate_spik_sens$ppt^2),
  endo       = demography_climate_spik_sens$Endo,
  herb       = demography_climate_spik_sens$Herbivory,
  size       = demography_climate_spik_sens$log_size_t0,
  y          = demography_climate_spik_sens$spikelet_t1,
  N          = nrow(demography_climate_spik_sens)
)

# ══════════════════════════════════════════════════════════════════════════
# 3. REFIT THE SAME (LINEAR / "*_l.stan") MODELS ON THE SENSITIVITY DATA
#    — identical stan files, sim_pars, and seed as the primary analysis
# ══════════════════════════════════════════════════════════════════════════

# fit_surv_sens <- stan(
#   file    = file.path(stan_dir, "survival_l.stan"),
#   data    = demography_surv_ppt_sens,
#   warmup  = sim_pars$warmup, iter = sim_pars$iter,
#   chains  = sim_pars$chains, control = sim_pars$control, seed = 13
# )
# 
# fit_grow_sens <- stan(
#   file    = file.path(stan_dir, "growth_l.stan"),
#   data    = demography_grow_ppt_sens,
#   warmup  = sim_pars$warmup, iter = sim_pars$iter,
#   chains  = sim_pars$chains, control = sim_pars$control, seed = 13
# )
# 
# fit_inf_sens <- stan(
#   file    = file.path(stan_dir, "inflorescence_l.stan"),
#   data    = demography_inf_ppt_sens,
#   warmup  = sim_pars$warmup, iter = sim_pars$iter,
#   chains  = sim_pars$chains, control = sim_pars$control, seed = 13
# )
# 
# fit_spik_sens <- stan(
#   file    = file.path(stan_dir, "spikelet_l.stan"),
#   data    = demography_spik_ppt_sens,
#   warmup  = sim_pars$warmup, iter = sim_pars$iter,
#   chains  = sim_pars$chains, control = sim_pars$control, seed = 13
# )

# Save the fits so you don't have to refit if the session dies
# saveRDS(fit_surv_sens, file.path(out_dir, "fit_surv_sens_noKERSON.rds"))
# saveRDS(fit_grow_sens, file.path(out_dir, "fit_grow_sens_noKERSON.rds"))
# saveRDS(fit_inf_sens,  file.path(out_dir, "fit_inf_sens_noKERSON.rds"))
# saveRDS(fit_spik_sens, file.path(out_dir, "fit_spik_sens_noKERSON.rds"))

fit_surv_sens <- readRDS(url("https://www.dropbox.com/scl/fi/nu2azmsok9mztqqdjjj99/fit_surv_sens_noKERSON.rds?rlkey=xbmxiz0j7woo71fv954jybz2u&dl=1"))
fit_grow_sens <- readRDS(url("https://www.dropbox.com/scl/fi/e0verxpbqgq0l4ovb5a63/fit_grow_sens_noKERSON.rds?rlkey=0to4b7zkpsog6svnq86kdc4ds&dl=1"))
fit_inf_sens  <- readRDS(url("https://www.dropbox.com/scl/fi/oufjkibr5uk2x1y43ppzh/fit_inf_sens_noKERSON.rds?rlkey=yjrve44p1saofezvgsn29w46x&dl=1"))
fit_spik_sens <- readRDS(url("https://www.dropbox.com/scl/fi/16lhmj4t0x8ku5ut35nx8/fit_spik_sens_noKERSON.rds?rlkey=51x4r2rm9ya28wr328xy8c2by&dl=1"))

# ══════════════════════════════════════════════════════════════════════════
# 4. CONVERGENCE DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════════
check_convergence <- function(fit, label, pars = c("b0", "bendo", "bclim",
                                                     "bendoclim", "bherb",
                                                     "bendoherb", "bherbclim",
                                                     "bendoherbclim")) {
  s <- summary(fit)$summary
  s <- s[grepl(paste0("^(", paste(pars, collapse = "|"), ")\\["), rownames(s)), ]
  cat("\n===", label, "===\n")
  cat("Max Rhat:  ", round(max(s[, "Rhat"], na.rm = TRUE), 3), "\n")
  cat("Min n_eff: ", round(min(s[, "n_eff"], na.rm = TRUE), 0), "\n")
  divergent <- rstan::get_num_divergent(fit)
  cat("Divergent transitions:", divergent, "\n")
  invisible(s)
}

check_convergence(fit_surv_sens, "Survival (sensitivity)")
check_convergence(fit_grow_sens, "Growth (sensitivity)")
check_convergence(fit_inf_sens,  "Inflorescence (sensitivity)")
check_convergence(fit_spik_sens, "Spikelet (sensitivity)")

# Trace plots (same coefficients as 03_Model_building.R) — inspect visually
bayesplot::mcmc_trace(
  as.array(fit_surv_sens),
  pars = c("b0[1]", "b0[2]", "b0[3]",
           "bendo[1]", "bendo[2]", "bendo[3]",
           "bendoclim[1]", "bendoclim[2]", "bendoclim[3]")
) + theme_bw()

# ══════════════════════════════════════════════════════════════════════════
# 5. EXTRACT DELTA (S+ - S-), Pr(DELTA>0), 90% CrI — SAME LINK FUNCTIONS AS
#    compute_delta_surv / compute_delta_grow / compute_delta_inf /
#    compute_delta_spik IN 05_Plot_vital_rate.R
# ══════════════════════════════════════════════════════════════════════════

compute_delta_surv <- function(clim_val, ps, herb_values = c(0, 1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      pred_Splus <- plogis(
        ps$b0[,sp] + ps$bendo[,sp]*1 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
          ps$bendoclim[,sp]*1*clim_val + ps$bendoherb[,sp]*1*h +
          ps$bherbclim[,sp]*h*clim_val + ps$bendoherbclim[,sp]*1*h*clim_val
      )
      pred_Sminus <- plogis(
        ps$b0[,sp] + ps$bendo[,sp]*0 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
          ps$bendoclim[,sp]*0*clim_val + ps$bendoherb[,sp]*0*h +
          ps$bherbclim[,sp]*h*clim_val + ps$bendoherbclim[,sp]*0*h*clim_val
      )
      data.frame(species = sp, herb = h, clim = clim_val,
                 Posterior_Sample = 1:n_post, delta = pred_Splus - pred_Sminus)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

compute_delta_grow <- function(clim_val, ps, herb_values = c(0, 1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      pred_Eplus  <- ps$b0[,sp] + ps$bendo[,sp]*1 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
        ps$bendoclim[,sp]*clim_val*1 + ps$bendoherb[,sp]*1*h + ps$bendoherbclim[,sp]*1*h*clim_val
      pred_Eminus <- ps$b0[,sp] + ps$bendo[,sp]*0 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
        ps$bendoclim[,sp]*clim_val*0 + ps$bendoherb[,sp]*0*h + ps$bendoherbclim[,sp]*0*h*clim_val
      data.frame(species = sp, herb = h, clim = clim_val,
                 Posterior_Sample = 1:n_post, delta = pred_Eplus - pred_Eminus)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

compute_delta_inf <- function(clim, ps, herb_values = c(0, 1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      eta_p <- ps$b0[,sp] + ps$bendo[,sp]*1 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim +
        ps$bendoclim[,sp]*1*clim + ps$bendoherb[,sp]*1*h +
        ps$bherbclim[,sp]*h*clim + ps$bendoherbclim[,sp]*1*h*clim
      eta_m <- ps$b0[,sp] + ps$bendo[,sp]*0 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim +
        ps$bendoclim[,sp]*0*clim + ps$bendoherb[,sp]*0*h +
        ps$bherbclim[,sp]*h*clim + ps$bendoherbclim[,sp]*0*h*clim
      data.frame(species = sp, herb = h, clim = clim,
                 Posterior_Sample = 1:n_post, delta = exp(eta_p) - exp(eta_m))
    }) %>% bind_rows()
  }) %>% bind_rows()
}

compute_delta_spik <- function(clim, ps, herb_values = c(0, 1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      eta_p <- ps$b0[,sp] + ps$bendo[,sp]*1 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim +
        ps$bendoclim[,sp]*clim*1 + ps$bendoherb[,sp]*1*h +
        ps$bherbclim[,sp]*h*clim + ps$bendoherbclim[,sp]*1*h*clim
      eta_m <- ps$b0[,sp] + ps$bendo[,sp]*0 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim +
        ps$bendoclim[,sp]*clim*0 + ps$bendoherb[,sp]*0*h +
        ps$bherbclim[,sp]*h*clim + ps$bendoherbclim[,sp]*0*h*clim
      data.frame(species = sp, herb = h, clim = clim,
                 Posterior_Sample = 1:n_post, delta = exp(eta_p) - exp(eta_m))
    }) %>% bind_rows()
  }) %>% bind_rows()
}

# Helper: sweep across each species' OWN climate range (as in 05) and
# summarise median / 90% CrI / Pr(Delta>0)
summarise_delta_over_clim <- function(compute_fn, ps, clim_vec_by_species,
                                       species_ids, model_label) {
  lapply(species_ids, function(sp) {
    lapply(clim_vec_by_species[[as.character(sp)]], function(cl) {
      compute_fn(cl, ps) %>% filter(species == sp)
    }) %>% bind_rows()
  }) %>%
    bind_rows() %>%
    group_by(species, herb, clim) %>%
    summarise(
      median_delta   = median(delta),
      lower_90       = quantile(delta, 0.05),
      upper_90       = quantile(delta, 0.95),
      prob_delta_gt0 = mean(delta > 0),
      .groups = "drop"
    ) %>%
    mutate(model = model_label)
}

# ── Posteriors: sensitivity fits ────────────────────────────────────────────
posterior_surv_sens <- rstan::extract(fit_surv_sens)
posterior_grow_sens <- rstan::extract(fit_grow_sens)
posterior_inf_sens  <- rstan::extract(fit_inf_sens)
posterior_spik_sens <- rstan::extract(fit_spik_sens)

# ── Posteriors: FULL-data fits (primary analysis, loaded from Dropbox) ─────
fit_surv_full <- readRDS(url("https://www.dropbox.com/scl/fi/mh5es9xqo4t608h12zg4q/fit_surv_abio_bio_endo_linear.rds?rlkey=akzrlhtqbrx3sut9h58aidp0v&dl=1"))
fit_grow_full <- readRDS(url("https://www.dropbox.com/scl/fi/mu3gpry42ad7fkfbtlvne/fit_grow_abio_bio_endo_linear.rds?rlkey=pz8uiqevdm5ogy7qvm369qyvb&dl=1"))
fit_inf_full  <- readRDS(url("https://www.dropbox.com/scl/fi/5lfkgq6d5a2vzx5h2t9yr/fit_inf_abio_bio_endo_linear.rds?rlkey=pfqvm7sg8un5c14slypn1aqa9&dl=1"))
fit_spik_full <- readRDS(url("https://www.dropbox.com/scl/fi/7ivmicuigz1pahg4vxa7q/fit_spik_abio_bio_endo_linear.rds?rlkey=8h3js8dnaue95ojom8evrkmkl&dl=1"))

posterior_surv_full <- rstan::extract(fit_surv_full)
posterior_grow_full <- rstan::extract(fit_grow_full)
posterior_inf_full  <- rstan::extract(fit_inf_full)
posterior_spik_full <- rstan::extract(fit_spik_full)

# NOTE on climate scaling: ppt_scaled is standardized on the FULL dataset in
# 03_Model_building.R, and we deliberately keep that same scaling for the
# sensitivity fit (do NOT re-standardize on the subset) so Delta is compared
# on an identical climate axis between the two models.

species_ids_3 <- 1:3
species_ids_2 <- 1:2

clim_range_surv <- lapply(species_ids_3, function(sp)
  seq(min(demography_surv_ppt_sens$clim[demography_surv_ppt_sens$Spp == sp]),
      max(demography_surv_ppt_sens$clim[demography_surv_ppt_sens$Spp == sp]),
      length.out = 30)) %>% setNames(as.character(species_ids_3))

clim_range_grow <- lapply(species_ids_3, function(sp)
  seq(min(demography_grow_ppt_sens$clim[demography_grow_ppt_sens$Spp == sp]),
      max(demography_grow_ppt_sens$clim[demography_grow_ppt_sens$Spp == sp]),
      length.out = 30)) %>% setNames(as.character(species_ids_3))

clim_range_inf <- lapply(species_ids_3, function(sp)
  seq(min(demography_inf_ppt_sens$clim[demography_inf_ppt_sens$Spp == sp]),
      max(demography_inf_ppt_sens$clim[demography_inf_ppt_sens$Spp == sp]),
      length.out = 30)) %>% setNames(as.character(species_ids_3))

clim_range_spik <- lapply(species_ids_2, function(sp)
  seq(min(demography_spik_ppt_sens$clim[demography_spik_ppt_sens$Spp == sp]),
      max(demography_spik_ppt_sens$clim[demography_spik_ppt_sens$Spp == sp]),
      length.out = 30)) %>% setNames(as.character(species_ids_2))

# ── Delta summaries: full vs. sensitivity, per vital rate ──────────────────
delta_surv <- bind_rows(
  summarise_delta_over_clim(compute_delta_surv, posterior_surv_full, clim_range_surv, species_ids_3, "Full data"),
  summarise_delta_over_clim(compute_delta_surv, posterior_surv_sens, clim_range_surv, species_ids_3, "Excl. KER/SON")
) %>% mutate(trait = "Survival")

delta_grow <- bind_rows(
  summarise_delta_over_clim(compute_delta_grow, posterior_grow_full, clim_range_grow, species_ids_3, "Full data"),
  summarise_delta_over_clim(compute_delta_grow, posterior_grow_sens, clim_range_grow, species_ids_3, "Excl. KER/SON")
) %>% mutate(trait = "Growth")

delta_inf <- bind_rows(
  summarise_delta_over_clim(compute_delta_inf, posterior_inf_full, clim_range_inf, species_ids_3, "Full data"),
  summarise_delta_over_clim(compute_delta_inf, posterior_inf_sens, clim_range_inf, species_ids_3, "Excl. KER/SON")
) %>% mutate(trait = "Inflorescence")

delta_spik <- bind_rows(
  summarise_delta_over_clim(compute_delta_spik, posterior_spik_full, clim_range_spik, species_ids_2, "Full data"),
  summarise_delta_over_clim(compute_delta_spik, posterior_spik_sens, clim_range_spik, species_ids_2, "Excl. KER/SON")
) %>% mutate(trait = "Spikelet")

delta_all <- bind_rows(delta_surv, delta_grow, delta_inf, delta_spik)

# ══════════════════════════════════════════════════════════════════════════
# 6. COMPACT COMPARISON TABLE FOR THE REVIEWER RESPONSE / SUPPLEMENT
#    One row per trait x species x herbivory treatment, evaluated at
#    representative precipitation values (min / median / max of each
#    species' range), full-data vs. excl. KER/SON, side by side.
# ══════════════════════════════════════════════════════════════════════════
make_summary_table <- function(delta_df) {
  delta_df %>%
    group_by(trait, species, herb, model) %>%
    group_modify(~ {
      n <- nrow(.x)
      .x %>%
        arrange(clim) %>%
        slice(c(1, ceiling(n / 2), n)) %>%
        mutate(clim_point = c("Low ppt", "Median ppt", "High ppt"))
    }) %>%
    ungroup() %>%
    mutate(clim_mm = round(ppt_mm(clim), 0)) %>%
    dplyr::select(trait, species, herb, clim_point, clim_mm, model,
                  median_delta, lower_90, upper_90, prob_delta_gt0) %>%
    arrange(trait, species, herb, clim_point, model)
}

sensitivity_table <- make_summary_table(delta_all) %>%
  mutate(across(c(median_delta, lower_90, upper_90, prob_delta_gt0), ~ round(.x, 3)))

# print(sensitivity_table, n = Inf)

# write.csv(
#   sensitivity_table,
#   file.path(out_dir, "sensitivity_KER_SON_delta_comparison.csv"),
#   row.names = FALSE
# )
message("Saved: sensitivity_KER_SON_delta_comparison.csv")

# Also export a LaTeX-ready version for the manuscript supplement, mirroring
# the xtable blocks used for the primary Delta tables in 05_Plot_vital_rate.R
# if (requireNamespace("xtable", quietly = TRUE)) {
#   library(xtable)
#   sens_latex <- sensitivity_table %>%
#     mutate(
#       species = case_when(
#         trait %in% c("Survival", "Growth", "Inflorescence") ~ SPP3_LABELS[as.character(species)],
#         TRUE ~ SPP2_LABELS[as.character(species)]
#       ),
#       herb = HERB_LABELS[as.character(herb)]
#     ) %>%
#     rename(
#       Trait = trait, Species = species, Herbivory = herb,
#       `Precip. level` = clim_point, `Precip. (mm)` = clim_mm,
#       Model = model, `Median $\\Delta$` = median_delta,
#       `Lower 90\\%` = lower_90, `Upper 90\\%` = upper_90,
#       `P($\\Delta>0$)` = prob_delta_gt0
#     )
#   
#   sens_xt <- xtable(
#     sens_latex,
#     caption = paste0(
#       "Sensitivity analysis excluding sites KER and SON. Posterior median, ",
#       "90\\% credible interval, and posterior probability of a positive ",
#       "endophyte effect [$\\Delta = $Pr(S+) $-$ Pr(S-), or the analogous ",
#       "quantity on the response scale for growth/inflorescence/spikelet ",
#       "production] for the full dataset and for the dataset excluding KER ",
#       "and SON, evaluated at low, median, and high precipitation within ",
#       "each species' observed range."
#     ),
#     label = "tab:sensitivity_KER_SON"
#   )
#   print(sens_xt, file = file.path(out_dir, "TableS_sensitivity_KER_SON.tex"),
#         include.rownames = FALSE, sanitize.text.function = identity,
#         caption.placement = "top")
#   message("Saved: TableS_sensitivity_KER_SON.tex")
# }

# ══════════════════════════════════════════════════════════════════════════
# 7. PUBLICATION-READY FIGURE — full data vs. excl. KER/SON, per trait
#    - x-axis on natural precipitation scale (mm), not the standardized clim
#    - species names italicized (facet columns)
#    - herbivory access / exclusion shown as rows
#    - shaded ribbon = 90% CrI, line = posterior median
#    - single shared legend, consistent with 09_Plot_effect_size.R style
# ══════════════════════════════════════════════════════════════════════════

make_sensitivity_plot <- function(trait_name, y_label, spp_labels) {
  
  df <- delta_all %>%
    filter(trait == trait_name) %>%
    mutate(
      clim_mm = ppt_mm(clim),
      species_lab = spp_labels[as.character(species)],
      herb_lab    = factor(HERB_LABELS[as.character(herb)],
                           levels = c("Herbivory access", "Herbivory exclusion")),
      model = factor(model, levels = c("Full data", "Excl. KER/SON"))
    )
  
  spp_in_trait <- unique(df$species_lab)
  spp_labeller <- setNames(paste0("italic('", spp_in_trait, "')"), spp_in_trait)
  
  ggplot(df, aes(x = clim_mm, y = median_delta, color = model, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_grid(
      herb_lab ~ species_lab,
      labeller = labeller(species_lab = as_labeller(spp_labeller, label_parsed))
    ) +
    scale_color_manual(values = MODEL_COLORS, name = NULL) +
    scale_fill_manual(values = MODEL_COLORS, name = NULL) +
    labs(
      title = trait_name,
      x = "Precipitation (mm)",
      y = y_label
    ) +
    sens_theme +
    theme(
      strip.text.x = element_text(size = 12, face = "italic")   # species names
    )
}

Fig_sens_surv <- make_sensitivity_plot("Survival",      "\u0394 survival probability (S+ \u2212 S\u2212)", SPP3_LABELS)
Fig_sens_grow <- make_sensitivity_plot("Growth",         "\u0394 growth (size scale, S+ \u2212 S\u2212)",   SPP3_LABELS)
Fig_sens_inf  <- make_sensitivity_plot("Inflorescence",  "\u0394 inflorescence count (S+ \u2212 S\u2212)",   SPP3_LABELS)
Fig_sens_spik <- make_sensitivity_plot("Spikelet",       "\u0394 spikelet count (S+ \u2212 S\u2212)",        SPP2_LABELS)

# ── Individual panels, each its own PDF (in case only one trait is needed) ──
export_sens_plot <- function(plot_obj, filename, width = 8, height = 5) {
  Cairo::CairoPDF(file.path(out_dir, filename), width = width, height = height)
  print(plot_obj)
  dev.off()
  message("Saved: ", filename)
}

export_sens_plot(Fig_sens_surv, "FigS_sensitivity_KER_SON_survival.pdf",      width = 8, height = 5)
export_sens_plot(Fig_sens_grow, "FigS_sensitivity_KER_SON_growth.pdf",        width = 8, height = 5)
export_sens_plot(Fig_sens_inf,  "FigS_sensitivity_KER_SON_inflorescence.pdf", width = 8, height = 5)
export_sens_plot(Fig_sens_spik, "FigS_sensitivity_KER_SON_spikelet.pdf",      width = 6, height = 5)

# ── Combined 2x2 supplementary figure — matches the FigS_coef_caterpillar_
#    combined layout convention used in 09_Plot_effect_size.R ──────────────
bump_strip <- function(p) {
  p + theme(strip.text.x = element_text(size = 13, face = "bold.italic"))
}

sens_body <- ggarrange(
  bump_strip(Fig_sens_surv) + theme(legend.position = "none"),
  bump_strip(Fig_sens_grow) + theme(legend.position = "none"),
  bump_strip(Fig_sens_inf)  + theme(legend.position = "none"),
  bump_strip(Fig_sens_spik) + theme(legend.position = "none"),
  ncol = 2, nrow = 2,
  labels     = c("(a)", "(b)", "(c)", "(d)"),
  font.label = list(size = 10, face = "plain")
)

shared_legend_sens <- ggpubr::get_legend(
  Fig_sens_surv + theme(legend.position = "bottom", legend.box = "horizontal")
)

FigS_sensitivity_combined <- ggarrange(
  sens_body,
  ggpubr::as_ggplot(shared_legend_sens),
  ncol = 1, heights = c(1, 0.06)
)

Cairo::CairoPDF(
  file.path(out_dir, "FigS_sensitivity_KER_SON_combined.pdf"),
  width = 14, height = 11
)
print(FigS_sensitivity_combined)
dev.off()
message("Saved: FigS_sensitivity_KER_SON_combined.pdf")

message("Sensitivity analysis complete. See sensitivity_KER_SON_delta_comparison.csv, ",
        "TableS_sensitivity_KER_SON.tex, and FigS_sensitivity_KER_SON_*.pdf ",
        "for the reviewer response / manuscript supplement.")
