# Project: Aridity and herbivory shape demographic benefits of fungal endophytes in cool-season grasses
# Purpose: Model comparison based on epdl/looic 
# Authors: Jacob K. Moutouama
# Date last modified (Y-M-D):2026-11-17
rm(list = ls())
source(here::here("R", "config.R"))   # shared paths + data_path() helper
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
library(bayestestR)

# Helper function to compute loo
compute_loo <- function(fit) {
  log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(log_lik))
  loo_result <- loo(log_lik, r_eff = r_eff, cores = 4)
  plot(loo_result)
  return(loo_result)
}

# Survival models ----
fit_surv_abio_bio_endo <- readRDS(data_path("fit_surv_abio_bio_endo.rds"))
fit_surv_abio_bio_endo_linear <- readRDS(data_path("fit_surv_abio_bio_endo_linear.rds"))

loo_surv <- compute_loo(fit_surv_abio_bio_endo)
loo_surv_linear <- compute_loo(fit_surv_abio_bio_endo_linear)

comp_surv <- loo::loo_compare(loo_surv, loo_surv_linear)
comp_surv

# Growth models ----
fit_grow_abio_bio_endo <- readRDS(data_path("fit_grow_abio_bio_endo.rds"))
fit_grow_abio_bio_endo_linear <- readRDS(data_path("fit_grow_abio_bio_endo_linear.rds"))

loo_grow <- compute_loo(fit_grow_abio_bio_endo)
loo_grow_linear <- compute_loo(fit_grow_abio_bio_endo_linear)

comp_grow <- loo::loo_compare(loo_grow, loo_grow_linear)
comp_grow

# Inflorescence models ----
fit_inf_abio_bio_endo <- readRDS(data_path("fit_inf_abio_bio_endo.rds"))
fit_inf_abio_bio_endo_linear <- readRDS(data_path("fit_inf_abio_bio_endo_linear.rds"))
fit_inf_abio_bio_endo_hurdle_linear<-readRDS(data_path("fit_inf_abio_bio_endo_hurdle_linear.rds"))
loo_inf <- compute_loo(fit_inf_abio_bio_endo)
loo_inf_linear <- compute_loo(fit_inf_abio_bio_endo_linear)
loo_inf_hurdle_linear <- compute_loo(fit_inf_abio_bio_endo_hurdle_linear)

comp_inf_l <- loo::loo_compare(loo_inf_hurdle_linear, loo_inf_linear)
comp_inf <- loo::loo_compare(loo_inf, loo_inf_linear)
comp_inf

# Spikelet models ----
fit_spik_ppt_abiotic <- readRDS(data_path("fit_spik_abio_bio_endo.rds"))
fit_spik_ppt_abiotic_linear <- readRDS(data_path("fit_spik_abio_bio_endo_linear.rds"))

loo_spik <- compute_loo(fit_spik_ppt_abiotic)
loo_spik_linear <- compute_loo(fit_spik_ppt_abiotic_linear)

comp_spik <- loo::loo_compare(loo_spik, loo_spik_linear)
comp_spik

# Helper to extract loo_compare summary into a tidy table with model labels
extract_loo_table <- function(comp_result, model_names, vital_rate) {
  tibble(
    Vital_rate = vital_rate,
    Model_type = model_names,
    ELPD_difference = comp_result[, "elpd_diff"],
    SE_difference   = comp_result[, "se_diff"]
  )
}

# Combine all comparisons into one table
comp_table <- bind_rows(
  extract_loo_table(comp_surv, c("Linear","Quadratic"), "Survival"),
  extract_loo_table(comp_grow, c("Linear","Quadratic"), "Growth"),
  extract_loo_table(comp_inf, c("Quadratic","Linear"), "Inflorescence"),
  extract_loo_table(comp_spik, c("Linear","Quadratic"), "Spikelet")
)

# Print table
print(comp_table)

# Export to CSV
#write_csv(comp_table, file.path(DATA_DIR, "loo_comparison_summary_linear_vs_original.csv"))