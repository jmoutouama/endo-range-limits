# Purpose: Model comparison based on epdl/looic 
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
fit_surv_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/tsyih2vbg04zf9odu2goe/fit_surv_abio_bio_endo.rds?rlkey=tp4no6tfv5mb6f85jwtkqfuyg&dl=1"))
fit_surv_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/khyez2xegn8j5elkchaf6/fit_surv_abio_bio_endo_linear.rds?rlkey=zh4hx9czjov9aivlmaycfcuq7&dl=1"))

loo_surv <- compute_loo(fit_surv_abio_bio_endo)
loo_surv_linear <- compute_loo(fit_surv_abio_bio_endo_linear)

comp_surv <- loo::loo_compare(loo_surv, loo_surv_linear)
comp_surv

# Growth models ----
fit_grow_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/4r5062xbfc66gqh5l6xbz/fit_grow_abio_bio_endo.rds?rlkey=spnn4nj0zvzfsss1kn3qsnocj&dl=1"))
fit_grow_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/o62tvjf8aqqz15gjxnrjn/fit_grow_abio_bio_endo_linear.rds?rlkey=xg1s6u5ctsluampm1l2zy1wqn&dl=1"))

loo_grow <- compute_loo(fit_grow_abio_bio_endo)
loo_grow_linear <- compute_loo(fit_grow_abio_bio_endo_linear)

comp_grow <- loo::loo_compare(loo_grow, loo_grow_linear)
comp_grow

# Flowering models ----
fit_flow_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/5717xz8nt6sph3neq6jj9/fit_flow_abio_bio_endo.rds?rlkey=p4s7391sdqgepd89x53tbgw82&dl=1"))
fit_flow_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/1v4f4thyh826qcuiiyhub/fit_flow_abio_bio_endo_linear.rds?rlkey=raj4ls5dcqkeeexvcqj8b495m&dl=1"))

loo_flow <- compute_loo(fit_flow_abio_bio_endo)
loo_flow_linear <- compute_loo(fit_flow_abio_bio_endo_linear)

comp_flow <- loo::loo_compare(loo_flow, loo_flow_linear)
comp_flow

# Spikelet models ----
fit_spik_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/pebuc3ysvv9rrr2dvihl2/fit_spik_abio_bio_endo.rds?rlkey=f1l4q9ucvk4h4236600zoyjci&dl=1"))
fit_spik_ppt_abiotic_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6fcebl4lw8mu94fz62hnh/fit_spik_abio_bio_endo_linear.rds?rlkey=zy25y44zocugs6shh68lwpy1q&dl=1"))

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
  extract_loo_table(comp_surv, c("Quadratic", "Linear"), "Survival"),
  extract_loo_table(comp_grow, c("Quadratic", "Linear"), "Growth"),
  extract_loo_table(comp_flow, c("Quadratic", "Linear"), "Flowering"),
  extract_loo_table(comp_spik, c("Quadratic", "Linear"), "Spikelet")
)

# Print table
print(comp_table)

# Export to CSV
write_csv(comp_table, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Dataloo_comparison_summary_linear_vs_original.csv")

