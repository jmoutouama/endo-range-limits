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
fit_surv_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/vu750s78j6zjnmbohp8a0/fit_surv_abio_bio_endo.rds?rlkey=52n3g1naikyqzj0ix57ve6e23&dl=1"))
fit_surv_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/mh5es9xqo4t608h12zg4q/fit_surv_abio_bio_endo_linear.rds?rlkey=akzrlhtqbrx3sut9h58aidp0v&dl=1"))

loo_surv <- compute_loo(fit_surv_abio_bio_endo)
loo_surv_linear <- compute_loo(fit_surv_abio_bio_endo_linear)

comp_surv <- loo::loo_compare(loo_surv, loo_surv_linear)
comp_surv

# Growth models ----
fit_grow_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/wdxva2182fcd2i6jzsjcr/fit_grow_abio_bio_endo.rds?rlkey=v79tfa2br3l0sd6mblas8rto2&dl=1"))
fit_grow_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/mu3gpry42ad7fkfbtlvne/fit_grow_abio_bio_endo_linear.rds?rlkey=pz8uiqevdm5ogy7qvm369qyvb&dl=1"))

loo_grow <- compute_loo(fit_grow_abio_bio_endo)
loo_grow_linear <- compute_loo(fit_grow_abio_bio_endo_linear)

comp_grow <- loo::loo_compare(loo_grow, loo_grow_linear)
comp_grow

# Inflorescence models ----
fit_inf_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/857ar8r838n5wauqo3dtz/fit_inf_abio_bio_endo.rds?rlkey=2mot44dp52xlxj19bov5jndly&dl=1"))
fit_inf_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/5lfkgq6d5a2vzx5h2t9yr/fit_inf_abio_bio_endo_linear.rds?rlkey=pfqvm7sg8un5c14slypn1aqa9&dl=1"))
fit_inf_abio_bio_endo_hurdle_linear<-readRDS(url("https://www.dropbox.com/scl/fi/v5kcirie79uae1dex6l1a/fit_inf_abio_bio_endo_hurdle_linear.rds?rlkey=pxp0pdptn95j5ubkjs819ou3x&dl=1"))
loo_inf <- compute_loo(fit_inf_abio_bio_endo)
loo_inf_linear <- compute_loo(fit_inf_abio_bio_endo_linear)
loo_inf_hurdle_linear <- compute_loo(fit_inf_abio_bio_endo_hurdle_linear)

comp_inf_l <- loo::loo_compare(loo_inf_hurdle_linear, loo_inf_linear)
comp_inf <- loo::loo_compare(loo_inf, loo_inf_linear)
comp_inf

# Spikelet models ----
fit_spik_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/a6o21wspoofthisda1aye/fit_spik_abio_bio_endo.rds?rlkey=sxsxzzfjjqopqnh5km3s92xpa&dl=1"))
fit_spik_ppt_abiotic_linear <- readRDS(url("https://www.dropbox.com/scl/fi/7ivmicuigz1pahg4vxa7q/fit_spik_abio_bio_endo_linear.rds?rlkey=8h3js8dnaue95ojom8evrkmkl&dl=1"))

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
#write_csv(comp_table, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/loo_comparison_summary_linear_vs_original.csv")