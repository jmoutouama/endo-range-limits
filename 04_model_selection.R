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

# Inflorescence models ----
fit_inf_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/rnkijsri04jtrczshfmp6/fit_inf_abio_bio_endo.rds?rlkey=kyxj4f6mpwdp5m78wi0p6v64r&dl=1"))
fit_inf_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6ngnypika10yc0jrvr531/fit_inf_abio_bio_endo_linear.rds?rlkey=822f09t91dd2jw4r8svwmdh29&dl=1"))

loo_inf <- compute_loo(fit_inf_abio_bio_endo)
loo_inf_linear <- compute_loo(fit_inf_abio_bio_endo_linear)

comp_inf <- loo::loo_compare(loo_inf, loo_inf_linear)
comp_inf

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
  extract_loo_table(comp_inf, c("Quadratic", "Linear"), "Inflorescence"),
  extract_loo_table(comp_spik, c("Quadratic", "Linear"), "Spikelet")
)

# Print table
print(comp_table)

# Export to CSV
#write_csv(comp_table, "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Data/loo_comparison_summary_linear_vs_original.csv")

# What drive the outcome of symbiosis ? 
# Survival models ----
fit_surv_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/hfqoeo14h5fs4chl6dodb/fit_surv_endo_clim.rds?rlkey=qcw7rtpz31khpeplm5wuco6qg&dl=1"))
fit_surv_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/y288qm224kc1rz0jyfzhe/fit_surv_endo_herb.rds?rlkey=f95fswt5g0xthvoezw1mpyiao&dl=1"))
fit_surv_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/khyez2xegn8j5elkchaf6/fit_surv_abio_bio_endo_linear.rds?rlkey=zh4hx9czjov9aivlmaycfcuq7&dl=1"))

loo_surv_endo_clim <- compute_loo(fit_surv_endo_clim)
loo_surv_endo_herb <- compute_loo(fit_surv_endo_herb)
loo_surv_abio_bio_endo_linear <- compute_loo(fit_surv_abio_bio_endo_linear)


comp_surv_driver <- loo::loo_compare(loo_surv_endo_clim, loo_surv_endo_herb,loo_surv_abio_bio_endo_linear)
comp_surv_driver

# Growth models ----
fit_grow_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/ufdm43yzjg0dqjliuo2ri/fit_grow_endo_clim.rds?rlkey=kq9mswjtx554g114nr65v0i8p&dl=1"))
fit_grow_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/avidlmm6htcy59nvsg2m6/fit_grow_endo_herb.rds?rlkey=ilyuu0qjcwxr7ydudbe08krer&dl=1"))
fit_grow_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/o62tvjf8aqqz15gjxnrjn/fit_grow_abio_bio_endo_linear.rds?rlkey=xg1s6u5ctsluampm1l2zy1wqn&dl=1"))

loo_grow_endo_clim <- compute_loo(fit_grow_endo_clim)
loo_grow_endo_herb <- compute_loo(fit_grow_endo_herb)
loo_grow_abio_bio_endo_linear <- compute_loo(fit_grow_abio_bio_endo_linear)


comp_grow_driver <- loo::loo_compare(loo_grow_endo_clim, loo_grow_endo_herb,loo_grow_abio_bio_endo_linear)
comp_grow_driver

# Inflorescence models ----
fit_inf_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/k0drs7vhqz86rizmrveyj/fit_inf_endo_clim.rds?rlkey=cy9l8foe16y5g7j7kmjusy53j&dl=1"))
fit_inf_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/xv6sxf4q6e1b6xbgts7j3/fit_inf_endo_herb.rds?rlkey=zgmilm8daq9kvvi93afi4z3sw&dl=1"))
fit_inf_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6ngnypika10yc0jrvr531/fit_inf_abio_bio_endo_linear.rds?rlkey=822f09t91dd2jw4r8svwmdh29&dl=1"))

loo_inf_endo_clim <- compute_loo(fit_inf_endo_clim)
loo_inf_endo_herb <- compute_loo(fit_inf_endo_herb)
loo_inf_abio_bio_endo_linear <- compute_loo(fit_inf_abio_bio_endo_linear)


comp_inf_driver <- loo::loo_compare(loo_inf_endo_clim, loo_inf_endo_herb,loo_inf_abio_bio_endo_linear)
comp_inf_driver

# Spikelets models ----
fit_spik_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/0j6g437ycf5z3udw97naz/fit_spik_endo_clim.rds?rlkey=270oqdolmcqq7u8wl2cdgn35c&dl=1"))
fit_spik_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/f3jmc2u2sebt255fe342y/fit_spik_endo_herb.rds?rlkey=nbng4yyuvb09sjx10sua8g16m&dl=1"))
fit_spik_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6fcebl4lw8mu94fz62hnh/fit_spik_abio_bio_endo_linear.rds?rlkey=zy25y44zocugs6shh68lwpy1q&dl=1"))

loo_spik_endo_clim <- compute_loo(fit_spik_endo_clim)
loo_spik_endo_herb <- compute_loo(fit_spik_endo_herb)
loo_spik_abio_bio_endo_linear <- compute_loo(fit_spik_abio_bio_endo_linear)


comp_spik_driver <- loo::loo_compare(loo_spik_endo_clim, loo_spik_endo_herb,loo_spik_abio_bio_endo_linear)
comp_spik_driver

comp_table_driver <- bind_rows(
  extract_loo_table(comp_surv_driver, c("Endophye x climate", "Endophye x herd","Endophye x climate x herb"), "Survival"),
  extract_loo_table(comp_grow_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb"), "Growth"),
  extract_loo_table(comp_inf_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb"), "Inflorescence"),
  extract_loo_table(comp_spik_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb"), "Spikelet")
)
