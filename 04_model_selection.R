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
fit_surv_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/fjq0aa56baziu3g6im205/fit_surv_abio_bio_endo.rds?rlkey=cqz351n8gwnpmfmcypiolx8tc&dl=1"))
fit_surv_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6vz1nxc2310os0u83skn6/fit_surv_abio_bio_endo_linear.rds?rlkey=mmghgcfx1glfwzg0t56wa42d4&dl=1"))

loo_surv <- compute_loo(fit_surv_abio_bio_endo)
loo_surv_linear <- compute_loo(fit_surv_abio_bio_endo_linear)

comp_surv <- loo::loo_compare(loo_surv, loo_surv_linear)
comp_surv

# Growth models ----
fit_grow_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/x5gpxjlzwfr6a997ka4gn/fit_grow_abio_bio_endo.rds?rlkey=otlvalgg26h7y0qcy84u33rf2&dl=1"))
fit_grow_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/jxujzmf5rcgptbzqbyu8u/fit_grow_abio_bio_endo_linear.rds?rlkey=hxs21k5algqp86md8aqsgc59c&dl=1"))

loo_grow <- compute_loo(fit_grow_abio_bio_endo)
loo_grow_linear <- compute_loo(fit_grow_abio_bio_endo_linear)

comp_grow <- loo::loo_compare(loo_grow, loo_grow_linear)
comp_grow

# Inflorescence models ----
fit_inf_abio_bio_endo <- readRDS(url("https://www.dropbox.com/scl/fi/hlwxsi38dhf1qswjyhwsp/fit_inf_abio_bio_endo.rds?rlkey=p96ec0vxfg71z5i9w6sft7mxp&dl=1"))
fit_inf_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/vn7ms9zrp91svfr9o9wy8/fit_inf_abio_bio_endo_linear.rds?rlkey=accz66armitumo1y2n142q3tj&dl=1"))

loo_inf <- compute_loo(fit_inf_abio_bio_endo)
loo_inf_linear <- compute_loo(fit_inf_abio_bio_endo_linear)

comp_inf <- loo::loo_compare(loo_inf, loo_inf_linear)
comp_inf

# Spikelet models ----
fit_spik_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/aptez3xb4v5n60l0fcf71/fit_spik_abio_bio_endo.rds?rlkey=wdp4g5r8vs2jg6rwwr70tuskc&dl=1"))
fit_spik_ppt_abiotic_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6vnqhv5f8q5r0xkyiz8jw/fit_spik_abio_bio_endo_linear.rds?rlkey=0fhzukeb0294bskqv7boheokm&dl=1"))

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
fit_surv_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/2n6zta5vijb6xt3djhxw7/fit_surv_endo_clim.rds?rlkey=pvrg561oekfzbua2a2p9cdfbm&dl=1"))
fit_surv_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/1b41a0twrz3x1kkeu5937/fit_surv_endo_herb.rds?rlkey=w6tn6iyy37o3zedq5soulgnck&dl=1"))
fit_surv_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6vz1nxc2310os0u83skn6/fit_surv_abio_bio_endo_linear.rds?rlkey=mmghgcfx1glfwzg0t56wa42d4&dl=1"))
fit_surv_abio_bio_endo_linear_wi <- readRDS(url("https://www.dropbox.com/scl/fi/mkkl86jzneey71hrsqay3/fit_surv_abio_bio_endo_linear_wi.rds?rlkey=zaxr4h1d4arhrygtq8dog674d&dl=1"))

loo_surv_endo_clim <- compute_loo(fit_surv_endo_clim)
loo_surv_endo_herb <- compute_loo(fit_surv_endo_herb)
loo_surv_abio_bio_endo_linear <- compute_loo(fit_surv_abio_bio_endo_linear)
loo_surv_abio_bio_endo_linear_wi <- compute_loo(fit_surv_abio_bio_endo_linear_wi)

comp_surv_driver <- loo::loo_compare(loo_surv_endo_clim, loo_surv_endo_herb,loo_surv_abio_bio_endo_linear,loo_surv_abio_bio_endo_linear_wi)
comp_surv_driver

# Growth models ----
fit_grow_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/rznv7exh29hss75ar3k1m/fit_grow_endo_clim.rds?rlkey=lc6rtgvweux5g2a2r8x81cmu6&dl=1"))
fit_grow_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/fi738wqjpuqznfq3k2a03/fit_grow_endo_herb.rds?rlkey=3mpe134uxsdcz22j5sqitgfpo&dl=1"))
fit_grow_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/jxujzmf5rcgptbzqbyu8u/fit_grow_abio_bio_endo_linear.rds?rlkey=hxs21k5algqp86md8aqsgc59c&dl=1"))
fit_grow_abio_bio_endo_linear_wi <- readRDS(url("https://www.dropbox.com/scl/fi/xxr1pjaye1b2pgs54t4um/fit_grow_abio_bio_endo_linear_wi.rds?rlkey=mzv0vdhh5gik8r0cve12lo66o&dl=1"))

loo_grow_endo_clim <- compute_loo(fit_grow_endo_clim)
loo_grow_endo_herb <- compute_loo(fit_grow_endo_herb)
loo_grow_abio_bio_endo_linear <- compute_loo(fit_grow_abio_bio_endo_linear)
loo_grow_abio_bio_endo_linear_wi <- compute_loo(fit_grow_abio_bio_endo_linear_wi)

comp_grow_driver <- loo::loo_compare(loo_grow_endo_clim, loo_grow_endo_herb,loo_grow_abio_bio_endo_linear,loo_grow_abio_bio_endo_linear_wi)
comp_grow_driver

# Inflorescence models ----
fit_inf_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/t9t32qq2hvd892bxa13ql/fit_inf_endo_clim.rds?rlkey=6ygzhdc9k7oj2o8uohcy7wbji&dl=1"))
fit_inf_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/bjo9gkxxpvk2cms4weyej/fit_inf_endo_herb.rds?rlkey=tskx6axw8vl0qjfy0t17v0wnu&dl=1"))
fit_inf_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/vn7ms9zrp91svfr9o9wy8/fit_inf_abio_bio_endo_linear.rds?rlkey=accz66armitumo1y2n142q3tj&dl=1"))
fit_inf_abio_bio_endo_linear_wi <- readRDS(url("https://www.dropbox.com/scl/fi/vu3rwpapiuhgyop1cpcwm/fit_inf_abio_bio_endo_linear_wi.rds?rlkey=32cbzcgbjq5tx912jbm1zzgn6&dl=1"))

loo_inf_endo_clim <- compute_loo(fit_inf_endo_clim)
loo_inf_endo_herb <- compute_loo(fit_inf_endo_herb)
loo_inf_abio_bio_endo_linear <- compute_loo(fit_inf_abio_bio_endo_linear)
loo_inf_abio_bio_endo_linear_wi <- compute_loo(fit_inf_abio_bio_endo_linear_wi)

comp_inf_driver <- loo::loo_compare(loo_inf_endo_clim, loo_inf_endo_herb,loo_inf_abio_bio_endo_linear,loo_inf_abio_bio_endo_linear_wi)
comp_inf_driver

# Spikelets models ----
fit_spik_endo_clim <- readRDS(url("https://www.dropbox.com/scl/fi/k7i6md4owzod6fligb2gm/fit_spik_endo_clim.rds?rlkey=cpwvsiz7w0pt3j7ckxhwfrhmg&dl=1"))
fit_spik_endo_herb <- readRDS(url("https://www.dropbox.com/scl/fi/7dqrsbmqadgh7pp5bbalw/fit_spik_endo_herb.rds?rlkey=zbsa9emyp0rcrx7sxhrs4eq44&dl=1"))
fit_spik_abio_bio_endo_linear <- readRDS(url("https://www.dropbox.com/scl/fi/6vnqhv5f8q5r0xkyiz8jw/fit_spik_abio_bio_endo_linear.rds?rlkey=0fhzukeb0294bskqv7boheokm&dl=1"))
fit_spik_abio_bio_endo_linear_wi <- readRDS(url("https://www.dropbox.com/scl/fi/4oan1tf3zl0nz7m5nz6pc/fit_spik_abio_bio_endo_linear_wi.rds?rlkey=thr7eubz4bi75m6041r96zyn1&dl=1"))

loo_spik_endo_clim <- compute_loo(fit_spik_endo_clim)
loo_spik_endo_herb <- compute_loo(fit_spik_endo_herb)
loo_spik_abio_bio_endo_linear <- compute_loo(fit_spik_abio_bio_endo_linear)
loo_spik_abio_bio_endo_linear_wi <- compute_loo(fit_spik_abio_bio_endo_linear_wi)


comp_spik_driver <- loo::loo_compare(loo_spik_endo_clim, loo_spik_endo_herb,loo_spik_abio_bio_endo_linear,loo_spik_abio_bio_endo_linear_wi)
comp_spik_driver

comp_table_driver <- bind_rows(
  extract_loo_table(comp_surv_driver, c("Endophye x climate", "Endophye x herd","Endophye x climate x herb"), "Survival"),
  extract_loo_table(comp_grow_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb"), "Growth"),
  extract_loo_table(comp_inf_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb"), "Inflorescence"),
  extract_loo_table(comp_spik_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb"), "Spikelet")
)

comp_table_driver <- bind_rows(
  extract_loo_table(comp_surv_driver, c("Endophye x climate", "Endophye x herd","Endophye x climate x herb","Endophye + climate + herb"), "Survival"),
  extract_loo_table(comp_grow_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb","Endophye + climate + herb"), "Growth"),
  extract_loo_table(comp_inf_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb","Endophye + climate + herb"), "Inflorescence"),
  extract_loo_table(comp_spik_driver, c("Endophye x climate", "Endophye x herb","Endophye x climate x herb","Endophye + climate + herb"), "Spikelet")
)
