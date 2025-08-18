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
bayes_R2(fit_surv_ppt)

#options(timeout = 120)
# Survival----
fit_surv_ppt_intercept <- readRDS(url("https://www.dropbox.com/scl/fi/ezsb8geutk7lvhopeq75g/fit_surv_ppt_intercept.rds?rlkey=kcee3sggvo613yizunj1g23im&dl=1"))
fit_surv_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/ty2yavh70dopswv2wbrij/fit_surv_ppt_abiotic.rds?rlkey=sqd5etosja0pw41pocx6x9pyu&dl=1"))
fit_surv_ppt_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/m5d7ggxlhqvprl3m8dlij/fit_surv_ppt_biotic.rds?rlkey=c5kr5q6fqq7q0n18rmk7qjtu5&dl=1"))
fit_surv_ppt<- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_surv_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))

log_lik_surv_ppt_intercept <- loo::extract_log_lik(fit_surv_ppt_intercept, merge_chains = FALSE)
r_eff_surv_ppt_intercept <- loo::relative_eff(exp(log_lik_surv_ppt_intercept))
loo_surv_ppt_intercept <- loo(log_lik_surv_ppt_intercept, r_eff = r_eff_surv_ppt_intercept, cores = 4)
plot(loo_surv_ppt_intercept)

log_lik_surv_ppt_abiotic <- loo::extract_log_lik(fit_surv_ppt_abiotic, merge_chains = FALSE)
r_eff_surv_ppt_abiotic <- loo::relative_eff(exp(log_lik_surv_ppt_abiotic))
loo_surv_ppt_abiotic <- loo(log_lik_surv_ppt_abiotic, r_eff = r_eff_surv_ppt_abiotic, cores = 4)
plot(loo_surv_ppt_abiotic)

log_lik_surv_ppt_biotic <- loo::extract_log_lik(fit_surv_ppt_biotic, merge_chains = FALSE)
r_eff_surv_ppt_biotic <- loo::relative_eff(exp(log_lik_surv_ppt_biotic))
loo_surv_ppt_biotic <- loo(log_lik_surv_ppt_biotic, r_eff = r_eff_surv_ppt_biotic, cores = 4)
plot(loo_surv_ppt_biotic)

log_lik_surv_ppt <- loo::extract_log_lik(fit_surv_ppt, merge_chains = FALSE)
r_eff_surv_ppt <- loo::relative_eff(exp(log_lik_surv_ppt))
loo_surv_ppt <- loo(log_lik_surv_ppt, r_eff = r_eff_surv_ppt, cores = 4)
plot(loo_surv_ppt)


(comp_surv <- loo::loo_compare(loo_surv_ppt_intercept, loo_surv_ppt_abiotic,loo_surv_ppt_biotic,loo_surv_ppt))

# Growth----
fit_grow_ppt_intercept <- readRDS(url("https://www.dropbox.com/scl/fi/ezsb8geutk7lvhopeq75g/fit_grow_ppt_intercept.rds?rlkey=kcee3sggvo613yizunj1g23im&dl=1"))
fit_grow_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/ty2yavh70dopswv2wbrij/fit_grow_ppt_abiotic.rds?rlkey=sqd5etosja0pw41pocx6x9pyu&dl=1"))
fit_grow_ppt_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/m5d7ggxlhqvprl3m8dlij/fit_grow_ppt_biotic.rds?rlkey=c5kr5q6fqq7q0n18rmk7qjtu5&dl=1"))
fit_grow_ppt<- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_grow_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))

log_lik_grow_ppt_intercept <- loo::extract_log_lik(fit_grow_ppt_intercept, merge_chains = FALSE)
r_eff_grow_ppt_intercept <- loo::relative_eff(exp(log_lik_grow_ppt_intercept))
loo_grow_ppt_intercept <- loo(log_lik_grow_ppt_intercept, r_eff = r_eff_grow_ppt_intercept, cores = 4)
plot(loo_grow_ppt_intercept)

log_lik_grow_ppt_abiotic <- loo::extract_log_lik(fit_grow_ppt_abiotic, merge_chains = FALSE)
r_eff_grow_ppt_abiotic <- loo::relative_eff(exp(log_lik_grow_ppt_abiotic))
loo_grow_ppt_abiotic <- loo(log_lik_grow_ppt_abiotic, r_eff = r_eff_grow_ppt_abiotic, cores = 4)
plot(loo_grow_ppt_abiotic)

log_lik_grow_ppt_biotic <- loo::extract_log_lik(fit_grow_ppt_biotic, merge_chains = FALSE)
r_eff_grow_ppt_biotic <- loo::relative_eff(exp(log_lik_grow_ppt_biotic))
loo_grow_ppt_biotic <- loo(log_lik_grow_ppt_biotic, r_eff = r_eff_grow_ppt_biotic, cores = 4)
plot(loo_grow_ppt_biotic)

log_lik_grow_ppt <- loo::extract_log_lik(fit_grow_ppt, merge_chains = FALSE)
r_eff_grow_ppt <- loo::relative_eff(exp(log_lik_grow_ppt))
loo_grow_ppt <- loo(log_lik_grow_ppt, r_eff = r_eff_grow_ppt, cores = 4)
plot(loo_grow_ppt)

(comp_grow <- loo::loo_compare(loo_grow_ppt_intercept, loo_grow_ppt_abiotic,loo_grow_ppt_biotic,loo_grow_ppt))

# Flowering----
fit_flow_ppt_intercept <- readRDS(url("https://www.dropbox.com/scl/fi/ezsb8geutk7lvhopeq75g/fit_flow_ppt_intercept.rds?rlkey=kcee3sggvo613yizunj1g23im&dl=1"))
fit_flow_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/ty2yavh70dopswv2wbrij/fit_flow_ppt_abiotic.rds?rlkey=sqd5etosja0pw41pocx6x9pyu&dl=1"))
fit_flow_ppt_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/m5d7ggxlhqvprl3m8dlij/fit_flow_ppt_biotic.rds?rlkey=c5kr5q6fqq7q0n18rmk7qjtu5&dl=1"))
fit_flow_ppt<- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_flow_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))

log_lik_flow_ppt_intercept <- loo::extract_log_lik(fit_flow_ppt_intercept, merge_chains = FALSE)
r_eff_flow_ppt_intercept <- loo::relative_eff(exp(log_lik_flow_ppt_intercept))
loo_flow_ppt_intercept <- loo(log_lik_flow_ppt_intercept, r_eff = r_eff_flow_ppt_intercept, cores = 4)
plot(loo_flow_ppt_intercept)

log_lik_flow_ppt_abiotic <- loo::extract_log_lik(fit_flow_ppt_abiotic, merge_chains = FALSE)
r_eff_flow_ppt_abiotic <- loo::relative_eff(exp(log_lik_flow_ppt_abiotic))
loo_flow_ppt_abiotic <- loo(log_lik_flow_ppt_abiotic, r_eff = r_eff_flow_ppt_abiotic, cores = 4)
plot(loo_flow_ppt_abiotic)

log_lik_flow_ppt_biotic <- loo::extract_log_lik(fit_flow_ppt_biotic, merge_chains = FALSE)
r_eff_flow_ppt_biotic <- loo::relative_eff(exp(log_lik_flow_ppt_biotic))
loo_flow_ppt_biotic <- loo(log_lik_flow_ppt_biotic, r_eff = r_eff_flow_ppt_biotic, cores = 4)
plot(loo_flow_ppt_biotic)

log_lik_flow_ppt <- loo::extract_log_lik(fit_flow_ppt, merge_chains = FALSE)
r_eff_flow_ppt <- loo::relative_eff(exp(log_lik_flow_ppt))
loo_flow_ppt <- loo(log_lik_flow_ppt, r_eff = r_eff_flow_ppt, cores = 4)
plot(loo_flow_ppt)

(comp_flow <- loo::loo_compare(loo_flow_ppt_intercept, loo_flow_ppt_abiotic,loo_flow_ppt_biotic,loo_flow_ppt))

# Spikelet----
fit_spik_ppt_intercept <- readRDS(url("https://www.dropbox.com/scl/fi/ezsb8geutk7lvhopeq75g/fit_spik_ppt_intercept.rds?rlkey=kcee3sggvo613yizunj1g23im&dl=1"))
fit_spik_ppt_abiotic <- readRDS(url("https://www.dropbox.com/scl/fi/ty2yavh70dopswv2wbrij/fit_spik_ppt_abiotic.rds?rlkey=sqd5etosja0pw41pocx6x9pyu&dl=1"))
fit_spik_ppt_biotic <- readRDS(url("https://www.dropbox.com/scl/fi/m5d7ggxlhqvprl3m8dlij/fit_spik_ppt_biotic.rds?rlkey=c5kr5q6fqq7q0n18rmk7qjtu5&dl=1"))
fit_spik_ppt<- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_spik_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))

log_lik_spik_ppt_intercept <- loo::extract_log_lik(fit_spik_ppt_intercept, merge_chains = FALSE)
r_eff_spik_ppt_intercept <- loo::relative_eff(exp(log_lik_spik_ppt_intercept))
loo_spik_ppt_intercept <- loo(log_lik_spik_ppt_intercept, r_eff = r_eff_spik_ppt_intercept, cores = 4)
plot(loo_spik_ppt_intercept)

log_lik_spik_ppt_abiotic <- loo::extract_log_lik(fit_spik_ppt_abiotic, merge_chains = FALSE)
r_eff_spik_ppt_abiotic <- loo::relative_eff(exp(log_lik_spik_ppt_abiotic))
loo_spik_ppt_abiotic <- loo(log_lik_spik_ppt_abiotic, r_eff = r_eff_spik_ppt_abiotic, cores = 4)
plot(loo_spik_ppt_abiotic)

log_lik_spik_ppt_biotic <- loo::extract_log_lik(fit_spik_ppt_biotic, merge_chains = FALSE)
r_eff_spik_ppt_biotic <- loo::relative_eff(exp(log_lik_spik_ppt_biotic))
loo_spik_ppt_biotic <- loo(log_lik_spik_ppt_biotic, r_eff = r_eff_spik_ppt_biotic, cores = 4)
plot(loo_spik_ppt_biotic)

log_lik_spik_ppt <- loo::extract_log_lik(fit_spik_ppt, merge_chains = FALSE)
r_eff_spik_ppt <- loo::relative_eff(exp(log_lik_spik_ppt))
loo_spik_ppt <- loo(log_lik_spik_ppt, r_eff = r_eff_spik_ppt, cores = 4)
plot(loo_spik_ppt)

(comp_spik <- loo::loo_compare(loo_spik_ppt_intercept, loo_spik_ppt_abiotic,loo_spik_ppt_biotic,loo_spik_ppt))
