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
fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/d0x30lqqcxnatupsm2hej/fit_grow_ppt.rds?rlkey=er2is1le25trin73an23ztfgm&dl=1"))
fit_grow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/3ayyysw9k68lessw5hv56/fit_grow_distance.rds?rlkey=3cu65tyq7gal3be38nsk3ve6b&dl=1"))
fit_grow_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/5mbsdo6c591noq8os7hks/fit_grow_geo_distance.rds?rlkey=qupj9g4jz8ui5weo0i97fwfow&dl=1"))

log_lik_grow_ppt <- loo::extract_log_lik(fit_grow_ppt, merge_chains = FALSE)
r_eff_grow_ppt <- loo::relative_eff(exp(log_lik_grow_ppt))
loo_grow_ppt <- loo(log_lik_grow_ppt, r_eff = r_eff_grow_ppt, cores = 4)
plot(loo_grow_ppt)

log_lik_grow_distance <- loo::extract_log_lik(fit_grow_distance, merge_chains = FALSE)
r_eff_grow_distance <- loo::relative_eff(exp(log_lik_grow_distance))
loo_grow_distance <- loo(log_lik_grow_distance, r_eff = r_eff_grow_distance, cores = 4)
plot(loo_grow_distance)

log_lik_grow_geo_distance <- loo::extract_log_lik(fit_grow_geo_distance, merge_chains = FALSE)
r_eff_grow_geo_distance <- loo::relative_eff(exp(log_lik_grow_geo_distance))
loo_grow_geo_distance <- loo(log_lik_grow_geo_distance, r_eff = r_eff_grow_geo_distance, cores = 4)
plot(loo_grow_geo_distance)

(comp_grow <- loo::loo_compare(loo_grow_ppt, loo_grow_distance,loo_grow_geo_distance))

# Flowering----
fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/zra9rhooij33qgpznbse6/fit_flow_ppt.rds?rlkey=4pse2luz1aj08fqn95rt8m72y&dl=1"))
fit_flow_distance <- readRDS(url("https://www.dropbox.com/scl/fi/bccl31vszjauwpod6kyrr/fit_flow_distance.rds?rlkey=0dx711bk0jjdx3yjk409kfodq&dl=1"))
fit_flow_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/klp5vh3c2q2og2ej05rzt/fit_flow_geo_distance.rds?rlkey=kc5s9bddrvcmo2swla3dzbr98&dl=1"))

log_lik_flow_ppt <- loo::extract_log_lik(fit_flow_ppt, merge_chains = FALSE)
r_eff_flow_ppt <- loo::relative_eff(exp(log_lik_flow_ppt))
loo_flow_ppt <- loo(log_lik_flow_ppt, r_eff = r_eff_flow_ppt, cores = 4)
# plot(loo_flow_ppt)

log_lik_flow_distance <- loo::extract_log_lik(fit_flow_distance, merge_chains = FALSE)
r_eff_flow_distance <- loo::relative_eff(exp(log_lik_flow_distance))
loo_flow_distance <- loo(log_lik_flow_distance, r_eff = r_eff_flow_distance, cores = 4)
# plot(loo_flow_distance)

log_lik_flow_geo_distance <- loo::extract_log_lik(fit_flow_geo_distance, merge_chains = FALSE)
r_eff_flow_geo_distance <- loo::relative_eff(exp(log_lik_flow_geo_distance))
loo_flow_geo_distance <- loo(log_lik_flow_geo_distance, r_eff = r_eff_flow_geo_distance, cores = 4)
# plot(loo_flow_distance)

(comp_flow <- loo::loo_compare(loo_flow_ppt, loo_flow_distance,loo_flow_geo_distance))

# Spikelet----
fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/g61i9urje1i0e2234522e/fit_spik_ppt.rds?rlkey=5sqw0l4a23ozzg17zodklgyjp&dl=1"))
fit_spik_distance <- readRDS(url("https://www.dropbox.com/scl/fi/ibgrtngs6k0sz5bzin3kv/fit_spik_distance.rds?rlkey=ty1utrssm9tffmv3kxr4wqm8v&dl=1"))
fit_spik_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/lpj0sbo7mur0d2zc1hf0y/fit_spik_geo_distance.rds?rlkey=wyxw0g1pfxgk1ugqki1rwbb84&dl=1"))

log_lik_spi_ppt <- loo::extract_log_lik(fit_spik_ppt, merge_chains = FALSE)
r_eff_spi_ppt <- loo::relative_eff(exp(log_lik_spi_ppt))
loo_spi_ppt <- loo(log_lik_spi_ppt, r_eff = r_eff_spi_ppt, cores = 4)
#plot(loo_spi_ppt)

log_lik_spi_distance <- loo::extract_log_lik(fit_spik_distance, merge_chains = FALSE)
r_eff_spi_distance <- loo::relative_eff(exp(log_lik_spi_distance))
loo_spi_distance <- loo(log_lik_spi_distance, r_eff = r_eff_spi_distance, cores = 4)
#plot(loo_spi_distance)

log_lik_spi_geo_distance <- loo::extract_log_lik(fit_spik_geo_distance, merge_chains = FALSE)
r_eff_spi_geo_distance <- loo::relative_eff(exp(log_lik_spi_geo_distance))
loo_spi_geo_distance <- loo(log_lik_spi_geo_distance, r_eff = r_eff_spi_geo_distance, cores = 4)
#plot(loo_spi_geo_distance)

(comp_spi <- loo::loo_compare(loo_spi_ppt, loo_spi_distance,loo_spi_geo_distance))
