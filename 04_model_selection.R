# Project:
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
#options(timeout = 120)
# Survival----
fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/9xa46n5v7u1lddxaj69cs/fit_surv_ppt.rds?rlkey=1mbkby4394s04j7qej4kvzeo2&dl=1"))
fit_surv_distance <- readRDS(url("https://www.dropbox.com/scl/fi/jn2a8wzcezmceplrwd356/fit_surv_distance.rds?rlkey=ow5bw2g31ce7af0quxjjrzlfv&dl=1"))
fit_surv_geo_distance <- readRDS(url("https://www.dropbox.com/scl/fi/bt0087bzg6664s8gjnv8p/fit_surv_geo_distance.rds?rlkey=luyktp37f34tkzlz4ut38ueaw&dl=1"))

log_lik_surv_ppt <- loo::extract_log_lik(fit_surv_ppt, merge_chains = FALSE)
r_eff_surv_ppt <- loo::relative_eff(exp(log_lik_surv_ppt))
loo_surv_ppt <- loo(log_lik_surv_ppt, r_eff = r_eff_surv_ppt, cores = 4)
plot(loo_surv_ppt)

log_lik_surv_distance <- loo::extract_log_lik(fit_surv_distance, merge_chains = FALSE)
r_eff_surv_distance <- loo::relative_eff(exp(log_lik_surv_distance))
loo_surv_distance <- loo(log_lik_surv_distance, r_eff = r_eff_surv_distance, cores = 4)
plot(loo_surv_distance)

log_lik_surv_geo_distance <- loo::extract_log_lik(fit_surv_geo_distance, merge_chains = FALSE)
r_eff_surv_geo_distance <- loo::relative_eff(exp(log_lik_surv_geo_distance))
loo_surv_geo_distance <- loo(log_lik_surv_distance, r_eff = r_eff_surv_geo_distance, cores = 4)
plot(loo_surv_geo_distance)

(comp_surv <- loo::loo_compare(loo_surv_ppt, loo_surv_distance,loo_surv_geo_distance))

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
