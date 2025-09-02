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


(comp_surv <- loo::loo_compare(loo_surv_ppt_abiotic,loo_surv_ppt_biotic,loo_surv_ppt))

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

(comp_grow <- loo::loo_compare( loo_grow_ppt_abiotic,loo_grow_ppt_biotic,loo_grow_ppt))

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

(comp_flow <- loo::loo_compare( loo_flow_ppt_abiotic,loo_flow_ppt_biotic,loo_flow_ppt))

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

(comp_spik <- loo::loo_compare(loo_spik_ppt_abiotic,loo_spik_ppt_biotic,loo_spik_ppt))

# Radar 
# Example: your loo_compare results
comp_surv <- data.frame(
  model = c("model1","model3","model2"),
  elpd_diff = c(0.0,-2.1,-5.2),
  se_diff = c(0.0,3.7,2.0)
)

comp_grow <- data.frame(
  model = c("model1","model3","model2"),
  elpd_diff = c(0.0,-2.1,-5.2),
  se_diff = c(0.0,3.7,2.0)
)

comp_flow <- data.frame(
  model = c("model1","model3","model2"),
  elpd_diff = c(0.0,-2.1,-5.2),
  se_diff = c(0.0,3.7,2.0)
)

comp_spik <- data.frame(
  model = c("model1","model3","model2"),
  elpd_diff = c(0.0,-2.1,-5.2),
  se_diff = c(0.0,3.7,2.0)
)

# Rename models
rename_models <- function(df){
  df$model <- recode(df$model,
                     model1 = "abiotic",
                     model2 = "biotic",
                     model3 = "abiotic_biotic")
  df
}

comp_surv <- rename_models(comp_surv)
comp_grow <- rename_models(comp_grow)
comp_flow <- rename_models(comp_flow)
comp_spik <- rename_models(comp_spik)

# Add a column to identify the response
comp_surv$response <- "survival"
comp_grow$response <- "growth"
comp_flow$response <- "flowering"
comp_spik$response <- "spikelet"

# Combine all into one table
comp_all <- bind_rows(comp_surv, comp_grow, comp_flow, comp_spik)

library(ggplot2)
library(dplyr)

# Example data
comp_all <- data.frame(
  model = c("abiotic","abiotic_biotic","biotic",
            "abiotic","abiotic_biotic","biotic",
            "abiotic","abiotic_biotic","biotic",
            "abiotic","abiotic_biotic","biotic"),
  elpd_diff = c(0.0,-2.1,-5.2,
                0.0,-2.1,-5.2,
                0.0,-2.1,-5.2,
                0.0,-2.1,-5.2),
  se_diff = c(0.0,3.7,2.0,
              0.0,3.7,2.0,
              0.0,3.7,2.0,
              0.0,3.7,2.0),
  response = c("Survival","Survival","Survival",
               "Growth","Growth","Growth",
               "Inflorescence","Inflorescence","Inflorescence",
               "Spikelet","Spikelet","Spikelet")
)

# Mark the best model per response
comp_all <- comp_all %>%
  group_by(response) %>%
  mutate(best = elpd_diff == max(elpd_diff)) %>%
  ungroup()

# Color-blind friendly palette
cb_palette <- c("abiotic" = "#E69F00",
                "abiotic_biotic" = "#009E73",
                "biotic" = "#56B4E9")

# Plot with circles and stars, showing real elpd_diff
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/bestdriver.pdf",
  width = 8,
  height = 6
)

ggplot(comp_all, aes(x = model, y = elpd_diff, color = model)) +
  
  # Error bars
  geom_errorbar(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff),
                width = 0.25) +
  
  # Circles for actual elpd_diff values
  geom_point(aes(fill = model), shape = 21, size = 4, stroke = 1.5) +
  
  # Star for best model
  geom_text(data = comp_all %>% filter(best),
            aes(label = "*",x=0.91, y = max(elpd_diff) + se_diff + 0.9), 
            color = "black", size = 8, hjust = 0) +
  # Colors
  scale_fill_manual(values = cb_palette) +
  scale_color_manual(values = cb_palette) +
  
  # Facets
  facet_wrap(~response, scales = "free_y") +
  
  # Labels and theme
  labs(x = "Driver", y = "ELPD Difference",
       fill = "Model", color = "Model") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, color = "black")
  ) +
  
  # Flip axes
  coord_flip()
dev.off()
