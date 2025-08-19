# Purpose: Plot  vital rate models (survival, growth, flowering and spikelet) as function of climate or distance from niche center.
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
library(lubridate)
library(patchwork)
# Install cmdstanr from R
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# # Install CmdStan backend (once)
# cmdstanr::install_cmdstan()
#library(cmdstanr)
# Define some basic functions that we'll use later
quote_bare <- function(...) {
  substitute(alist(...)) %>%
    eval() %>%
    sapply(deparse)
}
set.seed(13)
# Demographic data -----
# Merge the demographic census
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1", stringsAsFactors = F)
dat23 <- read.csv("https://www.dropbox.com/scl/fi/fkwm0dan6nx2eaeyxjrjw/census2023.csv?rlkey=hy9209t53j9n7vxhta7axl5jk&dl=1", stringsAsFactors = F)
dat24 <- read.csv("https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1", stringsAsFactors = F)
dat25<-read.csv("https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1", stringsAsFactors = F)
datherbivory <- read.csv("https://www.dropbox.com/scl/fi/2gnlfozxpd2u9gprzp9oi/herbivory.csv?rlkey=sz2cloqxtbc6ou29j97l3t10f&dl=1", stringsAsFactors = F)

# Climatic data ----
climate_site <- readRDS(url("https://www.dropbox.com/scl/fi/yjxwq78v2iruk6sshf8if/climate_census_years.rds?rlkey=k2dg37ofts7u6xcqiwvqqnmf4&dl=1"))
# head(climate_site)

# calculate the average spikelet and inflorescence number for each census
dat23 %>%
  mutate(spikelet_23 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T)), digit = 0) -> dat23_spike
dat24 %>%
  mutate(spikelet_24 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_24 = round(rowSums(across(attachedInf_24:brokenInf_24), na.rm = T), digit = 0)) -> dat24_spike
dat25 %>%
  mutate(spikelet_25 = round(rowMeans(across(Spikelet_A:Spikelet_C), na.rm = T), digit = 0), Inf_25 = round(rowSums(across(attachedInf_25:brokenInf_25), na.rm = T), digit = 0)) -> dat25_spike

# Check for duplicate Tag_IDs
datini$Tag_ID <- as.character(datini$Tag_ID)
dat23_spike$Tag_ID <- as.character(dat23_spike$Tag_ID)
dat23_spike %>% count(Tag_ID) %>% filter(n > 1)
dat24_spike$Tag_ID <- as.character(dat24_spike$Tag_ID)
dat25_spike$Tag_ID <- as.character(dat25_spike$Tag_ID)
#view(dat24_spike)
#dat24_spike %>% count(Tag_ID) %>% filter(n > 1)
#dat25_spike %>% count(Tag_ID) %>% filter(n > 1)

## Merge the initial data to the 23 to get the Tag ID for each elements -----
datini23_spike <- dat23_spike %>%
  left_join(
    datini %>%
      dplyr::select(Site, Species, Plot, Position, Tag_ID, Population, 
                    GreenhouseID, Clone, Endo),
    by = "Tag_ID"
  ) %>%
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit")))

# view(datini23_spike)

#Combine datini and dat23
combined_data <- bind_rows(datini[,c("Site","Species","Plot","Tag_ID","Population","Endo")], datini23_spike[,c("Site","Species","Plot","Tag_ID","Population","Endo" )])%>% 
  distinct(Tag_ID, .keep_all = TRUE)

dat24_spike_sp_site_tag<-left_join(x = dat24_spike, y = combined_data, by = c("Tag_ID")) %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_24","brokenInf_24"))) %>% 
  filter(!is.na(Species))
#dat24_spike_sp_site_tag %>% count(Tag_ID) %>% filter(n > 1) # No duplicate 
# view(dat24_spike_sp_site_tag)

dat25_spike_sp_site_tag<- dat25_spike %>% 
  dplyr::select(-any_of(c("Spikelet_A", "Spikelet_B", "Spikelet_C", "digit","attachedInf_25","brokenInf_25"))) 
#dat25_spike_sp_site_tag %>% count(Tag_ID) %>% filter(n > 1) # No duplicate 

dat2324 <- datini23_spike %>%
  left_join(
    dat24_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Inf_24, Tiller_24, tiller_herb_24, date_24, stroma_24, spikelet_24),
    by = "Tag_ID"
  ) 

#dat2324 %>% count(Tag_ID) %>% filter(n > 1) # No duplicate 
# dat2324 %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)

dat2425 <- dat24_spike_sp_site_tag %>%
  left_join(
    dat25_spike_sp_site_tag %>%
      dplyr::select(Tag_ID, Inf_25, Tiller_25, tiller_herb_25, date_25, stroma_25, spikelet_25),
    by = "Tag_ID"
  )

# Change variable names
dat2324%>%
  mutate(
    tiller_t = Tiller_23,
    tiller_t1 = Tiller_24,
    inf_t = Inf_23,
    inf_t1 = Inf_24,
    spikelet_t = spikelet_23,
    spikelet_t1 = spikelet_24,
    tiller_Herb_t = tiller_Herb_23,
    tiller_Herb_t1 = tiller_herb_24,
    date_t=date_23,
    date_t1=date_24,
    census_year=rep(2024,nrow(datini23_spike))
  ) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Tag_ID,
    Population,
    Endo,
    tiller_t,
    tiller_t1,
    inf_t,
    inf_t1,
    spikelet_t,
    spikelet_t1,
    tiller_Herb_t,
    tiller_Herb_t1,
    census_year
  ) -> dat2324_t_t1

# dat2324_t_t1 %>% count(Tag_ID) %>% filter(n > 1) 
# dat2324_t_t1 %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)
dat2425 %>%
  mutate(
    tiller_t = Tiller_24,
    tiller_t1 = Tiller_25,
    inf_t = Inf_24,
    inf_t1 = Inf_25,
    spikelet_t = spikelet_24,
    spikelet_t1 = spikelet_25,
    tiller_Herb_t = tiller_herb_24,
    tiller_Herb_t1 = tiller_herb_25,
    date_t = date_24,
    date_t1 = date_25,
    census_year=rep(2025,nrow(dat24_spike_sp_site_tag))
  ) %>%
  dplyr::select(
    Site,
    Species,
    Plot,
    Tag_ID,
    Population,
    Endo,
    tiller_t,
    tiller_t1,
    inf_t,
    inf_t1,
    spikelet_t,
    spikelet_t1,
    tiller_Herb_t,
    tiller_Herb_t1,
    census_year
  ) -> dat2425_t_t1

# dat2425_t_t1 %>% count(Tag_ID) %>% filter(n > 1)
# dat2425_t_t1 %>%
#   group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)
dat_t_t1 <- rbind(dat2324_t_t1, dat2425_t_t1)
# Find duplicates by Tag ID within each census year
# dup_tags_per_year <- dat_t_t1 %>%
#   count(Tag_ID, census_year) %>%
#   filter(n > 1)
# dup_tags_per_year

## Merge the demographic data with the herbivory data -----
dat_t_t1_herb <- left_join(x = dat_t_t1, y = datherbivory, by = c("Site", "Plot", "Species")) # Merge the demographic data with the herbivory data
# head(dat_t_t1_herb)
# unique(dat_t_t1_herb$Species)
# view(dat_t_t1_herb)

## Merge the demographic data with the climatic data -----
climate_site_unique <- climate_site %>%
  distinct(site, Species, census_year, .keep_all = TRUE)

#view(climate_site_unique)
dat_t_t1_herb$census_year<-as.character(dat_t_t1_herb$census_year)
climate_site_unique$census_year<-as.character(climate_site_unique$census_year)

dat_t_t1_herb_clim <- dat_t_t1_herb %>%
  left_join(
    climate_site_unique,
    by = c("Site" = "site", "Species", "census_year")
  )

# dat_t_t1_herb_clim %>%
# group_by(Tag_ID,census_year) %>%
#   summarise(tag_rep = n()) %>%
#   filter(tag_rep>1)

dat_t_t1_herb_clim %>%
  mutate(
    site_year=interaction(Site,census_year),
    ##if the plant was dead or size was NA at the start of the transition year, survival is NA
    ##if the plant was alive at the start of the transition year, it survived if tillers_t1>0
    surv1 = ifelse(tiller_t>0,tiller_t1>0,NA),
    site_species_plot = interaction(Site, Species, Plot),
    ##if the plant was alive at the start of the transition year and it survived, growth is the log ratio of tiller counts, else NA
    ##note there that growth is conditional on survival, which I think is how it should be
    grow = ifelse(tiller_t>0 & tiller_t1>0,log(tiller_t1/tiller_t),NA_real_)
  ) -> demography_climate

# Survival----
## Read and format survival data to build the model
demography_climate_surv <- demography_climate %>%
  filter(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, site_species_plot, site_year, Endo, Herbivory,
    tiller_t, surv1, cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    # Grouping indices for random effects
    Site              = as.integer(factor(Site)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_year         = as.integer(factor(site_year)),
    site_species_plot = as.integer(factor(site_species_plot)),
    
    # Survival and predictors
    log_size_t0 = log(tiller_t),
    surv_t1     = as.integer(surv1),
    ppt         = log(cum_ppt)
  )

# summary(demography_climate_surv$cum_ppt)
# any(demography_climate_surv$cum_ppt <= 0)
# Convert into list for Stan
demography_surv_ppt <- list(
  nSpp       = n_distinct(demography_climate_surv$Species),
  nSite      = n_distinct(demography_climate_surv$Site),
  nsite_year = n_distinct(demography_climate_surv$site_year),
  nPop       = n_distinct(demography_climate_surv$Population),
  nPlot      = n_distinct(demography_climate_surv$site_species_plot),
  Spp        = demography_climate_surv$Species,
  site       = demography_climate_surv$Site,
  site_year  = demography_climate_surv$site_year,
  pop        = demography_climate_surv$Population,
  plot       = demography_climate_surv$site_species_plot,
  clim       = as.vector(demography_climate_surv$ppt),
  endo       = demography_climate_surv$Endo,      # already 0/1
  herb       = demography_climate_surv$Herbivory, # already 0/1
  size       = demography_climate_surv$log_size_t0,
  y          = demography_climate_surv$surv_t1,
  N          = nrow(demography_climate_surv)
)


fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/0g5pn2igdi65vr3ky7heh/fit_surv_ppt.rds?rlkey=pkdj56oi1s3mfdn4p8nkgvlpy&dl=1"))
# Create a new data frame for generating predictions
climate_range <- seq(min(demography_surv_ppt$clim),
                     max(demography_surv_ppt$clim),
                     length.out = 30)
endo_status <- c(0, 1) # Endophyte negative and positive
herb_status <- c(0, 1) # Herbivory no and yes
species <- 1:3 # Species 1, 2, 3

# Create a data frame with all combinations
predictions <- expand.grid(
  clim = climate_range,
  endo = endo_status,
  herb = herb_status,
  species = species
)

# Extract posterior samples
posterior_samples <- rstan::extract(fit_surv_ppt)

# Function to calculate predictions based on the posterior samples
get_predictions <- function(clim,
                            endo,
                            herb,
                            species_index,
                            posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted survival (logit scale)
  logit_preds <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  # Convert logit to probability using logistic function
  pred_probs <- 1 / (1 + exp(-logit_preds))
  return(pred_probs)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probs_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probs_matrix[i, ] <- get_predictions(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_preds <- get_predictions(0.5, 1, 0, 1, posterior_samples)
species_2_preds <- get_predictions(0.2, 0, 1, 2, posterior_samples)
species_3_preds <- get_predictions(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probs_df <- as.data.frame(pred_probs_matrix)
colnames(pred_probs_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probs_df <- cbind(predictions, pred_probs_df)

# Reshape the data frame so we have long format for ggplot
pred_probs_long_df <- gather(pred_probs_df,
                             key = "Posterior_Sample",
                             value = "Pred_Survival",
                             -clim,
                             -endo,
                             -herb,
                             -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervals <- pred_probs_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    lower_95 = quantile(Pred_Survival, 0.025),
    upper_95 = quantile(Pred_Survival, 0.975),
    median = quantile(Pred_Survival, 0.5),
    mean = mean(Pred_Survival) # Calculate the mean survival probability
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_data <- data.frame(
  clim = demography_surv_ppt$clim,
  # Your climate data
  endo = demography_surv_ppt$endo,
  # Your endophyte status data
  herb = demography_surv_ppt$herb,
  # Your herbivory status data
  species = demography_surv_ppt$Spp,
  # Your species data
  y = demography_surv_ppt$y # Observed survival
)


# --- Calculate differences E+ - E- ---
diff_df <- pred_probs_long_df %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Survival[endo == 1]) - mean(Pred_Survival[endo == 0]),
    .groups = "drop"
  )

# Credible intervals for differences
diff_ci <- diff_df %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    mean = mean(diff),
    median = median(diff),
    .groups = "drop"
  ) %>%
  mutate(panel = "Difference (E+ - E-)")

cred_intervals <- cred_intervals %>%
  mutate(panel = "Predicted survival probability")

plot_data <- bind_rows(cred_intervals, diff_ci)


library(patchwork)

# Split data for clarity
top_panel_data <- top_panel_data %>%
  mutate(panel = "Predicted survival probability")

lower_panel_data <- lower_panel_data %>%
  mutate(panel = "Δ (E+ - E-)")

# Combine so they share the same panel factor
plot_data <- bind_rows(top_panel_data, lower_panel_data) %>%
  mutate(panel = factor(panel,
                        levels = c("Predicted survival probability", "Δ (E+ - E-)")))


# Add a panel column to observed data if needed
observed_data <- observed_data %>%
  mutate(panel = "Predicted survival probability")
Cairo::CairoPDF("/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/PrSurival_ppt_diff.pdf",
width = 12, height = 6.5)
ggplot(plot_data) +
  # survival probability
  geom_line(
    data = subset(plot_data, panel == "Predicted survival probability"),
    aes(x = exp(clim), y = mean, color = factor(endo), group = endo),
    size = 1
  ) +
  geom_ribbon(
    data = subset(plot_data, panel == "Predicted survival probability"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90,
        fill = factor(endo), group = endo),
    alpha = 0.3, color = NA
  ) +
  geom_point(
    data = observed_data,
    aes(x = exp(clim), y = y, color = factor(endo)),
    size = 2, position = position_jitter(width = 0, height = 0.02)
  ) +
  
  # delta panel
  geom_line(
    data = subset(plot_data, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), y = mean),
    color = "black", size = 1
  ) +
  geom_ribbon(
    data = subset(plot_data, panel == "Δ (E+ - E-)"),
    aes(x = exp(clim), ymin = lower_90, ymax = upper_90),
    fill = "#9B6B96", alpha = 0.5
  ) +
  geom_hline(
    data = subset(plot_data, panel == "Δ (E+ - E-)"),
    aes(yintercept = 0),
    linetype = "dashed", color = "black"
  ) +
  
  facet_grid(panel ~ species + herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(x = "Precipitation (mm)", y = "", color = "Endophyte", fill = "Endophyte") +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) +
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) +
  theme_bw() +
  theme(
    legend.position = c(0.06, 0.85),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    text = element_text(family = "Arial"),
    strip.text = element_text(size = 13)
  )

dev.off()


# Growth----
## Read and format survival data to build the model
demography_climate%>%
  filter(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site, Plot,site_year, site_species_plot, Endo, Herbivory,
    tiller_t, grow,cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    site_year = as.integer(factor(site_year)),
    Endo = as.integer(Endo),
    Herbivory = as.integer(Herbivory)
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    grow = grow,
    ppt = log(cum_ppt),
  ) -> demography_climate_grow

## Separate each variable to use the same model stan
### Precipitation
demography_grow_ppt <- list(
  nSpp = demography_climate_grow$Species %>% n_distinct(),
  nSite = demography_climate_grow$Site %>% n_distinct(),
  nsite_year = demography_climate_grow$site_year %>% n_distinct(),
  nPop = demography_climate_grow$Population %>% n_distinct(),
  nPlot = demography_climate_grow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_grow$Species,
  site = demography_climate_grow$Site,
  site_year = demography_climate_grow$site_year,
  pop = demography_climate_grow$Population,
  plot = demography_climate_grow$site_species_plot,
  clim = as.vector(demography_climate_grow$ppt),
  endo = demography_climate_grow$Endo,
  herb = demography_climate_grow$Herbivory,
  size = demography_climate_grow$log_size_t0,
  y = demography_climate_grow$grow,
  N = nrow(demography_climate_grow)
)

fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/5oduhrkn3l0cu5b9soju5/fit_grow_ppt.rds?rlkey=mpxxl4aejowhdm29pij9jv8kk&dl=1"))
posterior_samples <- rstan::extract(fit_grow_ppt)
predictions <- expand.grid(
  clim = seq(
    min(demography_grow_ppt$clim),
    max(demography_grow_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:3
)
# Function to calculate predictions based on the posterior samples
get_predictions_grow <- function(clim,
                                 endo,
                                 herb,
                                 species_index,
                                 posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predg <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  # Keep predg
  pred_probg <- predg
  return(pred_probg)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probg_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probg_matrix[i, ] <- get_predictions_grow(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_predg <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predg <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)
species_3_predg <- get_predictions_grow(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probg_df <- as.data.frame(pred_probg_matrix)
colnames(pred_probg_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probg_df <- cbind(predictions, pred_probg_df)

# Reshape the data frame so we have long format for ggplot
pred_probg_long_df <- gather(pred_probg_df,
                             key = "Posterior_Sample",
                             value = "Pred_Growth",
                             -clim,
                             -endo,
                             -herb,
                             -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalg <- pred_probg_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Growth, 0.05),
    upper_90 = quantile(Pred_Growth, 0.95),
    lower_95 = quantile(Pred_Growth, 0.025),
    upper_95 = quantile(Pred_Growth, 0.975),
    median = quantile(Pred_Growth, 0.5),
    mean = mean(Pred_Growth) # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_grow <- data.frame(
  clim = demography_grow_ppt$clim,
  # Your climate data
  endo = demography_grow_ppt$endo,
  # Your endophyte status data
  herb = demography_grow_ppt$herb,
  # Your herbivory status data
  species = demography_grow_ppt$Spp,
  # Your species data
  y = demography_grow_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Growth_ppt.pdf",
  useDingbats = F,
  height = 9,
  width = 7
)
ggplot(cred_intervalg, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(data = observed_grow,
             aes(
               x = exp(clim),
               y = y,
               color = factor(endo)
             ),
             size = 3,
             position = position_jitter(width = 0, height = 0.02)  # jitter only y-direction
  ) + # Observed data points
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "Predicted realtive growth",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.12, 0.23),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

# Flowering----
demography_climate %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, inf_t1, cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    site_year=as.integer(factor(site_year)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(Endo),
    Herbivory = as.integer(Herbivory)
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    flow_t1 = inf_t1,
    ppt = log(cum_ppt)
  ) -> demography_climate_flow

## Separate each variable to use the same model stan
### Precipitation 
demography_flow_ppt <- list(
  nSpp = demography_climate_flow$Species %>% n_distinct(),
  nSite = demography_climate_flow$Site %>% n_distinct(),
  nsite_year = demography_climate_flow$site_year %>% n_distinct(),
  nPop = demography_climate_flow$Population %>% n_distinct(),
  nPlot = demography_climate_flow$site_species_plot %>% n_distinct(),
  Spp = demography_climate_flow$Species,
  site = demography_climate_flow$Site,
  site_year = demography_climate_flow$site_year,
  pop = demography_climate_flow$Population,
  plot = demography_climate_flow$site_species_plot,
  clim = as.vector(demography_climate_flow$ppt),
  endo = demography_climate_flow$Endo,
  herb = demography_climate_flow$Herbivory,
  size = demography_climate_flow$log_size_t0,
  y = demography_climate_flow$flow_t1,
  N = nrow(demography_climate_flow)
)

fit_flow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/1j3ln3jxk94s56c9j193q/fit_flow_ppt.rds?rlkey=ag5bdlhngtg2gsfbfbx15purr&dl=1"))
posterior_samples <- rstan::extract(fit_flow_ppt)

predictions <- expand.grid(
  clim = seq(
    min(demography_flow_ppt$clim),
    max(demography_flow_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:3
)
# Function to calculate predictions based on the posterior samples
get_predictions_flow <- function(clim,
                                 endo,
                                 herb,
                                 species_index,
                                 posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predf <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  #  predf
  pred_probf <- exp(predf)
  return(pred_probf)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probf_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probf_matrix[i, ] <- get_predictions_flow(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_predf <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predf <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)
species_3_predf <- get_predictions_grow(-0.3, 1, 1, 3, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probf_df <- as.data.frame(pred_probf_matrix)
colnames(pred_probf_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probf_df <- cbind(predictions, pred_probf_df)

# Reshape the data frame so we have long format for ggplot
pred_probf_long_df <- gather(pred_probf_df,
                             key = "Posterior_Sample",
                             value = "Pred_Flow",
                             -clim,
                             -endo,
                             -herb,
                             -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalf <- pred_probf_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Flow, 0.05),
    upper_90 = quantile(Pred_Flow, 0.95),
    lower_95 = quantile(Pred_Flow, 0.025),
    upper_95 = quantile(Pred_Flow, 0.975),
    median = quantile(Pred_Flow, 0.5),
    mean = mean(Pred_Flow) # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_flow <- data.frame(
  clim = demography_flow_ppt$clim,
  # Your climate data
  endo = demography_flow_ppt$endo,
  # Your endophyte status data
  herb = demography_flow_ppt$herb,
  # Your herbivory status data
  species = demography_flow_ppt$Spp,
  # Your species data
  y = demography_flow_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Flow_ppt.pdf",
  useDingbats = F,
  height = 9,
  width = 7
)
ggplot(cred_intervalf, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(data = observed_flow,
             aes(
               x = exp(clim),
               y = y,
               color = factor(endo)
             ),
             size = 3,
             position = position_jitter(width = 0, height = 0.02)) + # Observed data points
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "AGHY", "2" = "ELVI", "3" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "# Inflorescences",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.1, 0.2),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

# Spikelet----
demography_climate %>%
  filter(Species %in% c("ELVI", "POAU")) %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(
    Species, Population, Site,site_year, Plot, site_species_plot, Endo, Herbivory,
    tiller_t, spikelet_t1, cum_ppt
  ) %>%
  na.omit() %>%
  mutate(
    Site = as.integer(factor(Site)),
    site_year=as.integer(factor(site_year)),
    Species = as.integer(factor(Species)),
    Population = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo = as.integer(Endo),
    Herbivory = as.integer(Herbivory)
  ) %>%
  mutate(
    log_size_t0 = log(tiller_t),
    spi_t1 = spikelet_t1,
    ppt = log(cum_ppt)
  ) -> demography_climate_spik

### Precipitation
demography_spik_ppt <- list(
  nSpp = demography_climate_spik$Species %>% n_distinct(),
  nSite = demography_climate_spik$Site %>% n_distinct(),
  nsite_year =demography_climate_spik$site_year %>% n_distinct(),
  nPop = demography_climate_spik$Population %>% n_distinct(),
  nPlot = demography_climate_spik$site_species_plot %>% n_distinct(),
  Spp = demography_climate_spik$Species,
  site = demography_climate_spik$Site,
  site_year = demography_climate_spik$site_year,
  pop = demography_climate_spik$Population,
  plot = demography_climate_spik$site_species_plot,
  clim = as.vector(demography_climate_spik$ppt),
  endo = demography_climate_spik$Endo,
  herb = demography_climate_spik$Herbivory,
  size = demography_climate_spik$log_size_t0,
  y = demography_climate_spik$spikelet_t1,
  N = nrow(demography_climate_spik)
)


fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/fg562lkkl077j8qoegb8q/fit_spik_ppt.rds?rlkey=cbfyy7tq7tja3e9mxujv0scm6&dl=1"))
posterior_samples <- rstan::extract(fit_spik_ppt)
predictions <- expand.grid(
  clim = seq(
    min(demography_spik_ppt$clim),
    max(demography_spik_ppt$clim),
    length.out = 30
  ),
  endo = c(0, 1),
  herb = c(0, 1),
  species = 1:2
)
# Function to calculate predictions based on the posterior samples
get_predictions_spk <- function(clim,
                                endo,
                                herb,
                                species_index,
                                posterior_samples) {
  b0 <- posterior_samples$b0[, species_index]
  bendo <- posterior_samples$bendo[, species_index]
  bherb <- posterior_samples$bherb[, species_index]
  bclim <- posterior_samples$bclim[, species_index]
  bendoclim <- posterior_samples$bendoclim[, species_index]
  bendoherb <- posterior_samples$bendoherb[, species_index]
  bclim2 <- posterior_samples$bclim2[, species_index]
  bendoclim2 <- posterior_samples$bendoclim2[, species_index]
  # Predicted growth
  predspk <- b0 +
    bendo * endo +
    bclim * clim +
    bherb * herb +
    bendoclim * clim * endo +
    bendoherb * endo * herb +
    bclim2 * clim ^ 2 +
    bendoclim2 * endo * clim ^ 2
  #  predf
  pred_probspk <- exp(predspk)
  return(pred_probspk)
}

# Apply the function to generate predictions for all combinations
n_posterior_samples <- length(posterior_samples$b0) # Number of posterior samples
# Initialize a matrix to hold predictions for each posterior sample
pred_probspk_matrix <- matrix(NA, nrow = nrow(predictions), ncol = n_posterior_samples)

# Generate predictions for each combination of climate, endophyte, herbivory, and species
for (i in 1:nrow(predictions)) {
  pred_probspk_matrix[i, ] <- get_predictions_spk(
    predictions$clim[i],
    predictions$endo[i],
    predictions$herb[i],
    predictions$species[i],
    posterior_samples
  )
}
species_1_predspk <- get_predictions_grow(0.5, 1, 0, 1, posterior_samples)
species_2_predspk <- get_predictions_grow(0.2, 0, 1, 2, posterior_samples)

# Convert the matrix into a data frame with the correct structure
pred_probspk_df <- as.data.frame(pred_probspk_matrix)
colnames(pred_probspk_df) <- paste("Posterior_Sample", 1:n_posterior_samples)

# Add the `predictions` columns (clim_s, endo_s, herb_s, species)
pred_probspk_df <- cbind(predictions, pred_probspk_df)

# Reshape the data frame so we have long format for ggplot
pred_probspk_long_df <- gather(pred_probspk_df,
                               key = "Posterior_Sample",
                               value = "Pred_Spik",
                               -clim,
                               -endo,
                               -herb,
                               -species)

# Calculate credible intervals (90% and 95%) and mean survival probability
cred_intervalspk <- pred_probspk_long_df %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Spik, 0.05),
    upper_90 = quantile(Pred_Spik, 0.95),
    lower_95 = quantile(Pred_Spik, 0.025),
    upper_95 = quantile(Pred_Spik, 0.975),
    median = quantile(Pred_Spik, 0.5),
    mean = mean(Pred_Spik) # Calculate the mean growth
  ) %>%
  ungroup()

# observed_data should have columns: clim_s, endo_s, herb_s, species, y_s (observed survival)
observed_spk <- data.frame(
  clim = demography_spik_ppt$clim,
  # Your climate data
  endo = demography_spik_ppt$endo,
  # Your endophyte status data
  herb = demography_spik_ppt$herb,
  # Your herbivory status data
  species = demography_spik_ppt$Spp,
  # Your species data
  y = demography_spik_ppt$y # Observed survival
)

# Plot the results with credible intervals, mean survival, and observed points using ggplot2
pdf(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/Spik_ppt.pdf",
  useDingbats = F,
  height = 7,
  width = 7
)
ggplot(cred_intervalspk, aes(
  x = exp(clim),
  y = mean,
  color = factor(endo)
)) +
  # geom_line(aes(y = median), linetype = "solid", size = 1) +  # Plot the median survival probability
  geom_line(aes(y = mean), linetype = "solid", size = 1) + # Plot the mean survival probability (dashed line)
  geom_ribbon(
    aes(
      ymin = lower_90,
      ymax = upper_90,
      fill = factor(endo)
    ),
    alpha = 0.3,
    color = NA
  ) + # Credible interval
  geom_point(data = observed_spk,
             aes(
               x = exp(clim),
               y = y,
               color = factor(endo)
             ),
             size = 3,
             position = position_jitter(width = 0, height = 0.02)) + # Observed data points
  facet_grid(species ~ herb,
             scales = "free_y",
             labeller = labeller(
               species = c("1" = "ELVI", "2" = "POAU"),
               herb = c("0" = "Unfenced", "1" = "Fenced")
             )) +
  labs(
    x = "Precipitation (mm)",
    y = "# Spikelets",
    color = "Endophyte",
    fill = "Endophyte",
    title = ""
  ) +
  scale_color_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                     labels = c("E-", "E+")) + # Change endophyte labels
  scale_fill_manual(values = c("0" = "tomato", "1" = "cornflowerblue"),
                    labels = c("E-", "E+")) + # Change fill labels
  theme_classic() +
  theme(
    legend.position = c(0.1, 0.35),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_text(size = 10),
    # Reduce legend title size
    legend.text = element_text(size = 12),
    # Adjust legend text size
    axis.title = element_text(size = 13),
    # Increase axis title size
    axis.text = element_text(size = 10),
    # Increase axis label size
    strip.text = element_text(size = 13)
  )
dev.off()

