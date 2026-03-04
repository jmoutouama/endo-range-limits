# -------------------------------
# 0. Clear workspace and load libraries
# -------------------------------
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)

# -------------------------------
# 1. Load data & Stan fit
# -------------------------------
demography <- readRDS(
  url("https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1")
)

fit_model <- readRDS(
  url("https://www.dropbox.com/scl/fi/cnefmiwlubfsovix6mxmm/herbivory_endo_siteyear.rds?rlkey=u7xtiftzeqdn50k9qo2njv75i&dl=1")
)

# -------------------------------
# 2. Extract posterior predictions
# -------------------------------
posterior <- rstan::extract(fit_model)

# -------------------------------
# 3. Prepare observed plot-level herbivory
# -------------------------------
observed_herb <- demography %>%
  filter(!is.na(tiller_Herb_t1)) %>%
  mutate(
    Endo_label = ifelse(Endo == 1, "E+", "E-")
  ) %>%
  group_by(Species,Plot, Endo_label, ppt_scaled) %>%
  summarise(
    herb_plot_mean = mean(tiller_Herb_t1, na.rm = TRUE),
    n_obs = sum(!is.na(tiller_Herb_t1)),
    .groups = "drop"
  )

# -------------------------------
# 4. Generate species-specific climate grids
# -------------------------------
species_ppt <- demography %>%
  group_by(Species) %>%
  summarise(
    ppt_min = min(ppt_scaled, na.rm = TRUE),
    ppt_max = max(ppt_scaled, na.rm = TRUE),
    .groups = "drop"
  )

pred_grid <- species_ppt %>%
  rowwise() %>%
  mutate(ppt_seq = list(seq(from = ppt_min, to = ppt_max, length.out = 20))) %>%
  unnest(cols = c(ppt_seq)) %>%
  crossing(Endo = c(0,1)) %>%
  mutate(Endo_label = ifelse(Endo == 1, "E+", "E-"))

# -------------------------------
# 5. Compute posterior predictions
# -------------------------------
species_unique <- unique(pred_grid$Species)
n_draws <- nrow(posterior$b0)
pred_matrix <- matrix(NA, nrow = nrow(pred_grid), ncol = n_draws)

for(i in seq_len(nrow(pred_grid))){
  sp <- pred_grid$Species[i]
  sp_idx <- which(species_unique == sp)
  endo <- pred_grid$Endo[i]
  ppt  <- pred_grid$ppt_seq[i]
  
  lp <- posterior$b0[, sp_idx] +
    posterior$bendo[, sp_idx] * endo +
    posterior$bclim[, sp_idx] * ppt +
    posterior$bendo_clim[, sp_idx] * endo * ppt
  
  pred_matrix[i, ] <- exp(lp)  # NB expected mean
}

pred_summary_clim <- pred_grid %>%
  mutate(
    median = apply(pred_matrix, 1, median),
    lower_90 = apply(pred_matrix, 1, quantile, probs = 0.05),
    upper_90 = apply(pred_matrix, 1, quantile, probs = 0.95)
  )

# -------------------------------
# 6. Convert scaled precipitation to mm
# -------------------------------
ppt_mean <- mean(demography$ppt_log, na.rm = TRUE)
ppt_sd   <- sd(demography$ppt_log, na.rm = TRUE)

pred_summary_clim <- pred_summary_clim %>%
  mutate(ppt_mm = exp(ppt_seq * ppt_sd + ppt_mean))

observed_herb <- observed_herb %>%
  mutate(ppt_mm = exp(ppt_scaled * ppt_sd + ppt_mean))

# -------------------------------
# 7. Plot predictions + observed points
# -------------------------------
Cairo::CairoPDF(
  "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/herb_endo_diff.pdf",
  width = 8, height = 8
)
ggplot(pred_summary_clim, aes(x = ppt_mm, y = median, color = Endo_label, fill = Endo_label)) +
  
  # Posterior 90% CI ribbon
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = 0.2, color = NA) +
  
  # Posterior median line
  geom_line(size = 1.2) +
  
  # Observed plot-level herbivory points
  geom_point(
    data = observed_herb,
    aes(x = ppt_mm, y = herb_plot_mean, color = Endo_label),
    size = 2.5,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  
  # Facet by species
  facet_wrap(~ Species, ncol = 2, scales = "free") +
  
  # Color scales
  scale_color_manual(values = c("E-" = "tomato", "E+" = "cornflowerblue")) +
  scale_fill_manual(values = c("E-" = "tomato", "E+" = "cornflowerblue")) +
  
  # Axis labels
  labs(
    x = "Precipitation (mm)",
    y = "Herbivory",
    color = "Endophyte",
    fill  = "Endophyte"
  ) +
  
  theme_test(base_size = 14) +
  theme(
    legend.position = c(0.15, 0.85),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()
# Compute posterior difference (E+ - E-) per draw
# Identify unique species and ppt combinations
species_unique <- unique(pred_grid$Species)
ppt_unique     <- unique(pred_grid$ppt_seq)

pred_diff_list <- list()

for(sp in species_unique){
  for(ppt_val in ppt_unique){
    # indices for this species & ppt
    idx_Eminus <- which(pred_grid$Species == sp & pred_grid$Endo == 0 & pred_grid$ppt_seq == ppt_val)
    idx_Eplus  <- which(pred_grid$Species == sp & pred_grid$Endo == 1 & pred_grid$ppt_seq == ppt_val)
    
    # difference for all posterior draws
    diff_vec <- pred_matrix[idx_Eplus, ] - pred_matrix[idx_Eminus, ]
    
    pred_diff_list[[paste(sp, ppt_val)]] <- diff_vec
  }
}

# Convert to matrix: rows = species × ppt, cols = draws
pred_diff <- do.call(rbind, pred_diff_list)

# Summarize posterior differences
pred_diff_summary <- data.frame(
  median_diff = apply(pred_diff, 1, median, na.rm = TRUE),
  lower_90    = apply(pred_diff, 1, quantile, probs = 0.05, na.rm = TRUE),
  upper_90    = apply(pred_diff, 1, quantile, probs = 0.95, na.rm = TRUE)
)

# Add species and ppt info
species_ppt_grid <- expand.grid(Species = species_unique, ppt_seq = ppt_unique)
pred_diff_summary <- cbind(species_ppt_grid, pred_diff_summary)

summary_table <- pred_diff_summary %>%
  group_by(Species) %>%
  summarise(
    median_diff_mean = mean(median_diff, na.rm = TRUE),
    lower_90_min    = min(lower_90, na.rm = TRUE),
    upper_90_max    = max(upper_90, na.rm = TRUE),
    prop_positive   = mean(median_diff > 0),   # fraction of ppt points where E+ > E-
    .groups = "drop"
  )

summary_table

