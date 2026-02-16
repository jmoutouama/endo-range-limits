
# Publication-ready Bayesian comparison of E+ vs E-
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)

# Load data

demography_climate <- readRDS(
  url("https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1")
)

fit_her_endo_year <- readRDS(
  url("https://www.dropbox.com/scl/fi/8v28c6tdwvingtstgnzeg/fit_her_endo_year.rds?rlkey=04vp4iale1i8wvycsqvepujyz&dl=1")
)

# Extract posterior predictions

posterior <- rstan::extract(fit_her_endo_year)

# Transform into median and 90% CI
pred_summary <- t(apply(posterior$pred, 2, function(x) quantile(x, probs = c(0.05, 0.5, 0.95))))
colnames(pred_summary) <- c("lower_90", "median", "upper_90")
pred_summary <- as.data.frame(pred_summary)
pred_summary$obs <- 1:nrow(pred_summary)

# Join with original demography data
pred_data <- demography_climate %>%
  mutate(obs = 1:n()) %>%
  left_join(pred_summary, by = "obs") %>%
  mutate(
    pred_count_median = exp(median),
    pred_count_lower  = exp(lower_90),
    pred_count_upper  = exp(upper_90),
    Endo_label = ifelse(Endo == 1, "E+", "E-")
  )


# Summarize by species & endophyte
pred_summary_species <- pred_data %>%
  group_by(Species, Endo_label) %>%
  summarise(
    median_count = median(pred_count_median, na.rm = TRUE),
    lower_90     = quantile(pred_count_median, 0.05, na.rm = TRUE),
    upper_90     = quantile(pred_count_median, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Compute posterior difference E+ - E- for each obs
posterior_diff <- pred_data %>%
  group_by(Species, obs, Endo_label) %>%
  summarise(pred_count_median = median(pred_count_median, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Endo_label, values_from = pred_count_median) %>%
  mutate(diff = `E+` - `E-`)

diff_summary <- posterior_diff %>%
  group_by(Species) %>%
  summarise(
    median_diff = median(diff, na.rm = TRUE),
    lower_90    = quantile(diff, 0.05, na.rm = TRUE),
    upper_90    = quantile(diff, 0.95, na.rm = TRUE),
    prob_pos    = mean(diff > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(diff_summary)

# -------------------------
# Plot: posterior median + CI + jittered posterior
# -------------------------
dodge_width <- 0.4

ggplot() +
  # Jittered posterior draws
  geom_jitter(
    data = pred_data,
    aes(y = pred_count_median, x = Species, color = Endo_label),
    height = 0, width = 0.1, alpha = 0.15, size = 1
  ) +
  # Median points with dodge
  geom_point(
    data = pred_summary_species,
    aes(y = median_count, x = Species, fill = Endo_label, color = Endo_label),
    shape = 21, size = 5, stroke = 1.2,
    position = position_dodge(width = dodge_width)
  ) +
  # 90% credible intervals
  geom_errorbar(
    data = pred_summary_species,
    aes(x = Species, ymin = lower_90, ymax = upper_90, color = Endo_label),
    width = 0.3, size = 1, position = position_dodge(width = dodge_width)
  ) +
  # Colors
  scale_color_manual(values = c("E-" = "tomato", "E+" = "cornflowerblue")) +
  scale_fill_manual(values = c("E-" = "tomato", "E+" = "cornflowerblue")) +
  # Labels
  labs(
    y = "Predicted herbivory counts",
    x = NULL,
    color = "Endophyte",
    fill = "Endophyte"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.85),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Extract posterior predictions: N x n_samples
pred_mat <- posterior$pred  # rows = observations, columns = posterior draws
n_samples <- ncol(pred_mat)

# Add observation info
pred_info <- demography_climate %>%
  mutate(obs = 1:n(), Endo_label = ifelse(Endo == 1, "E+", "E-"))

# Compute species-level posterior difference: E+ - E-
posterior_diff_species <- lapply(unique(pred_info$Species), function(sp) {
  
  idx <- which(pred_info$Species == sp)
  
  # Split indices by Endo status
  idx_Eplus  <- idx[pred_info$Endo_label[idx] == "E+"]
  idx_Eminus <- idx[pred_info$Endo_label[idx] == "E-"]
  
  # Mean across observations within species & Endo, for each posterior draw
  mean_Eplus  <- colMeans(pred_mat[idx_Eplus, , drop = FALSE])
  mean_Eminus <- colMeans(pred_mat[idx_Eminus, , drop = FALSE])
  
  # Difference E+ - E-
  diff <- mean_Eplus - mean_Eminus
  
  data.frame(
    Species = sp,
    draw = 1:n_samples,
    diff = diff
  )
}) %>% bind_rows()

# Summarize posterior differences per species
diff_summary <- posterior_diff_species %>%
  group_by(Species) %>%
  summarise(
    median_diff = median(diff),
    lower_90    = quantile(diff, 0.05),
    upper_90    = quantile(diff, 0.95),
    prob_pos    = mean(diff > 0),
    .groups = "drop"
  )

diff_summary
