# =============================================================================
# DOES SOIL INTERACT WITH SYMBIONT STATUS TO INFLUENCE DEMOGRAPHY?
# Full Bayesian Analysis via rstan
# Author: Jacob Moutouama
# =============================================================================
rm(list = ls())

# ── Packages ──────────────────────────────────────────────────────────────────
library(tidyverse)
library(rstan)
library(posterior)
library(bayesplot)
library(patchwork)
library(ggdist)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ── Colour palette ────────────────────────────────────────────────────────────
ENDO_COLORS <- c("S-" = "#E2492D", "S+" = "#3E6DC9")
# Plotmath expressions (not markdown) for the S+/S- superscript labels.
# Base R's plotmath renders these directly with no extra package, and it's
# the same mechanism already used for PANEL_LABELS/FOCAL_LABS elsewhere in
# this script, so it always draws correctly regardless of the graphics
# device (screen, pdf, tiff, patchwork-assembled figures, etc.).
# NOTE: this is a *list* of expressions, matched by position to
# breaks = endo_levels (i.e. c("S-", "S+")) wherever it's used below --
# unlike ENDO_COLORS this can't be a named character vector.
ENDO_LABELS <- list(
  bquote(italic(S)^"\u2212"),
  bquote(italic(S)^"+")
)
endo_levels <- c("S-", "S+")
site_order  <- c("SON", "KER", "BFL", "BAS", "COL", "HUN", "LAF")

# ── Shared theme (vr_theme) ───────────────────────────────────────────────────
vr_theme <- function() {
  theme_classic() +
    theme(
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line         = element_line(color = "black", linewidth = 0.1),
      axis.title        = element_text(size = 8),
      axis.text         = element_text(size = 6),
      axis.ticks.x      = element_line(color = "black", linewidth = 0.2),
      axis.ticks.y      = element_line(color = "black", linewidth = 0.2),
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 8),
      panel.spacing.y   = unit(0.2, "cm"),
      text              = element_text(family = "Arial"),
      strip.text.x      = element_text(size = 8, color = "black"),
      strip.text.y      = element_text(size = 8, color = "black"),
      strip.background  = element_rect(color = "black", fill = "grey80",
                                       linewidth = 0.2),
      plot.tag          = element_text(size = 9, face = "bold",
                                       margin = margin(r = 0))
    )
}

# ── Output directories ────────────────────────────────────────────────────────
output_dir <- "/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/"
fig_dir    <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Manuscript/Ecology letters/Manuscript/EcoletterR1"

save_fig <- function(p, filename, width_mm = 200, height_mm = 110) {
  ggsave(file.path(fig_dir, filename), plot = p,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, device = "pdf")
  message("Saved: ", filename)
}

TERM_COLORS <- c(
  "Main effect" = "#4D4D4D",
  "Interaction" = "#4D4D4D"
)

# =============================================================================
# PART 1 — DATA PREPARATION & EXPLORATORY DIAGNOSTICS
# =============================================================================

# ── 1.1  Column definitions ───────────────────────────────────────────────────
spike_indiv <- c("spike_inflo_a", "spike_inflo_b", "spike_inflo_c", "spike_inflo_d",
                 "spike_inflo_e", "spike_inflo_f", "spike_inflo_g", "spike_inflo_h",
                 "spike_inflo_i", "spike_infl_j")
spike_all   <- c(spike_indiv, "spike_inflos_unk")

# ── Birthday parser ───────────────────────────────────────────────────────────
parse_birthday <- function(x) {
  x[x == ""] <- NA
  x <- sub("/(2\\d)$", "/20\\1", x)
  as.Date(x, "%m/%d/%Y")
}

# ── n_Inflo decoder ───────────────────────────────────────────────────────────
decode_n_inflo <- function(x) {
  excel_epoch <- as.Date("1899-12-30")
  result      <- suppressWarnings(as.integer(x))
  date_vals   <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  is_date     <- !is.na(date_vals) & is.na(result)
  result[is_date] <- as.integer(date_vals[is_date] - excel_epoch)
  result
}

# ── Load & clean ──────────────────────────────────────────────────────────────
soils <- read.csv(
  "https://www.dropbox.com/scl/fi/tujihmutv87jufadtvbwj/reproduction_and_biomass-endo-soil.csv?rlkey=pa04j3jfae5cldlbc880w3xr3&dl=1",
  stringsAsFactors = FALSE,
  check.names      = FALSE
) %>%
  rename(Symbiont    = "Endo",
         seed_save   = "seed save",
         seed_squash = "seed squash") %>%
  mutate(
    across(all_of(c(spike_all,
                    "abg_mass_sans_inflo", "seed_save", "seed_squash",
                    "Inflo_mass", "total_inflo", "tot_spikelet",
                    "avg_spikelet", "total_spik_calc", "abg_mass_tot")),
           ~ suppressWarnings(as.numeric(.))),

    calc_n_inflo      = decode_n_inflo(n_Inflo),
    calc_n_measured   = rowSums(!is.na(across(all_of(spike_indiv)))),
    calc_tot_spikelet = if_else(
      rowSums(!is.na(across(all_of(spike_all)))) == 0,
      NA_real_,
      rowSums(across(all_of(spike_all)), na.rm = TRUE)
    ),
    calc_avg_spikelet = if_else(
      is.na(calc_n_inflo) | calc_n_inflo == 0,
      NA_real_,
      calc_tot_spikelet / calc_n_inflo
    ),
    calc_total_inflo  = total_inflo,
    calc_total_spik   = calc_avg_spikelet * calc_total_inflo,

    seed_save_clean   = if_else(is.na(seed_save),   0, seed_save),
    seed_squash_clean = if_else(is.na(seed_squash), 0, seed_squash),
    calc_Inflo_mass   = seed_save_clean + seed_squash_clean,
    calc_abg_mass_tot = case_when(
      !is.na(abg_mass_sans_inflo) ~ abg_mass_sans_inflo + calc_Inflo_mass,
      !is.na(abg_mass_tot)        ~ abg_mass_tot,
      TRUE                        ~ NA_real_
    ),

    Birthday = parse_birthday(Birthday),
    Birthday_final = if_else(
      !is.na(`Replacement_(NEW BIRTHDAY)`) & `Replacement_(NEW BIRTHDAY)` != "",
      parse_birthday(`Replacement_(NEW BIRTHDAY)`),
      Birthday
    ),
    age_days = as.numeric(as.Date("2025-01-08") - Birthday_final),
    age_std  = as.numeric(scale(age_days)),

    Site     = factor(Site, levels = site_order),
    Symbiont = factor(Symbiont, levels = c("E-", "E+")),
    Symbiont = factor(recode(as.character(Symbiont),
                             "E-" = "S-", "E+" = "S+"),
                      levels = endo_levels),
    Pop      = as.factor(Pop),
    Tray     = as.factor(Tray)
  )

# ── 1.2  Load and join 30-year precipitation normals ─────────────────────────
# Precipitation CSVs contain 30-year mean annual precipitation (mm/yr) for the
# collection location of each soil-origin site. The 'Precipitation' column is
# a site index (1–7); we map it to site names using the coordinates, which
# match the site_order locations confirmed by longitude:
#   1 = KER (-99.14°), 2 = BAS (-97.17°), 3 = SON (-100.56°),
#   4 = BFL (-97.78°), 5 = LAF (-92.01°), 6 = COL (-96.37°), 7 = HUN (-95.48°)
site_precip_map <- c("1" = "KER", "2" = "BAS", "3" = "SON", "4" = "BFL",
                     "5" = "LAF", "6" = "COL", "7" = "HUN")

Aghy_ppt <- read.csv(
  "https://www.dropbox.com/scl/fi/bst61ph31v0ewz6dyn79k/Agrostis_hyemalis_precipitation.csv?rlkey=xp1hq45eu7h8u6wnuznb9zlh0&dl=1",
  stringsAsFactors = FALSE
) %>%
  mutate(Site = site_precip_map[as.character(Precipitation)]) %>%
  dplyr::select(Site, precip_mean = mean)

poa_ppt <- read.csv(
  "https://www.dropbox.com/scl/fi/f6kkq559n795a3pgbdvsi/Poa_autumnalis_precipitation.csv?rlkey=adf5hwslj7uo74ui8c663y9hq&dl=1",
  stringsAsFactors = FALSE
) %>%
  mutate(Site = site_precip_map[as.character(Precipitation)]) %>%
  dplyr::select(Site, precip_mean = mean)

elvi_ppt <- read.csv(
  "https://www.dropbox.com/scl/fi/71ospl61twndy42ffz21m/Elymus_virginicus_precipitation.csv?rlkey=qihct3x5kcit7584sp78evcdy&dl=1",
  stringsAsFactors = FALSE
) %>%
  mutate(Site = site_precip_map[as.character(Precipitation)]) %>%
  dplyr::select(Site, precip_mean = mean)

# Join AGHY precipitation to soil data and standardise.
# precip_std is the covariate used in all models — standardised so that the
# Symbiont main effect is interpretable at mean precipitation.
aghysoils <- soils %>%
  filter(Species == "AGHY") %>%
  droplevels() %>%
  left_join(Aghy_ppt, by = "Site") %>%
  mutate(
    precip_std = as.numeric(scale(precip_mean)),
    Site       = factor(Site, levels = site_order),
    Site       = relevel(Site, ref = "LAF")
  )

# ── N clarification ───────────────────────────────────────────────────────────
# 196 pots were planted. Realised sample sizes for each response variable are
# lower because some plants died before the harvest census or failed to
# germinate after transplanting. We do not impute missing values; each model
# uses only the plants with observed data for that response.
N_planted <- 196L
N_total   <- nrow(aghysoils)
N_Sm      <- sum(aghysoils$Symbiont == "S-")
N_Sp      <- sum(aghysoils$Symbiont == "S+")

message(sprintf(
  "AGHY: %d pots planted | %d rows in data | S-: %d | S+: %d",
  N_planted, N_total, N_Sm, N_Sp
))
message("Realised N per response (excludes plants that died or failed to germinate):")
message(sprintf("  Biomass    : %d", sum(!is.na(aghysoils$calc_abg_mass_tot) &
                                          aghysoils$calc_abg_mass_tot > 0)))
message(sprintf("  Inflo      : %d", sum(!is.na(aghysoils$calc_total_inflo))))
message(sprintf("  Avg spklt  : %d", sum(!is.na(aghysoils$calc_avg_spikelet) &
                                           aghysoils$calc_avg_spikelet > 0)))

# ── Sample-size tables per response ───────────────────────────────────────────
aghysoils %>%
  filter(!is.na(calc_abg_mass_tot)) %>%
  count(Site, Symbiont) %>%
  complete(Site, Symbiont, fill = list(n = 0)) %>%
  arrange(Site)

aghysoils %>%
  filter(!is.na(calc_avg_spikelet)) %>%
  count(Site, Symbiont) %>%
  complete(Site, Symbiont, fill = list(n = 0)) %>%
  arrange(Site)

aghysoils %>%
  filter(!is.na(calc_total_inflo)) %>%
  count(Site, Symbiont) %>%
  complete(Site, Symbiont, fill = list(n = 0)) %>%
  arrange(Site)

# ── 1.3  Exploratory plots ────────────────────────────────────────────────────
resp_long <- aghysoils %>%
  dplyr::select(Site, Symbiont, precip_mean,
                calc_abg_mass_tot, calc_total_inflo,
                calc_total_spik, calc_avg_spikelet) %>%
  pivot_longer(-c(Site, Symbiont, precip_mean),
               names_to  = "response",
               values_to = "value") %>%
  dplyr::mutate(response = dplyr::recode(response,
                                         "calc_abg_mass_tot" = "Aboveground biomass (g)",
                                         "calc_total_inflo"  = "Total inflorescences",
                                         "calc_total_spik"   = "Total spikelets (projected)",
                                         "calc_avg_spikelet" = "Avg spikelets / inflo"))

p_precip_raw <- ggplot(resp_long %>% filter(!is.na(value), value > 0),
                       aes(x = precip_mean, y = value,
                           colour = Symbiont, fill = Symbiont)) +
  geom_point(alpha = 0.45, size = 1, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, formula = y ~ x) +
  scale_colour_manual(values = ENDO_COLORS, breaks = endo_levels,
                      labels = ENDO_LABELS, name = "Symbiont") +
  scale_fill_manual(values   = ENDO_COLORS, breaks = endo_levels,
                    labels = ENDO_LABELS, name = "Symbiont") +
  facet_wrap(~ response, scales = "free_y", ncol = 2) +
  labs(x = "30-yr mean annual precipitation of soil origin (mm/yr)",
       y = "Observed value") +
  vr_theme() +
  theme(legend.position = "top")

print(p_precip_raw)

# =============================================================================
# PART 2 — STAN MODEL BUILDING & SAMPLING
# =============================================================================

# ── 2.1  Design matrices ──────────────────────────────────────────────────────
# KEY CHANGE: Site factor replaced by precip_std as a continuous fixed covariate.
# This matches the field analysis structure (Symbiont × Precip) and allows
# direct comparison between greenhouse and field results.
# The null expectation for greenhouse data is that the Symbiont × precip_std
# interaction coefficient is near zero — i.e., soil-origin precipitation does
# not modify the endophyte effect when plants are grown in a common environment.
make_stan_data <- function(df, response_col, family = "negbinom") {
  d <- df %>%
    filter(!is.na(.data[[response_col]]),
           !is.na(Symbiont), !is.na(precip_std), !is.na(age_std)) %>%
    droplevels()

  if (family == "gaussian") {
    d <- d %>% filter(.data[[response_col]] > 0) %>% droplevels()
  }

  # Symbiont × precip_std mirrors the field model structure.
  # age_std included as a covariate to account for variation in plant age
  # at harvest (some plants were replaced after dying, giving them shorter
  # growing periods).
  X <- model.matrix(~ Symbiont * precip_std + age_std, data = d)

  y_vals <- d[[response_col]]
  if (family == "gaussian") y_vals <- log(y_vals)

  list(
    N      = nrow(d),
    K      = ncol(X),
    X      = X,
    y      = y_vals,
    y_cont = y_vals,
    n_pop  = nlevels(d$Pop),
    pop_id = as.integer(d$Pop),
    family = family,
    d_used = d
  )
}

sd_biomass <- make_stan_data(aghysoils, "calc_abg_mass_tot", "gaussian")
sd_inflo   <- make_stan_data(aghysoils, "calc_total_inflo",  "negbinom")
sd_totspi  <- make_stan_data(aghysoils, "calc_total_spik",   "negbinom")

# KER excluded from avg-spikelet: all KER S+ individuals had total_inflo == 0.
sd_avgspi  <- make_stan_data(
  aghysoils %>%
    filter(!is.na(calc_avg_spikelet),
           calc_avg_spikelet > 0,
           Site != "KER"),
  "calc_avg_spikelet", "negbinom"
)

# ── 2.2  Stan model: Gaussian / log-Normal (biomass) ─────────────────────────
# Response: log(aboveground biomass).
# Likelihood: Student-t with df = 3.
# Rationale for Student-t: log-biomass in small greenhouse experiments
# frequently shows heavier tails than a Gaussian due to a small number of
# exceptionally large or suppressed individuals. The Student-t with df = 3
# gives a robust likelihood that down-weights these outliers relative to a
# normal, while still being approximately normal for the bulk of the data.
# With only ~100–150 observations per model, a single outlier can
# substantially inflate sigma under a Gaussian; df = 3 guards against this.
# (A reviewer asking about this can be referred to Gelman et al. 2013,
# Bayesian Data Analysis, Ch. 17, on robust regression.)
stan_gaussian_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N]    y_cont;
  int<lower=1> n_pop;
  array[N] int<lower=1, upper=n_pop> pop_id;
}
parameters {
  vector[K]     beta;
  real<lower=0> sigma;
  vector[n_pop] z_pop;
  real<lower=0> sigma_pop;
}
transformed parameters {
  vector[n_pop] u_pop;
  u_pop = sigma_pop * z_pop;
}
model {
  beta      ~ normal(0, 1);
  sigma     ~ exponential(1);
  sigma_pop ~ exponential(1);
  z_pop     ~ normal(0, 1);
  y_cont    ~ student_t(3, X * beta + u_pop[pop_id], sigma);
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    real mu_n  = X[n] * beta + u_pop[pop_id[n]];
    y_rep[n]   = normal_rng(mu_n, sigma);
    log_lik[n] = student_t_lpdf(y_cont[n] | 3, mu_n, sigma);
  }
}
"

# ── 2.3  Stan model: Negative-Binomial (count responses) ─────────────────────
# Likelihood: neg_binomial_2_log.
# A Poisson likelihood was considered but overdispersion is expected for
# inflorescence and spikelet counts (variance > mean in all site × symbiont
# cells). The NB phi parameter absorbs this overdispersion.
stan_negbinom_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=0> y;
  int<lower=0>           y_max;
  int<lower=1>           n_pop;
  array[N] int<lower=1, upper=n_pop> pop_id;
}
parameters {
  vector[K]     beta;
  real          log_phi;
  vector[n_pop] z_pop;
  real<lower=0> sigma_pop;
}
transformed parameters {
  real<lower=0> phi;
  vector[n_pop] u_pop;
  phi   = exp(log_phi);
  u_pop = sigma_pop * z_pop;
}
model {
  beta      ~ normal(0, 1);
  log_phi   ~ normal(2, 1);
  sigma_pop ~ exponential(2);
  z_pop     ~ normal(0, 1);
  vector[N] log_mu;
  for (n in 1:N)
    log_mu[n] = X[n] * beta + u_pop[pop_id[n]];
  y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  array[N] int y_rep;
  vector[N]    log_lik;
  for (n in 1:N) {
    real log_mu_n    = X[n] * beta + u_pop[pop_id[n]];
    real safe_log_mu = fmin(log_mu_n, 10.0);
    real safe_phi    = fmin(phi, 1e6);
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_mu_n, phi);
    if (is_nan(safe_log_mu) || is_inf(safe_log_mu)) {
      y_rep[n] = -1;
    } else {
      real draw = neg_binomial_2_log_rng(safe_log_mu, safe_phi);
      y_rep[n]  = to_int(fmin(round(draw), y_max));
    }
  }
}
"

# ── 2.4  Compile models ───────────────────────────────────────────────────────
mod_gaussian <- stan_model(model_code = stan_gaussian_code,
                           model_name = "gaussian_sym_precip")
mod_negbinom <- stan_model(model_code = stan_negbinom_code,
                           model_name = "negbinom_sym_precip")

# ── 2.5  Sampling helper ──────────────────────────────────────────────────────
# MCMC settings: iter=2000, warmup=1000, chains=4.
# Reduced from iter=6000/warmup=2000 — the model has only 4 fixed-effect
# parameters (intercept, Symbiont, precip_std, interaction) plus age_std
# and a 2-level population random effect. This is a simple structure and
# 4000 post-warmup draws (4 chains × 1000) are more than sufficient for
# stable posterior summaries and LOO-CV.
run_stan <- function(stan_mod, data_list, family, seed = 13,
                     chains = 4, iter = 2000, warmup = 1000,
                     adapt_delta = 0.99) {
  if (family == "gaussian") {
    stan_data <- list(N      = data_list$N,
                      K      = data_list$K,
                      X      = data_list$X,
                      y_cont = data_list$y_cont,
                      n_pop  = data_list$n_pop,
                      pop_id = data_list$pop_id)
  } else {
    y_int <- as.integer(data_list$y)
    y_max <- as.integer(max(y_int, na.rm = TRUE) * 10L)
    stan_data <- list(N      = data_list$N,
                      K      = data_list$K,
                      X      = data_list$X,
                      y      = y_int,
                      y_max  = y_max,
                      n_pop  = data_list$n_pop,
                      pop_id = data_list$pop_id)
  }

  sampling(stan_mod, data = stan_data,
           chains  = chains, iter = iter, warmup = warmup,
           seed    = seed,
           control = list(adapt_delta   = adapt_delta,
                          max_treedepth = 12),
           refresh = 500)
}

# ── 2.6  Fit all four models ──────────────────────────────────────────────────
fit_biomass <- run_stan(mod_gaussian, sd_biomass, "gaussian")
fit_inflo   <- run_stan(mod_negbinom, sd_inflo,   "negbinom")
fit_totspi  <- run_stan(mod_negbinom, sd_totspi,  "negbinom")
fit_avgspi  <- run_stan(mod_negbinom, sd_avgspi,  "negbinom")

# ── 2.7  MCMC diagnostics ─────────────────────────────────────────────────────
check_hmc <- function(fit, label) {
  cat("\n===", label, "===\n")
  cat("Rhat summary (all should be < 1.01):\n")
  rhat_vals <- summary(fit)$summary[, "Rhat"]
  print(summary(rhat_vals))
  cat(sprintf("Divergent transitions  : %d\n", get_num_divergent(fit)))
  cat(sprintf("Max treedepth exceeded : %d\n", get_num_max_treedepth(fit)))
}

check_hmc(fit_biomass, "M1 Biomass")
check_hmc(fit_inflo,   "M2 Inflo")
check_hmc(fit_totspi,  "M3 Total spikelets (projected)")
check_hmc(fit_avgspi,  "M4 Avg spikelets")

p_trace_biomass <- mcmc_trace(fit_biomass, regex_pars = "beta\\[",
                              facet_args = list(ncol = 2)) +
  labs(title = "M1 \u2014 Biomass: trace plots (beta)") + vr_theme()
p_trace_inflo   <- mcmc_trace(fit_inflo, regex_pars = "beta\\[",
                              facet_args = list(ncol = 2)) +
  labs(title = "M2 \u2014 Inflo: trace plots (beta)") + vr_theme()

print(p_trace_biomass)
print(p_trace_inflo)

# ── 2.8  Posterior predictive checks ─────────────────────────────────────────
ppc_check <- function(fit, y_obs, label, xlim = NULL) {
  y_rep <- as.matrix(fit, pars = "y_rep")
  y_rep_sub <- y_rep[sample(nrow(y_rep), 500), ]
  
  p <- ppc_dens_overlay(y_obs, y_rep_sub) +
    labs(
      title = label,
      x = label,
      y = "Density"
    ) +
    vr_theme()
  
  if (!is.null(xlim))
    p <- p + coord_cartesian(xlim = xlim)
  
  p
}

p_ppc_biomass <- ppc_check(
  fit_biomass,
  sd_biomass$y,
  "",
  xlim = c(-8,5)
) +
  labs(x = "Log biomass", y = "Density")

p_ppc_inflo <- ppc_check(
  fit_inflo,
  sd_inflo$y,
  "",
  xlim = c(-1, 50)
) +
  labs(x = "Inflorescence count", y = "Density")

p_ppc_avgspi <- ppc_check(
  fit_avgspi,
  sd_avgspi$y,
  ""
) +
  labs(x = "Average spikelet production", y = "Density")
ppc_panel <- (p_ppc_biomass | p_ppc_inflo ) +
  plot_annotation(tag_levels = "A")
print(ppc_panel)
# ggsave(
#   file.path(fig_dir, "PPC_models_soilendo.pdf"),
#   plot   = ppc_panel,
#   width  = 170,
#   height = 100,
#   units  = "mm",
#   dpi    = 600,
#   device = cairo_pdf
# )

# =============================================================================
# PART 3 — POSTERIOR PREDICTIONS
# =============================================================================

# ── 3.1  Extract posterior fitted means across precipitation gradient ─────────
# FIX: Added ref_df argument. The back-transformation from precip_std to
# precip_mean (mm/yr) now uses mean/SD from the FULL aghysoils dataset (ref_df)
# rather than the filtered Stan dataset (d_used). This ensures that the
# prediction line's x-range matches the observed points exactly.
#
# Previously, using d_used for the back-transformation shifted/compressed the
# line because rows excluded by the Stan filter (those with NA in calc_total_spik
# or age_std) often include the most extreme precipitation values — causing the
# line to appear truncated relative to the raw observations.
extract_posterior_precip <- function(fit, sd_obj, link = "identity",
                                     n_grid = 50, ref_df = NULL) {
  draws_beta <- as.matrix(fit, pars = "beta")
  draws_upop <- as.matrix(fit, pars = "u_pop")
  d          <- sd_obj$d_used
  # Use full reference dataset for back-transform if supplied; fall back to d.
  ref        <- if (!is.null(ref_df)) ref_df else d
  n_iter     <- nrow(draws_beta)

  precip_seq <- seq(min(d$precip_std), max(d$precip_std), length.out = n_grid)

  map_dfr(endo_levels, function(sym) {
    sym_val <- if (sym == "S+") 1 else 0
    map_dfr(seq_along(precip_seq), function(i) {
      p <- precip_seq[i]
      # Design vector: intercept, Symbiont, precip_std, age_std=0, Symbiont:precip
      xvec <- c(1, sym_val, p, 0, sym_val * p)
      eta  <- draws_beta %*% xvec +
        rowMeans(draws_upop)   # marginalise over populations
      mu <- switch(link, log = exp(eta), identity = eta)
      tibble(
        Symbiont    = sym,
        precip_std  = p,
        # FIX: use ref$precip_mean (full dataset) not d$precip_mean (filtered)
        precip_mean = p * sd(ref$precip_mean) + mean(ref$precip_mean),
        draw        = seq_len(n_iter),
        mu          = as.numeric(mu)
      )
    })
  })
}

# Pass aghysoils as ref_df so back-transform uses the full dataset's mean/SD.
post_precip_biomass <- extract_posterior_precip(fit_biomass, sd_biomass, "identity",
                                                ref_df = aghysoils) %>%
  mutate(mu = exp(mu))
post_precip_inflo   <- extract_posterior_precip(fit_inflo,   sd_inflo,   "log",
                                                ref_df = aghysoils)
post_precip_totspi  <- extract_posterior_precip(fit_totspi,  sd_totspi,  "log",
                                                ref_df = aghysoils)
post_precip_avgspi  <- extract_posterior_precip(fit_avgspi,  sd_avgspi,  "log",
                                                ref_df = aghysoils)

# =============================================================================
# PART 4 — FIGURES
# =============================================================================
# Fig_Greenhouse_Combined: 2-column figure, panels (a)–(d).
#   Col 1 — Biomass: upper (predicted) / lower (Δ contrast), 2:1 height
#   Col 2 — Inflo:   upper (predicted) / lower (Δ contrast), 2:1 height
# Fig_Greenhouse_Coefficients: caterpillar for focal betas (separate figure).
#
# Panel tags (a–d) placed OUTSIDE plots via plot_annotation(tag_levels).
# Each sub-plot carries tag = "" so patchwork assigns the letter externally.
# =============================================================================

# ── 4.0  Display settings for biomass ─────────────────────────────────────────
# Biomass is modelled as a Gaussian on log(biomass) (see make_stan_data) and
# back-transformed with exp(). In grams on a linear y axis the fitted line is
# therefore an exponential curve.
#   BIOMASS_LOG_Y = TRUE          : log10 y axis (labelled in g) -> straight line.
#                                   Caption: "Biomass is shown on a logarithmic axis."
#   BIOMASS_DELTA_LOG_RATIO = TRUE: Δ panel shows log(S+ / S-) instead of the
#                                   difference in grams -> straight line.
#                                   FALSE keeps Δ in grams (curved but in g).
BIOMASS_LOG_Y           <- TRUE
BIOMASS_DELTA_LOG_RATIO <- FALSE

# ── 4.1  Data-building helpers ────────────────────────────────────────────────

build_upper_data <- function(post_df) {
  post_df %>%
    group_by(Symbiont, precip_mean) %>%
    summarise(
      med  = median(mu),
      lo95 = quantile(mu, 0.025),
      hi95 = quantile(mu, 0.975),
      .groups = "drop"
    ) %>%
    mutate(Symbiont = factor(Symbiont, levels = endo_levels))
}

# log_ratio = TRUE gives log(S+ / S-) per posterior draw (0 = no difference);
# FALSE gives the difference on the response scale (S+ - S-).
build_lower_data <- function(post_df, log_ratio = FALSE) {
  post_df %>%
    group_by(precip_mean, draw) %>%
    pivot_wider(names_from  = Symbiont,
                values_from = mu,
                id_cols     = c(precip_mean, draw)) %>%
    mutate(contrast = if (log_ratio) log(`S+`) - log(`S-`) else `S+` - `S-`) %>%
    group_by(precip_mean) %>%
    summarise(
      med  = median(contrast),
      lo95 = quantile(contrast, 0.025),
      hi95 = quantile(contrast, 0.975),
      .groups = "drop"
    )
}

build_obs_data <- function(raw_df, response_col) {
  raw_df %>%
    filter(!is.na(.data[[response_col]]), .data[[response_col]] > 0) %>%
    mutate(
      Symbiont = factor(Symbiont, levels = endo_levels),
      y_obs    = .data[[response_col]]
    ) %>%
    dplyr::select(precip_mean, Symbiont, y_obs)
}

# ── 4.2  Upper plot ───────────────────────────────────────────────────────────
make_upper_plot <- function(post_df, raw_df, response_col, y_label,
                            show_legend = TRUE, log_y = FALSE) {
  pd  <- build_upper_data(post_df)
  obs <- build_obs_data(raw_df, response_col)

  p <- ggplot() +
    geom_ribbon(
      data = pd,
      aes(x = precip_mean, ymin = lo95, ymax = hi95,
          fill = Symbiont, group = Symbiont),
      alpha = 0.2, color = NA
    ) +
    geom_line(
      data = pd,
      aes(x = precip_mean, y = med,
          color = Symbiont, group = Symbiont),
      linewidth = 0.6
    ) +
    geom_point(
      data = obs,
      aes(x = precip_mean, y = y_obs, color = Symbiont),
      alpha = 0.45, size = 2.5, shape = 16,
      show.legend = FALSE
    ) +
    scale_color_manual(values = ENDO_COLORS, breaks = endo_levels,
                       labels = ENDO_LABELS, name = "Symbiont") +
    scale_fill_manual(values  = ENDO_COLORS, breaks = endo_levels,
                      labels = ENDO_LABELS, name = "Symbiont") +
    labs(x = NULL, y = y_label) +
    vr_theme() +
    theme(
      legend.position = if (show_legend) c(0.2, 0.7) else "none",
      legend.key.size = unit(8, "pt"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank()
    )

  # Log y axis: a model linear in log(response) is drawn as a straight line.
  # Observed zeros are already removed in build_obs_data().
  if (log_y) {
    p <- p + scale_y_log10(
      labels = function(x) format(x, scientific = FALSE,
                                  drop0trailing = TRUE, trim = TRUE)
    )
  }
  p
}

# ── 4.3  Lower plot ───────────────────────────────────────────────────────────
make_lower_plot <- function(post_df,
                            x_label = "Soil origin MAP (mm/yr)",
                            log_ratio = FALSE) {
  pd <- build_lower_data(post_df, log_ratio = log_ratio)

  y_lab <- if (log_ratio) {
    expression(log ~ (italic(S)^"+" / italic(S)^"\u2212"))
  } else {
    expression(Delta ~ (italic(S)^"+" - italic(S)^"\u2212"))
  }

  ggplot(pd, aes(x = precip_mean)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95),
                fill = "#D9BFD6", alpha = 0.45, color = NA) +
    geom_line(aes(y = med), color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.4, color = "black") +
    labs(x = x_label, y = y_lab) +
    vr_theme() +
    theme(legend.position = "none")
}

# ── 4.4  Reusable column builder (kept for Part 4b / spikelets figure) ───────
make_vr_column <- function(post_df, raw_df, response_col, y_label,
                           show_legend = TRUE) {
  upper <- make_upper_plot(post_df, raw_df, response_col, y_label,
                           show_legend = show_legend)
  lower <- make_lower_plot(post_df)
  upper / lower + plot_layout(heights = c(2, 1))
}

# ── 4.5  Build predicted/Δ leaf plots ─────────────────────────────────────────
# Kept as flat (non-nested) plot objects, rather than pre-combined columns, so
# the final 2-row assembly in 4.7 can control every row's height explicitly
# instead of two different internal height ratios (2:1 vs 1:1) fighting when
# nested inside one outer plot_layout().
upper_biomass <- make_upper_plot(post_precip_biomass, aghysoils,
                                 "calc_abg_mass_tot", "Aboveground biomass (g)",
                                 show_legend = TRUE, log_y = BIOMASS_LOG_Y)
lower_biomass <- make_lower_plot(post_precip_biomass,
                                 log_ratio = BIOMASS_DELTA_LOG_RATIO)

upper_inflo <- make_upper_plot(post_precip_inflo, aghysoils,
                               "calc_total_inflo", "Total inflorescences",
                               show_legend = FALSE)
lower_inflo <- make_lower_plot(post_precip_inflo)

# ── 4.6  Coefficient caterpillar (saved as its own figure in 4.8) ────────────
# Focal betas from fit_biomass and fit_inflo:
#   beta[2] = SymbiontS+                 (main effect at mean precip)
#   beta[3] = precip_std                 (main effect of soil-origin precip)
#   beta[5] = SymbiontS+ × precip_std   (interaction; null expectation ~ 0)
# beta[1] (Intercept) and beta[4] (age_std) are nuisance — excluded.
FOCAL_IDX  <- c(2L, 3L, 5L)
FOCAL_LABS <- c(
  "beta[2]" = "Symbiont (S+)",
  "beta[3]" = "Precipitation (std)",
  "beta[5]" = "Symbiont \u00d7 Precipitation"
)

extract_focal_betas <- function(fit, model_label) {
  as.matrix(fit, pars = "beta")[, FOCAL_IDX, drop = FALSE] %>%
    as.data.frame() %>%
    setNames(paste0("beta[", FOCAL_IDX, "]")) %>%
    pivot_longer(everything(),
                 names_to  = "parameter",
                 values_to = "estimate") %>%
    mutate(model = model_label)
}

coef_long <- bind_rows(
  extract_focal_betas(fit_biomass, "Aboveground\nbiomass (g)"),
  extract_focal_betas(fit_inflo,   "Total\ninflorescences")
  # extract_focal_betas(fit_totspi,  "Total\nspikelets (projected)")
)

coef_summary <- coef_long %>%
  group_by(model, parameter) %>%
  summarise(
    median_est = median(estimate),
    lo90       = quantile(estimate, 0.05),
    hi90       = quantile(estimate, 0.95),
    lo95       = quantile(estimate, 0.025),
    hi95       = quantile(estimate, 0.975),
    prob_gt0   = mean(estimate > 0),
    prob_lt0   = mean(estimate < 0),
    .groups    = "drop"
  ) %>%
  mutate(
    strong_effect = prob_gt0 > 0.9 | prob_lt0 > 0.9,
    parameter = factor(parameter,
                       levels = rev(paste0("beta[", FOCAL_IDX, "]"))),
    model     = factor(model,
                       levels = c("Aboveground\nbiomass (g)",
                                  "Total\ninflorescences",
                                  "Total\nspikelets (projected)")),
    term_type = if_else(parameter == "beta[5]", "Interaction", "Main effect"),
    term_type = factor(term_type, levels = names(TERM_COLORS)),
    pt_fill   = if_else(strong_effect, as.character(term_type), "white")
  )

n_coef      <- nlevels(coef_summary$parameter)
shade_rects <- Filter(Negate(is.null), lapply(seq_len(n_coef), function(i) {
  if (i %% 2 == 1)
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = i - 0.49, ymax = i + 0.49,
             fill = "grey93", color = NA, alpha = 0.55)
}))

fig_coef <- ggplot(coef_summary, aes(y = parameter, color = term_type)) +
  shade_rects +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                linewidth = 0.35, width = 0) +
  geom_errorbar(aes(xmin = lo90, xmax = hi90),
                linewidth = 1.2, width = 0) +
  geom_point(aes(x = median_est, fill = pt_fill),
             shape = 21, size = 3.0, stroke = 0.9) +
  scale_fill_manual(values = c(TERM_COLORS, "white" = "white"),
                    guide  = "none") +
  scale_color_manual(
    values = TERM_COLORS,
    name   = NULL,
    guide  = guide_legend(
      override.aes = list(
        fill   = unname(TERM_COLORS),
        shape  = 21, size = 3, stroke = 0.9
      )
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey25", linewidth = 0.45) +
  facet_wrap(~ model, ncol = 1, scales = "free_x",
             strip.position = "right") +
  scale_y_discrete(labels = FOCAL_LABS, expand = expansion(add = 0.6)) +
  labs(x = "Posterior coefficient", y = NULL) +
  vr_theme() +
  theme(
    legend.position      = "none",
    legend.key.size      = unit(8, "pt"),
    panel.grid.major.y   = element_blank(),
    strip.text.y.right   = element_text(size = 8, color = "black", angle = 90,
                                        lineheight = 0.85),
    axis.text            = element_text(size = 8),
    axis.title           = element_text(size = 9)
  )

# ── 4.7  Assemble 2-row predicted/Δ figure with external panel tags ─────────
# Row 1 (upper, predicted):  upper_biomass | upper_inflo   (= a, b)
# Row 2 (lower, Δ contrast): lower_biomass | lower_inflo   (= c, d)
# The coefficient plot (previously panel e) is now a separate figure (4.8).
# Panel letters drawn exactly like panel_tag() in 05_Plot_vital_rate.R

# Panel letters drawn exactly like panel_tag() in 05_Plot_vital_rate.R
add_tag <- function(p, lab) {
  p + annotate("text", x = -Inf, y = Inf, label = lab,
               hjust = -0.6, vjust = 1.3, size = 4, fontface = "plain")
}

fig_greenhouse_combined <-
  (add_tag(upper_biomass +
             theme(legend.position      = c(0.98, 0.02),   # bottom-right of panel (a)
                   legend.justification = c(1, 0)), "(a)") |
     add_tag(upper_inflo, "(b)")) /
  (lower_biomass | lower_inflo) +
  plot_layout(heights = c(2, 1)) &
  theme(
    axis.title = element_text(size = 12),  # match vital-rate figures
    axis.text  = element_text(size = 8)
  )

print(fig_greenhouse_combined)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Combined.pdf"),
  plot   = fig_greenhouse_combined,
  width  = 160,
  height = 100,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Combined.tiff"),
  plot        = fig_greenhouse_combined,
  width       = 160,
  height      = 100,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# ── 4.8  Coefficient figure (standalone) ─────────────────────────────────────
# Same caterpillar plot as the former panel (e), saved on its own. The legend
# (main effect vs interaction) is shown here because there is no shared legend.
fig_coef_standalone <- fig_coef +
  theme(
    legend.position  = "none",
    legend.direction = "horizontal"
  )
print(fig_coef_standalone)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Coefficients.pdf"),
  plot   = fig_coef_standalone,
  width  = 110,
  height = 110,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Coefficients.tiff"),
  plot        = fig_coef_standalone,
  width       = 110,
  height      = 110,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# =============================================================================
# PART 4b — STANDALONE SPIKELETS FIGURE
# =============================================================================
# Two-panel figure for total spikelets (projected):
#   Left  — col_spike: upper (posterior predicted lines + observed points) /
#                      lower (Δ contrast S+ − S−), stacked 2:1
#   Right — fig_coef_spike: caterpillar plot for the three focal betas from
#                            fit_totspi (beta[2], beta[3], beta[5])
#
# This is the spikelet-only analogue of Fig_Greenhouse_Combined and is saved
# separately as Fig_Spikelets_Combined.pdf.
# =============================================================================

# ── 4b.1  Spikelet predicted column ──────────────────────────────────────────
col_spike <- make_vr_column(
  post_df      = post_precip_totspi,
  raw_df       = aghysoils,
  response_col = "calc_total_spik",
  y_label      = "Total spikelets",
  show_legend  = TRUE
)

# ── 4b.2  Spikelet coefficient caterpillar ────────────────────────────────────
# Extract focal betas from fit_totspi only.
coef_spike <- extract_focal_betas(fit_totspi, "Total\nspikelets (projected)")

coef_summary_spike <- coef_spike %>%
  group_by(model, parameter) %>%
  summarise(
    median_est = median(estimate),
    lo90       = quantile(estimate, 0.05),
    hi90       = quantile(estimate, 0.95),
    lo95       = quantile(estimate, 0.025),
    hi95       = quantile(estimate, 0.975),
    prob_gt0   = mean(estimate > 0),
    prob_lt0   = mean(estimate < 0),
    .groups    = "drop"
  ) %>%
  mutate(
    strong_effect = prob_gt0 > 0.9 | prob_lt0 > 0.9,
    parameter = factor(parameter,
                       levels = rev(paste0("beta[", FOCAL_IDX, "]"))),
    term_type = if_else(parameter == "beta[5]", "Interaction", "Main effect"),
    term_type = factor(term_type, levels = names(TERM_COLORS)),
    pt_fill   = if_else(strong_effect, as.character(term_type), "white")
  )

n_coef_spike  <- nlevels(coef_summary_spike$parameter)
shade_rects_spike <- Filter(Negate(is.null), lapply(seq_len(n_coef_spike), function(i) {
  if (i %% 2 == 1)
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = i - 0.49, ymax = i + 0.49,
             fill = "grey93", color = NA, alpha = 0.55)
}))

fig_coef_spike <- ggplot(coef_summary_spike, aes(y = parameter, color = term_type)) +
  shade_rects_spike +
  geom_errorbar(aes(xmin = lo95, xmax = hi95), linewidth = 0.35, width = 0) +
  geom_errorbar(aes(xmin = lo90, xmax = hi90), linewidth = 1.2,  width = 0) +
  geom_point(aes(x = median_est, fill = pt_fill),
             shape = 21, size = 3.0, stroke = 0.9) +
  scale_fill_manual(values = c(TERM_COLORS, "white" = "white"), guide = "none") +
  scale_color_manual(
    values = TERM_COLORS,
    name   = NULL,
    guide  = guide_legend(
      override.aes = list(fill = unname(TERM_COLORS),
                          shape = 21, size = 3, stroke = 0.9)
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey25", linewidth = 0.45) +
  scale_y_discrete(labels = FOCAL_LABS, expand = expansion(add = 0.6)) +
  labs(x = "Posterior coefficient (spikelets)", y = NULL) +
  vr_theme() +
  theme(
    legend.position    = "right",
    legend.key.size    = unit(8, "pt"),
    panel.grid.major.y = element_blank()
  )

# ── 4b.3  Assemble and save spikelet figure ───────────────────────────────────
fig_spikelets_combined <-
  (col_spike | fig_coef_spike) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    legend.position = "none",
    axis.text         = element_text(size = 10),
    plot.tag          = element_text(size = 12, face = "plain"),
    plot.tag.position = "topleft"
  )

print(fig_spikelets_combined)

ggsave(
  file.path(fig_dir, "Fig_Spikelets_Combined.pdf"),
  plot   = fig_spikelets_combined,
  width  = 173,
  height = 82.8,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "Fig_Spikelets_Combined.tiff"),
  plot        = fig_spikelets_combined,
  width       = 173,
  height      = 82.8,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# =============================================================================
# PART 5 — STATISTICAL SUMMARY TABLES
# =============================================================================

posterior_table <- function(fit, sd_obj, model_label) {
  coef_names <- colnames(sd_obj$X)
  coef_names <- gsub("\\(Intercept\\)",      "Intercept",              coef_names)
  coef_names <- gsub("SymbiontS\\+",         "Sym[S+]",                coef_names)
  coef_names <- gsub("precip_std",           "Precip (std)",           coef_names)
  coef_names <- gsub("Sym\\[S\\+\\]:Precip", "Sym[S+] \u00d7 Precip", coef_names)
  coef_names <- gsub("age_std",              "Age (std)",              coef_names)

  draws_beta <- as.data.frame(as.matrix(fit, pars = "beta"))
  colnames(draws_beta) <- coef_names

  draws_beta %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "draw") %>%
    group_by(Parameter) %>%
    summarise(
      Median  = round(median(draw), 3),
      Mean    = round(mean(draw),   3),
      SD      = round(sd(draw),     3),
      `Q2.5`  = round(quantile(draw, 0.025), 3),
      `Q97.5` = round(quantile(draw, 0.975), 3),
      `P(>0)` = round(mean(draw > 0), 3),
      .groups = "drop"
    ) %>%
    mutate(Model = model_label, .before = 1)
}

tab_biomass <- posterior_table(fit_biomass, sd_biomass, "Biomass (Gaussian)")
tab_inflo   <- posterior_table(fit_inflo,   sd_inflo,   "Inflo count (NB)")
tab_totspi  <- posterior_table(fit_totspi,  sd_totspi,  "Total spikelets projected (NB)")
tab_avgspi  <- posterior_table(fit_avgspi,  sd_avgspi,  "Avg spikelets (NB)")

all_tabs <- bind_rows(tab_biomass, tab_inflo, tab_totspi, tab_avgspi)

cat("\n\n")
cat("=============================================================\n")
cat(" BAYESIAN POSTERIOR SUMMARY — All models\n")
cat(" Key parameter: Sym[S+] x Precip\n")
cat(" Null expectation for greenhouse: this coefficient ~ 0\n")
cat(" (soil-origin precip does not modify endophyte effect\n")
cat("  when plants grow in a common environment)\n")
cat("=============================================================\n\n")
for (mod_label in unique(all_tabs$Model)) {
  cat("--- ", mod_label, " ---\n", sep = "")
  all_tabs %>% filter(Model == mod_label) %>% dplyr::select(-Model) %>% print(n = Inf)
  cat("\n")
}

write.csv(all_tabs,
          file      = file.path(output_dir, "Table1_PosteriorSummary_precip.csv"),
          row.names = FALSE)


# ── Full breakdown (for supplement) ───────────────────────────────────────────

ss_biomass <- aghysoils %>%
  filter(!is.na(calc_abg_mass_tot)) %>%
  count(Site, Symbiont) %>%
  mutate(Response = "Biomass")

ss_inflo <- aghysoils %>%
  filter(!is.na(calc_total_inflo)) %>%
  count(Site, Symbiont) %>%
  mutate(Response = "Inflorescence")

ss_spike <- aghysoils %>%
  filter(!is.na(calc_avg_spikelet)) %>%
  count(Site, Symbiont) %>%
  mutate(Response = "Avg spikelet")

sample_size_full <- bind_rows(ss_biomass, ss_inflo, ss_spike) %>%
  tidyr::complete(Response, Site, Symbiont, fill = list(n = 0)) %>%
  arrange(Response, Site)

write.csv(sample_size_full,
          file = file.path(output_dir, "Table_Sx_SampleSizes.csv"),
          row.names = FALSE)

# =============================================================================
# PART 6 — SENSITIVITY ANALYSIS: DROP KER & SON
# =============================================================================
# Tom's question on the soil-origin precipitation result (plants from wetter-
# origin soils had more biomass/inflorescences even under common greenhouse
# conditions): is that signal being driven by KER and SON specifically (e.g.
# their limestone-derived soils), or does a precipitation/climate signal
# persist even without those two sites?
#
# We refit each response model (same formula, same priors, same compiled Stan
# programs as PART 2) on the subset of the data with KER and SON removed, and
# compare the posterior for the precip_std main effect (beta[3]) and the
# Symbiont x precip_std interaction (beta[5]) to the full-data models fit
# earlier. If the precip_std effect holds up (median away from 0, similar
# sign/magnitude) with KER/SON dropped, that supports a genuine climate
# signal rather than an artifact of those two sites' soils.
# =============================================================================

aghysoils_sens <- aghysoils %>%
  filter(!Site %in% c("KER", "SON")) %>%
  droplevels()

message(sprintf(
  "Sensitivity subset (KER & SON dropped): %d of %d rows retained",
  nrow(aghysoils_sens), nrow(aghysoils)
))

# ── 6.1  Stan data + refit (reuses the same compiled models from PART 2) ─────
sd_biomass_sens <- make_stan_data(aghysoils_sens, "calc_abg_mass_tot", "gaussian")
sd_inflo_sens   <- make_stan_data(aghysoils_sens, "calc_total_inflo",  "negbinom")
sd_totspi_sens  <- make_stan_data(aghysoils_sens, "calc_total_spik",   "negbinom")
sd_avgspi_sens  <- make_stan_data(
  aghysoils_sens %>% filter(!is.na(calc_avg_spikelet), calc_avg_spikelet > 0),
  "calc_avg_spikelet", "negbinom"
)

fit_biomass_sens <- run_stan(mod_gaussian, sd_biomass_sens, "gaussian")
fit_inflo_sens   <- run_stan(mod_negbinom, sd_inflo_sens,   "negbinom")
fit_totspi_sens  <- run_stan(mod_negbinom, sd_totspi_sens,  "negbinom")
fit_avgspi_sens  <- run_stan(mod_negbinom, sd_avgspi_sens,  "negbinom")

check_hmc(fit_biomass_sens, "M1s Biomass (no KER/SON)")
check_hmc(fit_inflo_sens,   "M2s Inflo (no KER/SON)")
check_hmc(fit_totspi_sens,  "M3s Total spikelets (no KER/SON)")
check_hmc(fit_avgspi_sens,  "M4s Avg spikelets (no KER/SON)")

# ── 6.2  Posterior summary table (same layout as Table 1) ────────────────────
tab_biomass_sens <- posterior_table(fit_biomass_sens, sd_biomass_sens,
                                    "Biomass, no KER/SON (Gaussian)")
tab_inflo_sens   <- posterior_table(fit_inflo_sens,   sd_inflo_sens,
                                    "Inflo count, no KER/SON (NB)")
tab_totspi_sens  <- posterior_table(fit_totspi_sens,  sd_totspi_sens,
                                    "Total spikelets, no KER/SON (NB)")
tab_avgspi_sens  <- posterior_table(fit_avgspi_sens,  sd_avgspi_sens,
                                    "Avg spikelets, no KER/SON (NB)")

all_tabs_sens <- bind_rows(tab_biomass_sens, tab_inflo_sens,
                           tab_totspi_sens, tab_avgspi_sens)

cat("\n\n")
cat("=============================================================\n")
cat(" SENSITIVITY ANALYSIS \u2014 KER & SON dropped\n")
cat(" Question: does the soil-origin precip signal survive without\n")
cat(" KER/SON (i.e. is it just their limestone-derived soils), or is\n")
cat(" there a climate signal even excluding those two sites?\n")
cat("=============================================================\n\n")
for (mod_label in unique(all_tabs_sens$Model)) {
  cat("--- ", mod_label, " ---\n", sep = "")
  all_tabs_sens %>% filter(Model == mod_label) %>% dplyr::select(-Model) %>% print(n = Inf)
  cat("\n")
}

write.csv(all_tabs_sens,
          file      = file.path(output_dir, "TableS_PosteriorSummary_precip_noKERSON.csv"),
          row.names = FALSE)

# ── 6.3  Full-data vs. sensitivity coefficient comparison figure ────────────
# Side-by-side caterpillar plot (biomass + inflorescence, the two responses
# Tom asked about) contrasting "All sites" vs "No KER/SON" posteriors for
# beta[2] (Symbiont), beta[3] (Precip), and beta[5] (Symbiont x Precip).
# Uses one flat color for all terms (see TERM_COLORS note above) so the
# comparison is read purely off dataset, not term type.
coef_full_vs_sens <- bind_rows(
  extract_focal_betas(fit_biomass,      "Aboveground\nbiomass (g)")  %>% mutate(dataset = "All sites"),
  extract_focal_betas(fit_inflo,        "Total\ninflorescences")    %>% mutate(dataset = "All sites"),
  extract_focal_betas(fit_biomass_sens, "Aboveground\nbiomass (g)")  %>% mutate(dataset = "No KER/SON"),
  extract_focal_betas(fit_inflo_sens,   "Total\ninflorescences")    %>% mutate(dataset = "No KER/SON")
)

coef_summary_sens <- coef_full_vs_sens %>%
  group_by(dataset, model, parameter) %>%
  summarise(
    median_est = median(estimate),
    lo90       = quantile(estimate, 0.05),
    hi90       = quantile(estimate, 0.95),
    lo95       = quantile(estimate, 0.025),
    hi95       = quantile(estimate, 0.975),
    .groups    = "drop"
  ) %>%
  mutate(
    parameter = factor(parameter, levels = rev(paste0("beta[", FOCAL_IDX, "]"))),
    dataset   = factor(dataset, levels = c("All sites", "No KER/SON"))
  )

fig_sensitivity <- ggplot(coef_summary_sens,
                          aes(y = parameter, x = median_est, color = dataset)) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                position = position_dodge(width = 0.5),
                linewidth = 0.35, width = 0) +
  geom_errorbar(aes(xmin = lo90, xmax = hi90),
                position = position_dodge(width = 0.5),
                linewidth = 1.2, width = 0) +
  geom_point(position = position_dodge(width = 0.5), size = 2.6) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey25", linewidth = 0.45) +
  scale_color_manual(values = c("All sites" = "#4D4D4D", "No KER/SON" = "#0072B2"),
                     name = NULL) +
  facet_wrap(~ model, ncol = 1, scales = "free_x", strip.position = "right") +
  scale_y_discrete(labels = FOCAL_LABS, expand = expansion(add = 0.6)) +
  labs(x = "Posterior coefficient", y = NULL,
       title = "Sensitivity: all sites vs. KER & SON dropped") +
  vr_theme() +
  theme(
    legend.position     = "top",
    legend.key.size     = unit(8, "pt"),
    panel.grid.major.y  = element_blank(),
    strip.text.y.right  = element_text(size = 8, color = "black", angle = 90,
                                       lineheight = 0.85)
  )

print(fig_sensitivity)

ggsave(
  file.path(fig_dir, "FigS_Sensitivity_noKERSON.pdf"),
  plot   = fig_sensitivity,
  width  = 140,
  height = 120,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "FigS_Sensitivity_noKERSON.tiff"),
  plot        = fig_sensitivity,
  width       = 140,
  height      = 120,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# =============================================================================
# End of script
# =============================================================================
