# =============================================================================
# DOES SOIL INTERACT WITH SYMBIONT STATUS TO INFLUENCE DEMOGRAPHY?
# Full Bayesian Analysis via rstan
# Author: Jacob Moutouama
#
# REVISION NOTES (addressing coauthor comments):
# 1. Added 30-year precipitation normals (mm/yr) for each soil-origin site as a
#    fixed-effect covariate, replacing the categorical Site factor. This mirrors
#    the field analysis structure and allows direct comparison between greenhouse
#    and field results. The null expectation is that the Symbiont × Precip
#    interaction disappears in the greenhouse (soil origin drives demography
#    independently of endophyte status, because plants are not experiencing
#    that climate directly).
# 2. Clarified that realised N < 196 planted because some plants died or failed
#    to germinate before harvest.
# 3. Justified Student-t likelihood for biomass (see Section 2.2 comment).
# 4. Reduced MCMC iterations from iter=6000/warmup=2000 to iter=2000/warmup=1000
#    — adequate for this sample size and model complexity.
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
ENDO_COLORS <- c("S-" = "tomato", "S+" = "cornflowerblue")
ENDO_LABELS <- c("S-" = "S\u2212", "S+" = "S+")
endo_levels <- c("S-", "S+")
site_order  <- c("SON", "KER", "BFL", "BAS", "COL", "HUN", "LAF")

# ── Shared theme (vr_theme) ───────────────────────────────────────────────────
vr_theme <- function() {
  theme_classic() +
    theme(
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line         = element_line(color = "black", linewidth = 0.1),
      axis.title        = element_text(size = 14),
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
                                       margin = margin(r = 4))
    )
}

# ── Output directories ────────────────────────────────────────────────────────
output_dir <- "/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/"
fig_dir    <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure"

save_fig <- function(p, filename, width_mm = 200, height_mm = 110) {
  ggsave(file.path(fig_dir, filename), plot = p,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, device = "pdf")
  message("Saved: ", filename)
}

TERM_COLORS <- c(
  "Main effect" = "#4D4D4D",
  "Interaction" = "#D55E00"
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
  scale_colour_manual(values = ENDO_COLORS, labels = ENDO_LABELS,
                      name = "Symbiont") +
  scale_fill_manual(values   = ENDO_COLORS, labels = ENDO_LABELS,
                    name = "Symbiont") +
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
  y_rep     <- as.matrix(fit, pars = "y_rep")
  y_rep_sub <- y_rep[sample(nrow(y_rep), 200), ]
  p <- ppc_dens_overlay(y_obs, y_rep_sub) +
    labs(title = paste0("PPC \u2014 ", label)) +
    vr_theme()
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}

p_ppc_biomass <- ppc_check(fit_biomass, sd_biomass$y_cont, "Biomass (log scale)")
p_ppc_inflo   <- ppc_check(fit_inflo,   sd_inflo$y,        "Inflo count", xlim = c(-1, 50))
p_ppc_avgspi  <- ppc_check(fit_avgspi,  sd_avgspi$y,       "Avg spikelets")

ppc_panel <- (p_ppc_biomass | p_ppc_inflo | p_ppc_avgspi) +
  plot_annotation(tag_levels = "A")
print(ppc_panel)

# =============================================================================
# PART 3 — POSTERIOR PREDICTIONS
# =============================================================================

# ── 3.1  Extract posterior fitted means across precipitation gradient ─────────
# For each Symbiont level, predict across the observed range of precip_std,
# marginalising over age_std (set to 0 = mean age) and population random effects.
extract_posterior_precip <- function(fit, sd_obj, link = "identity",
                                     n_grid = 50) {
  draws_beta <- as.matrix(fit, pars = "beta")
  draws_upop <- as.matrix(fit, pars = "u_pop")
  d          <- sd_obj$d_used
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
        precip_mean = p * sd(d$precip_mean) + mean(d$precip_mean),
        draw        = seq_len(n_iter),
        mu          = as.numeric(mu)
      )
    })
  })
}

post_precip_biomass <- extract_posterior_precip(fit_biomass, sd_biomass, "identity") %>%
  mutate(mu = exp(mu))
post_precip_inflo   <- extract_posterior_precip(fit_inflo,   sd_inflo,   "log")
post_precip_totspi  <- extract_posterior_precip(fit_totspi,  sd_totspi,  "log")
post_precip_avgspi  <- extract_posterior_precip(fit_avgspi,  sd_avgspi,  "log")

# =============================================================================
# PART 4 — FIGURES
# =============================================================================
# Layout: 3-column combined figure.
#   Col 1 — Biomass:      upper (predicted) / lower (Δ contrast), 2:1 height
#   Col 2 — Inflo:        upper (predicted) / lower (Δ contrast), 2:1 height
#   Col 3 — Coefficients: caterpillar for focal betas, faceted by response
#
# Panel tags (A–E) placed OUTSIDE plots via plot_annotation(tag_levels).
# Each sub-plot carries tag = "" so patchwork assigns the letter externally.
# =============================================================================

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

build_lower_data <- function(post_df) {
  post_df %>%
    group_by(precip_mean, draw) %>%
    pivot_wider(names_from  = Symbiont,
                values_from = mu,
                id_cols     = c(precip_mean, draw)) %>%
    mutate(contrast = `S+` - `S-`) %>%
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
                            show_legend = TRUE) {
  pd  <- build_upper_data(post_df)
  obs <- build_obs_data(raw_df, response_col)

  ggplot() +
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
      alpha = 0.45, size = 2, shape = 16,
      show.legend = FALSE
    ) +
    scale_color_manual(values = ENDO_COLORS, labels = ENDO_LABELS,
                       name = "Symbiont") +
    scale_fill_manual(values  = ENDO_COLORS, labels = ENDO_LABELS,
                      name = "Symbiont") +
    labs(x = NULL, y = y_label) +
    vr_theme() +
    theme(
      legend.position = if (show_legend) c(0.4, 0.8) else "none",,
      legend.key.size = unit(8, "pt"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank()
    )
}

# ── 4.3  Lower plot ───────────────────────────────────────────────────────────
make_lower_plot <- function(post_df,
                            x_label = "Soil origin MAP (mm/yr)") {
  pd <- build_lower_data(post_df)

  ggplot(pd, aes(x = precip_mean)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95),
                fill = "#9B6B96", alpha = 0.45, color = NA) +
    geom_line(aes(y = med), color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.4, color = "black") +
    labs(x = x_label, y = "\u0394 (S+ \u2212 S\u2212)") +
    vr_theme() +
    theme(legend.position = "none")
}

# ── 4.4  Stack upper + lower into one labelled column ────────────────────────
# Tags are assigned via plot_annotation() on the outer combined figure,
# so each individual plot gets tag = "" here to avoid double-labelling.
# The two sub-plots within each column receive consecutive letters via
# the nested tagging in the final assembly below.
make_vr_column <- function(post_df, raw_df, response_col, y_label,
                           show_legend = TRUE) {
  upper <- make_upper_plot(post_df, raw_df, response_col, y_label,
                           show_legend = show_legend)
  lower <- make_lower_plot(post_df)
  upper / lower + plot_layout(heights = c(2, 1))
}

# ── 4.5  Build predicted/Δ columns ───────────────────────────────────────────
col_biomass <- make_vr_column(
  post_df      = post_precip_biomass,
  raw_df       = aghysoils,
  response_col = "calc_abg_mass_tot",
  y_label      = "Aboveground biomass (g)",
  show_legend  = TRUE
)

col_inflo <- make_vr_column(
  post_df      = post_precip_inflo,
  raw_df       = aghysoils,
  response_col = "calc_total_inflo",
  y_label      = "Total inflorescences",
  show_legend  = FALSE
)

# ── 4.6  Coefficient caterpillar (third column) ───────────────────────────────
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
                                  "Total\ninflorescences")),
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
    strip.text.y.right   = element_text(size = 14, color = "black", angle = 90)
  )

# ── 4.7  Assemble 3-column figure with external panel tags ───────────────────
# col_biomass contains 2 sub-plots (upper = A, lower = B)
# col_inflo   contains 2 sub-plots (upper = C, lower = D)
# fig_coef    is a single plot           (        = E)
#
# plot_annotation(tag_levels = "A") walks every leaf plot in patchwork order
# and assigns A, B, C, D, E from left to right, top to bottom.
# tag_prefix / tag_suffix can be added if brackets are preferred, e.g. "(A)".
# plot.tag.position = "topleft" puts each letter outside and above the panel.

col_biomass <- col_biomass + plot_annotation(tag_levels = list("A"))
col_inflo   <- col_inflo   + plot_annotation(tag_levels = list("B"))
fig_coef    <- fig_coef    + plot_annotation(tag_levels = list("C"))

fig_greenhouse_combined <-
  (col_biomass | col_inflo | fig_coef) +
  plot_layout(widths = c(1, 1, 1.4)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    axis.text         = element_text(size = 10),
    plot.tag = element_text(size = 12, face = "plain"),
    plot.tag.position = "topleft"
  )

print(fig_greenhouse_combined)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Combined.pdf"),
  plot   = fig_greenhouse_combined,
  width  = 277,
  height = 120,
  units  = "mm",
  device = cairo_pdf
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

# =============================================================================
# End of script
# =============================================================================
