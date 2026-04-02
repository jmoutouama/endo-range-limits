# =============================================================================
# DOES SOIL INTERACT WITH ENDOPHYTE STATUS TO INFLUENCE DEMOGRAPHY?
# Full Bayesian Analysis via rstan
# Author: Jacob Moutouama
# Structure:
#   PART 1 — Data preparation & exploratory diagnostics (no statistics)
#   PART 2 — Stan model building & sampling
#   PART 3 — Publication-quality figures (Ecology Letters standard)
#   PART 4 — Statistical summary tables
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────
library(tidyverse)
library(rstan)
library(posterior)      # tidy draws from stanfit
library(bayesplot)      # MCMC diagnostics
library(patchwork)      # combine ggplots
library(ggdist)         # halfeye / rain-cloud distributions

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ── Colour palette ────────────────────────────────────────────────────────────
# Only endo_col kept (no endo_est_col per request)
endo_col    <- c("E-" = "tomato", "E+" = "cornflowerblue")
endo_levels <- c("E-", "E+")
site_order  <- c("SON","KER","BFL","BAS","COL","HUN","LAF")

# ── ggplot theme for Ecology Letters ─────────────────────────────────────────
# EL: single column = 84 mm, double = 174 mm. Font ≥ 7 pt, preferably Arial/Helvetica.
theme_EL <- function(base_size = 9, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # axes
      axis.line        = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      # panel
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      # legend
      legend.position  = "right",
      legend.key.size  = unit(8, "pt"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      # strip (facets)
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold"),
      # plot tag (A, B, C …)
      plot.tag         = element_text(size = base_size + 1, face = "bold"),
      plot.title       = element_text(size = base_size, face = "bold",
                                       hjust = 0),
      plot.margin      = margin(4, 4, 4, 4, "pt")
    )
}

# ── Output directory ──────────────────────────────────────────────────────────
output_dir <- "/Users/jacobmoutouama/Desktop/Soil/output/"
fig_dir    <- file.path(output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# helper: save at EL double-column width (174 mm × 130 mm)
save_EL <- function(p, filename, width_mm = 174, height_mm = 130) {
  ggsave(file.path(fig_dir, filename), plot = p,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, device = "pdf")
  message("Saved: ", filename)
}


# =============================================================================
# PART 1 — DATA PREPARATION & EXPLORATORY DIAGNOSTICS
# =============================================================================
# Goal: reveal patterns in the raw data WITHOUT any inferential statistics.
# These plots are for the analyst — they guide model choice, not the paper.

# ── 1.1  Load & clean ─────────────────────────────────────────────────────────
soils <- read.csv(
  "/Users/jacobmoutouama/Desktop/Soil/reproduction_and_biomass - endo&soil.csv",
  stringsAsFactors = FALSE
) %>%
  rename(Endo = "Endo") %>%               # strip leading space
  mutate(
    across(c(abg_mass_tot, total_inflo,
             tot_spikelet, avg_spikelet), ~ suppressWarnings(as.numeric(.))),
    Birthday  = as.Date(Birthday, "%m/%d/%Y"),
    abg_mass_tot=log(abg_mass_tot),
    age_days  = as.numeric(as.Date("2025-01-01") - Birthday),
    age_std   = as.numeric(scale(age_days)),
    Site      = factor(Site, levels = site_order),
    Endo      = factor(Endo, levels = endo_levels),
    Pop       = as.factor(Pop),
    Tray      = as.factor(Tray)
  )

aghysoils <- soils %>% filter(Species == "AGHY")

message("AGHY n = ", nrow(aghysoils),
        "  |  E-: ", sum(aghysoils$Endo == "E-"),
        "  |  E+: ", sum(aghysoils$Endo == "E+"))

# ── 1.2  Sample sizes per Site × Endo ────────────────────────────────────────
sample_grid <- aghysoils %>%
  count(Site, Endo, name = "n") %>%
  complete(Site, Endo, fill = list(n = 0))

p_n <- ggplot(sample_grid, aes(x = Site, y = n, fill = Endo)) +
  geom_col(position = position_dodge(0.7), width = 0.65,
           colour = "white", linewidth = 0.25) +
  geom_text(aes(label = n), position = position_dodge(0.7),
            vjust = -0.4, size = 2.2, family = "sans") +
  scale_fill_manual(values = endo_col, name = "Endophyte") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Soil origin (site)", y = "Sample size",
       title = "Sample sizes per Site × Endophyte") +
  theme_EL()

# ── 1.3  Distribution of each response (log-scale where skewed) ───────────────
resp_long <- aghysoils %>%
  dplyr::select(Site, Endo, abg_mass_tot, total_inflo,
                tot_spikelet, avg_spikelet) %>%
  pivot_longer(-c(Site, Endo),
               names_to  = "response",
               values_to = "value") %>%
  dplyr::mutate(response = dplyr::recode(response,
                                         "abg_mass_tot" = "Aboveground biomass (g)",
                                         "total_inflo"  = "Total inflorescences",
                                         "tot_spikelet" = "Total spikelets",
                                         "avg_spikelet" = "Avg spikelets / inflo"))
p_dist <- ggplot(resp_long %>% filter(!is.na(value)),
                 aes(x = value, fill = Endo, colour = Endo)) +
  geom_density(alpha = 0.45, linewidth = 0.3, adjust = 1.2) +
  scale_fill_manual(values  = endo_col, name = "Endophyte") +
  scale_colour_manual(values = endo_col, name = "Endophyte") +
  facet_wrap(~ response, scales = "free", ncol = 2) +
  labs(x = "Observed value", y = "Density",
       title = "Response distributions by endophyte status") +
  theme_EL() +
  theme(legend.position = "top")

# ── 1.4  Raw means ± SE across sites (rain-cloud style) ──────────────────────
# One panel per response, showing the raw Site × Endo pattern
p_raw <- ggplot(resp_long %>% filter(!is.na(value)),
                aes(x = Site, y = value, fill = Endo, colour = Endo)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "pointrange",
               position = position_dodge(0.55),
               size = 0.35, linewidth = 0.5, shape = 21,
               colour = "black",
               aes(fill = Endo)) +
  scale_fill_manual(values  = endo_col, name = "Endophyte") +
  scale_colour_manual(values = endo_col, name = "Endophyte") +
  facet_wrap(~ response, scales = "free_y", ncol = 2) +
  labs(x = "Soil origin (site)", y = "Observed value (mean ± SE)",
       title = "Raw means ± SE — no statistics, patterns only") +
  theme_EL() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 35, hjust = 1))

# ── 1.5  Age covariate vs responses ──────────────────────────────────────────
resp_long <- resp_long %>%
  left_join(aghysoils %>% dplyr::select(Site, Endo, age_days),
            by = c("Site", "Endo"))

p_age <- ggplot(resp_long %>% filter(!is.na(value)),
                aes(x = age_days,
                    y = value,
                    colour = Endo)) +
  geom_point(alpha = 0.35, size = 0.7, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, formula = y ~ x) +
  scale_colour_manual(values = endo_col, name = "Endophyte") +
  facet_wrap(~ response, scales = "free_y", ncol = 2) +
  labs(x = "Age (days from transplant to 1 Jan 2025)",
       y = "Observed value",
       title = "Age covariate — checking for confounding") +
  theme_EL() +
  theme(legend.position = "top")

p_age

# ── Print Part 1 plots ────────────────────────────────────────────────────────
# (these are diagnostic; not saved at publication quality)
print(p_n)
print(p_dist)
print(p_raw)


# =============================================================================
# PART 2 — STAN MODEL BUILDING & SAMPLING
# =============================================================================
# Four models, one per response:
#   M1 — abg_mass_tot  : Gaussian (log-normal if right-skewed after check)
#   M2 — total_inflo   : Negative-Binomial
#   M3 — tot_spikelet  : Negative-Binomial
#   M4 — avg_spikelet  : Negative-Binomial
#
# Fixed effects: Endo × Site + age_std
# Random effect: Pop (4 seed-source populations; partial pooling)
# Priors: weakly informative — Normal(0,1) on standardised predictors,
#         Normal(0,0.5) on SD hyperparameters, Exponential(1) on phi (NB)

# ── 2.1  Design matrices ──────────────────────────────────────────────────────
# Remove rows missing both predictor and response
make_stan_data <- function(df, response_col, family = "negbinom") {
  d <- df %>%
    filter(!is.na(.data[[response_col]]),
           !is.na(Endo), !is.na(Site), !is.na(age_std)) %>%
    droplevels()

  # Dummy coding: reference = E-, SON
  X <- model.matrix(~ Endo * Site + age_std, data = d)

  list(
    N          = nrow(d),
    K          = ncol(X),
    X          = X,
    y          = d[[response_col]],
    y_cont     = d[[response_col]],   # used by Gaussian model
    n_pop      = nlevels(d$Pop),
    pop_id     = as.integer(d$Pop),
    family     = family,
    d_used     = d                    # keep data aligned with model
  )
}

sd_biomass <- make_stan_data(aghysoils, "abg_mass_tot",  "gaussian")
sd_inflo   <- make_stan_data(aghysoils, "total_inflo",   "negbinom")
sd_totspi  <- make_stan_data(aghysoils, "tot_spikelet",  "negbinom")
sd_avgspi  <- make_stan_data(aghysoils, "avg_spikelet",  "negbinom")


# ── 2.2  Stan model: Gaussian (biomass) ──────────────────────────────────────
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
  vector[K] beta;
  real<lower=0> sigma;
  vector[n_pop] z_pop;
  real<lower=0> sigma_pop;
}

transformed parameters {
  vector[n_pop] u_pop;
  u_pop = sigma_pop * z_pop;
}

model {
  beta ~ normal(0,1);
  sigma ~ exponential(1);
  sigma_pop ~ exponential(1);
  z_pop ~ normal(0,1);

  y_cont ~ student_t(3,X * beta + u_pop[pop_id], sigma);
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    real mu_n = X[n] * beta + u_pop[pop_id[n]];
    y_rep[n]    = normal_rng(mu_n, sigma);
    log_lik[n]  = normal_lpdf(y_cont[n] | mu_n, sigma);
  }
}
"

# ── 2.3  Stan model: Negative-Binomial (count responses) ─────────────────────
#
# FIX 1 — Divergences (fit_inflo):
#   Root cause: phi ~ exponential(2) allows phi → 0, creating a funnel in the
#   (phi, beta) joint posterior that HMC cannot traverse without divergences.
#   This is especially severe for inflorescence count which has many zeros
#   (KER site) and high variance.
#   Solution:
#     • Add lower=0.1 bound on phi in the parameters block — keeps the sampler
#       away from the degenerate region without strongly biasing the posterior.
#     • Tighten phi prior to gamma(2, 0.1), which has a mode at 10 and a long
#       right tail — more realistic for ecology count data than exponential(2).
#     • Use log(phi) parameterisation (log_phi unconstrained, phi = exp(log_phi))
#       so HMC sees a flat-ish geometry near phi ≈ 0 rather than a hard wall.
#
# FIX 2 — y_rep undefined (fit_totspi):
#   Root cause: neg_binomial_2_log_rng can sample arbitrarily large integers.
#   When mu is large (total spikelets can exceed 1000) and phi is small, the
#   RNG occasionally returns a value that overflows Stan's 32-bit signed int
#   (max 2^31 - 1 ≈ 2.1e9), producing an undefined/NA result.
#   Solution:
#     • Pass y_max (= max(y) * 10) as a data item — a generous safety ceiling.
#     • In generated quantities, clamp y_rep to y_max if the RNG draw would
#       overflow. This keeps y_rep valid for PPC without biasing the summary.
#     • Use real arithmetic for the RNG draw then cast, so overflow is caught
#       before assignment to int.

stan_negbinom_code <- "
data {
  int<lower=1> N;                       // number of observations
  int<lower=1> K;                       // number of predictors
  matrix[N, K] X;                       // design matrix
  array[N] int<lower=0> y;              // response variable
  int<lower=0> y_max;                   // FIX 2: safety ceiling for y_rep
  int<lower=1> n_pop;                   // number of populations
  array[N] int<lower=1, upper=n_pop> pop_id;  // population index
}
parameters {
  vector[K] beta;                       // fixed effects
  real log_phi;                         // FIX 1: unconstrained log-overdispersion
  vector[n_pop] z_pop;                  // non-centered RE (implies u_pop ~ N(0,sigma_pop))
  real<lower=0> sigma_pop;              // RE SD
}
transformed parameters {
  real<lower=0> phi;
  vector[n_pop] u_pop;
  phi   = exp(log_phi);                 // FIX 1: always positive, flat geometry
  u_pop = sigma_pop * z_pop;
}
model {
  // PRIORS
  beta      ~ normal(0, 1);
  log_phi   ~ normal(2, 1);             // FIX 1: prior on log scale; exp(2)≈7 is
                                        // a sensible centre for ecology count data
  sigma_pop ~ exponential(2);
  z_pop     ~ normal(0, 1);
  // LINEAR PREDICTOR
  vector[N] log_mu;
  for (n in 1:N)
    log_mu[n] = X[n] * beta + u_pop[pop_id[n]];
  // LIKELIHOOD
  y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  array[N] int y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    real log_mu_n = X[n] * beta + u_pop[pop_id[n]];
    // FIX 2: clamp log_mu before RNG to avoid integer overflow on large counts
    real safe_log_mu = fmin(log_mu_n, log(1e8));
    real draw_real   = neg_binomial_2_log_rng(safe_log_mu, phi);
    // Guard: clamp before casting — Stan has no C-style (int) cast
    y_rep[n]   = to_int(round(fmin(draw_real, y_max)));
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_mu_n, phi);
  }
}
"

# ── 2.4  Compile models ───────────────────────────────────────────────────────
mod_gaussian <- stan_model(model_code = stan_gaussian_code,
                           model_name = "gaussian_endo_soil")
mod_negbinom <- stan_model(model_code = stan_negbinom_code,
                           model_name = "negbinom_endo_soil")

# ── 2.5  Sampling helper ──────────────────────────────────────────────────────
# adapt_delta raised to 0.99 for NB models: log_phi reparameterisation helps
# geometry but the interaction-heavy design matrix benefits from a smaller step
# size. max_treedepth raised to 12 to give HMC room to find good proposals.
run_stan <- function(stan_mod, data_list, family, seed = 13,
                     chains = 4, iter = 6000, warmup = 2000,
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
    # FIX 2: generous ceiling — 10x observed max — passed to Stan to cap y_rep
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
           chains = chains, iter = iter, warmup = warmup,
           seed = seed,
           control = list(adapt_delta   = adapt_delta,
                          max_treedepth = 12),
           refresh = 500)
}

# ── 2.6  Fit all four models ──────────────────────────────────────────────────
message("\nFitting M1: Biomass (Gaussian) …")
fit_biomass <- run_stan(mod_gaussian, sd_biomass, "gaussian")

message("\nFitting M2: Inflorescence count (NB) …")
fit_inflo   <- run_stan(mod_negbinom, sd_inflo,  "negbinom")

message("\nFitting M3: Total spikelets (NB) …")
fit_totspi  <- run_stan(mod_negbinom, sd_totspi, "negbinom")

message("\nFitting M4: Average spikelets (NB) …")
fit_avgspi  <- run_stan(mod_negbinom, sd_avgspi, "negbinom")

message("\nAll models fitted.")

# ── 2.7  MCMC diagnostics ────────────────────────────────────────────────────
check_hmc <- function(fit, label) {
  cat("\n===", label, "===\n")
  cat("Rhat summary (should be < 1.01):\n")
  rhat_vals <- summary(fit)$summary[, "Rhat"]
  print(summary(rhat_vals))
  n_div <- get_num_divergent(fit)
  n_treedepth <- get_num_max_treedepth(fit)
  cat(sprintf("Divergent transitions: %d\n", n_div))
  cat(sprintf("Max treedepth exceeded: %d\n", n_treedepth))
}

check_hmc(fit_biomass, "M1 Biomass")
check_hmc(fit_inflo,   "M2 Inflo")
check_hmc(fit_totspi,  "M3 Total spikelets")
check_hmc(fit_avgspi,  "M4 Avg spikelets")

# Trace plots for beta parameters (first 6)
p_trace_biomass <- mcmc_trace(fit_biomass, regex_pars = "beta\\[",
                               facet_args = list(ncol = 2)) +
  labs(title = "M1 — Biomass: trace plots (beta)") + theme_EL()
p_trace_inflo   <- mcmc_trace(fit_inflo, regex_pars = "beta\\[",
                               facet_args = list(ncol = 2)) +
  labs(title = "M2 — Inflo: trace plots (beta)") + theme_EL()

print(p_trace_biomass)
print(p_trace_inflo)

# Posterior predictive checks
ppc_check <- function(fit, y_obs, label) {
  y_rep <- as.matrix(fit, pars = "y_rep")
  y_rep_sub <- y_rep[sample(nrow(y_rep), 200), ]
  p <- ppc_dens_overlay(y_obs, y_rep_sub) +
    labs(title = paste0("PPC — ", label)) +
    theme_EL()
  print(p)
}

ppc_check(fit_biomass, sd_biomass$y_cont, "Biomass")
ppc_check(fit_inflo,   sd_inflo$y,        "Inflo count")
ppc_check(fit_totspi,  sd_totspi$y,       "Total spikelets")
ppc_check(fit_avgspi,  sd_avgspi$y,       "Avg spikelets")

message("\nPart 2 complete — check diagnostics before Part 3.\n")


# =============================================================================
# PART 3 — PUBLICATION-QUALITY FIGURES (Ecology Letters standard)
# =============================================================================
# Strategy:
#   • Extract posterior draws for each Site × Endo combination
#   • Plot median + 50% and 95% credible intervals
#   • Use ggdist::stat_pointinterval for clean interval visualisation
#   • One multi-panel figure per response (panel A = distribution per site,
#     panel B = Endo contrast per site showing the interaction of interest)
#   • Combine with patchwork; tag panels A/B automatically
#   • Save at 600 dpi PDF (EL requirement)

# ── 3.1  Helper: extract posterior fitted means per Site × Endo ───────────────
# Returns a long data frame of posterior draws of the linear predictor
# (back-transformed to response scale) for each Site × Endo cell.

extract_posterior_means <- function(fit, stan_data_obj, link = "identity") {
  draws_beta <- as.matrix(fit, pars = "beta")           # iter × K
  draws_upop <- as.matrix(fit, pars = "u_pop")          # iter × n_pop

  d   <- stan_data_obj$d_used
  X   <- stan_data_obj$X
  pop <- stan_data_obj$pop_id
  n_iter <- nrow(draws_beta)

  # Compute posterior eta (linear predictor) for every observed row
  # then summarise at Site × Endo level by averaging over individuals
  grid <- d %>%
    dplyr::select(Site, Endo, Pop) %>%
    mutate(pop_int = as.integer(Pop)) %>%
    distinct()

  result <- map_dfr(seq_len(nrow(grid)), function(i) {
    rows <- which(d$Site == grid$Site[i] & d$Endo == grid$Endo[i])
    if (length(rows) == 0) return(NULL)
    # Average eta over individuals in this cell
    eta_mat <- matrix(NA_real_, nrow = n_iter, ncol = length(rows))
    for (j in seq_along(rows)) {
      r <- rows[j]
      eta_mat[, j] <- draws_beta %*% X[r, ] +
        draws_upop[, pop[r]]
    }
    eta_mean <- rowMeans(eta_mat)   # mean over individuals per draw
    mu <- switch(link,
      log      = exp(eta_mean),
      identity = eta_mean
    )
    tibble(
      Site  = grid$Site[i],
      Endo  = grid$Endo[i],
      draw  = seq_len(n_iter),
      mu    = mu
    )
  })
  result
}

message("Extracting posterior draws — this may take a minute …")
post_biomass <- extract_posterior_means(fit_biomass, sd_biomass, "identity")
post_inflo   <- extract_posterior_means(fit_inflo,   sd_inflo,   "log")
post_totspi  <- extract_posterior_means(fit_totspi,  sd_totspi,  "log")
post_avgspi  <- extract_posterior_means(fit_avgspi,  sd_avgspi,  "log")
message("Done.")

# ── 3.2  Helper: posterior contrast E+ minus E- per site ─────────────────────
endo_contrast <- function(post_df) {
  post_df %>%
    dplyr::group_by(Site, draw, Endo) %>%
    summarise(mu = mean(mu), .groups = "drop") %>%  # collapse duplicates
    pivot_wider(names_from = Endo, values_from = mu) %>%
    mutate(contrast = `E+` - `E-`)   # numeric subtraction now works
}

cont_biomass <- endo_contrast(post_biomass)
cont_inflo   <- endo_contrast(post_inflo)
cont_totspi  <- endo_contrast(post_totspi)
cont_avgspi  <- endo_contrast(post_avgspi)

# ── 3.3  Figure builder ───────────────────────────────────────────────────────
make_EL_figure <- function(post_df, raw_df, response_col,
                            y_label, contrast_df, tag_prefix = "A") {

  site_lev <- site_order[site_order %in% levels(factor(post_df$Site))]
  post_df  <- post_df %>% mutate(Site = factor(Site, levels = site_lev),
                                  Endo = factor(Endo, levels = endo_levels))
  raw_sub  <- raw_df  %>% filter(!is.na(.data[[response_col]])) %>%
    mutate(Site = factor(Site, levels = site_lev),
           Endo = factor(Endo, levels = endo_levels))
  contrast_df <- contrast_df %>%
    mutate(Site = factor(Site, levels = site_lev))

  # ── Panel A: raw jitter + posterior median & 95 % CrI per Site × Endo ──────
  pA <- ggplot() +
    geom_jitter(data  = raw_sub,
                aes(x = Site, y = .data[[response_col]], colour = Endo),
                width = 0.18, height = 0, alpha = 0.22, size = 0.7,
                shape = 16) +
    stat_pointinterval(
      data = post_df,
      aes(x = Site, y = mu, colour = Endo, fill = Endo),
      .width        = c(0.50, 0.95),
      position      = position_dodge(0.55),
      point_size    = 2.0,
      interval_size_range = c(0.4, 1.0),
      shape         = 21,
      point_colour  = "black",
      point_alpha   = 1
    ) +
    scale_colour_manual(values = endo_col, name = "Endophyte") +
    scale_fill_manual(values = endo_col) +  # keep fill mapping
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 21,
          fill = endo_col,
          colour = "black"  # match point borders
        )
      ),
      fill = "none"  # hide the extra fill legend
    ) +
    scale_x_discrete(name = "Soil origin (site)") +
    scale_y_continuous(name = y_label,
                       expand = expansion(mult = c(0.02, 0.08))) +
    theme_EL() +
    theme(legend.position = "top") +
    labs(tag = tag_prefix)

  # ── Panel B: E+ − E- contrast per site ──────────────────────────────────────
  # zero line = no difference; positive = E+ higher
  pB <- ggplot(contrast_df,
               aes(x = Site, y = contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.35, colour = "grey50") +
    stat_pointinterval(
      .width = c(0.50, 0.95),
      point_size    = 2.0,
      interval_size_range = c(0.4, 1.0),
      shape = 21,
      colour = "grey20",
      fill   = "grey60",
      point_colour = "black"
    ) +
    scale_x_discrete(name = "Soil origin (site)") +
    scale_y_continuous(
      name   = paste0("\u0394 ", y_label, "\n(E+ \u2212 E\u2212)"),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    theme_EL() +
    labs(tag = paste0(
      # next letter after tag_prefix
      intToUtf8(utf8ToInt(tag_prefix) + 1L)
    ))

  # ── Combine with patchwork ──────────────────────────────────────────────────
  pA / pB + plot_layout(heights = c(2, 1))
}

# ── 3.4  Generate and save figures ───────────────────────────────────────────
message("Building publication figures …")

fig_biomass <- make_EL_figure(
  post_biomass, aghysoils, "abg_mass_tot",
  "Aboveground biomass (g)", cont_biomass, "A"
)

fig_inflo <- make_EL_figure(
  post_inflo, aghysoils, "total_inflo",
  "Total inflorescences", cont_inflo, "A"
)

fig_totspi <- make_EL_figure(
  post_totspi, aghysoils, "tot_spikelet",
  "Total spikelets", cont_totspi, "A"
)

fig_avgspi <- make_EL_figure(
  post_avgspi, aghysoils, "avg_spikelet",
  "Avg spikelets per inflorescence", cont_avgspi, "A"
)

save_EL(fig_biomass, "Fig1_Biomass.pdf",     174, 150)
save_EL(fig_inflo,   "Fig2_Inflo.pdf",       174, 150)
save_EL(fig_totspi,  "Fig3_TotalSpikelets.pdf", 174, 150)
save_EL(fig_avgspi,  "Fig4_AvgSpikelets.pdf",   174, 150)

# ── 3.5  Optional: combined 4-panel figure (for supplementary or overview) ───
# Re-extract panel A only from each, relabel, combine 2 × 2
extract_panelA <- function(fig) fig[[1]]   # patchwork: first sub-plot

fig_all <- (extract_panelA(fig_biomass) + labs(tag = "A",
              title = "Aboveground biomass")) +
           (extract_panelA(fig_inflo)   + labs(tag = "B",
              title = "Inflorescences")) +
           (extract_panelA(fig_totspi)  + labs(tag = "C",
              title = "Total spikelets")) +
           (extract_panelA(fig_avgspi)  + labs(tag = "D",
              title = "Avg spikelets/inflo")) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")


library(patchwork)
library(ggplot2)

extract_panelA <- function(fig) fig[[1]]  # extract first sub-plot

fig_all <- (extract_panelA(fig_biomass) + labs(tag = "A", title = "Aboveground biomass")) +
  (extract_panelA(fig_inflo)   + labs(tag = "B", title = "Inflorescences")) +
  (extract_panelA(fig_totspi)  + labs(tag = "C", title = "Total spikelets")) +
  guide_area() +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position        = "bottom",
    legend.title.position  = "top",          # <-- title on its own line above keys
    legend.title           = element_text(face = "bold", hjust = 0.5),  # centred
    legend.text            = element_text(size = 10),
    legend.box.margin      = margin(t = 10, r = 0, b = 0, l = 0)
  )
save_EL(fig_all, "Fig_AllResponses_combined.pdf", 174, 200)

message("\nPart 3 complete — figures saved to: ", fig_dir, "\n")


# =============================================================================
# PART 4 — STATISTICAL SUMMARY TABLES
# =============================================================================
# Posterior summaries for all beta coefficients across models.
# Format ready for LaTeX (xtable) or Word (write.csv).

# ── 4.1  Helper: clean posterior summary for one model ───────────────────────
posterior_table <- function(fit, stan_data_obj, model_label) {
  d <- stan_data_obj$d_used

  # coefficient names from model matrix
  coef_names <- colnames(stan_data_obj$X)
  coef_names <- gsub("\\(Intercept\\)", "Intercept", coef_names)
  coef_names <- gsub("EndoE\\+",        "Endo[E+]",  coef_names)
  coef_names <- gsub("Site",            "Site:",      coef_names)
  coef_names <- gsub("EndoE\\+:Site",   "Endo[E+] × Site:", coef_names)
  coef_names <- gsub("age_std",         "Age (standardised)", coef_names)

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
      `P(>0)` = round(mean(draw > 0), 3),   # posterior prob > 0
      .groups = "drop"
    ) %>%
    mutate(Model = model_label, .before = 1)
}

tab_biomass <- posterior_table(fit_biomass, sd_biomass, "Biomass (Gaussian)")
tab_inflo   <- posterior_table(fit_inflo,   sd_inflo,   "Inflo count (NB)")
tab_totspi  <- posterior_table(fit_totspi,  sd_totspi,  "Total spikelets (NB)")
tab_avgspi  <- posterior_table(fit_avgspi,  sd_avgspi,  "Avg spikelets (NB)")

all_tabs <- bind_rows(tab_biomass, tab_inflo, tab_totspi, tab_avgspi)

# ── 4.2  Print to console ─────────────────────────────────────────────────────
cat("\n\n")
cat("=============================================================\n")
cat(" BAYESIAN POSTERIOR SUMMARY — All models\n")
cat(" Median, Mean, SD, 95% CrI, P(effect > 0)\n")
cat("=============================================================\n\n")
for (mod_label in unique(all_tabs$Model)) {
  cat("--- ", mod_label, " ---\n", sep = "")
  all_tabs %>%
    filter(Model == mod_label) %>%
    dplyr::select(-Model) %>%
    print(n = Inf)
  cat("\n")
}

# ── 4.3  Save as CSV (import into Word / Overleaf) ───────────────────────────
write.csv(all_tabs,
          file      = file.path(output_dir, "Table1_PosteriorSummary.csv"),
          row.names = FALSE)
message("Posterior summary table saved.")

# ── 4.4  Posterior probability of interaction (Site × Endo) ──────────────────
# Key ecological question: does soil modulate the endophyte effect?
# Extract all Endo:Site interaction betas and report P(beta ≠ direction)

interaction_summary <- function(fit, stan_data_obj, model_label) {
  coef_names <- colnames(stan_data_obj$X)
  interact_idx <- grep("EndoE\\+:Site", coef_names)

  if (length(interact_idx) == 0) {
    message(model_label, ": no interaction terms found.")
    return(NULL)
  }

  draws_beta <- as.matrix(fit, pars = "beta")

  map_dfr(interact_idx, function(k) {
    d <- draws_beta[, k]
    tibble(
      Model       = model_label,
      Term        = coef_names[k],
      Median      = round(median(d), 3),
      `95% CrI`   = sprintf("[%.3f, %.3f]",
                             quantile(d, 0.025), quantile(d, 0.975)),
      `P(beta>0)` = round(mean(d > 0), 3),
      `P(beta<0)` = round(mean(d < 0), 3),
      Evidence    = case_when(
        mean(d > 0) > 0.95 | mean(d < 0) > 0.95 ~ "Strong",
        mean(d > 0) > 0.80 | mean(d < 0) > 0.80 ~ "Moderate",
        TRUE                                      ~ "Weak"
      )
    )
  })
}

int_tabs <- bind_rows(
  interaction_summary(fit_biomass, sd_biomass, "Biomass (Gaussian)"),
  interaction_summary(fit_inflo,   sd_inflo,   "Inflo count (NB)"),
  interaction_summary(fit_totspi,  sd_totspi,  "Total spikelets (NB)"),
  interaction_summary(fit_avgspi,  sd_avgspi,  "Avg spikelets (NB)")
)

cat("\n=============================================================\n")
cat(" SITE × ENDOPHYTE INTERACTION SUMMARY\n")
cat(" (Core ecological question)\n")
cat("=============================================================\n\n")
print(int_tabs, n = Inf)

write.csv(int_tabs,
          file      = file.path(output_dir, "Table2_Interactions.csv"),
          row.names = FALSE)
message("Interaction table saved.")

message("\n\n=== ALL PARTS COMPLETE ===\n")
message("Figures: ", fig_dir)
message("Tables:  ", output_dir)

# =============================================================================
# End of script
# =============================================================================
