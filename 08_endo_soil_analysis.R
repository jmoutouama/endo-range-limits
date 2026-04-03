# =============================================================================
# DOES SOIL INTERACT WITH ENDOPHYTE STATUS TO INFLUENCE DEMOGRAPHY?
# Full Bayesian Analysis via rstan
# Author: Jacob Moutouama
# =============================================================================
rm(list = ls())

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
endo_col    <- c("E-" = "tomato", "E+" = "cornflowerblue")
endo_levels <- c("E-", "E+")
site_order  <- c("SON","KER","BFL","BAS","COL","HUN","LAF")

# ── ggplot theme for Ecology Letters ─────────────────────────────────────────
theme_EL <- function(base_size = 9, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(size = base_size, face = "plain"),
      axis.text         = element_text(size = base_size - 1, colour = "black"),
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      legend.position   = "right",
      legend.key.size   = unit(8, "pt"),
      legend.text       = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size - 1, face = "plain"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "plain"),
      plot.tag          = element_text(size = base_size + 1, face = "plain"),
      plot.title        = element_text(size = base_size, face = "plain", hjust = 0),
      plot.margin       = margin(4, 4, 4, 4, "pt")
    )
}

# ── Output directories ────────────────────────────────────────────────────────
output_dir <- "/Users/jacobmoutouama/Dropbox/Miller Lab/range limits model output/"
fig_dir    <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure"

save_EL <- function(p, filename, width_mm = 174, height_mm = 110) {
  ggsave(file.path(fig_dir, filename), plot = p,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, device = "pdf")
  message("Saved: ", filename)
}


# =============================================================================
# PART 1 — DATA PREPARATION & EXPLORATORY DIAGNOSTICS
# =============================================================================

# ── 1.1  Column definitions (from metadata) ───────────────────────────────────
# spike_inflo_a..j   : spikelets counted from each individually measured inflo
# spike_inflos_unk   : spikelets from inflos that cannot be traced individually
# tot_spikelet       : Excel sum of spikelets from counted inflos (QC only)
# avg_spikelet       : Excel mean spikelets among measured inflos (QC only)
# total_inflo        : total inflorescences COLLECTED per plant (field count)
# total_spik_calc    : Excel projected total = avg_spikelet × total_inflo
# Inflo_mass         : seed_save + seed_squash
# abg_mass_tot       : abg_mass_sans_inflo + Inflo_mass
#
# PRIMARY MODEL RESPONSES (recalculated from raw columns):
#   calc_abg_mass_tot  : total aboveground biomass (g)
#   calc_total_inflo   : total inflorescences collected
#   calc_total_spik    : projected total spikelets = calc_avg_spikelet × calc_total_inflo
#   calc_avg_spikelet  : mean spikelets among individually measured inflos (a-j)

spike_indiv <- c("spike_inflo_a", "spike_inflo_b", "spike_inflo_c", "spike_inflo_d",
                 "spike_inflo_e", "spike_inflo_f", "spike_inflo_g", "spike_inflo_h",
                 "spike_inflo_i", "spike_infl_j")
spike_all   <- c(spike_indiv, "spike_inflos_unk")

# ── Birthday parser (handles mixed US date formats & two-digit years) ─────────
parse_birthday <- function(x) {
  x[x == ""] <- NA
  x <- sub("/(\\d{2})$", "/20\\1", x)   # m/d/YY → m/d/20YY before parsing
  as.Date(x, "%m/%d/%Y")
}

# ── Load & clean ──────────────────────────────────────────────────────────────
soils <- read.csv(
  "https://www.dropbox.com/scl/fi/tujihmutv87jufadtvbwj/reproduction_and_biomass-endo-soil.csv?rlkey=pa04j3jfae5cldlbc880w3xr3&dl=1",
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  rename(Endo        = "Endo",
         seed_save   = "seed save",
         seed_squash = "seed squash") %>%
  mutate(
    # ── Force raw measurement columns to numeric ──────────────────────────────
    across(all_of(c(spike_all,
                    "abg_mass_sans_inflo", "seed_save", "seed_squash",
                    "Inflo_mass", "total_inflo", "tot_spikelet",
                    "avg_spikelet", "total_spik_calc", "abg_mass_tot")),
           ~ suppressWarnings(as.numeric(.))),

    # ── Recalculate response variables (order matters: each uses the previous) ─

    # (1) calc_n_measured: inflorescences individually measured (a-j only)
    #     unk cannot contribute to the mean — used as denominator only
    calc_n_measured   = rowSums(!is.na(across(all_of(spike_indiv)))),

    # (2) calc_avg_spikelet: mean spikelets among measured inflos only
    #     per metadata: "mean number of spikelets among inflorescences measured"
    calc_avg_spikelet = rowSums(across(all_of(spike_indiv)), na.rm = TRUE) /
                          calc_n_measured,

    # (3) calc_tot_spikelet: sum across counted inflos only — QC USE ONLY
    #     NOT a model response: underestimates plants where total_inflo > n_Inflo
    calc_tot_spikelet = ifelse(
      rowSums(!is.na(across(all_of(spike_all)))) == 0,
      NA_real_,
      rowSums(across(all_of(spike_all)), na.rm = TRUE)
    ),

    # (4) calc_total_inflo: total inflorescences collected (field count)
    #     per metadata: "total number of inflorescences collected from each plant"
    calc_total_inflo  = total_inflo,

    # (5) calc_total_spik: projected total spikelets = avg × total_inflo
    #     PRIMARY RESPONSE for total reproductive output — accounts for
    #     inflos collected but not individually counted
    #     per metadata: "multiplying average spikelet count and inflorescence count"
    calc_total_spik   = calc_avg_spikelet * calc_total_inflo,

    # (6) calc_Inflo_mass: seed_save + seed_squash
    calc_Inflo_mass   = seed_save + seed_squash,

    # (7) calc_abg_mass_tot: total aboveground mass
    #     log-transform applied inside make_stan_data() for Gaussian model only
    calc_abg_mass_tot = abg_mass_sans_inflo + calc_Inflo_mass,

    # ── Discrepancy flags (QC only — not used in models) ─────────────────────
    discrep_tot_spikelet = !is.na(tot_spikelet) &
                             tot_spikelet != calc_tot_spikelet,
    discrep_avg_spikelet = !is.na(avg_spikelet) &
                             round(avg_spikelet, 4) != round(calc_avg_spikelet, 4),
    discrep_total_inflo  = !is.na(total_inflo) &
                             total_inflo != calc_total_inflo,
    discrep_total_spik   = !is.na(total_spik_calc) &
                             round(total_spik_calc, 2) != round(calc_total_spik, 2),
    discrep_Inflo_mass   = !is.na(Inflo_mass) & !is.na(calc_Inflo_mass) &
                             round(Inflo_mass, 5) != round(calc_Inflo_mass, 5),
    discrep_abg_mass_tot = !is.na(abg_mass_tot) & !is.na(calc_abg_mass_tot) &
                             round(abg_mass_tot, 5) != round(calc_abg_mass_tot, 5),

    # ── Dates, age, factors ───────────────────────────────────────────────────
    Birthday = parse_birthday(Birthday),
    age_days = as.numeric(as.Date("2025-01-01") - Birthday),
    age_std  = as.numeric(scale(age_days)),
    Site     = factor(Site, levels = site_order),
    Endo     = factor(Endo, levels = endo_levels),
    Pop      = as.factor(Pop),
    Tray     = as.factor(Tray)
  )

# ── Sanity check: no implausible birthdays ────────────────────────────────────
bad_dates <- soils %>%
  filter(!is.na(Birthday),
         Birthday < as.Date("2024-01-01") | Birthday > as.Date("2025-12-31")) %>%
  select(Pop, Site, Endo, Birthday)
if (nrow(bad_dates) > 0) warning("Implausible Birthday values:\n", print(bad_dates))

# ── QC report ─────────────────────────────────────────────────────────────────
cat("=== Discrepancy report (Excel vs recalculated from raw columns) ===\n")
cat("  tot_spikelet  :", sum(soils$discrep_tot_spikelet,  na.rm = TRUE), "rows differ\n")
cat("  avg_spikelet  :", sum(soils$discrep_avg_spikelet,  na.rm = TRUE), "rows differ\n")
cat("  total_inflo   :", sum(soils$discrep_total_inflo,   na.rm = TRUE), "rows differ\n")
cat("  total_spik    :", sum(soils$discrep_total_spik,    na.rm = TRUE), "rows differ\n")
cat("  Inflo_mass    :", sum(soils$discrep_Inflo_mass,    na.rm = TRUE), "rows differ\n")
cat("  abg_mass_tot  :", sum(soils$discrep_abg_mass_tot,  na.rm = TRUE), "rows differ\n")
cat("All response variables use the recalculated (calc_*) versions.\n\n")

# NOTE: 39 rows differ from Excel avg_spikelet because:
# (1) Excel incorrectly included spike_inflos_unk in the mean numerator
# (2) Some rows used total_inflo as denominator instead of n measured inflos
# calc_avg_spikelet uses only spike_indiv columns / calc_n_measured — correct per metadata
summary(soils)

# ── Filter to AGHY only ───────────────────────────────────────────────────────
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

# ── 1.3  Distribution of each response ───────────────────────────────────────
# NOTE: calc_total_spik (projected total = avg × total_inflo) is the primary
#       response for total reproductive output — replaces calc_tot_spikelet
resp_long <- aghysoils %>%
  dplyr::select(Site, Endo,
                calc_abg_mass_tot, calc_total_inflo,
                calc_total_spik, calc_avg_spikelet) %>%
  pivot_longer(-c(Site, Endo),
               names_to  = "response",
               values_to = "value") %>%
  dplyr::mutate(response = dplyr::recode(response,
    "calc_abg_mass_tot" = "Aboveground biomass (g)",
    "calc_total_inflo"  = "Total inflorescences",
    "calc_total_spik"   = "Total spikelets (projected)",
    "calc_avg_spikelet" = "Avg spikelets / inflo"))

p_dist <- ggplot(resp_long %>% filter(!is.na(value)),
                 aes(x = value, fill = Endo, colour = Endo)) +
  geom_density(alpha = 0.45, linewidth = 0.3, adjust = 1.2) +
  scale_fill_manual(values  = endo_col, name = "Endophyte") +
  scale_colour_manual(values = endo_col, name = "Endophyte") +
  facet_wrap(~ response, scales = "free", ncol = 2) +
  labs(x = "Observed value", y = "Density") +
  theme_EL() +
  theme(legend.position = "top")

# ── 1.4  Raw means ± SE across sites ─────────────────────────────────────────
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
  labs(x = "Soil origin (site)", y = "Observed value (mean ± SE)") +
  theme_EL() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 35, hjust = 1))

# ── 1.5  Age covariate vs responses ──────────────────────────────────────────
age_summary <- aghysoils %>%
  group_by(Site, Endo) %>%
  summarise(age_days = mean(age_days, na.rm = TRUE), .groups = "drop")

resp_long <- resp_long %>%
  left_join(age_summary, by = c("Site", "Endo"))

p_age <- ggplot(resp_long %>% filter(!is.na(value)),
                aes(x = age_days, y = value, colour = Endo)) +
  geom_point(alpha = 0.35, size = 0.7, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, formula = y ~ x) +
  scale_colour_manual(values = endo_col, name = "Endophyte") +
  facet_wrap(~ response, scales = "free_y", ncol = 2) +
  labs(x = "Age (days from transplant to 1 Jan 2025)",
       y = "Observed value") +
  theme_EL() +
  theme(legend.position = "top")

print(p_n)
print(p_dist)
print(p_raw)
print(p_age)


# =============================================================================
# PART 2 — STAN MODEL BUILDING & SAMPLING
# =============================================================================

# ── 2.1  Design matrices ──────────────────────────────────────────────────────
make_stan_data <- function(df, response_col, family = "negbinom") {
  d <- df %>%
    filter(!is.na(.data[[response_col]]),
           !is.na(Endo), !is.na(Site), !is.na(age_std)) %>%
    droplevels()

  X <- model.matrix(~ Endo * Site + age_std, data = d)

  # Gaussian (biomass): log-transform here so calc_abg_mass_tot stays on gram
  # scale in plots; only the modelling step uses log
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

# ── Response variables per metadata definitions ───────────────────────────────
# calc_abg_mass_tot : total aboveground biomass (g)           → Gaussian
# calc_total_inflo  : total inflorescences collected          → NegBinom
# calc_total_spik   : projected total spikelets (avg×count)   → NegBinom  ← FIX
# calc_avg_spikelet : mean spikelets per measured inflo       → NegBinom
sd_biomass <- make_stan_data(aghysoils, "calc_abg_mass_tot",  "gaussian")
sd_inflo   <- make_stan_data(aghysoils, "calc_total_inflo",   "negbinom")
sd_totspi  <- make_stan_data(aghysoils, "calc_total_spik",    "negbinom")  # ← was calc_tot_spikelet
sd_avgspi  <- make_stan_data(aghysoils, "calc_avg_spikelet",  "negbinom")

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
  y_cont ~ student_t(3, X * beta + u_pop[pop_id], sigma);
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    real mu_n = X[n] * beta + u_pop[pop_id[n]];
    y_rep[n]   = normal_rng(mu_n, sigma);
    log_lik[n] = normal_lpdf(y_cont[n] | mu_n, sigma);
  }
}
"

# ── 2.3  Stan model: Negative-Binomial (count responses) ─────────────────────
stan_negbinom_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=0> y;
  int<lower=0> y_max;
  int<lower=1> n_pop;
  array[N] int<lower=1, upper=n_pop> pop_id;
}
parameters {
  vector[K] beta;
  real log_phi;
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
  vector[N] log_lik;
  for (n in 1:N) {
    real log_mu_n    = X[n] * beta + u_pop[pop_id[n]];
    real safe_log_mu = fmin(log_mu_n, 10.0);   // cap at ~22,000 counts
    real safe_phi    = fmin(phi, 1e6);          // cap phi
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
                           model_name = "gaussian_endo_soil")
mod_negbinom <- stan_model(model_code = stan_negbinom_code,
                           model_name = "negbinom_endo_soil")

# ── 2.5  Sampling helper ──────────────────────────────────────────────────────
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
fit_biomass <- run_stan(mod_gaussian, sd_biomass, "gaussian")
fit_inflo   <- run_stan(mod_negbinom, sd_inflo,   "negbinom")
fit_totspi  <- run_stan(mod_negbinom, sd_totspi,  "negbinom")
fit_avgspi  <- run_stan(mod_negbinom, sd_avgspi,  "negbinom")

# ── 2.7  MCMC diagnostics ────────────────────────────────────────────────────
check_hmc <- function(fit, label) {
  cat("\n===", label, "===\n")
  cat("Rhat summary (should be < 1.01):\n")
  rhat_vals <- summary(fit)$summary[, "Rhat"]
  print(summary(rhat_vals))
  cat(sprintf("Divergent transitions: %d\n",   get_num_divergent(fit)))
  cat(sprintf("Max treedepth exceeded: %d\n", get_num_max_treedepth(fit)))
}

check_hmc(fit_biomass, "M1 Biomass")
check_hmc(fit_inflo,   "M2 Inflo")
check_hmc(fit_totspi,  "M3 Total spikelets (projected)")
check_hmc(fit_avgspi,  "M4 Avg spikelets")

p_trace_biomass <- mcmc_trace(fit_biomass, regex_pars = "beta\\[",
                               facet_args = list(ncol = 2)) +
  labs(title = "M1 — Biomass: trace plots (beta)") + theme_EL()
p_trace_inflo   <- mcmc_trace(fit_inflo, regex_pars = "beta\\[",
                               facet_args = list(ncol = 2)) +
  labs(title = "M2 — Inflo: trace plots (beta)") + theme_EL()

print(p_trace_biomass)
print(p_trace_inflo)

ppc_check <- function(fit, y_obs, label, xlim = NULL) {
  y_rep     <- as.matrix(fit, pars = "y_rep")
  y_rep_sub <- y_rep[sample(nrow(y_rep), 200), ]
  p <- ppc_dens_overlay(y_obs, y_rep_sub) +
    labs(title = paste0("PPC — ", label)) +
    theme_EL()
  if(!is.null(xlim)){
    p <- p + coord_cartesian(xlim = xlim) +
      scale_x_continuous(expand = c(0, 0))
  }
  return(p)
}

p_biomass <- ppc_check(fit_biomass, sd_biomass$y_cont, "Biomass")
p_inflo <- ppc_check(
  fit_inflo,
  sd_inflo$y,
  "Inflo count",
  xlim = c(-1, 50)   # adjust as needed
)
p_inflo <- p_inflo + scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0))

p_avgspi <- ppc_check(
  fit_avgspi,
  sd_avgspi$y,
  "Avg spikelets"
)
ppc_panel <- (p_biomass | p_inflo | p_avgspi) +
  plot_annotation(tag_levels = "A")
# ggsave(
#   filename = file.path(fig_dir, "PPC_models_soilendo.pdf"),
#   plot = ppc_panel,
#   width = 12,   # adjust width
#   height = 4,   # adjust height
#   device = cairo_pdf   # high-quality PDF
# )
# =============================================================================
# PART 3 FIGURES
# =============================================================================

# ── 3.1  Extract posterior fitted means per Site × Endo ──────────────────────
extract_posterior_means <- function(fit, sd_obj, link = "identity") {
  draws_beta <- as.matrix(fit, pars = "beta")    # iter × K
  draws_upop <- as.matrix(fit, pars = "u_pop")   # iter × n_pop

  d      <- sd_obj$d_used
  X      <- sd_obj$X
  pop    <- sd_obj$pop_id
  n_iter <- nrow(draws_beta)

  grid <- d %>%
    dplyr::select(Site, Endo, Pop) %>%
    mutate(pop_int = as.integer(Pop)) %>%
    distinct()

  result <- map_dfr(seq_len(nrow(grid)), function(i) {
    rows <- which(d$Site == grid$Site[i] & d$Endo == grid$Endo[i])
    if (length(rows) == 0) return(NULL)
    eta_mat <- matrix(NA_real_, nrow = n_iter, ncol = length(rows))
    for (j in seq_along(rows)) {
      r <- rows[j]
      eta_mat[, j] <- draws_beta %*% X[r, ] + draws_upop[, pop[r]]
    }
    eta_mean <- rowMeans(eta_mat)
    mu <- switch(link,
      log      = exp(eta_mean),
      identity = eta_mean
    )
    tibble(Site = grid$Site[i], Endo = grid$Endo[i],
           draw = seq_len(n_iter), mu = mu)
  })
  result
}

post_biomass <- extract_posterior_means(fit_biomass, sd_biomass, "identity")
post_inflo   <- extract_posterior_means(fit_inflo,   sd_inflo,   "log")
post_totspi  <- extract_posterior_means(fit_totspi,  sd_totspi,  "log")
post_avgspi  <- extract_posterior_means(fit_avgspi,  sd_avgspi,  "log")

# ── 3.2  Posterior contrast E+ minus E- per site ─────────────────────────────
endo_contrast <- function(post_df) {
  post_df %>%
    dplyr::group_by(Site, draw, Endo) %>%
    summarise(mu = mean(mu), .groups = "drop") %>%
    pivot_wider(names_from = Endo, values_from = mu) %>%
    mutate(contrast = `E+` - `E-`)
}

cont_biomass <- endo_contrast(post_biomass)
cont_inflo   <- endo_contrast(post_inflo)
cont_totspi  <- endo_contrast(post_totspi)
cont_avgspi  <- endo_contrast(post_avgspi)

# ── 3.3  Figure builder ───────────────────────────────────────────────────────
make_EL_figure <- function(post_df, raw_df, response_col,
                           y_label, contrast_df, tag_prefix = "A") {
  site_lev    <- site_order[site_order %in% levels(factor(post_df$Site))]
  post_df     <- post_df %>% mutate(Site = factor(Site, levels = site_lev),
                                    Endo = factor(Endo, levels = endo_levels))
  raw_sub     <- raw_df  %>% filter(!is.na(.data[[response_col]]),
                                    .data[[response_col]] > 0) %>%   # drop zeros
    mutate(Site = factor(Site, levels = site_lev),
           Endo = factor(Endo, levels = endo_levels))
  contrast_df <- contrast_df %>% mutate(Site = factor(Site, levels = site_lev))
  
  pA <- ggplot() +
    geom_jitter(data = raw_sub,
                aes(x = Site, y = .data[[response_col]], colour = Endo),
                position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width  = 0.7),
                alpha = 0.55, size = 1, shape = 16) +
    stat_pointinterval(
      data          = post_df,
      aes(x = Site, y = mu, colour = Endo, fill = Endo),
      .width        = c(0.50, 0.95),
      position      = position_dodge(0.7),
      point_size    = 2.5,
      interval_size = 0.8,
      shape         = 21,
      point_alpha   = 1
    ) +
    scale_colour_manual(values = endo_col, name = "Endophyte") +
    scale_fill_manual(values   = endo_col) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape     = 21,
          fill      = endo_col,
          colour    = endo_col,
          size      = 3,
          linewidth = 0.8,
          alpha     = 1
        )
      ),
      fill = "none"
    ) +
    scale_x_discrete(name = "Soil origin (site)") +
    scale_y_continuous(name   = y_label,
                       expand = expansion(mult = c(0.02, 0.08))) +
    theme_EL() +
    theme(legend.position = "top") +
    labs(tag = tag_prefix)
  
  pB <- ggplot(contrast_df, aes(x = Site, y = contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.35, colour = "grey50") +
    stat_pointinterval(
      .width        = c(0.50, 0.95),
      point_size    = 2.5,
      interval_size = 0.8,
      shape         = 21,
      colour        = "grey20",
      fill          = "grey60",
      point_colour  = "black"
    ) +
    scale_x_discrete(name = "Soil origin (site)") +
    scale_y_continuous(
      name   = paste0("\u0394 ", y_label, "\n(E+ \u2212 E\u2212)"),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    theme_EL() +
    labs(tag = intToUtf8(utf8ToInt(tag_prefix) + 1L))
  
  pA / pB + plot_layout(heights = c(2, 1))
}

# Re-extract with back-transformation
post_biomass <- extract_posterior_means(fit_biomass, sd_biomass, "identity") %>%
  mutate(mu = exp(mu))   # grams scale

# Rebuild figures
fig_biomass <- make_EL_figure(post_biomass, aghysoils, "calc_abg_mass_tot",
                              "Aboveground biomass (g)",         cont_biomass, "A")
fig_inflo   <- make_EL_figure(post_inflo,   aghysoils, "calc_total_inflo",
                              "Total inflorescences",            cont_inflo,   "A")
fig_totspi  <- make_EL_figure(post_totspi,  aghysoils, "calc_total_spik",
                              "Total spikelets (projected)",     cont_totspi,  "A")
fig_avgspi  <- make_EL_figure(post_avgspi,  aghysoils, "calc_avg_spikelet",
                              "Avg spikelets per inflorescence", cont_avgspi,  "A")
# save_EL(fig_biomass, "Fig1_Biomass.pdf",              174, 150)
# save_EL(fig_inflo,   "Fig2_Inflo.pdf",                174, 150)
# save_EL(fig_totspi,  "Fig3_TotalSpikelets.pdf",       174, 150)
# save_EL(fig_avgspi,  "Fig4_AvgSpikelets.pdf",         174, 150)

# ── 3.5  Combined 3-panel figure + legend in D slot ──────────────────────────
extract_panelA <- function(fig) fig[[1]]

fig_all <- ((extract_panelA(fig_biomass) + labs(tag = "A")) +
              (extract_panelA(fig_inflo)   + labs(tag = "B")) +
              (extract_panelA(fig_avgspi)  + labs(tag = "C")) +
              guide_area() +
              plot_layout(ncol = 2, guides = "collect") &
              theme(
                legend.position       = "inside",
                legend.justification  = c(0.5, 0.5),
                legend.direction      = "vertical",
                legend.title.position = "top",
                legend.title          = element_text(face = "plain", hjust = 0.5),
                legend.text           = element_text(size = 10),
                legend.key.size       = unit(1.2, "cm"),
                legend.spacing.x      = unit(0.1, "cm"),   
                legend.spacing.y      = unit(0.2, "cm"),   
                legend.box.margin     = margin(t = 10, r = 0, b = 0, l = 0)
              ))

save_EL(fig_all, "Fig_AllResponses_combined.pdf", 174, 120)

message("\nPart 3 complete — figures saved to: ", fig_dir, "\n")


# =============================================================================
# PART 4 — STATISTICAL SUMMARY TABLES
# =============================================================================

# ── 4.1  Posterior summary for one model ─────────────────────────────────────
posterior_table <- function(fit, sd_obj, model_label) {
  coef_names <- colnames(sd_obj$X)
  coef_names <- gsub("\\(Intercept\\)", "Intercept",          coef_names)
  coef_names <- gsub("EndoE\\+",        "Endo[E+]",           coef_names)
  coef_names <- gsub("Site",            "Site:",              coef_names)
  coef_names <- gsub("EndoE\\+:Site",   "Endo[E+] × Site:",  coef_names)
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
      `P(>0)` = round(mean(draw > 0), 3),
      .groups = "drop"
    ) %>%
    mutate(Model = model_label, .before = 1)
}

tab_biomass <- posterior_table(fit_biomass, sd_biomass, "Biomass (Gaussian)")
tab_inflo   <- posterior_table(fit_inflo,   sd_inflo,   "Inflo count (NB)")
tab_totspi  <- posterior_table(fit_totspi,  sd_totspi,  "Total spikelets projected (NB)")  # ← label updated
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
  all_tabs %>% filter(Model == mod_label) %>% dplyr::select(-Model) %>% print(n = Inf)
  cat("\n")
}

# ── 4.3  Save as CSV ──────────────────────────────────────────────────────────
write.csv(all_tabs,
          file      = file.path(output_dir, "Table1_PosteriorSummary.csv"),
          row.names = FALSE)
message("Posterior summary table saved.")

# ── 4.4  Interaction summary ──────────────────────────────────────────────────
interaction_summary <- function(fit, sd_obj, model_label) {
  coef_names   <- colnames(sd_obj$X)
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
      `95% CrI`   = sprintf("[%.3f, %.3f]", quantile(d, 0.025), quantile(d, 0.975)),
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
  interaction_summary(fit_totspi,  sd_totspi,  "Total spikelets projected (NB)"),
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
