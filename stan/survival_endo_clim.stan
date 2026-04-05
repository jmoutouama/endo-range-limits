// Hierarchical Bayesian survival model with partial pooling across species
//
// Response:
//   Survival outcome modeled using a Bernoulli distribution with logit link.
//
// Fixed effects:
//   Species-specific coefficients are estimated for:
//     - Climate (precipitation; scaled)
//     - Endophyte presence
//     - Endophyte × Climate interaction
//
//   Species-level coefficients are drawn from global hyperpriors,
//   allowing partial pooling and information sharing among species.
//
// Random effects:
//   - Site-year effects (species-specific deviations with shared variance)
//   - Plot effects (shared variance)
//   - Source/population effects (shared variance)
//
// This specification tests whether endophyte effects on survival
// depend on climatic conditions.

data {
  int<lower=1> N;
  int<lower=1> nSpp;
  int<lower=1> nsite_year;
  int<lower=1> nPop;
  int<lower=1> nPlot;

  array[N] int<lower=1, upper=nSpp>       Spp;
  array[N] int<lower=1, upper=nsite_year> site_year;
  array[N] int<lower=1, upper=nPop>       pop;
  array[N] int<lower=1, upper=nPlot>      plot;

  array[N] int<lower=0, upper=1> y;       // survival outcome
  array[N] int<lower=0, upper=1> endo;    // endophyte presence
  vector[N] clim;                         // precipitation (scaled)
}

parameters {

  // Global means for species-level fixed effects
  real mu_b0;
  real mu_bclim;
  real mu_bendo;
  real mu_bendoclim;

  // Standard deviations describing variation among species
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bclim;
  real<lower=0> sigma_bendo;
  real<lower=0> sigma_bendoclim;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bclim;
  vector[nSpp] z_bendo;
  vector[nSpp] z_bendoclim;

  // Random effect standard deviations
  real<lower=0> sigma_site_year;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_pop;

  // Non-centered random effects
  matrix[nSpp, nsite_year] z_site_year;
  vector[nPlot] z_plot;
  vector[nPop] z_pop;
}

transformed parameters {

  // Species-specific regression coefficients via partial pooling
  vector[nSpp] b0        = mu_b0        + sigma_b0        * z_b0;
  vector[nSpp] bclim     = mu_bclim     + sigma_bclim     * z_bclim;
  vector[nSpp] bendo     = mu_bendo     + sigma_bendo     * z_bendo;
  vector[nSpp] bendoclim = mu_bendoclim + sigma_bendoclim * z_bendoclim;

  // Random effects (non-centered parameterization)
  matrix[nSpp, nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop]  pop_rfx  = sigma_pop  * z_pop;

  // Linear predictor
  vector[N] predS;

  for (i in 1:N) {
    predS[i] =
      b0[Spp[i]]
      + bclim[Spp[i]]     * clim[i]
      + bendo[Spp[i]]     * endo[i]
      + bendoclim[Spp[i]] * endo[i] * clim[i]
      + site_year_rfx[Spp[i], site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {

  // Priors for global regression coefficients
  mu_b0        ~ normal(0, 2);
  mu_bclim     ~ normal(0, 2);
  mu_bendo     ~ normal(0, 2);
  mu_bendoclim ~ normal(0, 2);

  // Priors for species-level variation
  sigma_b0        ~ normal(0, 1);
  sigma_bclim     ~ normal(0, 1);
  sigma_bendo     ~ normal(0, 1);
  sigma_bendoclim ~ normal(0, 1);

  // Non-centered species deviations
  z_b0        ~ normal(0, 1);
  z_bclim     ~ normal(0, 1);
  z_bendo     ~ normal(0, 1);
  z_bendoclim ~ normal(0, 1);

  // Priors for random-effect standard deviations
  sigma_site_year ~ normal(0, 1);
  sigma_plot      ~ normal(0, 1);
  sigma_pop       ~ normal(0, 1);

  // Non-centered random effects
  to_vector(z_site_year) ~ normal(0, 1);
  z_plot                 ~ normal(0, 1);
  z_pop                  ~ normal(0, 1);

  // Likelihood
  y ~ bernoulli_logit(predS);
}

generated quantities {

  // Pointwise log-likelihood for model comparison (LOO / WAIC)
  vector[N] log_lik;

  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);
}
