// Hierarchical Bayesian spikelet count model with partial pooling across species
// Response: spike count modeled with Negative Binomial
// Species-specific fixed effects with global hyperpriors
// Random effects for site-year, plot, and source population
//
// Fixed predictors:
//   - Climate (precipitation; scaled)
//   - Endophyte presence
//   - Herbivory
//   - Endophyte × Climate   (symbiont:climate)
//   - Endophyte × Herbivory (symbiont:herbivory)
//   NOTE: Herbivory × Climate and three-way interaction are excluded

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

  array[N] int<lower=0> y;
  array[N] int<lower=0, upper=1> endo;
  array[N] int<lower=0, upper=1> herb;

  vector[N] clim;
}

parameters {

  // Global means for species-level effects
  real mu_b0;
  real mu_bclim;
  real mu_bendo;
  real mu_bherb;
  real mu_bendoclim;    // symbiont:climate
  real mu_bendoherb;    // symbiont:herbivory

  // Species-level SDs
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bclim;
  real<lower=0> sigma_bendo;
  real<lower=0> sigma_bherb;
  real<lower=0> sigma_bendoclim;
  real<lower=0> sigma_bendoherb;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bclim;
  vector[nSpp] z_bendo;
  vector[nSpp] z_bherb;
  vector[nSpp] z_bendoclim;
  vector[nSpp] z_bendoherb;

  // Random-effect SDs
  real<lower=0> sigma_site_year;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_pop;

  // Non-centered random effects
  matrix[nSpp, nsite_year] z_site_year;
  vector[nPlot] z_plot;
  vector[nPop] z_pop;

  // Dispersion parameter for Negative Binomial
  real<lower=0> phi;
}

transformed parameters {

  // Species-specific coefficients
  vector[nSpp] b0        = mu_b0        + sigma_b0        * z_b0;
  vector[nSpp] bclim     = mu_bclim     + sigma_bclim     * z_bclim;
  vector[nSpp] bendo     = mu_bendo     + sigma_bendo     * z_bendo;
  vector[nSpp] bherb     = mu_bherb     + sigma_bherb     * z_bherb;
  vector[nSpp] bendoclim = mu_bendoclim + sigma_bendoclim * z_bendoclim;
  vector[nSpp] bendoherb = mu_bendoherb + sigma_bendoherb * z_bendoherb;

  // Random effects
  matrix[nSpp, nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop]  pop_rfx  = sigma_pop  * z_pop;

  // Linear predictor
  vector[N] predF;
  for (i in 1:N) {
    predF[i] =
      b0[Spp[i]]
      + bclim[Spp[i]]     * clim[i]
      + bendo[Spp[i]]     * endo[i]
      + bherb[Spp[i]]     * herb[i]
      + bendoclim[Spp[i]] * endo[i] * clim[i]   // symbiont:climate
      + bendoherb[Spp[i]] * endo[i] * herb[i]   // symbiont:herbivory
      + site_year_rfx[Spp[i], site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {

  // Priors for global means
  mu_b0        ~ normal(0, 2);
  mu_bclim     ~ normal(0, 2);
  mu_bendo     ~ normal(0, 2);
  mu_bherb     ~ normal(0, 2);
  mu_bendoclim ~ normal(0, 2);
  mu_bendoherb ~ normal(0, 2);

  // Priors for species SDs
  sigma_b0        ~ normal(0, 1);
  sigma_bclim     ~ normal(0, 1);
  sigma_bendo     ~ normal(0, 1);
  sigma_bherb     ~ normal(0, 1);
  sigma_bendoclim ~ normal(0, 1);
  sigma_bendoherb ~ normal(0, 1);

  // Non-centered species deviations
  z_b0        ~ normal(0, 1);
  z_bclim     ~ normal(0, 1);
  z_bendo     ~ normal(0, 1);
  z_bherb     ~ normal(0, 1);
  z_bendoclim ~ normal(0, 1);
  z_bendoherb ~ normal(0, 1);

  // Random-effect SD priors
  sigma_site_year ~ normal(0, 1);
  sigma_plot      ~ normal(0, 1);
  sigma_pop       ~ normal(0, 1);

  to_vector(z_site_year) ~ normal(0, 1);
  z_plot                 ~ normal(0, 1);
  z_pop                  ~ normal(0, 1);

  // Dispersion prior for Negative Binomial
  phi ~ exponential(1);

  // Likelihood
  y ~ neg_binomial_2_log(predF, phi);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | predF[i], phi);
}
