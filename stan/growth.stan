// Hierarchical Bayesian growth model with partial pooling across species
// Response: continuous growth (size t -> t+1) modeled with a Normal distribution
//
// Fixed predictors:
//   - Climate (scaled)
//   - Quadratic climate term (climate²)
//   - Endophyte presence
//   - Herbivory
//   - All two-way and three-way interactions
//
// Species-level coefficients are drawn from global hyperpriors,
// allowing information sharing among species.
//
// Random effects:
//   - Site-year effects
//   - Plot effects
//   - Source population effects

data {
  int<lower=1> N;
  int<lower=1> nSpp;
  int<lower=1> nsite_year;
  int<lower=1> nPop;
  int<lower=1> nPlot;

  array[N] int<lower=1, upper=nSpp> Spp;
  array[N] int<lower=1, upper=nsite_year> site_year;
  array[N] int<lower=1, upper=nPop> pop;
  array[N] int<lower=1, upper=nPlot> plot;

  vector[N] y;

  array[N] int<lower=0,upper=1> endo;
  array[N] int<lower=0,upper=1> herb;

  vector[N] clim;
  vector[N] clim2;
}

parameters {

  // Global means
  real mu_b0;
  real mu_bclim;
  real mu_bclim2;
  real mu_bendo;
  real mu_bherb;
  real mu_bendoclim;
  real mu_bherbclim;
  real mu_bendoherb;
  real mu_bendoherbclim;

  // Species-level variation
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bclim;
  real<lower=0> sigma_bclim2;
  real<lower=0> sigma_bendo;
  real<lower=0> sigma_bherb;
  real<lower=0> sigma_bendoclim;
  real<lower=0> sigma_bherbclim;
  real<lower=0> sigma_bendoherb;
  real<lower=0> sigma_bendoherbclim;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bclim;
  vector[nSpp] z_bclim2;
  vector[nSpp] z_bendo;
  vector[nSpp] z_bherb;
  vector[nSpp] z_bendoclim;
  vector[nSpp] z_bherbclim;
  vector[nSpp] z_bendoherb;
  vector[nSpp] z_bendoherbclim;

  // Random effect SDs
  real<lower=0> sigma_site_year;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_pop;

  // Non-centered random effects
  matrix[nSpp, nsite_year] z_site_year;
  vector[nPlot] z_plot;
  vector[nPop] z_pop;

  // Residual SD
  real<lower=0> sigma;
}

transformed parameters {

  vector[nSpp] b0            = mu_b0            + sigma_b0            * z_b0;
  vector[nSpp] bclim         = mu_bclim         + sigma_bclim         * z_bclim;
  vector[nSpp] bclim2        = mu_bclim2        + sigma_bclim2        * z_bclim2;
  vector[nSpp] bendo         = mu_bendo         + sigma_bendo         * z_bendo;
  vector[nSpp] bherb         = mu_bherb         + sigma_bherb         * z_bherb;
  vector[nSpp] bendoclim     = mu_bendoclim     + sigma_bendoclim     * z_bendoclim;
  vector[nSpp] bherbclim     = mu_bherbclim     + sigma_bherbclim     * z_bherbclim;
  vector[nSpp] bendoherb     = mu_bendoherb     + sigma_bendoherb     * z_bendoherb;
  vector[nSpp] bendoherbclim = mu_bendoherbclim + sigma_bendoherbclim * z_bendoherbclim;

  matrix[nSpp,nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop] pop_rfx = sigma_pop * z_pop;

  vector[N] predG;

  for (i in 1:N) {
    predG[i] =
      b0[Spp[i]]
      + bclim[Spp[i]] * clim[i]
      + bclim2[Spp[i]] * clim2[i]
      + bendo[Spp[i]] * endo[i]
      + bherb[Spp[i]] * herb[i]
      + bendoclim[Spp[i]] * endo[i] * clim[i]
      + bherbclim[Spp[i]] * herb[i] * clim[i]
      + bendoherb[Spp[i]] * endo[i] * herb[i]
      + bendoherbclim[Spp[i]] * endo[i] * herb[i] * clim[i]
      + site_year_rfx[Spp[i],site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {

  // Global priors
  mu_b0 ~ normal(0,2);
  mu_bclim ~ normal(0,2);
  mu_bclim2 ~ normal(0,2);
  mu_bendo ~ normal(0,2);
  mu_bherb ~ normal(0,2);
  mu_bendoclim ~ normal(0,2);
  mu_bherbclim ~ normal(0,2);
  mu_bendoherb ~ normal(0,2);
  mu_bendoherbclim ~ normal(0,2);

  // Species variation
  sigma_b0 ~ normal(0,1);
  sigma_bclim ~ normal(0,1);
  sigma_bclim2 ~ normal(0,1);
  sigma_bendo ~ normal(0,1);
  sigma_bherb ~ normal(0,1);
  sigma_bendoclim ~ normal(0,1);
  sigma_bherbclim ~ normal(0,1);
  sigma_bendoherb ~ normal(0,1);
  sigma_bendoherbclim ~ normal(0,1);

  // Non-centered
  z_b0 ~ normal(0,1);
  z_bclim ~ normal(0,1);
  z_bclim2 ~ normal(0,1);
  z_bendo ~ normal(0,1);
  z_bherb ~ normal(0,1);
  z_bendoclim ~ normal(0,1);
  z_bherbclim ~ normal(0,1);
  z_bendoherb ~ normal(0,1);
  z_bendoherbclim ~ normal(0,1);

  // Random effects
  sigma_site_year ~ normal(0,1);
  sigma_plot ~ normal(0,1);
  sigma_pop ~ normal(0,1);

  to_vector(z_site_year) ~ normal(0,1);
  z_plot ~ normal(0,1);
  z_pop ~ normal(0,1);

  sigma ~ normal(0,1);

  // Likelihood
  y ~ normal(predG, sigma);
}

generated quantities {

  vector[N] log_lik;

  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);
}
