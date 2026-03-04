// Hierarchical Bayesian growth model (endo × clim only)
// Response: continuous growth (size t -> t+1) modeled with Normal
// Fixed effects:
//   - Endophyte presence
//   - Climate
//   - Endophyte × Climate interaction
// Species-level coefficients drawn from global hyperpriors
// Random effects:
//   - Site-year (species-specific)
//   - Plot
//   - Source population

data {
  int<lower=1> N;                // number of observations
  int<lower=1> nSpp;             // number of species
  int<lower=1> nsite_year;       // number of site-year combinations
  int<lower=1> nPop;             // number of source populations
  int<lower=1> nPlot;            // number of plots

  array[N] int<lower=1, upper=nSpp> Spp;
  array[N] int<lower=1, upper=nsite_year> site_year;
  array[N] int<lower=1, upper=nPop> pop;
  array[N] int<lower=1, upper=nPlot> plot;

  vector[N] y;                   // continuous growth response
  array[N] int<lower=0, upper=1> endo;   // endophyte presence
  vector[N] clim;                         // climate covariate
}

parameters {
  // Global means
  real mu_b0;
  real mu_bendo;
  real mu_bclim;
  real mu_bendoclim;  // Endo × Clim

  // Species-level variation
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bendo;
  real<lower=0> sigma_bclim;
  real<lower=0> sigma_bendoclim;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bendo;
  vector[nSpp] z_bclim;
  vector[nSpp] z_bendoclim;

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
  vector[nSpp] b0        = mu_b0        + sigma_b0        * z_b0;
  vector[nSpp] bendo     = mu_bendo     + sigma_bendo     * z_bendo;
  vector[nSpp] bclim     = mu_bclim     + sigma_bclim     * z_bclim;
  vector[nSpp] bendoclim = mu_bendoclim + sigma_bendoclim * z_bendoclim;

  matrix[nSpp, nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop] pop_rfx = sigma_pop * z_pop;

  vector[N] predG;

  for (i in 1:N) {
    predG[i] =
      b0[Spp[i]]
      + bendo[Spp[i]] * endo[i]
      + bclim[Spp[i]] * clim[i]
      + bendoclim[Spp[i]] * endo[i] * clim[i]
      + site_year_rfx[Spp[i], site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {
  // Global priors
  mu_b0 ~ normal(0,2);
  mu_bendo ~ normal(0,2);
  mu_bclim ~ normal(0,2);
  mu_bendoclim ~ normal(0,2);

  // Species variation
  sigma_b0 ~ normal(0,1);
  sigma_bendo ~ normal(0,1);
  sigma_bclim ~ normal(0,1);
  sigma_bendoclim ~ normal(0,1);

  // Non-centered species deviations
  z_b0 ~ normal(0,1);
  z_bendo ~ normal(0,1);
  z_bclim ~ normal(0,1);
  z_bendoclim ~ normal(0,1);

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
