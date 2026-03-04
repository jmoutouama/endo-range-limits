// Hierarchical Bayesian survival model (main effects only)
data {
  int<lower=1> N;                // number of observations
  int<lower=1> nSpp;             // number of species
  int<lower=1> nsite_year;       // number of site-year combinations
  int<lower=1> nPop;             // number of source populations
  int<lower=1> nPlot;            // number of plots

  array[N] int<lower=1, upper=nSpp>       Spp;
  array[N] int<lower=1, upper=nsite_year> site_year;
  array[N] int<lower=1, upper=nPop>       pop;
  array[N] int<lower=1, upper=nPlot>      plot;

  array[N] int<lower=0, upper=1> y;       // survival outcome
  array[N] int<lower=0, upper=1> endo;    // endophyte presence
  array[N] int<lower=0, upper=1> herb;    // herbivory treatment
  vector[N] clim;                         // precipitation (scaled)
}

parameters {

  // Global means for species-level main effects
  real mu_b0;
  real mu_bclim;
  real mu_bendo;
  real mu_bherb;

  // Species-level standard deviations
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bclim;
  real<lower=0> sigma_bendo;
  real<lower=0> sigma_bherb;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bclim;
  vector[nSpp] z_bendo;
  vector[nSpp] z_bherb;

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

  // Species-specific main-effect coefficients
  vector[nSpp] b0    = mu_b0    + sigma_b0    * z_b0;
  vector[nSpp] bclim = mu_bclim + sigma_bclim * z_bclim;
  vector[nSpp] bendo = mu_bendo + sigma_bendo * z_bendo;
  vector[nSpp] bherb = mu_bherb + sigma_bherb * z_bherb;

  // Random effects
  matrix[nSpp, nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop]  pop_rfx  = sigma_pop  * z_pop;

  // Linear predictor
  vector[N] predS;
  for (i in 1:N) {
    predS[i] =
      b0[Spp[i]]
      + bclim[Spp[i]] * clim[i]
      + bendo[Spp[i]] * endo[i]
      + bherb[Spp[i]] * herb[i]
      + site_year_rfx[Spp[i], site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {
  // Priors for global means
  mu_b0    ~ normal(0, 2);
  mu_bclim ~ normal(0, 2);
  mu_bendo ~ normal(0, 2);
  mu_bherb ~ normal(0, 2);

  // Priors for species-level variation
  sigma_b0    ~ normal(0, 1);
  sigma_bclim ~ normal(0, 1);
  sigma_bendo ~ normal(0, 1);
  sigma_bherb ~ normal(0, 1);

  // Non-centered species deviations
  z_b0    ~ normal(0, 1);
  z_bclim ~ normal(0, 1);
  z_bendo ~ normal(0, 1);
  z_bherb ~ normal(0, 1);

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
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);
}
