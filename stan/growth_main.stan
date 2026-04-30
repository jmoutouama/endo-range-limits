// Hierarchical Bayesian growth model (endophyte effect only)
data {
  int<lower=1> N;
  int<lower=1> nSpp;
  int<lower=1> nsite_year;
  int<lower=1> nPop;
  int<lower=1> nPlot;

  array[N] int<lower=1, upper=nSpp>  Spp;
  array[N] int<lower=1, upper=nsite_year> site_year;
  array[N] int<lower=1, upper=nPop>  pop;
  array[N] int<lower=1, upper=nPlot> plot;

  vector[N] y;               // growth response
  array[N] int<lower=0,upper=1> endo;
}

parameters {
  // Global means for main effects
  real mu_b0;
  real mu_bendo;

  // Species-level variation
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bendo;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bendo;

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
  // Species-level main-effect coefficients
  vector[nSpp] b0    = mu_b0    + sigma_b0    * z_b0;
  vector[nSpp] bendo = mu_bendo + sigma_bendo * z_bendo;

  // Random effects
  matrix[nSpp, nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop]  pop_rfx  = sigma_pop  * z_pop;

  // Linear predictor
  vector[N] predG;
  for (i in 1:N) {
    predG[i] =
      b0[Spp[i]]
      + bendo[Spp[i]] * endo[i]
      + site_year_rfx[Spp[i], site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {
  // Global priors
  mu_b0    ~ normal(0,2);
  mu_bendo ~ normal(0,2);

  // Species-level variation
  sigma_b0    ~ normal(0,1);
  sigma_bendo ~ normal(0,1);

  // Non-centered species deviations
  z_b0    ~ normal(0,1);
  z_bendo ~ normal(0,1);

  // Random effects
  sigma_site_year ~ normal(0,1);
  sigma_plot      ~ normal(0,1);
  sigma_pop       ~ normal(0,1);

  to_vector(z_site_year) ~ normal(0,1);
  z_plot ~ normal(0,1);
  z_pop  ~ normal(0,1);

  // Residual SD
  sigma ~ normal(0,1);

  // Likelihood
  y ~ normal(predG, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);
}
