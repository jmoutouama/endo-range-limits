// Hierarchical Bayesian survival/flowering model (Bernoulli-logit)
// Response: spike count (integer count, not 0/1 anymore)
// Species-specific fixed effects with global hyperpriors
// Random effects for site-year, plot, and source population

data {
  int<lower=1> N;            // Number of observations
  int<lower=1> nSpp;         // Number of species
  int<lower=1> nsite_year;   // Number of site_year combinations
  int<lower=1> nPop;         // Number of source populations
  int<lower=1> nPlot;        // Number of plots

  array[N] int<lower=1, upper=nSpp> Spp;
  array[N] int<lower=1, upper=nsite_year> site_year;
  array[N] int<lower=1, upper=nPop> pop;
  array[N] int<lower=1, upper=nPlot> plot;

  array[N] int<lower=0> y;       // Spike count
  array[N] int<lower=0, upper=1> endo; // Endophyte status
  array[N] int<lower=0, upper=1> herb; // Herbivory status

  vector[N] clim;   // Climate covariate
}

parameters {
  // Global means for species-level effects
  real mu_b0;
  real mu_bclim;
  real mu_bendo;
  real mu_bherb;
  real mu_bendoclim;
  real mu_bherbclim;
  real mu_bendoherb;
  real mu_bendoherbclim;

  // Species-level SDs
  real<lower=0> sigma_b0;
  real<lower=0> sigma_bclim;
  real<lower=0> sigma_bendo;
  real<lower=0> sigma_bherb;
  real<lower=0> sigma_bendoclim;
  real<lower=0> sigma_bherbclim;
  real<lower=0> sigma_bendoherb;
  real<lower=0> sigma_bendoherbclim;

  // Non-centered species deviations
  vector[nSpp] z_b0;
  vector[nSpp] z_bclim;
  vector[nSpp] z_bendo;
  vector[nSpp] z_bherb;
  vector[nSpp] z_bendoclim;
  vector[nSpp] z_bherbclim;
  vector[nSpp] z_bendoherb;
  vector[nSpp] z_bendoherbclim;

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
  vector[nSpp] b0            = mu_b0            + sigma_b0            * z_b0;
  vector[nSpp] bclim         = mu_bclim         + sigma_bclim         * z_bclim;
  vector[nSpp] bendo         = mu_bendo         + sigma_bendo         * z_bendo;
  vector[nSpp] bherb         = mu_bherb         + sigma_bherb         * z_bherb;
  vector[nSpp] bendoclim     = mu_bendoclim     + sigma_bendoclim     * z_bendoclim;
  vector[nSpp] bherbclim     = mu_bherbclim     + sigma_bherbclim     * z_bherbclim;
  vector[nSpp] bendoherb     = mu_bendoherb     + sigma_bendoherb     * z_bendoherb;
  vector[nSpp] bendoherbclim = mu_bendoherbclim + sigma_bendoherbclim * z_bendoherbclim;

  // Random effects
  matrix[nSpp, nsite_year] site_year_rfx = sigma_site_year * z_site_year;
  vector[nPlot] plot_rfx = sigma_plot * z_plot;
  vector[nPop]  pop_rfx  = sigma_pop  * z_pop;

  // Linear predictor
  vector[N] predF;
  for (i in 1:N) {
    predF[i] =
      b0[Spp[i]]
      + bclim[Spp[i]]         * clim[i]
      + bendo[Spp[i]]         * endo[i]
      + bherb[Spp[i]]         * herb[i]
      + bendoclim[Spp[i]]     * endo[i] * clim[i]
      + bherbclim[Spp[i]]     * herb[i] * clim[i]
      + bendoherb[Spp[i]]     * endo[i] * herb[i]
      + bendoherbclim[Spp[i]] * endo[i] * herb[i] * clim[i]
      + site_year_rfx[Spp[i], site_year[i]]
      + plot_rfx[plot[i]]
      + pop_rfx[pop[i]];
  }
}

model {
  // Priors for global means
  mu_b0            ~ normal(0, 2);
  mu_bclim         ~ normal(0, 2);
  mu_bendo         ~ normal(0, 2);
  mu_bherb         ~ normal(0, 2);
  mu_bendoclim     ~ normal(0, 2);
  mu_bherbclim     ~ normal(0, 2);
  mu_bendoherb     ~ normal(0, 2);
  mu_bendoherbclim ~ normal(0, 2);

  // Priors for species SDs
  sigma_b0            ~ normal(0, 1);
  sigma_bclim         ~ normal(0, 1);
  sigma_bendo         ~ normal(0, 1);
  sigma_bherb         ~ normal(0, 1);
  sigma_bendoclim     ~ normal(0, 1);
  sigma_bherbclim     ~ normal(0, 1);
  sigma_bendoherb     ~ normal(0, 1);
  sigma_bendoherbclim ~ normal(0, 1);

  // Non-centered species deviations
  z_b0            ~ normal(0, 1);
  z_bclim         ~ normal(0, 1);
  z_bendo         ~ normal(0, 1);
  z_bherb         ~ normal(0, 1);
  z_bendoclim     ~ normal(0, 1);
  z_bherbclim     ~ normal(0, 1);
  z_bendoherb     ~ normal(0, 1);
  z_bendoherbclim ~ normal(0, 1);

  // Random-effect SD priors
  sigma_site_year ~ normal(0, 1);
  sigma_plot      ~ normal(0, 1);
  sigma_pop       ~ normal(0, 1);

  // Non-centered random effects
  to_vector(z_site_year) ~ normal(0, 1);
  z_plot                 ~ normal(0, 1);
  z_pop                  ~ normal(0, 1);

  // Dispersion for Negative Binomial
  phi ~ exponential(1);

  // Likelihood
  y ~ neg_binomial_2_log(predF, phi);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | predF[i], phi);
}
