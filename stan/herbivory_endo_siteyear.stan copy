// File: herbivory_endo_siteyear_ZINB.stan
data {
  int<lower=1> N;           // number of observations
  int<lower=1> nSpp;        // number of species
  int<lower=1> nSite;       // number of sites
  int<lower=1> nPop;        // number of populations
  int<lower=1> nPlot;       // number of plots
  int<lower=1> nsite_year;  // number of site-year combinations

  int<lower=0> y[N];        // observed herbivory counts
  int<lower=1> Spp[N];      // species index
  int<lower=1> site[N];     // site index
  int<lower=1> pop[N];      // population index
  int<lower=1> plot[N];     // plot index
  int<lower=1,upper=nsite_year> site_year[N]; // site-year index
  int<lower=0,upper=1> endo[N]; // Endophyte status
  vector[N] clim;           // optional climate covariate (ppt_scaled)
}

parameters {
  // Fixed effects
  vector[nSpp] b0;         // species intercepts
  vector[nSpp] bendo;      // Endo effect
  vector[nSpp] bclim;      // climate effect

  // Random effects
  real<lower=0> site_tau;
  vector[nSite] site_rfx;

  real<lower=0> pop_tau;
  vector[nPop] pop_rfx;

  real<lower=0> plot_tau;
  vector[nPlot] plot_rfx;

  vector<lower=0>[nSpp] site_year_tau;
  matrix[nSpp, nsite_year] site_year_rfx;

  real<lower=0> phi;        // overdispersion for NB
  real<lower=0, upper=1> zi; // probability of structural zero
}

transformed parameters {
  vector[N] pred; // predicted log counts
  for (i in 1:N){
    pred[i] = b0[Spp[i]] 
              + bendo[Spp[i]] * endo[i] 
              + bclim[Spp[i]] * clim[i]
              + site_rfx[site[i]]
              + pop_rfx[pop[i]]
              + plot_rfx[plot[i]]
              + site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  // Priors for fixed effects
  b0 ~ normal(0, 2);
  bendo ~ normal(0, 1);
  bclim ~ normal(0, 1);

  // Random effects
  site_tau ~ exponential(1);
  site_rfx ~ normal(0, site_tau);

  pop_tau ~ exponential(1);
  pop_rfx ~ normal(0, pop_tau);

  plot_tau ~ exponential(1);
  plot_rfx ~ normal(0, plot_tau);

  site_year_tau ~ inv_gamma(0.1, 0.1);
  for (s in 1:nSpp)
    for (sy in 1:nsite_year)
      site_year_rfx[s, sy] ~ normal(0, site_year_tau[s]);

  phi ~ gamma(2, 1);
  zi ~ beta(1, 1); // weakly informative prior for zero-inflation

  // Zero-inflated negative binomial likelihood
  for (i in 1:N){
    if (y[i] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | zi),
                            bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | pred[i], phi));
    else
      target += bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | pred[i], phi);
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    if (y[i] == 0)
      log_lik[i] = log_sum_exp(bernoulli_lpmf(1 | zi),
                               bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | pred[i], phi));
    else
      log_lik[i] = bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | pred[i], phi);
  }
}
