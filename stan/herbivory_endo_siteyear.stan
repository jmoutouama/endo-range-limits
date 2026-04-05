// File: herbivory_demography_NB_trunc.stan

data {
  int<lower=1> N;           // total data points
  int<lower=1> nSpp;        // number of species
  int<lower=1> nSite;       // number of sites
  int<lower=1> nPop;        // number of populations
  int<lower=1> nPlot;       // number of plots
  int<lower=1> nsite_year;  // number of site-year combos

  // Grouping indices
  int<lower=1, upper=nSpp> Spp[N];
  int<lower=1, upper=nSite> site[N];
  int<lower=1, upper=nPop> pop[N];
  int<lower=1, upper=nPlot> plot[N];
  int<lower=1, upper=nsite_year> site_year[N];

  // Predictors
  vector[N] endo;   // 0/1
  vector[N] herb;   // 0/1
  vector[N] clim;   // scaled ppt

  // Response
  int<lower=0> y[N];
}

parameters {
  // ------------------
  // Fixed effects
  // ------------------
  vector[nSpp] b0;         // species intercept
  vector[nSpp] bendo;      // endophyte effect per species
  vector[nSpp] bclim;      // climate effect per species
  vector[nSpp] bendo_clim; // endophyte × climate per species

  // ------------------
  // Random effects
  // ------------------
  real<lower=0> site_tau;
  vector[nSite] site_rfx;

  real<lower=0> pop_tau;
  vector[nPop] pop_rfx;

  real<lower=0> plot_tau;
  vector[nPlot] plot_rfx;

  vector<lower=0>[nSpp] siteyear_tau;
  matrix[nSpp, nsite_year] siteyear_rfx;

  // ------------------
  // NB overdispersion
  // ------------------
  real<lower=0> phi;
}

transformed parameters {
  vector[N] pred; // linear predictor

  for (i in 1:N){
    pred[i] = b0[Spp[i]]
              + bendo[Spp[i]] * endo[i]
              + bclim[Spp[i]] * clim[i]
              + bendo_clim[Spp[i]] * endo[i] * clim[i]
              + site_rfx[site[i]]
              + pop_rfx[pop[i]]
              + plot_rfx[plot[i]]
              + siteyear_rfx[Spp[i], site_year[i]];
  }
}

model {
  // ------------------
  // Priors for fixed effects
  // ------------------
  b0 ~ normal(0, 500);
  bendo ~ normal(0, 100);
  bclim ~ normal(0, 100);
  bendo_clim ~ normal(0, 100);

  // ------------------
  // Random effects priors
  // ------------------
  site_tau ~ inv_gamma(0.2, 0.2);
  site_rfx ~ normal(0, site_tau);

  pop_tau ~ inv_gamma(0.2, 0.2);
  pop_rfx ~ normal(0, pop_tau);

  plot_tau ~ inv_gamma(0.2, 0.2);
  plot_rfx ~ normal(0, plot_tau);

  siteyear_tau ~ inv_gamma(0.2, 0.2);
  for (s in 1:nSpp)
    for (sy in 1:nsite_year)
      siteyear_rfx[s, sy] ~ normal(0, siteyear_tau[s]);

  phi ~ gamma(2, 1);

  // ------------------
  // Sampling: zero-truncated NB
  // ------------------
  for (i in 1:N){
    y[i] ~ neg_binomial_2_log(pred[i], phi);
    target += - log1m(neg_binomial_2_log_lpmf(0 | pred[i], phi));
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | pred[i], phi);
  }
}
