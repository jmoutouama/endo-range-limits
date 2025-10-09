data {
  int<lower=1> nSpp;
  int<lower=1> nsite_year;
  int<lower=1> nPop;
  int<lower=1> N;
  int<lower=1> nPlot;

  int<lower=1> Spp[N];
  int<lower=1> site_year[N];
  int<lower=1> plot[N];
  int<lower=1> pop[N];

  int<lower=0> y[N];
  vector[N] clim;
}

parameters {
  vector[nSpp] b0;
  vector[nSpp] bclim;

  real<lower=0> plot_tau;
  vector[nPlot] plot_rfx;

  real<lower=0> pop_tau;
  vector[nPop] pop_rfx;

  vector<lower=0>[nSpp] site_year_tau;
  matrix[nSpp, nsite_year] site_year_rfx;

  real<lower=0> phi;
}

transformed parameters {
  vector[N] pred;

  for (i in 1:N){
    pred[i] = b0[Spp[i]] +
              bclim[Spp[i]] * clim[i] +
              plot_rfx[plot[i]] +
              pop_rfx[pop[i]] +
              site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  b0 ~ normal(0,5);
  bclim ~ normal(0,5);
  phi ~ normal(0,5);

  plot_tau ~ normal(0,1);
  plot_rfx ~ normal(0, plot_tau);

  pop_tau ~ normal(0,1);
  pop_rfx ~ normal(0, pop_tau);

  site_year_tau ~ normal(0,1);
  for (i in 1:nSpp)
    for (j in 1:nsite_year)
      site_year_rfx[i,j] ~ normal(0, site_year_tau[i]);

  y ~ neg_binomial_2_log(pred, phi);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | pred[i], phi);
}
