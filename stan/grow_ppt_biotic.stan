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

  vector[N] y;
  int<lower=0,upper=1> endo[N];
  int<lower=0,upper=1> herb[N];
}

parameters {
  vector[nSpp] b0;       // Intercept per species
  vector[nSpp] bendo;    // Endophyte effect per species
  vector[nSpp] bherb;    // Herbivory effect per species

  real<lower=0> plot_tau;
  vector[nPlot] plot_rfx;

  real<lower=0> pop_tau;
  vector[nPop] pop_rfx;

  real<lower=0> site_year_tau;
  matrix[nSpp, nsite_year] site_year_rfx;

  real<lower=0> sigma;
}

transformed parameters {
  vector[N] pred;

  for (i in 1:N){
    pred[i] = b0[Spp[i]] +
              bendo[Spp[i]] * endo[i] +
              bherb[Spp[i]] * herb[i] +
              plot_rfx[plot[i]] +
              pop_rfx[pop[i]] +
              site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  b0 ~ normal(0, 5);
  bendo ~ normal(0, 5);
  bherb ~ normal(0, 5);
  sigma ~ normal(0, 5);

  plot_tau ~ normal(0, 1);
  plot_rfx ~ normal(0, plot_tau);

  pop_tau ~ normal(0, 1);
  pop_rfx ~ normal(0, pop_tau);

  site_year_tau ~ normal(0, 1);
  for (i in 1:nSpp)
    for (j in 1:nsite_year)
      site_year_rfx[i,j] ~ normal(0, site_year_tau);

  y ~ normal(pred, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | pred[i], sigma);
}
