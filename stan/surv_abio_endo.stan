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
  int<lower=0,upper=1> y[N]; // Survival outcome
  int<lower=0,upper=1> endo[N];   // Endophyte presence
  vector[N] clim;            // Climate covariate (e.g., precipitation)
}

parameters {
  vector[nSpp] b0;            // Intercept per species
  vector[nSpp] bendo;         // Endophyte main effect
  vector[nSpp] bclim;         // Linear climate effect
  vector[nSpp] bendoclim;     // Endophyte × climate interaction
  vector[nSpp] bclim2;        // Quadratic climate effect

  real<lower=0> plot_tau;     
  vector[nPlot] plot_rfx;     
  real<lower=0> pop_tau;      
  vector[nPop] pop_rfx;       
  vector<lower=0>[nSpp] site_year_tau;
  matrix[nSpp, nsite_year] site_year_rfx;
}

transformed parameters {
  vector[N] predS;

  for (i in 1:N) {
    predS[i] = b0[Spp[i]] +
               // Main effects
               bendo[Spp[i]] * endo[i] +
               bclim[Spp[i]] * clim[i] +
               // Endophyte × climate interactions
               bendoclim[Spp[i]] * endo[i] * clim[i] +
               // Quadratic climate effects
               bclim2[Spp[i]] * square(clim[i]) +
               // Random effects
               plot_rfx[plot[i]] +
               pop_rfx[pop[i]] +
               site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  // Priors
  b0 ~ normal(0, 5);
  bendo ~ normal(0, 5);
  bclim ~ normal(0, 5);
  bendoclim ~ normal(0, 5);
  bclim2 ~ normal(0, 5);

  plot_tau ~ normal(0, 1);
  plot_rfx ~ normal(0, plot_tau);
  pop_tau ~ normal(0, 1);
  pop_rfx ~ normal(0, pop_tau);
  site_year_tau ~ normal(0, 1);
  for (s in 1:nSpp)
    site_year_rfx[s] ~ normal(0, site_year_tau[s]);

  // Likelihood
  y ~ bernoulli_logit(predS);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);
  }
}
