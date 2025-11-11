data {
  // Indices
  int<lower=1> nSpp;           // Number of species
  int<lower=1> nsite_year;     // Number of site-years
  int<lower=1> nPop;           // Number of source populations
  int<lower=1> N;              // Number of observations for growth
  int<lower=1> nPlot;          // Number of plots

  // Observation-level data
  int<lower=1> Spp[N];         // Species index
  int<lower=1> site_year[N];   // Site-year index
  int<lower=1> plot[N];        // Plot index
  int<lower=1> pop[N];         // Population index
  vector[N] y;                 // Continuous response (growth from t to t+1)

  // Covariates
  int<lower=0,upper=1> endo[N]; // Endophyte status (biotic symbiosis)
  int<lower=0,upper=1> herb[N]; // Herbivory status (biotic stress)
  vector[N] clim;               // Climate covariate (abiotic stress)
}

parameters {
  // Fixed effects
  vector[nSpp] b0;             // Species-specific intercepts
  vector[nSpp] bendo;          // Endophyte main effect
  vector[nSpp] bherb;          // Herbivory main effect
  vector[nSpp] bclim;          // Climate main effect
  vector[nSpp] bendoclim;      // Endophyte × Climate interaction
  vector[nSpp] bendoherb;      // Endophyte × Herbivory interaction
  vector[nSpp] bendoherbclim;  // Endophyte × Herbivory × Climate interaction
  vector[nSpp] bclim2;         // Quadratic climate effect

  // Random effects
  real<lower=0> plot_tau;      
  vector[nPlot] plot_rfx;      
  real<lower=0> pop_tau;       
  vector[nPop] pop_rfx;        
  vector<lower=0>[nSpp] site_year_tau;
  matrix[nSpp, nsite_year] site_year_rfx;

  // Residual standard deviation
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] predG;  // Predicted mean growth

  for (i in 1:N) {
    predG[i] = b0[Spp[i]] + 
               // Main effects for each covariate: baseline growth, endophyte, climate, herbivory
               bendo[Spp[i]] * endo[i] +
               bclim[Spp[i]] * clim[i] +
               bherb[Spp[i]] * herb[i] +
               // Two-way interactions between endophyte and stressors
               bendoclim[Spp[i]] * endo[i] * clim[i] +
               bendoherb[Spp[i]] * endo[i] * herb[i] +
               // Three-way interaction between endophyte, herbivory, and climate
               bendoherbclim[Spp[i]] * endo[i] * herb[i] * clim[i] +
               // Quadratic climate effect
               bclim2[Spp[i]] * square(clim[i]) +
               // Random effects to account for unobserved variation
               plot_rfx[plot[i]] +               // Plot-level variation
               pop_rfx[pop[i]] +                 // Population-level variation
               site_year_rfx[Spp[i], site_year[i]]; // Site-year-level variation
  }
}

model {
  // Priors for fixed effects
  b0 ~ normal(0, 1);
  bendo ~ normal(0, 1);
  bherb ~ normal(0, 1);
  bclim ~ normal(0, 1);
  bendoclim ~ normal(0, 1);
  bendoherb ~ normal(0, 1);
  bendoherbclim ~ normal(0, 1);
  bclim2 ~ normal(0, 1);
  sigma ~ normal(0, 1);

  // Priors for random effects to capture variation across plots, populations, and site-years
  plot_tau ~ inv_gamma(0.1,0.1);
  plot_rfx ~ normal(0, plot_tau);
  pop_tau ~ inv_gamma(0.1,0.1);
  pop_rfx ~ normal(0, pop_tau);
  site_year_tau ~ inv_gamma(0.1,0.1);
  for (s in 1:nSpp)
    site_year_rfx[s] ~ normal(0, site_year_tau[s]);

  // Likelihood: models growth as a normal function of all predictors
  // Interactions test whether endophyte benefits are modulated by climate, herbivory, or both
  y ~ normal(predG, sigma);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for model comparison (e.g., WAIC, LOO)
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);
  }
}
