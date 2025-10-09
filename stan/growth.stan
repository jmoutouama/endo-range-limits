data {
  // Indices
  int<lower=1> nSpp;           // Number of species
  int<lower=1> nsite_year;     // Number of site_years
  int<lower=1> nPop;           // Number of source populations
  int<lower=1> N;              // Number of observations for growth
  int<lower=1> nPlot;          // Number of plots
  
  // Observation-level data
  int<lower=1> Spp[N];         // Species index
  int<lower=1> site_year[N];   // site_year index
  int<lower=1> plot[N];        // Plot index
  int<lower=1> pop[N];         // Population index
  vector[N] y;                 // Continuous response (growth from t to t+1)
  
  // Binary covariates
  int<lower=0,upper=1> endo[N]; // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N]; // Herbivory status (1 = present, 0 = absent)
  
  // Continuous covariate
  vector[N] clim;               // Climate covariate (e.g., precipitation, PET, or Mahalanobis distance)
}

parameters {
  // Fixed effects
  vector[nSpp] b0;             // Species-specific intercepts
  vector[nSpp] bendo;          // Endophyte main effect
  vector[nSpp] bherb;          // Herbivory main effect
  vector[nSpp] bclim;          // Climate effect
  vector[nSpp] bendoclim;      // Endophyte × Climate interaction
  vector[nSpp] bendoherb;      // Endophyte × Herbivory interaction
  vector[nSpp] bclim2;         // Quadratic climate effect
  vector[nSpp] bendoclim2;     // Endophyte × Climate² interaction
  
  // Random effects
  real<lower=0> plot_tau;      // SD for plot-level random effects
  vector[nPlot] plot_rfx;      // Plot-level random effects
  
  real<lower=0> pop_tau;       // SD for population-level random effects
  vector[nPop] pop_rfx;        // Population-level random effects
  
  vector<lower=0>[nSpp] site_year_tau;     // SD for site_year random effects (per species)
  matrix[nSpp, nsite_year] site_year_rfx;  // site_year random effects per species
  
  // Residual variance
  real<lower=0> sigma;         // Residual SD
}

transformed parameters {
  vector[N] predG;  // Predicted mean growth
  
  for (i in 1:N) {
    predG[i] = b0[Spp[i]] +
               // Main effects
               bendo[Spp[i]] * endo[i] +
               bherb[Spp[i]] * herb[i] +
               bclim[Spp[i]] * clim[i] +
               // Interactions
               bendoclim[Spp[i]] * endo[i] * clim[i] +
               bendoherb[Spp[i]] * endo[i] * herb[i] +
               // Quadratic climate terms
               bclim2[Spp[i]] * square(clim[i]) +
               bendoclim2[Spp[i]] * endo[i] * square(clim[i]) +
               // Random effects
               plot_rfx[plot[i]] +
               pop_rfx[pop[i]] +
               site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  // Priors for fixed effects
  b0 ~ normal(0, 5);
  bendo ~ normal(0, 5);
  bherb ~ normal(0, 5);
  bclim ~ normal(0, 5);
  bendoclim ~ normal(0, 5);
  bendoherb ~ normal(0, 5);
  bclim2 ~ normal(0, 5);
  bendoclim2 ~ normal(0, 5);
  sigma ~ normal(0, 5);
  
  // Random effect priors
  plot_tau ~ normal(0, 1);
  plot_rfx ~ normal(0, plot_tau);
  
  pop_tau ~ normal(0, 1);
  pop_rfx ~ normal(0, pop_tau);
  
  site_year_tau ~ normal(0, 1);
  for (g in 1:nSpp)
    for (t in 1:nsite_year)
      site_year_rfx[g, t] ~ normal(0, site_year_tau[g]);
  
  // Likelihood
  y ~ normal(predG, sigma);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for model comparison
  
  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);
}
