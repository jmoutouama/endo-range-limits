data {
  // Indices
  int<lower=1> nSpp;         // Number of species
  int<lower=1> nsite_year;   // Number of site_years
  int<lower=1> nPop;         // Number of source populations
  int<lower=1> N;            // Number of observations for flowering model
  int<lower=1> nPlot;        // Number of plots

  // Observation-level data
  int<lower=1> Spp[N];       // Species index
  int<lower=1> site_year[N]; // Site-year index
  int<lower=1> plot[N];      // Plot index
  int<lower=1> pop[N];       // Population index
  int<lower=0> y[N];         // Flowering counts at t+1
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive)
  vector[N] clim;            // Climate covariate
}

parameters {
  // Fixed effects (species-specific)
  vector[nSpp] b0;              
  vector[nSpp] bendo;           
  vector[nSpp] bclim;           // Linear climate effect
  vector[nSpp] bendoclim;       // Interaction: endophyte x climate
  vector[nSpp] bclim2;          // Quadratic climate effect

  // Random effects
  real<lower=0> plot_tau;              
  vector[nPlot] plot_rfx;              
  real<lower=0> pop_tau;               
  vector[nPop] pop_rfx;                
  vector<lower=0>[nSpp] site_year_tau; 
  matrix[nSpp, nsite_year] site_year_rfx; 

  real<lower=0> phi; 
}

transformed parameters {
  vector[N] predF; 

  for (i in 1:N) {
    predF[i] =
      b0[Spp[i]] +                              
      // Main effects
      bendo[Spp[i]] * endo[i] +
      bclim[Spp[i]] * clim[i] +
      // Interaction
      bendoclim[Spp[i]] * endo[i] * clim[i] +
      // Quadratic climate effect
      bclim2[Spp[i]] * square(clim[i]) +
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
  bclim ~ normal(0, 5);  
  bendoclim ~ normal(0, 5);  
  bclim2 ~ normal(0, 5);  
  phi ~ normal(0, 5); 

  // Priors for random effects
  plot_tau ~ normal(0, 1);
  plot_rfx ~ normal(0, plot_tau);  

  pop_tau ~ normal(0, 1);
  pop_rfx ~ normal(0, pop_tau);    

  site_year_tau ~ normal(0, 1);     
  for (f in 1:nSpp)
    site_year_rfx[f] ~ normal(0, site_year_tau[f]);

  // Likelihood
  y ~ neg_binomial_2_log(predF, phi);
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | predF[i], phi);
  }
}
