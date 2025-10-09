data {
  // Indices
  int<lower=1> nSpp;         // Number of species
  int<lower=1> nsite_year;   // Number of site_years
  int<lower=1> nPop;         // Number of source populations
  int<lower=1> N;            // Number of observations for flowering model
  int<lower=1> nPlot;        // Number of plots

  // Inflorescence response data
  int<lower=1> Spp[N];       // Species index for each observation
  int<lower=1> site_year[N]; // site_year index for each observation
  int<lower=1> plot[N];      // Plot index for each observation
  int<lower=1> pop[N];       // Population index for each observation
  int<lower=0> y[N];         // Flowering counts at t+1
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N];  // Herbivory status (1 = present, 0 = absent)
  vector[N] clim;            // Climate covariate (e.g., precipitation, PET, Mahalanobis Distance)
}

parameters {
  // Fixed effects coefficients (species-specific)
  vector[nSpp] b0;          // Intercept for each species
  vector[nSpp] bendo;       // Effect of endophyte presence
  vector[nSpp] bherb;       // Effect of herbivory
  vector[nSpp] bclim;       // Effect of climate variable
  vector[nSpp] bendoclim;   // Interaction between endophyte and climate
  vector[nSpp] bendoherb;   // Interaction between endophyte and herbivory
  vector[nSpp] bclim2;      // Quadratic effect of climate
  vector[nSpp] bendoclim2;  // Interaction of endophyte and quadratic climate effect

  // Random effects
  real<lower=0> plot_tau;              // SD for plot-level random effect
  vector[nPlot] plot_rfx;              // Random effect for each plot
  real<lower=0> pop_tau;               // SD for population-level random effect
  vector[nPop] pop_rfx;                // Random effect for each population
  vector<lower=0>[nSpp] site_year_tau; // SD for site_year-level random effects
  matrix[nSpp, nsite_year] site_year_rfx; // Random effects for each species x site_year

  real<lower=0> phi; // Dispersion parameter for negative binomial
}

transformed parameters {
  vector[N] predF; // Predicted flowering counts for each observation

  // Loop over all observations to compute linear predictor
  for (i in 1:N) {
    predF[i] = 
      b0[Spp[i]] +                               // Species-specific intercept
      bendo[Spp[i]] * endo[i] +                  // Main effect of endophyte
      bclim[Spp[i]] * clim[i] +                  // Main effect of climate
      bherb[Spp[i]] * herb[i] +                  // Main effect of herbivory
      bendoclim[Spp[i]] * clim[i] * endo[i] +   // Endophyte x climate interaction
      bendoherb[Spp[i]] * endo[i] * herb[i] +   // Endophyte x herbivory interaction
      bclim2[Spp[i]] * square(clim[i]) +        // Quadratic climate effect
      bendoclim2[Spp[i]] * endo[i] * square(clim[i]) + // Endophyte x quadratic climate
      plot_rfx[plot[i]] +                        // Plot-level random effect
      pop_rfx[pop[i]] +                          // Population-level random effect
      site_year_rfx[Spp[i], site_year[i]];       // Site_year-level random effect for each species
  }
}

model {
  // Priors for fixed effects (weakly informative)
  b0 ~ normal(0, 5);    
  bendo ~ normal(0, 5);   
  bherb ~ normal(0, 5); 
  bclim ~ normal(0, 5);  
  bendoclim ~ normal(0, 5);  
  bendoherb ~ normal(0, 5); 
  bclim2 ~ normal(0, 5);  
  bendoclim2 ~ normal(0, 5);
  phi ~ normal(0, 5); // Prior for dispersion parameter

  // Priors for random effects
  plot_tau ~ normal(0, 1);
  plot_rfx ~ normal(0, plot_tau);   // Random effect for each plot

  pop_tau ~ normal(0, 1);
  pop_rfx ~ normal(0, pop_tau);     // Random effect for each population

  site_year_tau ~ normal(0, 1);     // SD for site_year random effects per species
  // Site_year-level random effects for each species (loop index f for flowering species)
  for (f in 1:nSpp) {
    for (sy in 1:nsite_year) {
      site_year_rfx[f, sy] ~ normal(0, site_year_tau[f]);
    }
  }

  // Likelihood: Negative binomial for flowering counts
  y ~ neg_binomial_2_log(predF, phi);
}

generated quantities {
  vector[N] log_lik; // Log-likelihood for each observation (for model comparison)

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | predF[i], phi);
  }
}
