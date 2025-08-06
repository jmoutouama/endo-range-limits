data {
  // indices
  int<lower=1> nSpp;         // Number of species in the model
  int<lower=1> nsite_year;           // Number of site_years
  int<lower=1> nPop;            // Number of source populations
  int<lower=1> N;               // Number of data points for the survival model
  int<lower=1> nPlot;          // Number of plots
  // Survival data
  int<lower=1> Spp[N];     // Species index for each survival data point
  int<lower=1> site_year[N];     // site_year_year index for each survival data point
  int<lower=1> plot[N];       // Plot index for each survival data point
  int<lower=1> pop[N];        // Population index for each survival data point
  int<lower=0,upper=1> y[N];  // Survival status at time t+1 (1 = survived, 0 = did not survive)
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N];  // Herbivory status (1 = herbivory present, 0 = no herbivory)
  vector[N] clim;             // Climate covariate (e.g., precipitation, PET, or Mahalanobis Distance)
}

parameters {
  // Fixed effects coefficients
  vector[nSpp] b0;        // Intercept term for each species
  vector[nSpp] bendo;     // Coefficient for endophyte status for each species
  vector[nSpp] bherb;     // Coefficient for herbivory for each species
  vector[nSpp] bclim;     // Coefficient for climate variable for each species
  vector[nSpp] bendoclim; // Coefficient for the interaction between endophyte status and climate
  vector[nSpp] bendoherb; // Coefficient for the interaction between endophyte status and herbivory
  vector[nSpp] bclim2;    // Coefficient for the quadratic effect of climate for each species
  vector[nSpp] bendoclim2; // Coefficient for the interaction between endophyte status and the quadratic effect of climate

  // Random effects variances
  real<lower=0> plot_tau;      // Variance for plot-level random effects
  vector[nPlot] plot_rfx;   // Random effects for each plot
  real<lower=0> pop_tau;       // Variance for population-level random effects
  vector[nPop] pop_rfx;      // Random effects for each population
  vector<lower=0>[nSpp] site_year_tau;      // Variance for site_year-level random effects
  matrix[nSpp, nsite_year] site_year_rfx; // Random effects for each species at each site_year
}

transformed parameters {
  vector[N] predS;  // Predicted survival probabilities
  
  // Loop over all survival data points and calculate the predicted survival probability
  for (isurv in 1:N) {
    predS[isurv] = b0[Spp[isurv]] + 
                // Main effects for each covariate
                bendo[Spp[isurv]] * endo[isurv] +
                bclim[Spp[isurv]] * clim[isurv] +
                bherb[Spp[isurv]] * herb[isurv] + 
                // 2-way interaction effects between covariates
                bendoclim[Spp[isurv]] * clim[isurv] * endo[isurv] + // Correct indexing
                bendoherb[Spp[isurv]] * endo[isurv] * herb[isurv] + // Correct indexing
                // Quadratic climate effects
                bclim2[Spp[isurv]] * square(clim[isurv]) +  
                bendoclim2[Spp[isurv]] * endo[isurv] * square(clim[isurv]) + 
                // Random effects
                plot_rfx[plot[isurv]] +  // Plot-level random effect
                pop_rfx[pop[isurv]] +    // Population-level random effect
                site_year_rfx[Spp[isurv], site_year[isurv]]; // site_year-level random effect for each species at each site_year
  }
}

model {
  // Priors for fixed effects (Normally distributed with mean 0 and standard deviation 10)
  b0 ~ normal(0, 5);    
  bendo ~ normal(0, 5);   
  bherb ~ normal(0, 5); 
  bclim ~ normal(0, 5);  
  bendoclim ~ normal(0, 5);  
  bendoherb ~ normal(0, 5); 
  bclim2 ~ normal(0, 5);  
  bendoclim2 ~ normal(0, 5);
  
  // Priors for random effects variances (inverse gamma distribution)
  plot_tau ~ normal(0, 1);
  for (i in 1:nPlot) {
    plot_rfx[i] ~ normal(0, plot_tau);  // Random effects for each plot
  }
  
  pop_tau ~ normal(0, 1);
  for (i in 1:nPop) {
    pop_rfx[i] ~ normal(0, pop_tau);   // Random effects for each population
  }
  
  // Priors for site_year-level random effects variances (separate for each species)
  site_year_tau ~ normal(0, 1);  // Separate variance for each species
  for (i in 1:nSpp) {
    for (j in 1:nsite_year) {
      site_year_rfx[i, j] ~ normal(0, site_year_tau[i]);  // Random effects for each site_year and species
    }
  }

  // Likelihood for survival (Bernoulli logistic regression)
  y ~ bernoulli_logit(predS);  // Survived or not, based on the predicted survival probabilities
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for each survival data point
  
  // Loop over all survival data points and calculate the log-likelihood
  for (i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);  // Log-likelihood for Bernoulli distribution
  }
}
