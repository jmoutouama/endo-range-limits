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
  int<lower=0,upper=1> endo[N];   // Endophyte presence (biotic symbiosis)
  int<lower=0,upper=1> herb[N];   // Herbivory presence (biotic stress)
  vector[N] clim;            // Climate covariate (abiotic stress)
}

parameters {
  vector[nSpp] b0;            
  vector[nSpp] bendo;         
  vector[nSpp] bherb;         
  vector[nSpp] bclim;         
  vector[nSpp] bendoclim;     
  vector[nSpp] bendoherb;     
  vector[nSpp] bendoherbclim; 
  vector[nSpp] bherbclim;     
  vector[nSpp] bclim2;        

  real<lower=0> plot_tau;     
  vector[nPlot] plot_rfx;     
  real<lower=0> pop_tau;      
  vector[nPop] pop_rfx;       
  vector<lower=0>[nSpp] site_year_tau;
  matrix[nSpp, nsite_year] site_year_rfx;
}

transformed parameters {
  vector[N] predS;

  for (isurv in 1:N) {
    predS[isurv] = b0[Spp[isurv]] + 
                // Main effects for each covariate: baseline survival, endophyte, climate, herbivory
                bendo[Spp[isurv]] * endo[isurv] +
                bclim[Spp[isurv]] * clim[isurv] +
                bherb[Spp[isurv]] * herb[isurv] + 
                // 2-way interactions between endophyte and stressors
                bendoclim[Spp[isurv]] * clim[isurv] * endo[isurv] + 
                bendoherb[Spp[isurv]] * endo[isurv] * herb[isurv] + 
                bherbclim[Spp[isurv]] * herb[isurv] * clim[isurv] +  
                // 3-way interaction between endophyte, herbivory, and climate
                bendoherbclim[Spp[isurv]] * endo[isurv] * herb[isurv] * clim[isurv] + 
                // Quadratic climate effects
                bclim2[Spp[isurv]] * square(clim[isurv]) +  
                // Random effects to account for unobserved variation
                plot_rfx[plot[isurv]] +        // Plot-level variation
                pop_rfx[pop[isurv]] +          // Population-level variation
                site_year_rfx[Spp[isurv], site_year[isurv]]; // Site-year-level variation
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
  bherbclim ~ normal(0, 5);       
  bendoherbclim ~ normal(0, 5);   
  bclim2 ~ normal(0, 5);          

  // Priors for random effects to capture variation across plots, populations, and site-years
  plot_tau ~ inv_gamma(0.1,0.1);
  plot_rfx ~ normal(0, plot_tau);
  pop_tau ~ inv_gamma(0.1,0.1);
  pop_rfx ~ normal(0, pop_tau);
  site_year_tau ~ inv_gamma(0.1,0.1);
  for (s in 1:nSpp)
    site_year_rfx[s] ~ normal(0, site_year_tau[s]);

  // Likelihood: models survival as a logistic function of all predictors
  // Interactions test whether endophyte benefits are modulated by climate, herbivory, or both
  y ~ bernoulli_logit(predS);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    // Log-likelihood for model comparison (e.g., WAIC, LOO)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);
  }
}