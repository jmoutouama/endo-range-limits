data {
  // Indices
  int<lower=1> nSpp;         
  int<lower=1> nsite_year;   
  int<lower=1> nPop;         
  int<lower=1> N;            
  int<lower=1> nPlot;        

  // Observation-level data
  int<lower=1> Spp[N];       
  int<lower=1> site_year[N]; 
  int<lower=1> plot[N];      
  int<lower=1> pop[N];       
  int<lower=0> y[N];         
  int<lower=0,upper=1> endo[N];  
  int<lower=0,upper=1> herb[N];  
  vector[N] clim;            
}

parameters {
  // Fixed effects (species-specific)
  vector[nSpp] b0;              
  vector[nSpp] bendo;           
  vector[nSpp] bherb;           
  vector[nSpp] bclim;           
  vector[nSpp] bendoclim;       
  vector[nSpp] bendoherb;       
  vector[nSpp] bendoherbclim;   

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
      bendo[Spp[i]] * endo[i] +
      bclim[Spp[i]] * clim[i] +
      bherb[Spp[i]] * herb[i] +
      bendoclim[Spp[i]] * clim[i] * endo[i] +
      bendoherb[Spp[i]] * endo[i] * herb[i] +
      bendoherbclim[Spp[i]] * endo[i] * herb[i] * clim[i] +
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
  bendoherbclim ~ normal(0, 5);
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
