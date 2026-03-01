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
  vector[N] y;                

  // Covariates
  int<lower=0,upper=1> endo[N]; 
  int<lower=0,upper=1> herb[N]; 
  vector[N] clim;               
}

parameters {
  // Fixed effects
  vector[nSpp] b0;            
  vector[nSpp] bendo;         
  vector[nSpp] bherb;         
  vector[nSpp] bclim;         

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
  vector[N] predG;

  for (i in 1:N) {
    predG[i] = b0[Spp[i]] + 
               bendo[Spp[i]] * endo[i] +
               bclim[Spp[i]] * clim[i] +
               bherb[Spp[i]] * herb[i] +
               plot_rfx[plot[i]] +
               pop_rfx[pop[i]] +
               site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  // Priors for fixed effects
  b0 ~ normal(0, 1);              
  bendo ~ normal(0, 1);           
  bherb ~ normal(0, 1);           
  bclim ~ normal(0, 1);               
  sigma ~ normal(0, 1);

  // Random effects
  plot_tau ~ inv_gamma(0.1,0.1);
  plot_rfx ~ normal(0, plot_tau);
  pop_tau ~ inv_gamma(0.1,0.1);
  pop_rfx ~ normal(0, pop_tau);
  site_year_tau ~ inv_gamma(0.1,0.1);
  for (s in 1:nSpp)
    site_year_rfx[s] ~ normal(0, site_year_tau[s]);

  // Likelihood
  y ~ normal(predG, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);
  }
}
