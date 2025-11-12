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
}

parameters {
  // Fixed effects (species-specific)
  vector[nSpp] b0;              
  vector[nSpp] bendo;           
  vector[nSpp] bherb;           
  vector[nSpp] bendoherb;       

  // Random effects
  real<lower=0> plot_tau;              
  vector[nPlot] plot_rfx;              
  real<lower=0> pop_tau;               
  vector[nPop] pop_rfx;                
  vector<lower=0>[nSpp] site_year_tau; 
  matrix[nSpp, nsite_year] site_year_rfx; 

  // NB overdispersion and zero-inflation
  real<lower=0> phi; 
  real<lower=0, upper=1> zi; // probability of structural zeros
}

transformed parameters {
  vector[N] predF; 

  for (i in 1:N) {
    predF[i] =
      b0[Spp[i]] +                              
      bendo[Spp[i]] * endo[i] +
      bherb[Spp[i]] * herb[i] +
      bendoherb[Spp[i]] * endo[i] * herb[i] +
      plot_rfx[plot[i]] +
      pop_rfx[pop[i]] +
      site_year_rfx[Spp[i], site_year[i]];
  }
}

model {
  // Priors
  b0 ~ normal(0, 1);    
  bendo ~ normal(0, 1);   
  bherb ~ normal(0, 1); 
  bendoherb ~ normal(0, 1);
  phi ~ gamma(2,0.1); 
  zi ~ beta(1, 1);  // weak prior for zero inflation

  // Random effects
  plot_tau ~ inv_gamma(0.1,0.1);
  plot_rfx ~ normal(0, plot_tau);  

  pop_tau ~ inv_gamma(0.1,0.1);
  pop_rfx ~ normal(0, pop_tau);    

  site_year_tau ~ inv_gamma(0.1,0.1);     
  for (f in 1:nSpp)
    site_year_rfx[f] ~ normal(0, site_year_tau[f]);

  // Zero-inflated negative binomial likelihood
  for (i in 1:N) {
    if (y[i] == 0)
      target += log_sum_exp(
        bernoulli_lpmf(1 | zi), // structural zero
        bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | predF[i], phi) // NB zero
      );
    else
      target += bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | predF[i], phi);
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    if (y[i] == 0)
      log_lik[i] = log_sum_exp(
        bernoulli_lpmf(1 | zi),
        bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | predF[i], phi)
      );
    else
      log_lik[i] = bernoulli_lpmf(0 | zi) + neg_binomial_2_log_lpmf(y[i] | predF[i], phi);
  }
}
