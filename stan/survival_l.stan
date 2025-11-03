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
  int<lower=0,upper=1> herb[N];   // Herbivory presence
  vector[N] clim;            // Climate covariate
}

parameters {
  vector[nSpp] b0;            
  vector[nSpp] bendo;         
  vector[nSpp] bherb;         
  vector[nSpp] bclim;         
  vector[nSpp] bendoclim;     
  vector[nSpp] bendoherb;     
  vector[nSpp] bendoherbclim; 

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
                    bendo[Spp[isurv]] * endo[isurv] +
                    bclim[Spp[isurv]] * clim[isurv] +
                    bherb[Spp[isurv]] * herb[isurv] + 
                    bendoclim[Spp[isurv]] * clim[isurv] * endo[isurv] + 
                    bendoherb[Spp[isurv]] * endo[isurv] * herb[isurv] + 
                    bendoherbclim[Spp[isurv]] * endo[isurv] * herb[isurv] * clim[isurv] + 
                    plot_rfx[plot[isurv]] +
                    pop_rfx[pop[isurv]] +
                    site_year_rfx[Spp[isurv], site_year[isurv]];
  }
}

model {
  // Priors
  b0 ~ normal(0, 5);              
  bendo ~ normal(0, 5);           
  bherb ~ normal(0, 5);           
  bclim ~ normal(0, 5);           
  bendoclim ~ normal(0, 5);       
  bendoherb ~ normal(0, 5);       
  bendoherbclim ~ normal(0, 5);   

  plot_tau ~ normal(0, 1);
  plot_rfx ~ normal(0, plot_tau);
  pop_tau ~ normal(0, 1);
  pop_rfx ~ normal(0, pop_tau);
  site_year_tau ~ normal(0, 1);
  for (s in 1:nSpp)
    site_year_rfx[s] ~ normal(0, site_year_tau[s]);

  y ~ bernoulli_logit(predS);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);
  }
}
