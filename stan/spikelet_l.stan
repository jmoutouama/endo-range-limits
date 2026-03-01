data {
  // indices
  int<lower=1> nSpp;       // Number of species        
  int<lower=1> nsite_year; // Number of site_years    
  int<lower=1> nPop;       // Number of source populations       
  int<lower=1> N;          // Number of data points for the inflorescence model
  int<lower=1> nPlot;      // Number of plots   

  // Flowering data
  int<lower=1> Spp[N];         // Species index
  int<lower=1> site_year[N];   // site_year index
  int<lower=1> plot[N];        // Plot index 
  int<lower=1> pop[N];         // Population index
  int<lower=0> y[N];           // Response variable: flowering at time t+1
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N];  // Herbivory status (1 = herbivory present, 0 = absent)
  vector[N] clim;               // Climate covariate
}

parameters {
  // Fixed effects
  vector[nSpp] b0;              
  vector[nSpp] bendo;           
  vector[nSpp] bherb;           
  vector[nSpp] bclim;           
  vector[nSpp] bendoclim;       
  vector[nSpp] bendoherb;       
  vector[nSpp] bendoherbclim;   
  vector[nSpp] bherbclim;       // climate × herbivory interaction

  // Random effects
  real<lower=0> plot_tau;         
  vector[nPlot] plot_rfx;        

  real<lower=0> pop_tau;          
  vector[nPop] pop_rfx;          

  vector<lower=0>[nSpp] site_year_tau;  
  matrix[nSpp, nsite_year] site_year_rfx;  

  real<lower=0> phi;  // Dispersion parameter for negative binomial
}

transformed parameters {
  vector[N] predF; // Predicted flowering response

  for (iflow in 1:N) {
    predF[iflow] = b0[Spp[iflow]] +
                    // Main effects
                    bendo[Spp[iflow]] * endo[iflow] +
                    bclim[Spp[iflow]] * clim[iflow] +
                    bherb[Spp[iflow]] * herb[iflow] +
                    // Two-way interactions
                    bendoclim[Spp[iflow]] * clim[iflow] * endo[iflow] +
                    bendoherb[Spp[iflow]] * endo[iflow] * herb[iflow] +
                    bherbclim[Spp[iflow]] * clim[iflow] * herb[iflow] +  // UPDATED NAME
                    // Three-way interaction
                    bendoherbclim[Spp[iflow]] * endo[iflow] * herb[iflow] * clim[iflow] +
                    // Random effects
                    plot_rfx[plot[iflow]] +
                    pop_rfx[pop[iflow]] +
                    site_year_rfx[Spp[iflow], site_year[iflow]];
  }
}

model {
  // Priors on fixed effects
  b0 ~ normal(0, 1);
  bendo ~ normal(0, 1);
  bherb ~ normal(0, 1);
  bclim ~ normal(0, 1);
  bendoclim ~ normal(0, 1);
  bendoherb ~ normal(0, 1);
  bendoherbclim ~ normal(0, 1);
  bherbclim ~ normal(0, 1); // UPDATED

  phi ~ gamma(2, 1);

  // Priors for random effects
  plot_tau ~ inv_gamma(0.1, 0.1);
  plot_rfx ~ normal(0, plot_tau);

  pop_tau ~ inv_gamma(0.1, 0.01);
  pop_rfx ~ normal(0, pop_tau);

  site_year_tau ~ inv_gamma(0.1, 0.1);
  for (i in 1:nSpp) {
    for (j in 1:nsite_year) {
      site_year_rfx[i, j] ~ normal(0, site_year_tau[i]);
    }
  }

  // Likelihood
  y ~ neg_binomial_2_log(predF, phi);
}

generated quantities {
  vector[N] log_lik;
  for (nfi in 1:N) {
    log_lik[nfi] = neg_binomial_2_log_lpmf(y[nfi] | predF[nfi], phi);
  }
}
