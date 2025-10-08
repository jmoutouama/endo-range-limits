data {
  // Indices
  int<lower=1> nSpp;         // Number of species
  int<lower=1> nsite_year;   // Number of site_years
  int<lower=1> nPop;         // Number of populations
  int<lower=1> N;            // Number of survival observations
  int<lower=1> nPlot;        // Number of plots

  // Survival data
  int<lower=1> Spp[N];       // Species index for each observation
  int<lower=1> site_year[N]; // Site_year index for each observation
  int<lower=1> plot[N];      // Plot index for each observation
  int<lower=1> pop[N];       // Population index for each observation
  int<lower=0,upper=1> y[N]; // Survival outcome (1 = survived, 0 = died)
  int<lower=0,upper=1> endo[N];  // Endophyte presence (1 = positive, 0 = negative)
  vector[N] clim;             // Climate covariate (e.g., precipitation)
}

parameters {
  // Fixed effects
  vector[nSpp] b0;            // Intercept for each species
  vector[nSpp] bendoclim;     // Coefficient for linear endophyte × climate interaction
  vector[nSpp] bendoclim2;    // Coefficient for quadratic endophyte × climate^2 interaction

  // Random effects
  real<lower=0> plot_tau;         // SD for plot-level random effects
  vector[nPlot] plot_rfx;         // Plot-level random effects
  real<lower=0> pop_tau;          // SD for population-level random effects
  vector[nPop] pop_rfx;           // Population-level random effects
  vector<lower=0>[nSpp] site_year_tau;        // SD for site_year-level random effects (per species)
  matrix[nSpp, nsite_year] site_year_rfx;    // Site_year-level random effects per species
}

transformed parameters {
  vector[N] predS;  // Predicted survival probabilities (logit scale)

  // Loop over all observations to compute predicted survival
  for (i in 1:N) {
    predS[i] = b0[Spp[i]] +                              // Species-specific intercept
               bendoclim[Spp[i]] * endo[i] * clim[i] +  // Linear endo × climate interaction
               bendoclim2[Spp[i]] * endo[i] * square(clim[i]) + // Quadratic endo × climate interaction
               plot_rfx[plot[i]] +                        // Plot-level random effect
               pop_rfx[pop[i]] +                          // Population-level random effect
               site_year_rfx[Spp[i], site_year[i]];      // Site_year-level random effect
  }
}

model {
  // Priors for fixed effects
  b0 ~ normal(0,5);
  bendoclim ~ normal(0,5);
  bendoclim2 ~ normal(0,5);

  // Priors and distributions for random effects
  plot_tau ~ normal(0,1);
  plot_rfx ~ normal(0, plot_tau);

  pop_tau ~ normal(0,1);
  pop_rfx ~ normal(0, pop_tau);

  site_year_tau ~ normal(0,1);  // Separate SD per species
  for (s in 1:nSpp)
    for (t in 1:nsite_year)
      site_year_rfx[s, t] ~ normal(0, site_year_tau[s]);

  // Likelihood: Bernoulli survival with logit link
  y ~ bernoulli_logit(predS);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for each observation

  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | predS[i]);  // For model comparison
}
