data {
  // Basic hierarchical structure
  int<lower=1> n_sites;       // Total number of sites
  int<lower=1> n_pops;        // Total number of source populations
  int<lower=1> n_plots;       // Total number of plots
  int<lower=1> n_species;     // Total number of species

  // Survival data
  int<lower=1> n_s;           // Number of survival observations
  int<lower=1> site_s[n_s];   // Site index for each observation
  int<lower=1> pop_s[n_s];    // Population index for each observation
  int<lower=1> plot_s[n_s];   // Plot index for each observation
  int<lower=1> species_s[n_s]; // Species index for each observation
  int<lower=0,upper=1> y_s[n_s]; // Survival outcome (0 = dead, 1 = alive)
  vector[n_s] endo_s;         // Endophyte status (1 = positive, 0 = negative)
  vector[n_s] herb_s;         // Herbivory (1 = full fence, 0 =  half fence)
  vector[n_s] temp_s;         // Temperature covariate
}

parameters {
  // Fixed effects (shared across species)
  real b0_s;             // Intercept
  real bendo_s;           // Effect of endophyte
  real bherb_s;           // Effect of herbivory
  real btemp_s;           // Effect of temperature
  real bendotemp_s;       // Interaction: endophyte x temperature
  real bherbtemp_s;       // Interaction: herbivory x temperature
  real bendoherb_s;       // Interaction: endophyte x herbivory

  // Random effect standard deviations (shared across species)
  real<lower=0> site_tau;   // Standard deviation of site random effects
  real<lower=0> pop_tau;    // Standard deviation of population random effects
  real<lower=0> plot_tau;   // Standard deviation of plot random effects

  // Species-specific random effects
  matrix[n_sites, n_species] site_rfx; // Site random effect for each species
  matrix[n_pops, n_species] pop_rfx;   // Population random effect for each species
  matrix[n_plots, n_species] plot_rfx; // Plot random effect for each species
}

transformed parameters {
  // Linear predictor for survival
  vector[n_s] predS;  // Linear predictor for each observation

  for (i in 1:n_s) {
    int sp = species_s[i]; // Species for this observation

    // Linear predictor: fixed effects + species-specific random effects
    predS[i] = b0_s +
               bendo_s * endo_s[i] +
               bherb_s * herb_s[i] +
               btemp_s * temp_s[i] +
               bendotemp_s * temp_s[i] * endo_s[i] +
               bherbtemp_s * temp_s[i] * herb_s[i] +
               bendoherb_s * endo_s[i] * herb_s[i] +
               plot_rfx[plot_s[i], sp] +   // species-specific plot effect
               pop_rfx[pop_s[i], sp] +     // species-specific population effect
               site_rfx[site_s[i], sp];    // species-specific site effect
  }
}

model {
  // Priors on fixed effects
  b0_s ~ normal(0, 3);
  bendo_s ~ normal(0, 3);
  bherb_s ~ normal(0, 3);
  btemp_s ~ normal(0, 3);
  bendotemp_s ~ normal(0, 3);
  bherbtemp_s ~ normal(0, 3);
  bendoherb_s ~ normal(0, 3);

  // Priors on random effect standard deviations
  site_tau ~ inv_gamma(0.1, 0.1);  // weakly informative prior
  pop_tau ~ inv_gamma(0.1, 0.1);
  plot_tau ~ inv_gamma(0.1, 0.1);

  // Priors on species-specific random effects
  for (sp in 1:n_species) {
    for (s in 1:n_sites) site_rfx[s, sp] ~ normal(0, site_tau); // species-specific site effect
    for (p in 1:n_pops)  pop_rfx[p, sp] ~ normal(0, pop_tau);   // species-specific pop effect
    for (pl in 1:n_plots) plot_rfx[pl, sp] ~ normal(0, plot_tau); // species-specific plot effect
  }

  // Likelihood
  y_s ~ bernoulli_logit(predS);  // Bernoulli survival outcome
}

generated quantities {
  // Log-likelihood for model comparison or WAIC/LOO
  vector[n_s] log_lik;
  for (i in 1:n_s) {
    log_lik[i] = bernoulli_logit_lpmf(y_s[i] | predS[i]);
  }
}
