data {
  // Data for growth sub-model (g)
  int<lower=1> nSpp;         // Number of species in the model
  int<lower=1> nsite_year;   // Number of site_years
  int<lower=1> nPop;         // Number of source populations
  int<lower=1> N;            // Number of data points for the growth model
  int<lower=1> nPlot;        // Number of plots

  // Indices for categorical variables
  int<lower=1> Spp[N];       // Species index for each data point
  int<lower=1> site_year[N]; // site_year index for each data point
  int<lower=1> plot[N];      // Plot index for each data point
  int<lower=1> pop[N];       // Population index for each data point

  // Response variable
  vector[N] y;               // Observed growth (continuous response variable)

  // Binary covariates
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N];  // Herbivory status (1 = present, 0 = absent)
}

parameters {
  // Fixed effects coefficients
  vector[nSpp] b0;           // Intercept term for each species
  vector[nSpp] bendo;        // Coefficient for endophyte status for each species
  vector[nSpp] bherb;        // Coefficient for herbivory for each species
  vector[nSpp] bendoherb;    // Coefficient for the interaction between endophyte and herbivory
  vector[nSpp] bherb2;       // Coefficient for the quadratic effect of herbivory
  vector[nSpp] bendoherb2;   // Coefficient for the quadratic interaction between endophyte and herbivory

  // Random effects variances
  real<lower=0> plot_tau;         // Variance for plot-level random effects
  vector[nPlot] plot_rfx;         // Random effects for each plot

  real<lower=0> pop_tau;          // Variance for population-level random effects
  vector[nPop] pop_rfx;           // Random effects for each population

  real<lower=0> site_year_tau;    // Variance for site_year-level random effects (shared across species)
  matrix[nSpp, nsite_year] site_year_rfx;  // Random effects for each species at each site_year

  real<lower=0> sigma;            // Residual standard deviation
}

transformed parameters {
  vector[N] predG;  // Predicted growth values

  // Loop over all growth data points and calculate the predicted growth value
  for (igrow in 1:N) {
    predG[igrow] = b0[Spp[igrow]] + 
                   // Main effects
                   bendo[Spp[igrow]] * endo[igrow] +
                   bherb[Spp[igrow]] * herb[igrow] +
                   // Two-way interaction
                   bendoherb[Spp[igrow]] * endo[igrow] * herb[igrow] +
                   // Quadratic herbivory effects
                   bherb2[Spp[igrow]] * square(herb[igrow]) +  
                   bendoherb2[Spp[igrow]] * endo[igrow] * square(herb[igrow]) + 
                   // Random effects
                   plot_rfx[plot[igrow]] +     // Plot-level random effect
                   pop_rfx[pop[igrow]] +       // Population-level random effect
                   site_year_rfx[Spp[igrow], site_year[igrow]];  // site_year-level random effect for each species at each site_year
  }
}

model {
  // Priors for fixed effects (Normally distributed with mean 0 and SD 5)
  b0 ~ normal(0, 5);    
  bendo ~ normal(0, 5);   
  bherb ~ normal(0, 5);  
  bendoherb ~ normal(0, 5);  
  bherb2 ~ normal(0, 5);  
  bendoherb2 ~ normal(0, 5);
  sigma ~ normal(0, 5);  // Prior for residual standard deviation

  // Priors for random effects variances
  plot_tau ~ normal(0, 1);
  for (i in 1:nPlot) {
    plot_rfx[i] ~ normal(0, plot_tau);  // Random effects for each plot
  }
  
  pop_tau ~ normal(0, 1);
  for (i in 1:nPop) {
    pop_rfx[i] ~ normal(0, pop_tau);   // Random effects for each population
  }

  site_year_tau ~ normal(0, 1);  // Shared variance for all species at site_year level
  for (i in 1:nSpp) {
    for (j in 1:nsite_year) {
      site_year_rfx[i, j] ~ normal(0, site_year_tau);  // Random effects for each species at each site_year
    }
  }

  // Likelihood for growth (Normal distribution)
  y ~ normal(predG, sigma);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for each growth data point

  // Loop over all growth data points and calculate the log-likelihood
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);  // Log-likelihood for Normal distribution
  }
}
