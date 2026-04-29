
data {
  int<lower=0> n;
  int<lower=0> n_species;
  int<lower=0> n_endo;
  int<lower=0> n_fence;
  int<lower=0> n_site;
  int<lower=0> n_year;
  int<lower=0> n_plot;
  array[n] int y_tillers;
  array[n] int y_damaged;
  array[n] int species;
  array[n] int endo;
  array[n] int fence;
  array[n] int site;
  array[n] int year;
  array[n] int plot;
}

parameters {
  vector[n_species] beta0;
  vector[n_species] beta_endo;
  vector[n_species] beta_fence;
  vector[n_site] eps_site;
  real<lower=0> sigma_site;
  vector[n_year] eps_year;
  real<lower=0> sigma_year;
  vector[n_plot] eps_plot;
  real<lower=0> sigma_plot;
}

transformed parameters{
  vector[n] eta;
  for(i in 1:n){
    eta[i]=beta0[species[i]] + beta_endo[species[i]]*endo[i] + beta_fence[species[i]]*fence[i] +
    eps_site[site[i]] +
    eps_year[year[i]] +
    eps_plot[plot[i]];
  }
}

model {
  to_vector(beta0)~normal(0,1);
  to_vector(beta_endo)~normal(0,1);
  to_vector(beta_fence)~normal(0,1);
  eps_site ~ normal(0, sigma_site);
  eps_year ~ normal(0, sigma_year);
  eps_plot ~ normal(0, sigma_plot);
  sigma_site~normal(0,1);
  sigma_year~normal(0,1);
  sigma_plot~normal(0,1);
  y_damaged~binomial_logit(y_tillers,eta);
}

generated quantities{
  array[n] int y_rep;
  for(i in 1:n){
    y_rep[i]=binomial_rng(y_tillers[i],inv_logit(eta[i]));
  }
}