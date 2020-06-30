// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  matrix[N,S] SpMatrix;
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values
  vector[N] phos;	// Phosphorous measurements for each site
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  vector[2] lambdas;
  vector[2] alphas;

  matrix[2,S] alpha_sp_tilde;
  matrix<lower = 0>[2,S] local_shrinkage;

  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Declare all the transformed parameters to be used
  real c2;
  real tau;
  matrix[2,S] alpha_sp;
  matrix<lower = 0>[2,S] local_shrinkage_tilde;
  matrix[N,S] alpha_terms;
  vector[N] interaction_effects;
  vector[N] lambda;

  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)
  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  for(i in 1:2){
    local_shrinkage_tilde[i,] = sqrt( c2 * square(local_shrinkage[i,]) ./ (c2 + square(tau) * square(local_shrinkage[i,])) );
    alpha_sp[i,] = tau * local_shrinkage_tilde[i,] .* alpha_sp_tilde[i,];	// alpha_sp ~ normal(0, tau*local_shrinkage_sp_tilde)
  }
  
  // Now calculate the matrix of alpha values based on the phosphorous values for each site and the vector of lambda values
  lambda = exp(lambdas[1] + lambdas[2] * phos);
  for(i in 1:N){
    alpha_terms[i,] = exp(alphas[1] + alpha_sp[1,] + (alphas[2] + alpha_sp[2,]) * phos[i]);
    interaction_effects[i] = sum(alpha_terms[i,] .* SpMatrix[i,]);
  }
}

model{
  // Make the vector of expected fecundities
  vector[N] F_hat;

  // set regular priors
  lambdas ~ normal(0, 100);
  alphas ~ normal(0, 100);

  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  for(i in 1:2){
    alpha_sp_tilde[1,] ~ normal(0,1);
    local_shrinkage[1,] ~ cauchy(0,1);
  }
  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  F_hat = lambda ./ (1 + interaction_effects);
  Fecundity ~ poisson(F_hat);
}

