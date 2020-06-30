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
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  real log_lambda;
  real alpha_hat;
  vector[S] alpha_sp_tilde;
  vector<lower = 0>[S] local_shrinkage;
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  vector[S] alpha_sp;
  vector[S] local_shrinkage_tilde;
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)
  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  local_shrinkage_tilde = sqrt( c2 * square(local_shrinkage) ./ (c2 + square(tau) * square(local_shrinkage)) );
  alpha_sp = tau * local_shrinkage_tilde .* alpha_sp_tilde;	// alpha_sp ~ normal(0, tau*local_shrinkage_tilde)
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  vector[S] alpha_terms;
  

  // set regular priors
  alpha_hat ~ normal(0, 100);
  log_lambda ~ normal(0, 100);

  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  alpha_sp_tilde ~ normal(0,1);
  local_shrinkage ~ cauchy(0,1);
  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  alpha_terms = exp(alpha_hat + alpha_sp);
  interaction_effects = SpMatrix * alpha_terms;
  for(i in 1:N){ 
    F_hat[i] = exp(log_lambda)/(1 + interaction_effects[i]);
  }
  Fecundity ~ poisson(F_hat);
}

