// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  int reserve[N];
  matrix[N,S] SpMatrix;
  vector[N] env;
  int<lower = 0> Intra[S];
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  matrix[2,2] lambdas;
  vector[2] alpha_generic_tilde;
  vector[2] alpha_intra_tilde;
  matrix[2,S] alpha_hat_ij_tilde;
  matrix[2,S] alpha_hat_eij_tilde;
  matrix<lower = 0>[2,S] local_shrinkage_ij;
  matrix<lower = 0>[2,S] local_shrinkage_eij;
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  matrix[2,S] alpha_hat_ij;
  matrix[2,S] local_shrinkage_ij_tilde;
  matrix[2,S] alpha_hat_eij;
  matrix[2,S] local_shrinkage_eij_tilde;
  vector[2] alpha_generic;
  vector[2] alpha_intra;

  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  for(i in 1:2){
    for(s in 1:S){
      local_shrinkage_ij_tilde[i,s] = sqrt( c2 * square(local_shrinkage_ij[i,s]) / (c2 + square(tau) * square(local_shrinkage_ij[i,s])) );
      alpha_hat_ij[i,s] = tau * local_shrinkage_ij_tilde[i,s] * alpha_hat_ij_tilde[i,s];

      local_shrinkage_eij_tilde[i,s] = sqrt( c2 * square(local_shrinkage_eij[i,s]) / (c2 + square(tau) * square(local_shrinkage_eij[i,s])) );
      alpha_hat_eij[i,s] = tau * local_shrinkage_eij_tilde[i,s] * alpha_hat_eij_tilde[i,s];
    }
  }

  // scale the lambdas and alphas values
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;
  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;
  alpha_generic[2] = 0.5 * alpha_generic_tilde[2];
  alpha_intra[2] = 0.5 * alpha_intra_tilde[2];
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  matrix[N,S] alpha_eij;
  vector[N] lambda_ei;

  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  for(i in 1:2){
    lambdas[i,] ~ normal(0, 1);
  }

  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  for(i in 1:2){
    alpha_hat_ij_tilde[i,] ~ normal(0,1);
    local_shrinkage_ij[i,] ~ cauchy(0,1);

    alpha_hat_eij_tilde[i,] ~ normal(0,1);
    local_shrinkage_eij[i,] ~ normal(0,1);
  }
  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambdas[reserve[i],1] + lambdas[reserve[i],2] * env[i]);
    for(s in 1:S){
        alpha_eij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) * alpha_hat_ij[reserve[i],s] + ((1-Intra[s]) * alpha_generic[2] + (1-Intra[s]) * alpha_hat_eij[reserve[i],s] + Intra[s] * alpha_intra[2]) * env[i]);
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
    
    F_hat[i] = lambda_ei[i] / (1 + interaction_effects[i]);
  }
  Fecundity ~ poisson(F_hat);
}
