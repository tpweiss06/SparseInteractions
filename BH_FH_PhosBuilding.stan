// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  matrix[N,S] SpMatrix;
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha species' specific phosphorous slope values
  real slab_df;		// effective degrees of freedom for significant species' specific phosphorous slope values
  vector[N] phos;	// Phosphorous measurements for each site

  // A vector of 1's and 0's indicating whether previous model fits support 
  // 	inclusion of the alpha species's specific intercept terms as well
  //	as the number of species' specific terms to be used.
  int InclusionIntercept[S];	
  int<lower = 0> S_intercept;
  real intercept_mean[S_intercept];
  real intercept_sd[S_intercept];
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  vector[2] lambdas;
  vector[2] alphas;
  vector[S_intercept] alpha_sp_intercept;

  vector[S] alpha_sp_phos_tilde;
  vector<lower = 0>[S] local_shrinkage;

  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Declare all the transformed parameters to be used
  real c2;
  real tau;
  vector[S] alpha_sp_phos;
  vector<lower = 0>[S] local_shrinkage_tilde;
  
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)
  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  local_shrinkage_tilde = sqrt( c2 * square(local_shrinkage) ./ (c2 + square(tau) * square(local_shrinkage)) );
  alpha_sp_phos = tau * local_shrinkage_tilde .* alpha_sp_phos_tilde;	// alpha_sp ~ normal(0, tau*local_shrinkage_sp_tilde)  
}

model{
  // Make the vector of expected fecundities, matrix of alpha terms, vector of summed effects, and the lambda vector
  vector[N] F_hat;
  matrix[N,S] alpha_terms;
  vector[N] interaction_effects;
  vector[N] lambda;
  int index;

  // set regular priors
  lambdas ~ normal(0, 100);
  alphas ~ normal(0, 100);
  
  // set alpha_sp_intercept informative priors from previous fit
  for(i in 1:S_intercept){
    alpha_sp_intercept[i] ~ normal(intercept_mean[i], intercept_sd[i]);
  }

  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  alpha_sp_phos_tilde ~ normal(0,1);
  local_shrinkage ~ cauchy(0,1);
  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  // First calculate the matrix of alpha values based on the phosphorous values for each site and the vector of lambda values
  lambda = exp(lambdas[1] + lambdas[2] * phos);
  for(i in 1:N){
    index = 1;
    for(s in 1:S){
      if(InclusionIntercept[s] == 1){
        alpha_terms[i,s] = exp(alphas[1] + alpha_sp_intercept[index] + (alphas[2] + alpha_sp_phos[s]) * phos[i]);
        index = index + 1;
      }else{
        alpha_terms[i,s] = exp(alphas[1] + (alphas[2] + alpha_sp_phos[s]) * phos[i]);
      }
    }
    interaction_effects[i] = sum(alpha_terms[i,] .* SpMatrix[i,]);
  }
  // Now use those values to calculate the vector of expected fecundities
  F_hat = lambda ./ (1 + interaction_effects);
  Fecundity ~ poisson(F_hat);
}

