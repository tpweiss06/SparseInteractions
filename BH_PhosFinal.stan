// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  matrix[N,S] SpMatrix;
  vector[N] phos;	// Phosphorous measurements for each site

  // A vector of 1's and 0's indicating whether previous model fits support 
  // 	inclusion of the alpha species's specific intercept terms as well
  //	as the number of species' specific terms to be used.
  int InclusionIntercept[S];	
  int<lower = 0> S_intercept;
  real intercept_mean[S_intercept];
  real intercept_sd[S_intercept];
  int InclusionSlope[S];
  int<lower = 0> S_slope;
}

parameters{
  vector[2] lambdas;
  vector[2] alphas;
  vector[S_intercept] alpha_sp_intercept;
  vector[S_slope] alpha_sp_slope;
}

model{
  // Make the vector of expected fecundities, matrix of alpha terms, vector of summed effects, and the lambda vector
  vector[N] F_hat;
  matrix[N,S] alpha_terms;
  vector[N] interaction_effects;
  vector[N] lambda;
  int index_intercept;
  int index_slope;

  // set regular priors
  lambdas ~ normal(0, 100);
  alphas ~ normal(0, 100);
  //alpha_sp_intercept ~ normal(0, 100);
  alpha_sp_slope ~ normal(0, 100);
  
  // set alpha_sp_intercept informative priors from previous fit
  for(i in 1:S_intercept){
    alpha_sp_intercept[i] ~ normal(intercept_mean[i], intercept_sd[i]);
  }

  // implement the biological model
  // First calculate the matrix of alpha values based on the phosphorous values for each site and the vector of lambda values
  lambda = exp(lambdas[1] + lambdas[2] * phos);
  for(i in 1:N){
    index_intercept = 1;
    index_slope = 1;
    for(s in 1:S){
      if(InclusionIntercept[s] == 1){
        if(InclusionSlope[s] == 1){
          alpha_terms[i,s] = exp(alphas[1] + alpha_sp_intercept[index_intercept] + (alphas[2] + alpha_sp_slope[index_slope]) * phos[i]);
          index_intercept = index_intercept + 1;
          index_slope = index_slope + 1;
        }else{
          alpha_terms[i,s] = exp(alphas[1] + alpha_sp_intercept[index_intercept] + alphas[2] * phos[i]);
          index_intercept = index_intercept + 1;
        }
      }else{
        if(InclusionSlope[s] == 1){
          alpha_terms[i,s] = exp(alphas[1] + (alphas[2] + alpha_sp_slope[index_slope]) * phos[i]);
          index_slope = index_slope + 1;
        }else{
          alpha_terms[i,s] = exp(alphas[1] + alphas[2] * phos[i]);
        }
      }
    }
    interaction_effects[i] = sum(alpha_terms[i,] .* SpMatrix[i,]);
  }
  // Now use those values to calculate the vector of expected fecundities
  F_hat = lambda ./ (1 + interaction_effects);
  Fecundity ~ poisson(F_hat);
}

