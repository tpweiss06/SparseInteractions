// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;  // Number of plots
  int<lower = 1> S;  // Number of species
  int Fecundity[N];  // Fecundity of the focal species in each plot
  int reserve[N];    // Indicator variable for the reserve each plot is located in
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  vector[N] env;  // Environmental values for each plot
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  int Inclusion_ij[2,S];  // Boolean indicator variables to identify the species x reserve intercept parameters identified for inclusion in the final model
  int Inclusion_eij[2,S]; // Boolean indicator variables to identify the species x reserve slope parameters identified for inclusion in the final model
}

parameters{
  matrix[2,2] lambdas;
  vector[2] alpha_generic_tilde;
  vector[2] alpha_intra_tilde;
  matrix[2,S] alpha_hat_ij;
  matrix[2,S] alpha_hat_eij;
}

transformed parameters{
  vector[2] alpha_generic;
  vector[2] alpha_intra;

  // scale the alpha values
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
    alpha_hat_ij[i,] ~ normal(0,1);
    alpha_hat_eij[i,] ~ normal(0,1);
  }

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambdas[reserve[i],1] + lambdas[reserve[i],2] * env[i]);
    for(s in 1:S){
      alpha_eij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + Inclusion_ij[reserve[i],s] * alpha_hat_ij[reserve[i],s] + ((1-Intra[s]) * alpha_generic[2] + Inclusion_eij[reserve[i],s] * alpha_hat_eij[reserve[i],s] + Intra[s] * alpha_intra[2]) * env[i]);
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
    F_hat[i] = lambda_ei[i] / (1 + interaction_effects[i]);
  }
  Fecundity ~ poisson(F_hat);
}
