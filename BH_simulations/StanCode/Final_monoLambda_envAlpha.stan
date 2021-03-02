data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Nt[N];
  int Ntp1[N];
  matrix[N,S] SpMatrix;
  vector[N] env;

  int Inclusion_ij[S];
  int Inclusion_eij[S];
}

parameters{
  vector[2] lambdas;    // 1: intercept, 2: slope
  vector[2] alpha_generic_tilde;
  vector[2] alpha_intra_tilde;
  vector[S] alpha_hat_ij_tilde;
  vector[S] alpha_hat_eij_tilde;
}

transformed parameters{
  vector[2] alpha_generic;
  vector[2] alpha_intra;
  vector[S] alpha_hat_ij;
  vector[S] alpha_hat_eij;

  // scale the lambdas and alphas values
  alpha_generic[1] = 0.75 * alpha_generic_tilde[1] - 2;
  alpha_intra[1] = 0.75 * alpha_intra_tilde[1] - 2;
  alpha_generic[2] = alpha_generic_tilde[2] * 0.5;
  alpha_intra[2] = alpha_intra_tilde[2] * 0.5;
  for(s in 1:S){
       alpha_hat_ij[s] = 0.75 * alpha_hat_ij_tilde[s] - 2;
       alpha_hat_eij[s] = alpha_hat_eij_tilde[s] * 0.5;
  }
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] Ntp1_hat;
  vector[N] interaction_effects;
  matrix[N,S] alpha_eij;
  vector[N] lambda_ei;

  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  lambdas ~ normal(0,1);
  alpha_hat_ij_tilde ~ normal(0,1);
  alpha_hat_eij_tilde ~ normal(0,1);

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambdas[1] + lambdas[2]*env[i]);
    for(s in 1:S){
      if(Inclusion_ij[s] == 1){
        if(Inclusion_eij[s] == 1){
          alpha_eij[i,s] = exp(alpha_generic[1] + alpha_hat_ij[s] + (alpha_generic[2] + alpha_hat_eij[s]) * env[i]);
        }else{
          alpha_eij[i,s] = exp(alpha_generic[1] + alpha_hat_ij[s] + alpha_generic[2] * env[i]);
        }
      }else{
        if(Inclusion_eij[s] == 1){
          alpha_eij[i,s] = exp(alpha_generic[1] + (alpha_generic[2] + alpha_hat_eij[s]) * env[i]);
        }else{
          alpha_eij[i,s] = exp(alpha_generic[1] + alpha_generic[2] * env[i]);
        }
      }
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) + exp(alpha_intra[1] + alpha_intra[2] * env[i]) * Nt[i];
    Ntp1_hat[i] = Nt[i] * lambda_ei[i] / (1 + interaction_effects[i]);
    if(Ntp1_hat[i] > 0){
      Ntp1[i] ~ poisson(Ntp1_hat[i]);
    }
  }
}

