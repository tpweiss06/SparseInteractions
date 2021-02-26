data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Nt[N];
  int Ntp1[N];
  matrix[N,S] SpMatrix;
  vector[N] env;

  int Inclusion_ij[S];
  int Inclusion_eij[S];
  
  // Include the data for the posterior predictive check
  int<lower = 1> N_ppc;
  int<lower = 0> Nt_ppc[N_ppc];
  matrix[N_ppc,S] SpMatrix_ppc;
  vector[N_ppc] env_ppc;
}

parameters{
  vector[2] lambdas;    // 1: intercept, 2: slope
  vector[2] alphas_tilde;
  vector[S] alpha_hat_ij;
  vector[S] alpha_hat_eij;
}

transformed parameters{
  vector[2] alphas;

  // scale the lambdas and alphas values
  alphas[1] = 1.75 * alphas_tilde[1] - 7;
  alphas[2] = alphas_tilde[2] * 0.5;
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
  alphas_tilde ~ normal(0,1);
  lambdas ~ normal(0,1);
  alpha_hat_ij ~ normal(0,1);
  alpha_hat_eij ~ normal(0,1);

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambdas[1] + lambdas[2]*env[i]);
    for(s in 1:S){
      if(Inclusion_ij[s] == 1){
        if(Inclusion_eij[s] == 1){
          alpha_eij[i,s] = exp(alphas[1] + alpha_hat_ij[s] + (alphas[2] + alpha_hat_eij[s]) * env[i]);
        }else{
          alpha_eij[i,s] = exp(alphas[1] + alpha_hat_ij[s] + alphas[2] * env[i]);
        }
      }else{
        if(Inclusion_eij[s] == 1){
          alpha_eij[i,s] = exp(alphas[1] + (alphas[2] + alpha_hat_eij[s]) * env[i]);
        }else{
          alpha_eij[i,s] = exp(alphas[1] + alphas[2] * env[i]);
        }
      }
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
    Ntp1_hat[i] = Nt[i] * lambda_ei[i] / (1 + interaction_effects[i]);
    if(Ntp1_hat[i] > 0){
      Ntp1[i] ~ poisson(Ntp1_hat[i]);
    }
  }
}

