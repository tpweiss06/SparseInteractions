// This script will fit the Beverton-Holt model to data with no environmental covariates or random effects

data{
  int<lower = 1> N;
  int<lower = 1> S_total;
  int<lower = 1> S_eff;
  int Fecundity[N];
  int Inclusion[S_total];
  matrix[N,S_total] SpMatrix;
}

parameters{
  real log_lambda;
  real alpha_hat;
  vector[S_eff] alpha_sp;
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  vector[S_total] alpha_terms;
  int index;

  // set priors
  alpha_hat ~ normal(0, 1);
  alpha_sp ~ normal(0, 1);
  log_lambda ~ normal(0, 10);

  // implement the biological model
  index = 1;
  for(s in 1:S_total){
    if(Inclusion[s] == 1){
      alpha_terms[s] = exp(alpha_hat + alpha_sp[index]);
      index = index + 1;
    }else{
      alpha_terms[s] = exp(alpha_hat);
    }
  }
  interaction_effects = SpMatrix * alpha_terms;
  for(i in 1:N){ 
    F_hat[i] = exp(log_lambda)/(1 + interaction_effects[i]);
  }
  Fecundity ~ poisson(F_hat);
}








