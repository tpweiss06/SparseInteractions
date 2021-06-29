data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Nt[N];
  int Ntp1[N];
  matrix[N,S] SpMatrix;
  vector[N] env;
  int<lower = 0> Intra[S];
  int Inclusion_ij[S];
}

parameters{
  real lambda_opt;
  real<lower = 0> lambda_max;
  real<lower = 0> lambda_width;
  real alpha_generic_tilde;
  real alpha_intra_tilde;
  vector[S] alpha_hat_ij;
}

transformed parameters{
  real alpha_generic;
  real alpha_intra;

  // scale the lambdas and alpha values
  alpha_generic = 0.75 * alpha_generic_tilde - 2;
  alpha_intra = 0.75 * alpha_intra_tilde - 2;
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] Ntp1_hat;
  vector[N] interaction_effects;
  row_vector[S] alpha_ij;
  vector[N] lambda_ei;

  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_hat_ij ~ normal(0,1);
  lambda_opt ~ normal(0, 1);
  lambda_max ~ normal(0, 7.5);
  lambda_width ~ normal(0, 1);

  

  // implement the biological model
  for(s in 1:S){
        alpha_ij[s] = exp((1-Intra[s]) * alpha_generic + Intra[s] * alpha_intra + Inclusion_ij[s] * alpha_hat_ij[s]);
  }
  for(i in 1:N){
    lambda_ei[i] = lambda_max * exp(-1*((lambda_opt - env[i])/(2*lambda_width))^2);
    interaction_effects[i] = sum(alpha_ij .* SpMatrix[i,]);
    Ntp1_hat[i] = Nt[i] * lambda_ei[i] / (1 + interaction_effects[i]);
    if(Ntp1_hat[i] > 0){
      Ntp1[i] ~ poisson(Ntp1_hat[i]);
    }
  }
}
