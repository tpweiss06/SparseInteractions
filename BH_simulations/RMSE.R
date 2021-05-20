# calculate the root mean squared error (RMSE) for our growth predictions and
#       lambda, intra, and generic terms

setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Get the data needed for the ppc calculations
Focal <- 1
FullSim <- read.csv("Simulations/simulation_perturb2.csv")
TrueVals <- read.csv("Simulations/parameters_perturb2.csv")
max_N <- 200
S <- 15
ppc_data <- subset(FullSim, (species == Focal) & (run > max_N) & (time == 0) & (thinned == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
     SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) & (time == 0) & (thinned == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

# load in the final fit
load("StanFits/monoLambda_envAlpha/N10_FinalFit.rdata")
Post10 <- Posteriors
load("StanFits/monoLambda_envAlpha/N50_FinalFit.rdata")
Post50 <- Posteriors
load("StanFits/monoLambda_envAlpha/N200_FinalFit.rdata")
Post200 <- Posteriors

# Calculate the posterior growth predictions
PostLength <- length(Post10$alpha_generic[,1])
# calculate the posterior distributions of the interaction coefficients and lambdas
alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S, 3))
lambda_ei <- array(NA, dim = c(PostLength, N_ppc, 3))
for(i in 1:N_ppc){
     lambda_ei[,i,1] <- exp(Post10$lambdas[,1] + Post10$lambdas[,2]*env_ppc[i])
     lambda_ei[,i,2] <- exp(Post50$lambdas[,1] + Post50$lambdas[,2]*env_ppc[i])
     lambda_ei[,i,3] <- exp(Post200$lambdas[,1] + Post200$lambdas[,2]*env_ppc[i])
     for(s in 1:S){
          if(s == Focal){
               alpha_eij[,i,s,1] <- exp(Post10$alpha_intra[,1] + Post10$alpha_intra[,2] * env_ppc[i])
               alpha_eij[,i,s,2] <- exp(Post50$alpha_intra[,1] + Post50$alpha_intra[,2] * env_ppc[i])
               alpha_eij[,i,s,3] <- exp(Post200$alpha_intra[,1] + Post200$alpha_intra[,2] * env_ppc[i])
          }else{
               alpha_eij[,i,s,1] <- exp(Post10$alpha_generic[1] + Inclusion_ij[s] * Post10$alpha_hat_ij[,s] +
                                           (Post10$alpha_generic[,2] + Inclusion_eij[s] * Post10$alpha_hat_eij[,s]) * env_ppc[i])
               alpha_eij[,i,s,2] <- exp(Post50$alpha_generic[1] + Inclusion_ij[s] * Post50$alpha_hat_ij[,s] +
                                             (Post50$alpha_generic[,2] + Inclusion_eij[s] * Post50$alpha_hat_eij[,s]) * env_ppc[i])
               alpha_eij[,i,s,3] <- exp(Post200$alpha_generic[1] + Inclusion_ij[s] * Post200$alpha_hat_ij[,s] +
                                             (Post200$alpha_generic[,2] + Inclusion_eij[s] * Post200$alpha_hat_eij[,s]) * env_ppc[i])
          }
     }
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- array(data = NA, dim = c(PostLength, N_ppc, 3))
for(n in 1:3){
     for(i in 1:PostLength){
          for(j in 1:N_ppc){
               SigmaTerm <- sum(alpha_eij[i,j,,n] * SpMatrix_ppc[j,])
               Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j,n] / (1 + SigmaTerm)
               Growth_pred[i,j,n] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          }
     }
}

# load in the true values from the graph stuff data
TrueGenericInterceptVals <- rep(NA, 3)
TrueGenericSlopeVals <- rep(NA, 3)
load("StanFits/monoLambda_envAlpha/N10_GraphStuff.rdata")
TrueGenericInterceptVals[1] <- TrueGenericSlope
TrueGenericSlopeVals[1] <- TrueGenericSlope
rm(LambdaEsts, PredVals, AlphaEsts, Inclusion_eij, Inclusion_ij)
load("StanFits/monoLambda_envAlpha/N50_GraphStuff.rdata")
TrueGenericInterceptVals[2] <- TrueGenericSlope
TrueGenericSlopeVals[2] <- TrueGenericSlope
rm(LambdaEsts, PredVals, AlphaEsts, Inclusion_eij, Inclusion_ij)
load("StanFits/monoLambda_envAlpha/N200_GraphStuff.rdata")
TrueGenericInterceptVals[3] <- TrueGenericSlope
TrueGenericSlopeVals[3] <- TrueGenericSlope
rm(LambdaEsts, PredVals, AlphaEsts, Inclusion_eij, Inclusion_ij)


# Now calculate the posterior RMSE for growth predictions
GrowthRMSE <- matrix(NA, nrow = PostLength, ncol = 3)
for(n in 1:3){
     for(i in 1:PostLength){
          deviation_sq <- (Growth_pred[i,,n] - Growth_ppc)^2
          GrowthRMSE[i,n] <- sqrt(sum(deviation_sq) / N_ppc)
     }
}

colMeans(GrowthRMSE)
HDInterval::hdi(GrowthRMSE[,3])   
# N = 10:  0.495 (0.353, 0.665)
# N = 50:  0.315 (0.211, 0.520)
# N = 200: 0.227 (0.195, 0.288)





