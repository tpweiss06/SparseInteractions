# This script will fit the model with a monotonic lambda and constant alphas
#    (with respect to the environment) to the corresponding data, using a series
#    of dataset sizes (10, 20, 50, 100, 200), and then examine the accuracy of
#    parameter estimates and posterior predictive checks on 300 out of sample
#    data points.

rm(list = ls())
library(rstan)
library(here)
library(HDInterval)
library(RColorBrewer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the current sample size and associated prefix for all graph and result
#    file names
N <- 50
max_N <- 200
FilePrefix <- paste("N", N, "_", sep = "")

# Now assign the focal species and the file paths for the stan models
Focal <- 8
PrelimStanPath <- here("BH_simulations/Box/StanCode/Prelim_optLambda_constAlpha.stan")
FinalStanPath <- here("BH_simulations/Box/StanCode/Final_optLambda_constAlpha.stan")

# Load in the appropriate data
FullSim <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/simulation_perturb_opt_const.csv"))
TrueVals <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/parameters_perturb_opt_const.csv"))
TrueAlphas <- exp(TrueVals$alpha.5)
TrueLambdaOpt <- TrueVals$z.env[Focal]
TrueLambdaMax <- TrueVals$lambda.max[Focal]
TrueLambdaWidth <- TrueVals$sigma.env[Focal]

# assign some universal values to be used across model fits and graphs
S <- 15
Intra <- rep(0, S)
Intra[Focal] <- 1
tau0 <- 1
slab_df <- 4
slab_scale <- sqrt(2)
AllSpecies <- 1:S
OtherSpecies <- AllSpecies[-Focal]

# Set initial values to avoid initial problems with the random number generator
ChainInitials <- list(lambda_max = 1, lambda_opt = 0, lambda_width = 0.5, 
                      alpha_generic_tilde = 0, alpha_hat_ij_tilde = rep(0, S),
                      local_shrinkage_ij = rep(5, S), c2_tilde = 1.25, tau_tilde = 15, 
                      alpha_intra_tilde = 0)
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# save the values for the posterior predictive checks
ppc_data <- subset(FullSim, (species == Focal) & (run > max_N) & (time == 0) & (thinned == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
        SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) &
                                           (time == 0) & (thinned == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) &
                           (time == 1) & (thinned == 0))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

# Create the data vectors to be passed to rstan for subsequent model fits
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij")

# Set the local values to pass to rstan
FullData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0) & (thinned == 0))
ThinData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0) & (thinned == 1))
Nt <- c(FullData$pop, ThinData$pop)
env <- c(FullData$run.env, ThinData$run.env)
SpMatrix <- matrix(data = NA, nrow = 2*N, ncol = S)
for(s in 1:S){
        SpMatrix[1:N,s] <- subset(FullSim, (species == s) & (run <= N) & (time == 0) & (thinned == 0))$pop
        SpMatrix[(N+1):(2*N),s] <- subset(FullSim, (species == s) & (run <= N) & (time == 0) & (thinned == 1))$pop
}
Ntp1 <- c(subset(FullSim, (species == Focal) & (run <= N) & (time == 1) & (thinned == 0))$pop,
          subset(FullSim, (species == Focal) & (run <= N) & (time == 1) & (thinned == 1))$pop)


# Now run the preliminary fit of the model to assess parameter shrinkage
N <- length(Nt)
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, init = InitVals, control = list(max_treedepth = 15, adapt_delta = 0.99))
PrelimPosteriors <- extract(PrelimFit)
FitFileName <- paste(here("BH_simulations/Box/StanFits/optLambda_constAlpha/"), FilePrefix, "PrelimFit.rdata", sep = "")
save(PrelimFit, PrelimPosteriors, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
pairs(PrelimFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
traceplot(PrelimFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
traceplot(PrelimFit, pars = "alpha_hat_ij")
traceplot(PrelimFit, pars = c("local_shrinkage_ij", "c2", "tau"))
acf(PrelimPosteriors$alpha_generic)
acf(PrelimPosteriors$alpha_intra)
acf(PrelimPosteriors$lambda_opt)
acf(PrelimPosteriors$lambda_max)
acf(PrelimPosteriors$lambda_width)

########### Posterior Predictive Check---------
PostLength <- length(PrelimPosteriors$alpha_generic)
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
        alpha_ij[,s] <- exp(PrelimPosteriors$alpha_generic + PrelimPosteriors$alpha_hat_ij[,s])
}
alpha_intra <- exp(PrelimPosteriors$alpha_intra)
# calculate the posterior distributions of lambda_ei
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        lambda_ei[,i] <- PrelimPosteriors$lambda_max * 
                exp(-1 * ( (PrelimPosteriors$lambda_opt - env_ppc[i]) / (2 * PrelimPosteriors$lambda_width) )^2)
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
        for(j in 1:N_ppc){
                SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,]) + alpha_intra[i]*Nt_ppc[j]
                Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
                Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
                Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
        }
}

# Calculate preliminary fit results for the ppc
PrelimPredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
PrelimCIwidths <- rep(NA, N_ppc)
PrelimDevVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
     PrelimPredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
     PrelimPredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
     PrelimCIwidths[i] <- PrelimPredVals[3,i] - PrelimPredVals[2,i]
     PrelimDevVals[1,i] <- mean(Growth_dev[,i], na.rm = TRUE)
     PrelimDevVals[2:3,i] <- HDInterval::hdi(Growth_dev[,i])
}
PrelimCIwidth <- mean(PrelimCIwidths)

# Calculate the deviation in the lambda estimates
PrelimLambdaDevs <- matrix(data = NA, nrow = PostLength, ncol = 3)
PrelimLambdaDevs[,1] <- PrelimPosteriors$lambda_opt - TrueLambdaOpt
PrelimLambdaDevs[,2] <- PrelimPosteriors$lambda_max - TrueLambdaMax
PrelimLambdaDevs[,3] <- PrelimPosteriors$lambda_width - TrueLambdaWidth

# Now calculate the alpha estimates
PrelimAlphaEsts <- matrix(data = NA, nrow = 3, ncol = S+1)
for(s in 1:S){
     AlphaDev <- alpha_ij[,s] - as.numeric(TrueAlphas[OtherSpecies[s]])
     PrelimAlphaEsts[1,s] <- mean(AlphaDev)
     PrelimAlphaEsts[2:3,s] <- HDInterval::hdi(AlphaDev)
}
IntraDev <- alpha_intra - TrueAlphas[Focal]
PrelimAlphaEsts[1,S+1] <- mean(IntraDev)
PrelimAlphaEsts[2:3,S+1] <- HDInterval::hdi(IntraDev)

## Final model ---------
# Determine the parameters that should be included and run the final model
plot(PrelimFit, pars = "alpha_hat_ij")

Inclusion_ij <- rep(0, S)
IntLevel <- 0.5
for(s in 1:S){
     Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
     if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
          Inclusion_ij[s] <- 1
     }
     if(s == Focal){
             Inclusion_ij[s] <- 0
     }
}
Inclusion_ij

# Reset initial conditions to allow faster fitting
ChainInitials <- list(lambda_max = mean(PrelimPosteriors$lambda_max), 
                      lambda_opt = mean(PrelimPosteriors$lambda_opt), 
                      lambda_width = mean(PrelimPosteriors$lambda_width), 
                      alpha_generic_tilde = mean(PrelimPosteriors$alpha_generic_tilde), 
                      alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), 
                      local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
                      c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
                      alpha_intra_tilde = mean(PrelimPosteriors$alpha_intra_tilde))
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals)
FinalPosteriors <- extract(FinalFit)
FitFileName <- paste(here("BH_simulations/Box/StanFits/optLambda_constAlpha/"), FilePrefix, "FinalFit.rdata", sep = "")
save(FinalFit, FinalPosteriors, Inclusion_ij, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
 pairs(FinalFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
 hist(summary(FinalFit)$summary[,"Rhat"])
 hist(summary(FinalFit)$summary[,"n_eff"])
 traceplot(FinalFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
 which(Inclusion_ij == 1)
 traceplot(FinalFit, pars = "alpha_hat_ij")
 
 # Double check the autocorrelation in a few potentially suspect traceplots
 acf(FinalPosteriors$alpha_intra)
 acf(FinalPosteriors$alpha_generic)
 acf(FinalPosteriors$lambda_opt)
 acf(FinalPosteriors$lambda_max)
 acf(FinalPosteriors$lambda_width)

########### Posterior Predictive Check
PostLength <- length(FinalPosteriors$alpha_generic)
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
     if(Inclusion_ij[s] == 1){
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_generic + FinalPosteriors$alpha_hat_ij[,s])
     }else{
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_generic)
     }
     
}
alpha_intra <- exp(FinalPosteriors$alpha_intra)

# calculate the posterior distributions of lambda_ei
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
     lambda_ei[,i] <- FinalPosteriors$lambda_max * 
             exp(-1 * ( (FinalPosteriors$lambda_opt - env_ppc[i]) / (2 * FinalPosteriors$lambda_width) )^2)
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,]) + alpha_intra[i] * Nt_ppc[j]
          Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
          Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
     }
}

# Calculate final fit results for the ppc
FinalPredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
FinalCIwidths <- rep(NA, N_ppc)
FinalDevVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
     FinalPredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
     FinalPredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
     FinalCIwidths[i] <- FinalPredVals[3,i] - FinalPredVals[2,i]
     FinalDevVals[1,i] <- mean(Growth_dev[,i], na.rm = TRUE)
     FinalDevVals[2:3,i] <- HDInterval::hdi(Growth_dev[,i])
}
FinalCIwidth <- mean(FinalCIwidths)

# Calculate the deviation in the lambda estimates
FinalLambdaDevs <- matrix(data = NA, nrow = PostLength, ncol = 3)
FinalLambdaDevs[,1] <- FinalPosteriors$lambda_opt - TrueLambdaOpt
FinalLambdaDevs[,2] <- FinalPosteriors$lambda_max - TrueLambdaMax
FinalLambdaDevs[,3] <- FinalPosteriors$lambda_width - TrueLambdaWidth

# Now calculate the alpha estimates
FinalAlphaEsts <- matrix(data = NA, nrow = 3, ncol = S+1)
for(s in 1:S){
     AlphaDev <- alpha_ij[,s] - as.numeric(TrueAlphas[OtherSpecies[s]])
     FinalAlphaEsts[1,s] <- median(AlphaDev)
     FinalAlphaEsts[2:3,s] <- HDInterval::hdi(AlphaDev)
}
IntraDev <- alpha_intra - TrueAlphas[Focal]
FinalAlphaEsts[1,S+1] <- mean(IntraDev)
FinalAlphaEsts[2:3,S+1] <- HDInterval::hdi(IntraDev)


