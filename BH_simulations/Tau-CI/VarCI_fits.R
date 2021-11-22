# This script will fit the model for varying fixed values of tau and using different
#    cut-offs for the CI's used to determine non-generic terms. It will save
#    parameter values for graphing, but not ppc vals.

rm(list = ls())
library(here)
library(HDInterval)
library(RColorBrewer)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Define the CI thresholds we will use
CIs <- c(0.5, 0.8, 0.95)

# Now assign the file paths for the stan models
PrelimStanPath <- here("BH_simulations/Main/StanCode/Prelim_monoLambda_envAlpha.stan")
FinalStanPath <- here("BH_simulations/Main/StanCode/Final_monoLambda_envAlpha.stan")

# Load in the appropriate data
load(here("BH_simulations/test_multiple_simulations.RData"))
ExampleSim <- simulations[[8]]
TrueVals <- ExampleSim[[1]]
FullSim <- ExampleSim[[2]]
Focal <- which(TrueVals$focal == 1)
TrueAlphaMeans <- TrueVals$alpha.8 #This simulation has species 8 as the focal
TrueAlphaSlopes <- TrueVals$alpha.env 

# assign some universal values to be used across model fits and graphs
N <- 50
S <- 15
Intra <- rep(0, S)
Intra[Focal] <- 1
slab_df <- 4 
slab_scale <- sqrt(2)
tau0 <- 1

# Create the data vectors to be passed to rstan for subsequent model fits
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij", "Inclusion_eij")

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
N <- 2*N
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 15))
PrelimPosteriors <- extract(PrelimFit)

AllInclusion_ij <- matrix(0, nrow = S, ncol = 3)
AllInclusion_eij <- matrix(0, nrow = S, ncol = 3)
for(i in 1:3){
     for(s in 1:S){
          Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = CIs[i])
          Ints_eij <- hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = CIs[i])
          if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
               AllInclusion_ij[s,i] <- 1
          }
          if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
               AllInclusion_eij[s,i] <- 1
          }
          if(s == Focal){
               AllInclusion_ij[s,i] <- 0
               AllInclusion_eij[s,i] <- 0
          }
     }
}

AllInclusion_ij
AllInclusion_eij
colSums(AllInclusion_eij)
colSums(AllInclusion_ij)

# Run the final fit of the model
FinalCIFits <- vector(mode = "list", length = 3)
Rhats <- vector(mode = "list", length = 3)
Neffs <- vector(mode = "list", length = 3)
for(i in 1:3){
     Inclusion_eij <- AllInclusion_eij[,i]
     Inclusion_ij <- AllInclusion_ij[,i]
     FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                      chains = 3, control = list(adapt_delta = 0.9))
     FinalCIFits[[i]] <- extract(FinalFit)
     Rhats[[i]] <- summary(FinalFit)$summary[,"Rhat"]
     Neffs[[i]] <- summary(FinalFit)$summary[,"n_eff"]
}

# Quickly examine the Rhat and Neff values for problems
quartz()
par(mfrow = c(3,2))
for(i in 1:3){
     hist(Rhats[[i]])
     hist(Neffs[[i]])
}

# If the fit looks good, safe the final output here
FitFileName <- here("BH_simulations/Tau-CI/CIFits.rdata")
save(FinalCIFits, AllInclusion_ij, AllInclusion_eij, file = FitFileName)

########### Now calculate the parameter deviations

# Now calculate the parameter deviations, starting with lambdas
LambdaEsts <- array(data = NA, dim = c(3,3,2))
for(i in 1:3){
     LambdaEsts[i,1,1] <- mean(FinalCIFits[[i]]$lambdas[,1] - TrueVals$lambda.mean[Focal])
     LambdaEsts[i,2:3,1] <- hdi(FinalCIFits[[i]]$lambdas[,1] - TrueVals$lambda.mean[Focal])
     LambdaEsts[i,1,2] <- mean(FinalCIFits[[i]]$lambdas[,2] - TrueVals$lambda.env[Focal])
     LambdaEsts[i,2:3,2] <- hdi(FinalCIFits[[i]]$lambdas[,2] - TrueVals$lambda.env[Focal])
}


# Now calculate the alpha values
AlphaEsts <- array(data = NA, dim = c(3,3,2,S))
for(i in 1:3){
     for(s in 1:S){
          if(s == Focal){
               # First the intercept
               Intercept <- FinalCIFits[[i]]$alpha_intra[,1]
               AlphaEsts[i,1,1,s] <- mean(Intercept - TrueAlphaMeans[s])
               AlphaEsts[i,2:3,1,s] <- hdi(Intercept - TrueAlphaMeans[s])
               # Now the slope
               Slope <- FinalCIFits[[i]]$alpha_intra[,2]
               AlphaEsts[i,1,2,s] <- mean(Slope - TrueAlphaSlopes[s])
               AlphaEsts[i,2:3,2,s] <- hdi(Slope - TrueAlphaSlopes[s])
          }else{
               # First the intercept
               Intercept <- FinalCIFits[[i]]$alpha_generic[,1] + FinalCIFits[[i]]$alpha_hat_ij[,s] * AllInclusion_ij[s,i]
               AlphaEsts[i,1,1,s] <- mean(Intercept - TrueAlphaMeans[s])
               AlphaEsts[i,2:3,1,s] <- hdi(Intercept - TrueAlphaMeans[s])
               # Now the slope
               Slope <- FinalCIFits[[i]]$alpha_generic[,2] + FinalCIFits[[i]]$alpha_hat_eij[,s] * AllInclusion_eij[s,i]
               AlphaEsts[i,1,2,s] <- mean(Slope - TrueAlphaSlopes[s])
               AlphaEsts[i,2:3,2,s] <- hdi(Slope - TrueAlphaSlopes[s])
          }
     }
}

# Finally, calculate the "true" generic alpha that the model is attempting to estimate
TrueGenericIntercept <- rep(NA, 3)
TrueGenericSlope <- rep(NA, 3)
for(i in 1:3){
     GenericIntercepts <- 1 - AllInclusion_ij[,i]
     GenericSlopes <- 1 - AllInclusion_eij[,i]
     GenericIntercepts[Focal] <- 0
     GenericSlopes[Focal] <- 0
     TrueGenericIntercept[i] <- log(sum(colSums(SpMatrix) * GenericIntercepts * exp(TrueAlphaMeans)) / sum(colSums(SpMatrix) * GenericIntercepts))
     Totals <- colSums(SpMatrix) * GenericSlopes
     TrueGenericSlope[i] <- mean(TrueAlphaSlopes*GenericSlopes*(Totals/max(Totals)))
}



# Finally, save all the necessary results for the figures
FileName <- here("BH_simulations/Tau-CI/CIFits_GraphStuff.rdata")
save(LambdaEsts, AlphaEsts, AllInclusion_eij, AllInclusion_ij, TrueGenericIntercept, TrueGenericSlope,
     CIs, file = FileName)
