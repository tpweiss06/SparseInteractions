# This script will fit the main formulation of the model from the manuscript
#    (monotonic alpha*env and lambda*env relationships) to many simulations.
#    For each simulation, the script will fit the model to 10, 50, and 200 data 
#    points in each "full community" and "thinned" treatment. From each model 
#    fit, the script will save the mean parameter deviation for each parameter,
#    the number of non-generic terms identified (both intercept and slope), the
#    total number of non-generic species involved (since one species could have
#    both or only one term), the RMSE to the out-of-sample data, and key model
#    diagnostics for both the preliminary and final model fits.

rm(list = ls())
setwd("/project/commbayes/SparseInteractions/BH_sims")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(HDInterval)

# Set the current number of simulations to run through
SimIDs <- 901:1000
OutputFile <- "MainSimResults10.csv"

# Load in the simulation data
load("BH_simulations_1200.RData")

# Create a data frame with columns for each of the values specified in the file
#    header
Nsims <- length(SimIDs)
SimID <- rep(SimIDs, each = 3)
Size <- rep(c(10, 50, 200), Nsims)
LambdaIntDev <- rep(NA, Nsims*3)
LambdaSlopeDev <- rep(NA, Nsims*3)
IntraIntDev <- rep(NA, Nsims*3)
IntraSlopeDev <- rep(NA, Nsims*3)
GenericIntDev <- rep(NA, Nsims*3)
GenericSlopeDev <- rep(NA, Nsims*3)
NonGenericIntDev <- rep(NA, Nsims*3)
NonGenericSlopeDev <- rep(NA, Nsims*3)
NumNonGenericInt <- rep(NA, Nsims*3)
NumNonGenericSlope <- rep(NA, Nsims*3)
NumNonGenericSpecies <- rep(NA, Nsims*3)
RMSE <- rep(NA, Nsims*3)
PrelimMaxRhat <- rep(NA, Nsims*3)
PrelimNumDiv <- rep(NA, Nsims*3)
PrelimMeanNeff <- rep(NA, Nsims*3)
FinalMaxRhat <- rep(NA, Nsims*3)
FinalNumDiv <- rep(NA, Nsims*3)
FinalMeanNeff <- rep(NA, Nsims*3)
MainSimResults <- data.frame(SimID, Size, LambdaIntDev, LambdaSlopeDev, IntraIntDev,
                             IntraSlopeDev, GenericIntDev, GenericSlopeDev, NonGenericIntDev,
                             NonGenericSlopeDev, NumNonGenericInt, NumNonGenericSlope,
                             NumNonGenericSpecies, RMSE, PrelimMaxRhat, PrelimNumDiv,
                             PrelimMeanNeff, FinalMaxRhat, FinalNumDiv, FinalMeanNeff)

# Establish thresholds to use for max(Rhat), number of divergent transitions,
#    and mean(n_eff)
RhatThresh <- 1.1
PrelimDivThresh <- 5
FinalDivThresh <- 0
NeffThresh <- 1500

# Set the credible interval threshold for inclusion as a non-generic term
IntLevel <- 0.5

# Set the file paths to be passed to the nodes
PrelimStanPath <- "/project/commbayes/SparseInteractions/BH_sims/Prelim_monoLambda_envAlpha.stan"
FinalStanPath <- "/project/commbayes/SparseInteractions/BH_sims/Final_monoLambda_envAlpha.stan"

# Create the data vectors to be passed to the nodes
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij", "Inclusion_eij")

# Make a function to pass to the nodes
ModelFit <- function(i){
     # Set the size of the dataset and the SimID to use
     N <- MainSimResults$Size[i]
     CurID <- MainSimResults$SimID[i]
     CurParams <- simulations[[CurID]][[1]]
     CurData <- simulations[[CurID]][[2]]
     
     # Determine the focal species for the current simulation
     Focal <- which(CurParams$focal == 1)
     
     # assign some universal values to be used across model fits and graphs
     S <- 15
     Intra <- rep(0, S)
     Intra[Focal] <- 1
     tau0 <- 1
     slab_df <- 4 
     slab_scale <- sqrt(2) 
     
     # Set the local values to pass to rstan
     FullData <- subset(CurData, (species == Focal) & (run <= N) & (time == 0) & (thinned == 0))
     ThinData <- subset(CurData, (species == Focal) & (run <= N) & (time == 0) & (thinned == 1))
     Nt <- c(FullData$pop, ThinData$pop)
     env <- c(FullData$run.env, ThinData$run.env)
     SpMatrix <- matrix(data = NA, nrow = 2*N, ncol = S)
     for(s in 1:S){
          SpMatrix[1:N,s] <- subset(CurData, (species == s) & (run <= N) & (time == 0) & (thinned == 0))$pop
          SpMatrix[(N+1):(2*N),s] <- subset(CurData, (species == s) & (run <= N) & (time == 0) & (thinned == 1))$pop
     }
     Ntp1 <- c(subset(CurData, (species == Focal) & (run <= N) & (time == 1) & (thinned == 0))$pop,
               subset(CurData, (species == Focal) & (run <= N) & (time == 1) & (thinned == 1))$pop)
     
     # Now run the preliminary fit of the model
     N <- 2*N
     PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                       chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 15))
     
     # Now assess model diagnostics
     PrelimMaxRhat <- max(summary(PrelimFit)$summary[,"Rhat"])
     PrelimMeanNeff <- mean(summary(PrelimFit)$summary[,"n_eff"])
     # Get the number of divergent transitions
     sp <- get_sampler_params(PrelimFit, inc_warmup = FALSE)
     PrelimNumDiv <- sum(sapply(sp, FUN = function(x) {
          if ("divergent__" %in% colnames(x)) return(sum(x[,"divergent__"]))
          else return(0)
     }))
     
     if((PrelimMaxRhat > RhatThresh) | (PrelimMeanNeff < NeffThresh) | (PrelimNumDiv > PrelimDivThresh)){
          CurResults <- list(LambdaIntDev = NA, LambdaSlopeDev = NA, IntraIntDev = NA,
                          IntraSlopeDev = NA, GenericIntDev = NA, GenericSlopeDev = NA, NonGenericIntDev = NA,
                          NonGenericSlopeDev = NA, NumNonGenericInt = NA, NumNonGenericSlope = NA,
                          NumNonGenericSpecies = NA, RMSE = NA, PrelimMaxRhat = PrelimMaxRhat, PrelimNumDiv = PrelimNumDiv,
                          PrelimMeanNeff = PrelimMeanNeff, FinalMaxRhat = NA, FinalNumDiv = NA, FinalMeanNeff = NA)
          
     }else{
          PrelimPosteriors <- extract(PrelimFit)
          Inclusion_ij <- rep(0, S)
          Inclusion_eij <- rep(0, S)
          for(s in 1:S){
               Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
               Ints_eij <- hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = IntLevel)
               if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
                    Inclusion_ij[s] <- 1
               }
               if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
                    Inclusion_eij[s] <- 1
               }
               if(s == Focal){
                    Inclusion_ij[s] <- 0
                    Inclusion_eij[s] <- 0
               }
          }
          NumNonGenericInt <- sum(Inclusion_ij)
          NumNonGenericSlope <- sum(Inclusion_eij)
          NumNonGenericSpecies <- sum((Inclusion_ij + Inclusion_eij) > 0)
          
          # Set initial conditions with values from the preliminary fit
          ChainInitials <- list(lambdas = colMeans(PrelimPosteriors$lambdas), 
                                alpha_generic_tilde = colMeans(PrelimPosteriors$alpha_generic_tilde), 
                                alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), 
                                local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
                                c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
                                alpha_hat_eij_tilde = colMeans(PrelimPosteriors$alpha_hat_eij_tilde), 
                                local_shrinkage_eij = colMeans(PrelimPosteriors$local_shrinkage_eij),
                                alpha_intra_tilde = colMeans(PrelimPosteriors$alpha_intra_tilde))
          InitVals <- list(ChainInitials, ChainInitials, ChainInitials)
          
          # Run the final fit of the model
          FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                           chains = 3, init = InitVals, control = list(adapt_delta = 0.99, max_treedepth = 15))
          
          # Now assess model diagnostics
          FinalMaxRhat <- max(summary(FinalFit)$summary[,"Rhat"])
          FinalMeanNeff <- mean(summary(FinalFit)$summary[,"n_eff"])
          # Get the number of divergent transitions
          sp <- get_sampler_params(FinalFit, inc_warmup = FALSE)
          FinalNumDiv <- sum(sapply(sp, FUN = function(x) {
               if ("divergent__" %in% colnames(x)) return(sum(x[,"divergent__"]))
               else return(0)
          }))
          
          if((FinalMaxRhat > RhatThresh) | (FinalMeanNeff < NeffThresh) | (FinalNumDiv > FinalDivThresh)){
               CurResults <- list(LambdaIntDev = NA, LambdaSlopeDev = NA, IntraIntDev = NA,
                               IntraSlopeDev = NA, GenericIntDev = NA, GenericSlopeDev = NA, NonGenericIntDev = NA,
                               NonGenericSlopeDev = NA, NumNonGenericInt = NumNonGenericInt, NumNonGenericSlope = NumNonGenericSlope,
                               NumNonGenericSpecies = NumNonGenericSpecies, RMSE = NA, PrelimMaxRhat = PrelimMaxRhat, PrelimNumDiv = PrelimNumDiv,
                               PrelimMeanNeff = PrelimMeanNeff, FinalMaxRhat = FinalMaxRhat, FinalNumDiv = FinalNumDiv, FinalMeanNeff = FinalMeanNeff)
               
          }else{
               FinalPosteriors <- extract(FinalFit)
               
               ##### Calculate the mean deviation of parameter estimates
               LambdaIntDev <- mean(FinalPosteriors$lambdas[,1]) - CurParams$lambda.mean[Focal]
               LambdaSlopeDev <- mean(FinalPosteriors$lambdas[,2]) - CurParams$lambda.env[Focal]
               IntraIntDev <- mean(FinalPosteriors$alpha_intra[,1]) - CurParams[Focal, 3+Focal]
               IntraSlopeDev <- mean(FinalPosteriors$alpha_intra[,2]) - CurParams$alpha.env[Focal]
               # To calculate the generic deviations, calculate what the model would consider
               #    the "true" generic alpha based on which species it includes
               GenericIntercepts <- 1 - Inclusion_ij
               GenericSlopes <- 1 - Inclusion_eij
               GenericIntercepts[Focal] <- 0
               GenericSlopes[Focal] <- 0
               TrueGenericIntercept <- log(sum(colSums(SpMatrix) * GenericIntercepts * exp(CurParams[, 3+Focal])) / sum(colSums(SpMatrix) * GenericIntercepts))
               Totals <- colSums(SpMatrix) * GenericSlopes
               TrueGenericSlope <- mean(CurParams$alpha.env*GenericSlopes*(Totals/max(Totals)))
               GenericIntDev <- mean(FinalPosteriors$alpha_generic[,1]) - TrueGenericIntercept
               GenericSlopeDev <- mean(FinalPosteriors$alpha_generic[,2]) - TrueGenericSlope
               # To calculate the average deviation of the non-generic terms, we need to 
               #    calculate each individual deviation first
               AllNonGenericIntDev <- NULL
               AllNonGenericSlopeDev <- NULL
               for(s in 1:S){
                    if(Inclusion_ij[s] == 1){
                         CurInt <- mean(FinalPosteriors$alpha_generic[,1] + FinalPosteriors$alpha_hat_ij[,s])
                         CurIntDev <- CurInt - CurParams[s,3+Focal]
                         AllNonGenericIntDev <- c(AllNonGenericIntDev, CurIntDev)
                    }
                    if(Inclusion_eij[s] == 1){
                         CurSlope <- mean(FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s])
                         CurSlopeDev <- CurSlope - CurParams$alpha.env[s]
                         AllNonGenericSlopeDev <- c(AllNonGenericSlopeDev, CurSlopeDev)
                    }
               }
               if(length(AllNonGenericIntDev) > 0){
                    NonGenericIntDev <- mean(AllNonGenericIntDev)
               }else{
                    NonGenericIntDev <- NA
               }
               if(length(AllNonGenericSlopeDev) > 0){
                    NonGenericSlopeDev <- mean(AllNonGenericSlopeDev)
               }else{
                    NonGenericSlopeDev <- NA
               }
               #### Calculate the Root Mean Squared Error
               # First, calculate the observed growth in the out-of-sample data for othe posterior predictive check
               max_N <- 200
               ppc_data <- subset(CurData, (species == Focal) & (run > max_N) & (time == 0) & (thinned == 0))
               ppc_points <- which(ppc_data$pop > 0)
               ppc_runs <- ppc_data$run[ppc_points]
               N_ppc <- length(ppc_points)
               Nt_ppc <- ppc_data$pop[ppc_points]
               env_ppc <- ppc_data$run.env[ppc_points]
               SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
               for(s in 1:S){
                    SpMatrix_ppc[,s] <- subset(CurData, (species == s) & (run %in% ppc_runs) & (time == 0) & (thinned == 0))$pop
               }
               Ntp1_ppc <- subset(CurData, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
               Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)
               
               # Now generate the predicted values from the model fit
               PostLength <- length(FinalPosteriors$alpha_generic[,1])
               # Calculate the posterior distributions of the interaction coefficients and lambdas
               alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
               lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
               for(n in 1:N_ppc){
                    lambda_ei[,n] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[n])
                    for(s in 1:S){
                         if(s == Focal){
                              alpha_eij[,n,s] <- exp(FinalPosteriors$alpha_intra[,1] + FinalPosteriors$alpha_intra[,2] * env_ppc[n])
                         }else{
                              alpha_eij[,n,s] <- exp(FinalPosteriors$alpha_generic[1] + Inclusion_ij[s] * FinalPosteriors$alpha_hat_ij[,s] +
                                                          (FinalPosteriors$alpha_generic[,2] + Inclusion_eij[s] * FinalPosteriors$alpha_hat_eij[,s]) * env_ppc[n])
                         }
                    }
               }
               
               # Use the above quantities to calculate the posterior prediction intervals for the new data
               Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
               for(p in 1:PostLength){
                    for(j in 1:N_ppc){
                         SigmaTerm <- sum(alpha_eij[p,j,] * SpMatrix_ppc[j,])
                         Ntp1_pred <- Nt_ppc[j] * lambda_ei[p,j] / (1 + SigmaTerm)
                         Growth_pred[p,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
                    }
               }
               
               # Now calculate the posterior RMSE for growth predictions
               GrowthRMSE <- rep(NA, PostLength)
               for(p in 1:PostLength){
                         deviation_sq <- (Growth_pred[p,] - Growth_ppc)^2
                         GrowthRMSE[p] <- sqrt(sum(deviation_sq) / N_ppc)
               }
               RMSE <- mean(GrowthRMSE)
          
               # Finally, save all the relevant output into a list to return
               CurResults <- list(LambdaIntDev = LambdaIntDev, LambdaSlopeDev = LambdaSlopeDev, IntraIntDev = IntraIntDev,
                               IntraSlopeDev = IntraSlopeDev, GenericIntDev = GenericIntDev, GenericSlopeDev = GenericSlopeDev, NonGenericIntDev = NonGenericIntDev,
                               NonGenericSlopeDev = NonGenericSlopeDev, NumNonGenericInt = NumNonGenericInt, NumNonGenericSlope = NumNonGenericSlope,
                               NumNonGenericSpecies = NumNonGenericSpecies, RMSE = RMSE, PrelimMaxRhat = PrelimMaxRhat, PrelimNumDiv = PrelimNumDiv,
                               PrelimMeanNeff = PrelimMeanNeff, FinalMaxRhat = FinalMaxRhat, FinalNumDiv = FinalNumDiv, FinalMeanNeff = FinalMeanNeff)
               
          }
     }
   return(CurResults)
}

# Extract the results, population the Results matrix, and save the output
for(i in 1:length(SimIDs)){
   CurFit <- ModelFit(i)
     MainSimResults$LambdaIntDev[i] <- CurFit$LambdaIntDev
     MainSimResults$LambdaSlopeDev[i] <- CurFit$LambdaSlopeDev
     MainSimResults$IntraIntDev[i] <- CurFit$IntraIntDev
     MainSimResults$IntraSlopeDev[i] <- CurFit$IntraSlopeDev
     MainSimResults$GenericIntDev[i] <- CurFit$GenericIntDev
     MainSimResults$GenericSlopeDev[i] <- CurFit$GenericSlopeDev
     MainSimResults$NonGenericIntDev[i] <- CurFit$NonGenericIntDev
     MainSimResults$NonGenericSlopeDev[i] <- CurFit$NonGenericSlopeDev
     MainSimResults$NumNonGenericInt[i] <- CurFit$NumNonGenericInt
     MainSimResults$NumNonGenericSlope[i] <- CurFit$NumNonGenericSlope
     MainSimResults$NumNonGenericSpecies[i] <- CurFit$NumNonGenericSpecies
     MainSimResults$RMSE[i] <- CurFit$RMSE
     MainSimResults$PrelimMaxRhat[i] <- CurFit$PrelimMaxRhat
     MainSimResults$PrelimMeanNeff[i] <- CurFit$PrelimMeanNeff
     MainSimResults$PrelimNumDiv[i] <- CurFit$PrelimNumDiv
     MainSimResults$FinalMaxRhat[i] <- CurFit$FinalMaxRhat
     MainSimResults$FinalMeanNeff[i] <- CurFit$FinalMeanNeff
     MainSimResults$FinalNumDiv[i] <- CurFit$FinalNumDiv
     write.csv(MainSimResults, file = OutputFile, row.names = FALSE, quote = FALSE)
}

# Save the results
# write.csv(AllSims, file = "MainSimResults_10.csv", row.names = FALSE, quote = FALSE)



