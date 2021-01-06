# This script will attempt to fit the sparse interactions model to scenarios 1
#       through 4 as laid out in "read_me_simulations.docx"
# For each scenario, the model will be fit to both simulations (a and b) and
#       the estimated parameters will be compared to the relevant parameter files
#       to assess model performance
# NOTE: For these tests, I'm simply going to use species 1 as the focal species
#       for simplicity. This will be run on Teton.

library(rstan)
library(HDInterval)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the working directory and load in necessary libraries
setwd("/project/commbayes/SparseInteractions/BH_sims/")

# Create some useful vectors to pass to the cluster
ScenarioNames <- c("none", "dem", "env", "env_dem")
SimFiles <- NULL
ParamFiles <- NULL
Prefixes <- NULL
for(i in 1:4){
     SimFiles <- c(SimFiles, paste("simulation_", i, "a_", ScenarioNames[i], ".csv", sep = ""))
     SimFiles <- c(SimFiles, paste("simulation_", i, "b_", ScenarioNames[i], ".csv", sep = ""))
     ParamFiles <- c(ParamFiles, paste("parameters_", i, "a_", ScenarioNames[i], ".csv", sep = ""))
     ParamFiles <- c(ParamFiles, paste("parameters_", i, "b_", ScenarioNames[i], ".csv", sep = ""))
     Prefixes <- c(Prefixes, paste(i, "a_", ScenarioNames[i], sep = ""))
     Prefixes <- c(Prefixes, paste(i, "b_", ScenarioNames[i], sep = ""))
}

################################# Create a function to send to the cluster
ModelFit <- function(CurScen){	
     # Set some of the universal input data for the stan model and a vector to send
     #    the objects to the cluster
     N <- 50
     S <- 10
     tau0 <- 1
     slab_scale <- log(2)
     slab_df <- 25
     PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "tau0", "slab_scale", "slab_df")
     FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Inclusion_ij", "Inclusion_eij")
     
     SimData <- read.csv(SimFiles[CurScen])
     
     # First get rid of any data with species 1 starting at 0 abundance
     Gen0Data <- subset(SimData, species == 1 & time == 0)
     Gen1Data <- subset(SimData, species == 1 & time == 1)
     
     # Extract the necessary data for the model
     env <- Gen0Data$run.env
     Nt <- round(Gen0Data$pop)
     Ntp1 <- round(Gen1Data$pop)
     
     # Need to divide by initial population to get per capita fecundity to match
     SpMatrix <- matrix(data = NA, nrow = N, ncol = S)
     for(s in 1:S){
          SpMatrix[,s] <- round(subset(SimData, species == s & time == 0)$pop)
     }

     # Now run the perliminary fit of the model to assess parameter shrinkage
     PrelimFit <- stan(file = "BH_FH_Preliminary_Sims.stan", data = PrelimDataVec, iter = 6000,
                       chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 20),
                       warmup = 5000)
     PrelimPosteriors <- extract(PrelimFit)

     # # Save some diagnostic values
     PrelimRhats <- summary(PrelimFit)$summary[,"Rhat"]
     PrelimNeffs <- summary(PrelimFit)$summary[,"n_eff"]

     # Determine which parameters warrant inclusion in the final model
     Inclusion_ij <- rep(0, S)
     Inclusion_eij <- rep(0, S)
     IntLevel <- 0.5
     for(s in 1:S){
             Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
             Ints_eij <- hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = IntLevel)
             if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
                     Inclusion_ij[s] <- 1
             }
             if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
                     Inclusion_eij[s] <- 1
             }
     }

     # Run the final fit of the model
     FinalFit <- stan(file = "BH_Final_Sims.stan", data = FinalDataVec, iter = 6000,
                      chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 20),
                      warmup = 5000)
     FinalPosteriors <- extract(PrelimFit)

     # Save some diagnostic values
     FinalRhats <- summary(FinalFit)$summary[,"Rhat"]
     FinalNeffs <- summary(FinalFit)$summary[,"n_eff"]

     # Calculate lambda and alpha posteriors across the environment
     EnvVals <- sort(unique(env))
     EnvLength <- length(EnvVals)
     PostLength <- dim(FinalPosteriors$lambdas)[1]
     Lambda <- matrix(data = NA, nrow = PostLength, ncol = EnvLength)
     Alphas <- array(data = NA, dim = c(S, PostLength, EnvLength))
     for(i in 1:PostLength){
          Lambda[i,] <- exp(FinalPosteriors$lambdas[i,1] + FinalPosteriors$lambdas[i,2]*EnvVals)
          for(s in 1:S){
               if(Inclusion_ij[s] == 1){
                    if(Inclusion_eij[s] == 1){
                         Alphas[s,i,] <- exp(FinalPosteriors$alphas[i,1] + FinalPosteriors$alpha_hat_ij[i,s] + (FinalPosteriors$alphas[i,2] + FinalPosteriors$alpha_hat_eij[s])*EnvVals)
                    }else{
                         Alphas[s,i,] <- exp(FinalPosteriors$alphas[i,1] + FinalPosteriors$alpha_hat_ij[i,s] + (FinalPosteriors$alphas[i,2])*EnvVals)
                    }
               }else{
                    if(Inclusion_eij[s] == 1){
                         Alphas[s,i,] <- exp(FinalPosteriors$alphas[i,1] + (FinalPosteriors$alphas[i,2] + FinalPosteriors$alpha_hat_eij[s])*EnvVals)
                    }else{
                         Alphas[s,i,] <- exp(FinalPosteriors$alphas[i,1] + FinalPosteriors$alphas[i,2]*EnvVals)
                    }
               }
          }
     }

     # Generate mean estimates and 95% credible intervals
     LambdaEsts <- matrix(data = NA, nrow = 3, ncol = EnvLength)
     AlphaEsts <- array(data = NA, dim = c(S, 3, EnvLength))
     for(e in 1:EnvLength){
          LambdaEsts[1,e] <- mean(Lambda[,e])
          LambdaEsts[2:3,e] <- hdi(Lambda[,e])
          for(s in 1:S){
               AlphaEsts[s,1,e] <- mean(Alphas[s,,e])
               AlphaEsts[s,2:3,e] <- hdi(Alphas[s,,e])
          }
     }

     # Load in the relevant parameters to compare on the graph
     TrueVals <- read.csv(ParamFiles[CurScen])

     # Plot parameter estimate distributions with vertical lines at true values
     # A 2 x 6 figure with the first plot being lambda across the environment and the
     #    next 10 plots being the alphas across the environment
     FigName <- paste("Results_", Prefixes[CurScen], ".pdf", sep = "")
     pdf(file = FigName, width = 10, height = 7, onefile = FALSE, paper = "special")
          par(mfrow = c(2,6))
          # First plot lambda
          plot(NA, NA, main = "", xlab = "", ylab = "Lambda", xlim = range(EnvVals), ylim = c(0, 5))
          abline(a = TrueVals$lambda.mean[1], b = TrueVals$lambda.env[1], lwd = 2)
          lines(x = EnvVals, y = LambdaEsts[1,], col = "darkred")
          lines(x = EnvVals, y = LambdaEsts[2,], col = "darkred", lty = 2)
          lines(x = EnvVals, y = LambdaEsts[3,], col = "darkred", lty = 2)
          # Now the alpha values
          for(s in 1:S){
               yLabel <- paste("alpha", s, sep = " ")
               plot(NA, NA, main = "", xlab = "", ylab = yLabel, xlim = range(EnvVals), ylim = c(0, 0.02))
               abline(h = TrueVals[1,4+s], lwd = 2)
               lines(x = EnvVals, y = AlphaEsts[s,1,], col = "darkred")
               lines(x = EnvVals, y = AlphaEsts[s,2,], col = "darkred", lty = 2)
               lines(x = EnvVals, y = AlphaEsts[s,3,], col = "darkred", lty = 2)
          }
          mtext("Environmental value", side = 1, outer = TRUE)
     dev.off()

     # Return the Rhats and things
     ModelDiagnostics <- list(PrelimRhats = PrelimRhats, PrelimNeffs = PrelimNeffs, 
                              FinalRhats = FinalRhats, FinalNeffs = FinalNeffs)
     return(ModelDiagnostics)
}

# Run the model for each scenario
PrelimRhats <- vector(mode = "list", length = 8)
PrelimNeffs <- vector(mode = "list", length = 8)
FinalRhats <- vector(mode = "list", length = 8)
FinalNeffs <- vector(mode = "list", length = 8)
for(i in 1:4){
        Results <- ModelFit(i)
        PrelimRhats[[i]] <- Results$PrelimRhats
        PrelimNeffs[[i]] <- Results$PrelimNeffs
        FinalRhats[[i]] <- Results$FinalRhats
        FinalNeffs[[i]] <- Results$FinalNeffs
}

save(PrelimRhats, PrelimNeffs, FinalRhats, FinalNeffs, file = "ModelDiagnostics.rdata")
