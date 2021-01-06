# This script will fit the sparse model code to a simulation from Chhaya with 1000
#    data points. It will fit the model to variously sized subsets of the data and
#    use the final data points (751-1000) for a posterior predictive check. This 
#    initial stage will run on Teton and save the model output and some diagnostics.
#    Then I will use a separate script to visualize the comparison between the 
#    posterior predicitve interval and the "true" data from the simulation.

library(rstan)
library(HDInterval)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the working directory and load in the simulation data
setwd("/project/commbayes/SparseInteractions/BH_sims/")
FullSim <- read.csv("")

# Set the size of the datasets we'll use
N_vec <- c(50, 100, 250, 500, 750)

################################# Create a function to fit the model with an 
#    argument determining the number of data points to use.
ModelFit <- function(NumData){	
     # Set some of the universal input data for the stan model and a vector to send
     #    the objects to the cluster
     N <- NumData
     S <- 10
     tau0 <- 1
     slab_scale <- log(2)
     slab_df <- 25
     PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "tau0", "slab_scale", "slab_df")
     FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Inclusion_ij", "Inclusion_eij")
     
     # subset the data to only consider species 1 and N data points
     Gen0Data <- subset(SimData, species == 1 & time == 0 & run <= N)
     Gen1Data <- subset(SimData, species == 1 & time == 1 & run <= N)
     
     # Extract the necessary data for the model
     env <- Gen0Data$run.env
     Nt <- round(Gen0Data$pop)
     Ntp1 <- round(Gen1Data$pop)
     
     # Need to divide by initial population to get per capita fecundity to match
     SpMatrix <- matrix(data = NA, nrow = N, ncol = S)
     for(s in 1:S){
          SpMatrix[,s] <- round(subset(SimData, species == s & time == 0 & run <= N)$pop)
     }

     # Now run the perliminary fit of the model to assess parameter shrinkage
     PrelimFit <- stan(file = "BH_FH_Preliminary_PPC.stan", data = PrelimDataVec, iter = 6000,
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
     FinalFit <- stan(file = "BH_Final_PPC.stan", data = FinalDataVec, iter = 6000,
                      chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 20),
                      warmup = 5000)
     FinalPosteriors <- extract(PrelimFit)

     # Save some diagnostic values
     FinalRhats <- summary(FinalFit)$summary[,"Rhat"]
     FinalNeffs <- summary(FinalFit)$summary[,"n_eff"]

     # Return the posteriors and diagnostics
     ModelResults <- list(PrelimRhats = PrelimRhats, PrelimNeffs = PrelimNeffs, PrelimPosteriors = PrelimPosteriors, 
                              FinalRhats = FinalRhats, FinalNeffs = FinalNeffs, FinalPosteriors = FinalPosteriors)
     return(ModelResults)
}

# Run the model for each number of data points
PostPredResults <- vector(mode = "list", length = length(N_vec))
for(i in 1:length(N_vec)){
     PostPredResults[[i]] <- ModelFit(i)
}

save(PostPredResults, file = "PostPredResults.rdata")
