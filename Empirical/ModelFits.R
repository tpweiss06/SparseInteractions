# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

rm(list = ls())
library(rstan)
library(here)

FocalLetter <- "W" # "W" or "A"
FocalPrefix <- "WAAC" # "Waitzia" or "ARCA"
FocalSpecies <- "Waitzia.acuminata" # "Waitzia.acuminata" or "Arctotheca.calendula"

EnvCov <- "Phos" # "Phos" or "Shade"
EnvCol <- 71  # 72 for Canopy or 71 for Phosphorous

# Load in the data and subset out the current focal species.
SpData <- read.csv(here("Empirical/water_full_env.csv"))
SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
SpData <- na.omit(SpData) 
FocalLetter
SpDataFocal <- subset(SpData, Focal.sp.x == FocalLetter)

# Next continue to extract the data needed to run the model. 
N <- as.integer(nrow(SpDataFocal))
Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
reserve <- as.integer(SpDataFocal$Reserve.x)
env <- as.vector(scale(SpDataFocal[,EnvCol]))

# Now calculate the total number of species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpAbunds <- SpDataFocal[,10:69]
AllSpNames <- names(SpDataFocal[10:69])
SpTotals <- colSums(AllSpAbunds)
SpToKeep <- SpTotals > 0
S <- sum(SpToKeep)
SpMatrix <- matrix(NA, nrow = N, ncol = S)
i <- 1
for(s in 1:ncol(AllSpAbunds)){
     if(SpToKeep[s] == 1){
          SpMatrix[,i] <- AllSpAbunds[,s]
          i <- i + 1
     }
}
SpNames <- AllSpNames[SpToKeep]
Intra <- ifelse(SpNames == FocalSpecies, 1, 0)

# load rstan, set the parameters for the Finnish Horseshoe, and create a vector
#       of all the data objects needed for the model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
tau0 <- 2
slab_scale <- 1
slab_df <- 2*S - 1
DataVec <- c("N", "S", "Fecundity", "reserve", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")

# Now run a perliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = here("Empirical/StanCode/BH_FH_Preliminary.stan"), data = DataVec, iter = 3000, 
                  chains = 3)
PrelimPosteriors <- extract(PrelimFit)

##### Diagnostic plots
# First check the distribution of Rhats and effective sample sizes
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
# Next check the correlation among key model parameters and identify any
#       divergent transitions
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
# Finally, check for autocorrelation in the posteriors of key model parameters
acf(PrelimPosteriors$lambdas[,1,1])
acf(PrelimPosteriors$lambdas[,1,2])
acf(PrelimPosteriors$lambdas[,2,1])
acf(PrelimPosteriors$lambdas[,2,2])
acf(PrelimPosteriors$alpha_generic[,1])
acf(PrelimPosteriors$alpha_generic[,2])
acf(PrelimPosteriors$alpha_intra[,1])
acf(PrelimPosteriors$alpha_intra[,2])

#### If the diagnostic plots don't reveal any problems wiht the model fit, now
#       move on to determining which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
Inclusion_ij <- matrix(data = 0, nrow = 2, ncol = S)
Inclusion_eij <- matrix(data = 0, nrow = 2, ncol = S)
IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
for(i in 1:2){
     for(s in 1:S){
          Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,i,s], credMass = IntLevel)
          Ints_eij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_eij[,i,s], credMass = IntLevel)
          if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
               Inclusion_ij[i,s] <- 1
          }
     }
}
sum(Inclusion_ij)
sum(Inclusion_eij)

DataVec <- c("N", "S", "Fecundity", "reserve", "SpMatrix", "env", "Intra",
                  "Inclusion_ij", "Inclusion_eij")
FinalFit <- stan(file = here("Empirical/StanCode/BH_Final.stan"), data = DataVec, iter = 3000, chains = 3)
FinalPosteriors <- extract(FinalFit)

# Diagnostic figures
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
acf(FinalPosteriors$lambdas[,1,1])
acf(FinalPosteriors$lambdas[,1,2])
acf(FinalPosteriors$lambdas[,2,1])
acf(FinalPosteriors$lambdas[,2,2])
acf(FinalPosteriors$alpha_generic[,1])
acf(FinalPosteriors$alpha_generic[,2])
acf(FinalPosteriors$alpha_intra[,1])
acf(FinalPosteriors$alpha_intra[,2])

FileName <- paste(here("Empirical/StanFits/"), FocalPrefix, "_", EnvCov, "_FinalFit.rdata", sep = "")
save(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
     Inclusion_eij, tau0, slab_scale, slab_df, Intra, file = FileName)

