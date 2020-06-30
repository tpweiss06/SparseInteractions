#### This script will run the Beverton-Holt model using the Finnish Horseshoe prior
#       to implement shrinkage of the species-specific contributions to alpha
#       (both in terms of the intercept and the slope with phosphorous)

setwd("/project/commbayes/SparseInteractions")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in the data and sort it
SpData <- read.csv("water_full_env.csv")
SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
SpData <- na.omit(SpData) 
SpDataFocal <- subset(SpData, Focal.sp.x == "W")

# From here we need to calculate and create the Intra vector, the
#   SpMatrix of abundances for all heterospecifics, and other necessary
#   objects like N, S, and Fecundity
shade <- scale(SpDataFocal$Canopy)
phos <- scale(SpDataFocal$Colwell.P)
Intra <- as.integer(SpDataFocal$Waitzia.acuminata)
N <- as.integer(dim(SpDataFocal)[1])
Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
Species <- names(SpDataFocal[10:69]) 
Species <- setdiff(Species, "Waitzia.acuminata")
TempSpMatrix <- subset(SpDataFocal, select = Species)
# Now discount any columns with 0 abundance
SpTotals <- colSums(TempSpMatrix)
SpToKeep <- SpTotals > 0
SpMatrix <- matrix(NA, nrow = nrow(TempSpMatrix), ncol = sum(SpToKeep)+1)
SpMatrix[,1] <- Intra
s <- 2
for(i in 1:ncol(TempSpMatrix)){
     if(SpToKeep[i]){
          SpMatrix[,s] <- TempSpMatrix[,i]
          s <- s + 1
     }
}
print(s)
S <- ncol(SpMatrix)
print(S)

# Do a preliminary fit to compile the stan model and check for convergence, mixing,
#    autocorrelation, etc.
tau0 <- 1
slab_scale <- log(2)
slab_df <- 25 # follows the online implementation. this might need tuning
DataVec <- c("N", "S", "Fecundity", "SpMatrix", "tau0", "slab_scale", "slab_df", "phos")
nIter <- 3000
nChains <- 3
PrelimFit <- stan(file = "BH_FH_Phos.stan", data = DataVec, iter = nIter, chains = nChains,
                  control = list(adapt_delta = 0.99, max_treedepth = 15))

# Calculate the shrinkage coefficients for the alpha_sp values ccording to 
#       Bahdra et al. (2019) and Piironen and Vehtari (2017):
#              kappa_i = 1/(1 + N*sigma_tilde^-2*lambda_i^2*tau^2*s^2)
# Equation 3.14 in Piironen and Vehtari (2017)
post <- extract(PrelimFit)
# First calculate the pseudo-variance (approximate to normal) from the poisson likelihood
PostLength <- nIter/2 * nChains
F_hat <- matrix(NA, nrow = N, ncol = PostLength)
SigmaTilde <- rep(NA, PostLength)
for(n in 1:PostLength){
        InteractionEffects <- rep(NA, N)
        lambda <- exp(post$lambdas[n,1] + post$lambdas[n,2] * phos)
        for(i in 1:N){
                alpha_terms<- exp(post$alphas[n,1] + post$alpha_sp[n,1,] + (post$alphas[n,2] + post$alpha_sp[n,2,]) * phos[i])
                InteractionEffects[i] <- sum(alpha_terms * SpMatrix[i,])
        }
        F_hat[,n] <- lambda / (1 + InteractionEffects)
        SigmaTilde[n] <- 1 / mean(F_hat[,n])
}

# Now calculate the pseudo-posteriors for the shrinkage coefficients
ShrinkCoef <- array(NA, dim = c(2,S,nIter))
for(i in 1:2){
        SpConstTerm <- N*SigmaTilde^2 * post$tau^2 * slab_scale^2
        for(s in 1:S){
                ShrinkCoef[1,s,] <- 1 / (1 + SpConstTerm * post$local_shrinkage[,1,s]^2)
                ShrinkCoef[2,s,] <- 1 / (1 + SpConstTerm * post$local_shrinkage[,2,s]^2)
        }
}

save(PrelimFit, ShrinkCoef, file = "Current_BH_FH_Phos_fit.rdata")

