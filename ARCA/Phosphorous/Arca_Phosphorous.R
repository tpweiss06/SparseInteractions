# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

setwd("~/Desktop/Wyoming/SparseInteractions/")

FocalLetter <- "A" # "W"
FocalPrefix <- "ARCA" # "Waitzia"
FocalSpecies <- "Arctotheca.calendula" # "Waitzia.acuminata"

EnvCov <- "Phos" # "Shade"
EnvCol <- 71  # 72 for Canopy

# Load in the data and subset out the current focal species.
SpData <- read.csv("water_full_env.csv")
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
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
tau0 <- 2
slab_scale <- 1
slab_df <- 2*S - 1
DataVec <- c("N", "S", "Fecundity", "reserve", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")

# Now run a perliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = "Model Code/BH_FH_Preliminary.stan", data = DataVec, iter = 3000, 
                  chains = 3)
save(PrelimFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, tau0, slab_scale, slab_df,
     file = paste("ARCA/Phosphorous/Model Fits/" FocalPrefix, EnvCov, "PrelimFit.rdata", sep = "_"))
     
# Now evaluate some of the model diagnostics and save them as figures
PrelimPosteriors <- extract(PrelimFit)
Rhats <- summary(PrelimFit)$summary[,"Rhat"]   # max < 1.1
Neffs <- summary(PrelimFit)$summary[,"n_eff"]  # min > 2000

# First make a figure with the Rhat and effective sample sizes as histograms
FigName <- paste(FocalPrefix, EnvCov, "PrelimDiagnostics.pdf", sep = "_")
pdf(file = FigName, width = 5, height = 4, onefile = FALSE, paper = "special")
     par(mfrow = c(1,2))
     hist(Rhats, main = "", xlab = expression(hat(R)))
     hist(Neffs, main = "", xlab = "Effective sample size")
dev.off()

# Now make a figure showing the autocorrelation for some key parameters
FigName <- paste(FocalPrefix, EnvCov, "PrelimAutocorrelation.pdf", sep = "_")
pdf(file = FigName, width = 7, height = 5, onefile = FALSE, paper = "special")
     par(mfrow = c(2,4))
     for(i in 1:2){
          acf(PrelimPosteriors$alphas[,i])
          for(j in 1:2){
               acf(PrelimPosteriors$lambdas[,i,j])
          }
     }
     acf(PrelimPosteriors$c2)
     acf(PrelimPosteriors$tau)
dev.off()

# Now create a visual of some key parameter estimates
FigName <- paste(FocalPrefix, EnvCov, "PrelimNonShrinkageEstimates.pdf", sep = "_")
pdf(file = FigName, width = 5, height = 3, onefile = FALSE, paper = "special")
     plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95,
          fill_color = "salmon", pars = c("lambdas", "alphas"))
dev.off()

# Finally, determine which parameters warrant inclusion in the final model and create
#    a visualization of those
Inclusion_ij <- matrix(data = 0, nrow = 2, ncol = S)
Inclusion_eij <- matrix(data = 0, nrow = 2, ncol = S)
IntLevel <- 0.5
for(i in 1:2){
     for(s in 1:S){
          Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,i,s], credMass = IntLevel)
          Ints_eij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_eij[,i,s], credMass = IntLevel)
          if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
               Inclusion_ij[i,s] <- 1
          }
     }
}
InclusionColVals <- c(Inclusion_ij[1,], Inclusion_ij[2,], Inclusion_eij[1,], Inclusion_eij[2,])
InclusionCols <- ifelse(InclusionColVals == 1, "green", "red")

FigName <- paste(FocalPrefix, EnvCov, "InclusionEstimates.pdf", sep = "_")
pdf(file = FigName, width = 4, height = 8, onefile = FALSE, paper = "special")
     plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95,
          fill_color = InclusionCols, pars = c("alpha_hat_ij", "alpha_hat_eij"))
dev.off()

# Now check if some of the basic model diagnostics warrant running a final model
if(max(Rhats) < 1.1 & min(Neffs) > 1000){
     DataVec <- c("N", "S", "Fecundity", "reserve", "SpMatrix", "env",
                  "Inclusion_ij", "Inclusion_eij")
     FinalFit <- stan(file = "BH_Final.stan", data = DataVec, iter = 6000, chains = 3, 
                      control = list(adapt_delta = 0.999, max_treedepth = 15))
     save(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
          Inclusion_eij, tau0, slab_scale, slab_df,
          file = paste(FocalPrefix, EnvCov, "FinalFit.rdata", sep = "_"))
     
     # Make similar visualizations of model diagnostics
     FinalPosteriors <- extract(FinalFit)
     Rhats <- summary(FinalFit)$summary[,"Rhat"]
     Neffs <- summary(FinalFit)$summary[,"n_eff"]
     
     FigName <- paste(FocalPrefix, EnvCov, "FinalDiagnostics.pdf", sep = "_")
     pdf(file = FigName, width = 5, height = 4, onefile = FALSE, paper = "special")
          par(mfrow = c(1,2))
          hist(Rhats, main = "", xlab = expression(hat(R)))
          hist(Neffs, main = "", xlab = "Effective sample size")
     dev.off()
     
     FigName <- paste(FocalPrefix, EnvCov, "FinalAutocorrelation.pdf", sep = "_")
     pdf(file = FigName, width = 7, height = 5, onefile = FALSE, paper = "special")
          par(mfrow = c(2,3))
          for(i in 1:2){
               acf(PrelimPosteriors$alphas[,i])
               for(j in 1:2){
                    acf(PrelimPosteriors$lambdas[,i,j])
               }
          }
     dev.off()
     
     FigName <- paste(FocalPrefix, EnvCov, "FinalNonShrinkageEstimates.pdf", sep = "_")
     pdf(file = FigName, width = 5, height = 3, onefile = FALSE, paper = "special")
          plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95,
               fill_color = "salmon", pars = c("lambdas", "alphas"))
     dev.off()
}












