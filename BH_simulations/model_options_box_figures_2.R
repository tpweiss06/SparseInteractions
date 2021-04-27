library(tidyverse)
library(patchwork)
library(HDInterval)
setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Create lists for the results from different sizes of datasets
FilePrefixes <- c("N10_", "N50_", "N200_")
n.samples <- length(FilePrefixes)
ppcPreds <- vector(mode = "list", length = n.samples)
Alphas <- vector(mode = "list", length = n.samples)
Lambdas <- vector(mode = "list", length = n.samples)
InterceptInclusion <- vector(mode = "list", length = n.samples)
SlopeInclusion <- vector(mode = "list", length = n.samples)
GenericDeviations <- vector(mode = "list", length = n.samples)

for(i in 1:3){
     FileName <- paste("StanFits/monoLambda_envAlpha/",
                       FilePrefixes[i], "GraphStuff.rdata", sep = "")
     load(FileName)
     ppcPreds[[i]] <- PredVals
     Alphas[[i]] <- AlphaEsts
     Lambdas[[i]] <- LambdaEsts
     InterceptInclusion[[i]] <- Inclusion_ij
     SlopeInclusion[[i]] <- Inclusion_eij
     FileName <- paste("StanFits/monoLambda_envAlpha/",
                       FilePrefixes[i], "FinalFit.rdata", sep = "")
     load(FileName)
     # Calculate the mean and hdi for the posterior of the deviations in the
     #      generic alpha terms
     InterceptDev <- Posteriors$alpha_generic[,1] - TrueGenericIntercept
     SlopeDev <- Posteriors$alpha_generic[,2] - TrueGenericSlope
     GenericDeviations[[i]] <- matrix(data = NA, nrow = 2, ncol = 3)
     GenericDeviations[[i]][1,1] <- mean(InterceptDev)
     GenericDeviations[[i]][2,1] <- mean(SlopeDev)
     GenericDeviations[[i]][1,2:3] <- HDInterval::hdi(InterceptDev)
     GenericDeviations[[i]][2,2:3] <- HDInterval::hdi(SlopeDev)
}
