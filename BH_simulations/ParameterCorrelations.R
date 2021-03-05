# This file will create graphs of the parameter correlations for the four different
#    scenarios

setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")
library(rstan)

# First, load in all the scenarios
load("StanFits/monoLambda_constAlpha/N200_FinalFit.rdata")
Scen1 <- FinalPosteriors
rm(FinalFit, FinalPosteriors, Inclusion_ij)
load("StanFits/optLambda_constAlpha/N200_FinalFit.rdata")
Scen2 <- FinalPosteriors
rm(FinalFit, FinalPosteriors, Inclusion_ij)
load("StanFits/monoLambda_envAlpha/N200_FinalFit.rdata")
Scen3 <- FinalPosteriors
rm(FinalFit, FinalPosteriors, Inclusion_ij, Inclusion_eij)
load("StanFits/optLambda_envAlpha/N200_FinalFit.rdata")
Scen4 <- FinalPosteriors
rm(FinalFit, FinalPosteriors, Inclusion_ij, Inclusion_eij)

# Scen1: show correlation between lambda intercept and generic alpha
# Scen2: show correlation between lambda_max and generic alpha
# Scen3: show correlation between lambda intercept and generic alpha intercept
# Scen4: show correlation between lambda_opt and generic slope;
pdf(file = "Results/ParameterCorrelations.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2))
     plot(x = Scen1$lambdas[,1], y = Scen1$alpha_generic, xlab = "Lambda intercept",
          ylab = "Generic alpha", main = "Scenario 1", las = 1)
     plot(x = Scen2$lambda_max, y = Scen2$alpha_generic, xlab = "Lambda maximum",
          ylab = "Generic alpha", main = "Scenario 2", las = 1)
     plot(x = Scen3$lambdas[,1], y = Scen3$alpha_generic[,1], xlab = "Lambda intercept",
          ylab = "Generic alpha intercept", main = "Scenario 3", las = 1)
     plot(x = Scen4$lambda_opt, y = Scen4$alpha_generic[,2], xlab = "Lambda optimum",
          ylab = "Generic alpha slope", main = "Scenario 4", las = 1)
dev.off()


