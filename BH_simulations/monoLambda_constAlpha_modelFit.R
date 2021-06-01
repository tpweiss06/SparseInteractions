# This script will fit the model with a monotonic lambda and constant alphas
#    (with respect to the environment) to the corresponding data, using a series
#    of dataset sizes (10, 20, 50, 100, 200), and then examine the accuracy of
#    parameter estimates and posterior predictive checks on 300 out of sample
#    data points.

#setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")
setwd("~/Documents/Work/Current Papers/SparseInteractions/BH_simulations/")

# Set the current sample size and associated prefix for all graph and result
#    file names
N <- 80
max_N <- 200
FilePrefix <- paste("N", N, "_", sep = "")

# Now assign the focal species and the file paths for the stan models
Focal <- 1
PrelimStanPath <- "StanCode/Prelim_monoLambda_constAlpha.stan"
FinalStanPath <- "StanCode/Final_monoLambda_constAlpha.stan"

# Load in the appropriate data
FullSim <- read.csv("Simulations/simulation_perturb2_const.csv")
TrueVals <- read.csv("Simulations/parameters_perturb2_const.csv")
TrueAlphas <- TrueVals$alpha.1

# Load necessary libraries
library(rstan)
library(HDInterval)
library(RColorBrewer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# assign some universal values to be used across model fits and graphs
S <- 15
Intra <- rep(0, S)
Intra[Focal] <- 1
tau0 <- 1
slab_df <- 4
slab_scale <- sqrt(2)

# Set initial values to avoid initial problems with the random number generator
ChainInitials <- list(lambdas = c(TrueVals$lambda.mean[Focal], TrueVals$lambda.env[Focal]), 
                      alpha_generic_tilde = mean(TrueAlphas)/0.75 + 1,  
                      alpha_hat_ij_tilde = rep(0, S), c2_tilde = 1.25,
                      local_shrinkage_ij = rep(5, S), tau_tilde = 15,
                      alpha_intra_tilde = TrueAlphas[Focal])
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
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
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

# Now run the perliminary fit of the model to assess parameter shrinkage
N <- length(Nt)
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, init = InitVals, control = list(max_treedepth = 15, adapt_delta = 0.995))
PrelimPosteriors <- extract(PrelimFit)
FitFileName <- paste("StanFits/monoLambda_constAlpha/", FilePrefix, "PrelimFit_b.rdata", sep = "")
save(PrelimFit, PrelimPosteriors, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
traceplot(PrelimFit, pars = "lambdas")
traceplot(PrelimFit, pars = "alpha_generic")
traceplot(PrelimFit, pars = "alpha_intra")
traceplot(PrelimFit, pars = "alpha_hat_ij")
acf(PrelimPosteriors$alpha_generic[,1])
acf(PrelimPosteriors$alpha_generic[,2])
acf(PrelimPosteriors$alpha_intra[,1])
acf(PrelimPosteriors$alpha_intra[,2])
PlotSamples <- sample(1:S, size = 4, replace = FALSE)
quartz()
for(i in 1:4){
        acf(PrelimPosteriors$alpha_hat_ij[,PlotSamples[i]])
}

### PPC of prelim fits--not used anymore-----------
########### Posterior Predictive Check
# PostLength <- length(PrelimPosteriors$alpha_generic)
# # calculate the posterior distributions of the interaction coefficients
# alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
# for(s in 1:S){
#         alpha_ij[,s] <- exp(PrelimPosteriors$alpha_generic + PrelimPosteriors$alpha_hat_ij[,s])
# }
# alpha_intra <- exp(PrelimPosteriors$alpha_intra)
# # calculate the posterior distributions of lambda_ei
# lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
# for(i in 1:N_ppc){
#         lambda_ei[,i] <- exp(PrelimPosteriors$lambdas[,1] + PrelimPosteriors$lambdas[,2]*env_ppc[i])
# }
# 
# # use the above quantities to calculate the posterior prediction intervals for the new data
# Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
# Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
# for(i in 1:PostLength){
#         for(j in 1:N_ppc){
#                 SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,]) + alpha_intra[i]*Nt_ppc[j]
#                 Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
#                 Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
#                 Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
#         }
# }
# 
# # Calculate preliminary fit results for the ppc
# PrelimPredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
# PrelimCIwidths <- rep(NA, N_ppc)
# PrelimDevVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
# for(i in 1:N_ppc){
#      PrelimPredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
#      PrelimPredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
#      PrelimCIwidths[i] <- PrelimPredVals[3,i] - PrelimPredVals[2,i]
#      PrelimDevVals[1,i] <- mean(Growth_dev[,i], na.rm = TRUE)
#      PrelimDevVals[2:3,i] <- HDInterval::hdi(Growth_dev[,i])
# }
# PrelimCIwidth <- mean(PrelimCIwidths)
# 
# # Calculate the accuracy of parameter estimates from the preliminary fits
# PrelimLambdaEsts <- matrix(data = NA, nrow = 3, ncol = 2)
# PrelimLambdaEsts[1,1] <- mean(PrelimPosteriors$lambdas[,1])
# PrelimLambdaEsts[1,2] <- mean(PrelimPosteriors$lambdas[,2])
# PrelimLambdaEsts[2:3,1] <- HDInterval::hdi(PrelimPosteriors$lambdas[,1])
# PrelimLambdaEsts[2:3,2] <- HDInterval::hdi(PrelimPosteriors$lambdas[,2])
# 
# # Now calculate the alpha estimates
# PrelimAlphaEsts <- matrix(data = NA, nrow = 3, ncol = S+1)
# for(s in 1:S){
#      AlphaDev <- alpha_ij[,s] - as.numeric(TrueAlphas[OtherSpecies[s]])
#      PrelimAlphaEsts[1,s] <- mean(AlphaDev)
#      PrelimAlphaEsts[2:3,s] <- HDInterval::hdi(AlphaDev)
# }
# IntraDev <- alpha_intra - TrueAlphas[Focal]
# PrelimAlphaEsts[1,S+1] <- mean(IntraDev)
# PrelimAlphaEsts[2:3,S+1] <- HDInterval::hdi(IntraDev)

### Final model -----------
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
ChainInitials <- list(lambdas = colMeans(PrelimPosteriors$lambdas), 
                      alpha_generic_tilde = mean(PrelimPosteriors$alpha_generic_tilde), 
                      alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), 
                      local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
                      c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
                      alpha_intra_tilde = mean(PrelimPosteriors$alpha_intra_tilde))
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals, control = list(max_treedepth = 15))
FinalPosteriors <- rstan::extract(FinalFit)
FitFileName <- paste("StanFits/monoLambda_constAlpha/", FilePrefix, "FinalFit_b.rdata", sep = "")
save(FinalFit, FinalPosteriors, Inclusion_ij, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = "lambdas")
traceplot(FinalFit, pars = "alpha_generic")
traceplot(FinalFit, pars = "alpha_intra")
which(Inclusion_ij == 1)
traceplot(FinalFit, pars = "alpha_hat_ij")

# Double check the autocorrelation
acf(FinalPosteriors$lambdas[,1])
acf(FinalPosteriors$lambdas[,2])
acf(FinalPosteriors$alpha_generic[,1])
acf(FinalPosteriors$alpha_generic[,2])
acf(FinalPosteriors$alpha_intra[,1])
acf(FinalPosteriors$alpha_intra[,2])
for(s in 1:S){
        if(Inclusion_ij[s] == 1){
                quartz()
                acf(FinalPosteriors$alpha_hat_ij[,s])
        }
}


########### Posterior Predictive Check
PostLength <- length(FinalPosteriors$alpha_generic)
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
     if(s == Focal){
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_intra)
     }else{
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_generic + 
                                      Inclusion_ij[s] * FinalPosteriors$alpha_hat_ij[,s])
     }
     
}
alpha_intra <- exp(FinalPosteriors$alpha_intra)
# calculate the posterior distributions of lambda_ei
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
     lambda_ei[,i] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,])
          Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
          Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
     }
}

# Calculate final fit results for the ppc
PredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
        PredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
        PredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
}

# Calculate the accuracy of parameter estimates from the preliminary fits
LambdaEsts <- matrix(data = NA, nrow = 3, ncol = 2)
LambdaEsts[1,1] <- mean(FinalPosteriors$lambdas[,1] - TrueVals$lambda.mean[Focal])
LambdaEsts[2:3,1] <- hdi(FinalPosteriors$lambdas[,1] - TrueVals$lambda.mean[Focal])
LambdaEsts[1,2] <- mean(FinalPosteriors$lambdas[,2] - TrueVals$lambda.env[Focal])
LambdaEsts[2:3,2] <- hdi(FinalPosteriors$lambdas[,2] - TrueVals$lambda.env[Focal])

# Now calculate the alpha estimates
AlphaEsts <- matrix(data = NA, nrow = 3, ncol = S)
for(s in 1:S){
        if(s == Focal){
                # First the intercept
                Intercept <- FinalPosteriors$alpha_intra
                AlphaEsts[1,s] <- mean(Intercept - TrueAlphas[s])
                AlphaEsts[2:3,s] <- hdi(Intercept - TrueAlphas[s])
        }else{
                # First the intercept
                Intercept <- FinalPosteriors$alpha_generic + FinalPosteriors$alpha_hat_ij[,s] * Inclusion_ij[s]
                AlphaEsts[1,s] <- mean(Intercept - TrueAlphas[s])
                AlphaEsts[2:3,s] <- hdi(Intercept - TrueAlphas[s])
        }
}

# Finally, calculate the "true" generic alpha that the model is attempting to estimate
# alpha_generic * sigma(N) = sigma(N*alpha_ij)
GenericIntercepts <- 1 - Inclusion_ij
GenericIntercepts[Focal] <- 0
TrueGenericIntercept <- log(sum(colSums(SpMatrix) * GenericIntercepts * exp(TrueAlphas)) / sum(colSums(SpMatrix) * GenericIntercepts))

# Finally, save all the necessary results for the figures
FileName <- paste("StanFits/monoLambda_constAlpha/", FilePrefix, "GraphStuff.rdata", sep = "")
save(PredVals, Growth_ppc, LambdaEsts, AlphaEsts, Inclusion_ij,
     TrueGenericIntercept, 
     file = FileName)

### Old figure code------------------
# Plot the results from the ppc
FigName <- paste("Results/monoLambda_constAlpha_test/", FilePrefix, "ppc.pdf", sep = "")
Pred_yRange <- c(-2.5,2)
Dev_yRange <- c(-1.2,1.2)
xRange <- range(Growth_ppc)
ppcCol <- "purple"
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     #par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,2))
     # Upper left: Prelim ppc estimates
     # plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
     #      ylab = "Predicted growth", las = 1)
     # points(x = Growth_ppc, y = PrelimPredVals[1,], pch = 1, col = ppcCol)
     # segments(x0 = Growth_ppc, y0 = PrelimPredVals[2,], x1 = Growth_ppc, y1 = PrelimPredVals[3,], col = ppcCol)
     # abline(a = 0, b = 1, lty = 2)
     # Message <- paste("CI width: ", round(PrelimCIwidth, digits = 2), sep = "")
     # mtext(Message, side = 3, adj = 0.1, line = -1.5, col = ppcCol)
     # mtext("Preliminary model fit", side = 3, line = 0.5)
     # Upper right: Final ppc estimates
     plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = PredVals[1,], pch = 1, col = ppcCol)
     segments(x0 = Growth_ppc, y0 = PredVals[2,], x1 = Growth_ppc, y1 = PredVals[3,], col = ppcCol)
     abline(a = 0, b = 1, lty = 2)
     #Message <- paste("CI width: ", round(FinalCIwidth, digits = 2), sep = "")
     #mtext(Message, side = 3, adj = 0.1, line = -1.5, col = ppcCol)
     mtext("Final model fit", side = 3, line = 0.5)
dev.off()
# 
# # Plot the accuracy of lambda estimates
# # Now create the lambda graph
PrelimLambdaIntercept <- PrelimPosteriors$lambdas[,1]
PrelimLambdaSlope <- PrelimPosteriors$lambdas[,2]
FinalLambdaIntercept <- FinalPosteriors$lambdas[,1]
FinalLambdaSlope <- FinalPosteriors$lambdas[,2]
TrueLambdaIntercept <- TrueVals$lambda.mean[Focal]
TrueLambdaSlope <- TrueVals$lambda.env[Focal]
PrelimInterceptDeviation <- PrelimLambdaIntercept - TrueLambdaIntercept
PrelimSlopeDeviation <- PrelimLambdaSlope - TrueLambdaSlope
FinalInterceptDeviation <- FinalLambdaIntercept - TrueLambdaIntercept
FinalSlopeDeviation <- FinalLambdaSlope - TrueLambdaSlope

LambdaCol <- "forestgreen"
InterceptRange <- c(-2, 4)
SlopeRange <- c(-0.5, 0.5)
FigName <- paste("Results/monoLambda_constAlpha/", FilePrefix, "lambdas.pdf", sep = "")
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), mar = c(5,4,2,2) + 0.1)
     # Upper left: Prelim lambda intercept
     plot(density(PrelimInterceptDeviation), main = "", xlab = "Intercept deviation",
          col = LambdaCol, xlim = InterceptRange, las = 1)
     abline(v = 0, lty = 2)
     abline(v = hdi(PrelimInterceptDeviation), lty = 3, col = LambdaCol)
     abline(v = mean(PrelimInterceptDeviation), lty = 1, col = LambdaCol)
     mtext("Preliminary model fit", side = 3, line = 1)
     # Upper right: Final lambda intercept
     plot(density(FinalInterceptDeviation), main = "", xlab = "Intercept deviation",
          col = LambdaCol, xlim = InterceptRange, las = 1)
     abline(v = 0, lty = 2)
     abline(v = hdi(FinalInterceptDeviation), lty = 3, col = LambdaCol)
     abline(v = mean(FinalInterceptDeviation), lty = 1, col = LambdaCol)
     mtext("Final model fit", side = 3, line = 1)
     # Lower left: Prelim lambda slope
     plot(density(PrelimSlopeDeviation), main = "", xlab = "Slope deviation",
          col = LambdaCol, xlim = SlopeRange, las = 1)
     abline(v = 0, lty = 2)
     abline(v = hdi(PrelimSlopeDeviation), lty = 3, col = LambdaCol)
     abline(v = mean(PrelimSlopeDeviation), lty = 1, col = LambdaCol)
     # Lower right: Final lambda slope
     plot(density(FinalSlopeDeviation), main = "", xlab = "Slope deviation",
          col = LambdaCol, xlim = SlopeRange, las = 1)
     abline(v = 0, lty = 2)
     abline(v = hdi(FinalSlopeDeviation), lty = 3, col = LambdaCol)
     abline(v = mean(FinalSlopeDeviation), lty = 1, col = LambdaCol)
# dev.off()
# 
# # Now plot the alpha estimates versus the true values
# xSeq <- 1:(S+1)
# pchSeq <- c(rep(1,S), 16)
# AlphaRange <- c(-0.02, 0.25)
# Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
# estCols <- Dark2Cols[1:2]
# FigName <- paste("Results/monoLambda_constAlpha/", FilePrefix, "alphas.pdf", sep = "")
# pdf(file = FigName, width = 10, height = 3, onefile = FALSE, paper = "special")
#      par(mar = c(5,4,2,2) + 0.1, mfrow = c(1,2), oma = c(0,2,0,0))
#      # Preliminary fit results
#      plot(x = NA, y = NA, xlim = c(1,S+1), ylim = AlphaRange,
#           xlab = "", ylab = "", las = 1)
#      points(x = xSeq, y = PrelimAlphaEsts[1,], col = estCols[1], pch = pchSeq)
#      segments(x0 = xSeq, y0 = PrelimAlphaEsts[2,], x1 = xSeq, y1 = PrelimAlphaEsts[3,], col = estCols[1])
#      abline(h = 0, lty = 2)
#      mtext("Preliminary model fit", side = 3, line = 1)
#      mtext("Alpha deviations", side = 2, line = 0.5, outer = TRUE, adj = 0.6)
#      # Final fit results
#      plot(x = NA, y = NA, xlim = c(1,S+1), ylim = AlphaRange,
#           xlab = "", ylab = "", las = 1)
#      points(x = xSeq, y = FinalAlphaEsts[1,], col = estCols[1], pch = pchSeq)
#      segments(x0 = xSeq, y0 = FinalAlphaEsts[2,], x1 = xSeq, y1 = FinalAlphaEsts[3,], col = estCols[1])
#      abline(h = 0, lty = 2)
#      mtext("Final model fit", side = 3, line = 1)
# 
#      mtext("Species", side = 1, line = -2, outer = TRUE)     
# dev.off()
# 
# # Finally, make a plot of the fixed parameter priors
# LambdaIntercept <- rnorm(n = 10000, mean = 0, sd = 1)
# LambdaSlope <- rnorm(n = 10000, mean = 0, sd = 1)
# TransformedLambdaIntercept <- exp(LambdaIntercept)
# 
# AlphaIntercept <- rnorm(n = 10000, mean = -2, sd = 0.75)
# AlphaSlope <- NULL
# TransformedAlphaIntercept <- exp(AlphaIntercept)
# 
# FigName <- paste("Results/monoLambda_constAlpha/priors.pdf", sep = "")
# pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
#      par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,3))
#      plot(density(LambdaIntercept), xlab = "Lambda intercept", ylab = "", main = "")
#      plot(density(LambdaSlope), xlab = "Lambda slope", ylab = "", main = "")
#      plot(density(TransformedLambdaIntercept), xlab = "exp(Lambda intercept)", ylab = "", main = "")
#      
#      plot(density(AlphaIntercept), xlab = "Generic alpha intercept", ylab = "", main = "")
#      plot(NA, NA, xlab = "", ylab = "", main = "", xlim = c(0,1), ylim = c(0,1), axes = FALSE)
#      plot(density(TransformedAlphaIntercept), xlab = "exp(Generic alpha intercept)", ylab = "", main = "")
# dev.off()
# 
# which(Inclusion_ij == 1)
# # Predicted non-generic competitors:
# # Focal = 8
# #    3, 6, 9, 15

