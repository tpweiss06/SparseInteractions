# This script will fit the model with a monotonic lambda and constant alphas
#    (with respect to the environment) to the corresponding data, using a series
#    of dataset sizes (10, 20, 50, 100, 200), and then examine the accuracy of
#    parameter estimates and posterior predictive checks on 300 out of sample
#    data points.

setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Set the current sample size and associated prefix for all graph and result
#    file names
N <- 200
max_N <- 200
FilePrefix <- paste("N", N, "_", sep = "")

# Now assign the focal species and the file paths for the stan models
Focal <- 11
PrelimStanPath <- "StanCode/Prelim_monoLambda_envAlpha.stan"
FinalStanPath <- "StanCode/Final_monoLambda_envAlpha.stan"

# Load in the appropriate data
FullSim <- read.csv("Simulations/simulation_3.csv")
TrueVals <- read.csv("Simulations/parameters_3.csv")
TrueAlphaMeans <- TrueVals$alpha.11
TrueAlphaSlopes <- TrueVals$alpha.env.gen + TrueVals$alpha.env.spec

# Load necessary libraries
library(rstan)
library(HDInterval)
library(RColorBrewer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# assign some universal values to be used across model fits and graphs
S <- 14
OtherSpecies <- setdiff(1:15, Focal)
tau0 <- 1
slab_scale <- log(2)
slab_df <- 25

# Set initial values to avoid initial problems with the random number generator
ChainInitials <- list(lambdas_tilde = c(1,0), alpha_generic_tilde = c(0,0), alpha_hat_ij_tilde = rep(0, S),
                      local_shrinkage_ij = rep(5, S), c2_tilde = 1.25, tau_tilde = 15,
                      alpha_hat_eij_tilde = rep(0, S), local_shrinkage_eij = rep(5, S),
                      alpha_intra_tilde = c(0,0))
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# save the values for the posterior predictive checks
ppc_data <- subset(FullSim, (species == Focal) & (run > max_N) & (time == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
     SpMatrix_ppc[,s] <- subset(FullSim, (species == OtherSpecies[s]) & (run %in% ppc_runs) & (time == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

# Create the data vectors to be passed to rstan for subsequent model fits
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "tau0", "slab_scale", "slab_df")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Inclusion_ij", "Inclusion_eij")

# Set the local values to pass to rstan
CurData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0))
Nt <- CurData$pop
env <- CurData$run.env
SpMatrix <- matrix(data = NA, nrow = N, ncol = S)
for(s in 1:S){
     SpMatrix[,s] <- subset(FullSim, (species == OtherSpecies[s]) & (run <= N) & (time == 0))$pop
}
Ntp1 <- subset(FullSim, (species == Focal) & (run <= N) & (time == 1))$pop

# Now run the preliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, init = InitVals, 
                  control = list(adapt_delta = 0.99, max_treedepth = 15))
PrelimPosteriors <- extract(PrelimFit)
FitFileName <- paste("StanFits/monoLambda_envAlpha/", FilePrefix, "PrelimFit.rdata", sep = "")
save(PrelimFit, PrelimPosteriors, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
quartz()
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
traceplot(PrelimFit, pars = "lambdas")
traceplot(PrelimFit, pars = c("alpha_generic", "alpha_intra"))
traceplot(PrelimFit, pars = "alpha_hat_ij")
traceplot(PrelimFit, pars = "alpha_hat_eij")
par(mfrow = c(2,2))
acf(PrelimPosteriors$alpha_generic[,1])
acf(PrelimPosteriors$alpha_generic[,2])
acf(PrelimPosteriors$alpha_intra[,1])
acf(PrelimPosteriors$alpha_intra[,2])
PlotSamples <- sample(1:S, size = 4, replace = FALSE)
par(mfrow = c(2,2))
for(i in 1:4){
        acf(PrelimPosteriors$alpha_hat_ij[,PlotSamples[i]])
}
for(i in 1:4){
        acf(PrelimPosteriors$alpha_hat_eij[,PlotSamples[i]])
}

########### Posterior Predictive Check
PostLength <- length(PrelimPosteriors$alpha_generic[,1])
# calculate the posterior distributions of the interaction coefficients and lambda
Prelim_alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
Prelim_alpha_intra <- matrix(NA, nrow = PostLength, ncol = N_ppc)
Prelim_lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        Prelim_lambda_ei[,i] <- exp(PrelimPosteriors$lambdas[,1] + PrelimPosteriors$lambdas[,2]*env_ppc[i])
        for(s in 1:S){
                Prelim_alpha_eij[,i,s] <- exp(PrelimPosteriors$alpha_generic[,1] + PrelimPosteriors$alpha_hat_ij[,s] +
                                 (PrelimPosteriors$alpha_generic[,2] + PrelimPosteriors$alpha_hat_eij[,s])*env_ppc[i])
        }
        Prelim_alpha_intra[,i] <- exp(PrelimPosteriors$alpha_intra[,1] + PrelimPosteriors$alpha_intra[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
        for(j in 1:N_ppc){
                SigmaTerm <- sum(Prelim_alpha_eij[i,j,] * SpMatrix_ppc[j,]) + Prelim_alpha_intra[i,j] * Nt_ppc[j]
                Ntp1_pred <- Nt_ppc[j] * Prelim_lambda_ei[i,j] / (1 + SigmaTerm)
                Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
                Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
        }
}

# Calculate preliminary fit results for the ppc
PrelimPredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
PrelimCIwidths <- rep(NA, N_ppc)
PrelimDevVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
     PrelimPredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
     PrelimPredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
     PrelimCIwidths[i] <- PrelimPredVals[3,i] - PrelimPredVals[2,i]
     PrelimDevVals[1,i] <- mean(Growth_dev[,i], na.rm = TRUE)
     PrelimDevVals[2:3,i] <- HDInterval::hdi(Growth_dev[,i])
}
PrelimCIwidth <- mean(PrelimCIwidths)

# Determine the parameters that should be included and run the final model
plot(PrelimFit, pars = "alpha_hat_ij")
plot(PrelimFit, pars = "alpha_hat_eij")

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
Inclusion_ij
Inclusion_eij

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals, 
                 control = list(adapt_delta = 0.99, max_treedepth = 15))
FinalPosteriors <- extract(FinalFit)
FitFileName <- paste("StanFits/monoLambda_envAlpha/", FilePrefix, "FinalFit.rdata", sep = "")
save(FinalFit, FinalPosteriors, Inclusion_ij, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
quartz()
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = "lambdas")
traceplot(FinalFit, pars = c("alpha_generic", "alpha_intra"))
which(Inclusion_ij == 1)
traceplot(FinalFit, pars = "alpha_hat_ij")
which(Inclusion_eij == 1)
traceplot(FinalFit, pars = "alpha_hat_eij")

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
        if(Inclusion_eij[s] == 1){
                quartz()
                acf(FinalPosteriors$alpha_hat_eij[,s])
        }
}

########### Posterior Predictive Check
PostLength <- length(FinalPosteriors$alpha_generic[,1])
# calculate the posterior distributions of the interaction coefficients anad lambdas
Final_alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
Final_alpha_intra <- matrix(NA, nrow = PostLength, ncol = N_ppc)
Final_lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        Final_lambda_ei[,i] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[i])
        for(s in 1:S){
                if(Inclusion_ij[s] == 1){
                        if(Inclusion_eij[s] == 1){
                                Final_alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        FinalPosteriors$alpha_hat_ij[,s] +
                                        (FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s])*env_ppc[i])
                        }else{
                                Final_alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        FinalPosteriors$alpha_hat_ij[,s] +
                                        FinalPosteriors$alpha_generic[,2]*env_ppc[i])
                        }
                }else{
                        if(Inclusion_eij[s] == 1){
                                Final_alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        (FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s])*env_ppc[i])
                        }else{
                                Final_alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        FinalPosteriors$alpha_generic[,2]*env_ppc[i])
                        }
                }
        }
        Final_alpha_intra[,i] <- exp(FinalPosteriors$alpha_intra[,1] + FinalPosteriors$alpha_intra[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(Final_alpha_eij[i,j,] * SpMatrix_ppc[j,]) + Final_alpha_intra[i,j] * Nt_ppc[j]
          Ntp1_pred <- Nt_ppc[j] * Final_lambda_ei[i,j] / (1 + SigmaTerm)
          Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
     }
}

# Calculate final fit results for the ppc
FinalPredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
FinalCIwidths <- rep(NA, N_ppc)
FinalDevVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
     FinalPredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
     FinalPredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
     FinalCIwidths[i] <- FinalPredVals[3,i] - FinalPredVals[2,i]
     FinalDevVals[1,i] <- mean(Growth_dev[,i], na.rm = TRUE)
     FinalDevVals[2:3,i] <- HDInterval::hdi(Growth_dev[,i])
}
FinalCIwidth <- mean(FinalCIwidths)

# Plot the results from the ppc
FigName <- paste("Results/monoLambda_envAlpha/", FilePrefix, "ppc.pdf", sep = "")
Pred_yRange <- c(-3.5,2.75)
Dev_yRange <- c(-1.2,1.2)
xRange <- range(Growth_ppc)
ppcCol <- "purple"
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,2))
     # Upper left: Prelim ppc estimates
     plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
          ylab = "Predicted growth", las = 1)
     points(x = Growth_ppc, y = PrelimPredVals[1,], pch = 1, col = ppcCol)
     segments(x0 = Growth_ppc, y0 = PrelimPredVals[2,], x1 = Growth_ppc, y1 = PrelimPredVals[3,], col = ppcCol)
     abline(a = 0, b = 1, lty = 2)
     Message <- paste("CI width: ", round(PrelimCIwidth, digits = 2), sep = "")
     mtext(Message, side = 3, adj = 0.1, line = -3.2, col = ppcCol)
     mtext("Preliminary model fit", side = 3, line = 0.5)
     # Upper right: Final ppc estimates
     plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = FinalPredVals[1,], pch = 1, col = ppcCol)
     segments(x0 = Growth_ppc, y0 = FinalPredVals[2,], x1 = Growth_ppc, y1 = FinalPredVals[3,], col = ppcCol)
     abline(a = 0, b = 1, lty = 2)
     Message <- paste("CI width: ", round(FinalCIwidth, digits = 2), sep = "")
     mtext(Message, side = 3, adj = 0.1, line = -3.2, col = ppcCol)
     mtext("Final model fit", side = 3, line = 0.5)
     # Lower left: Prelim ppc deviations
     plot(x = NA, y = NA, xlim = xRange, ylim = Dev_yRange, xlab = "",
          ylab = "Deviation", las = 1)
     points(x = Growth_ppc, y = PrelimDevVals[1,], pch = 1, col = ppcCol)
     segments(x0 = Growth_ppc, y0 = PrelimDevVals[2,], x1 = Growth_ppc, y1 = PrelimDevVals[3,], col = ppcCol)
     abline(h = 0, lty = 2)
     # Lower right: Final ppc deviations
     plot(x = NA, y = NA, xlim = xRange, ylim = Dev_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = FinalDevVals[1,], pch = 1, col = ppcCol)
     segments(x0 = Growth_ppc, y0 = FinalDevVals[2,], x1 = Growth_ppc, y1 = FinalDevVals[3,], col = ppcCol)
     abline(h = 0, lty = 2)
     
     mtext("True growth", side = 1, line = -2, outer = TRUE)
dev.off()

# Now create the lambda graph
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

InterceptRange <- c(-2, 4)
SlopeRange <- c(-1.5, 1.5)
LambdaCol <- "forestgreen"
FigName <- paste("Results/monoLambda_envAlpha/", FilePrefix, "lambdas.pdf", sep = "")
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
dev.off()

# Now create the alpha graph
PrelimInterceptDevs <- matrix(data = NA, nrow = 3, ncol = S+1)
PrelimSlopeDevs <- matrix(data = NA, nrow = 3, ncol = S+1)
FinalInterceptDevs <- matrix(data = NA, nrow = 3, ncol = S+1)
FinalSlopeDevs <- matrix(data = NA, nrow = 3, ncol = S+1)
for(s in 1:S){
     PrelimIntercepts <- PrelimPosteriors$alpha_generic[,1] + PrelimPosteriors$alpha_hat_ij[,s]
     PrelimSlopes <- PrelimPosteriors$alpha_generic[,2] + PrelimPosteriors$alpha_hat_ij[,s]
     if(Inclusion_ij[s] == 1){
          FinalIntercepts <- FinalPosteriors$alpha_generic[,1] + FinalPosteriors$alpha_hat_ij[,s]
     }else{
          FinalIntercepts <- FinalPosteriors$alpha_generic[,1]
     }
     if(Inclusion_eij[s] == 1){
          FinalSlopes <- FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s]
     }else{
          FinalSlopes <- FinalPosteriors$alpha_generic[,2]
     }
     PrelimInterceptDev <- PrelimIntercepts - TrueAlphaMeans[OtherSpecies[s]]
     PrelimSlopeDev <- PrelimSlopes - TrueAlphaSlopes[OtherSpecies[s]]
     FinalInterceptDev <- FinalIntercepts - TrueAlphaMeans[OtherSpecies[s]]
     FinalSlopeDev <- FinalSlopes - TrueAlphaSlopes[OtherSpecies[s]]
     
     PrelimInterceptDevs[1,s] <- mean(PrelimInterceptDev)
     PrelimInterceptDevs[2:3,s] <- hdi(PrelimInterceptDev)
     PrelimSlopeDevs[1,s] <- mean(PrelimSlopeDev)
     PrelimSlopeDevs[2:3,s] <- hdi(PrelimSlopeDev)
     FinalInterceptDevs[1,s] <- mean(FinalInterceptDev)
     FinalInterceptDevs[2:3,s] <- hdi(FinalInterceptDev)
     FinalSlopeDevs[1,s] <- mean(FinalSlopeDev)
     FinalSlopeDevs[2:3,s] <- hdi(FinalSlopeDev)
}
# Now calculate the intra values
PrelimInterceptDev <- PrelimPosteriors$alpha_intra[,1] - TrueAlphaMeans[Focal]
PrelimSlopeDev <- PrelimPosteriors$alpha_intra[,2] - TrueAlphaSlopes[Focal]
FinalInterceptDev <- FinalPosteriors$alpha_intra[,1] - TrueAlphaMeans[Focal]
FinalSlopeDev <- FinalPosteriors$alpha_intra[,2] - TrueAlphaSlopes[Focal]
PrelimInterceptDevs[1,S+1] <- mean(PrelimInterceptDev)
PrelimInterceptDevs[2:3,S+1] <- hdi(PrelimInterceptDev)
PrelimSlopeDevs[1,S+1] <- mean(PrelimSlopeDev)
PrelimSlopeDevs[2:3,S+1] <- hdi(PrelimSlopeDev)
FinalInterceptDevs[1,S+1] <- mean(FinalInterceptDev)
FinalInterceptDevs[2:3,S+1] <- hdi(FinalInterceptDev)
FinalSlopeDevs[1,S+1] <- mean(FinalSlopeDev)
FinalSlopeDevs[2:3,S+1] <- hdi(FinalSlopeDev)

InterceptSeq <- 1:(S+1) - 0.15
SlopeSeq <- 1:(S+1) + 0.15
pchSeq <- c(rep(1,S), 16)
xRange <- c(0.5, S+1.5)
yRange <- c(-2, 5.5)
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
estCols <- Dark2Cols[1:2]
FigName <- paste("Results/monoLambda_envAlpha/", FilePrefix, "alphas.pdf", sep = "")
pdf(file = FigName, width = 10, height = 3, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(1,2), oma = c(0,2,0,0))
     # Preliminary fit results
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange,
          xlab = "", ylab = "", las = 1, xaxt = "n")
     axis(1, at = 1:S, labels = FALSE, tcl = -0.25)
     axis(1, at = seq(3, 15, by = 3))
     points(x = InterceptSeq, y = PrelimInterceptDevs[1,], col = estCols[1], pch = pchSeq)
     points(x = SlopeSeq, y = PrelimSlopeDevs[1,], col = estCols[2], pch = pchSeq)
     segments(x0 = InterceptSeq, y0 = PrelimInterceptDevs[2,], x1 = InterceptSeq,
              y1 = PrelimInterceptDevs[3,], col = estCols[1])
     segments(x0 = SlopeSeq, y0 = PrelimSlopeDevs[2,], x1 = SlopeSeq,
              y1 = PrelimSlopeDevs[3,], col = estCols[2])
     abline(h = 0, lty = 2)
     mtext("Preliminary model fit", side = 3, line = 1)
     mtext("Alpha deviations", side = 2, line = -1, outer = TRUE, adj = 0.6)
     # Final fit results
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange,
          xlab = "", ylab = "", las = 1, xaxt = "n")
     axis(1, at = 1:S, labels = FALSE, tcl = -0.25)
     axis(1, at = seq(3, 15, by = 3))
     points(x = InterceptSeq, y = FinalInterceptDevs[1,], col = estCols[1], pch = pchSeq)
     points(x = SlopeSeq, y = FinalSlopeDevs[1,], col = estCols[2], pch = pchSeq)
     segments(x0 = InterceptSeq, y0 = FinalInterceptDevs[2,], x1 = InterceptSeq,
              y1 = FinalInterceptDevs[3,], col = estCols[1])
     segments(x0 = SlopeSeq, y0 = FinalSlopeDevs[2,], x1 = SlopeSeq,
              y1 = FinalSlopeDevs[3,], col = estCols[2])
     abline(h = 0, lty = 2)
     legend(x = "bottomleft", legend = c("Intercept", "Slope"), pch = 1, col = estCols,
            bty = "n")
     mtext("Final model fit", side = 3, line = 1)

     mtext("Species", side = 1, line = -2, outer = TRUE)     
dev.off()

# Finally, make a plot of the fixed parameter priors
LambdaIntercept <- rnorm(n = 10000, mean = 0, sd = 1)
LambdaSlope <- rnorm(n = 10000, mean = 0, sd = 1)
TransformedLambdaIntercept <- exp(LambdaIntercept)

AlphaIntercept <- rnorm(n = 10000, mean = -2, sd = 0.75)
AlphaSlope <- rnorm(n = 10000, mean = 0, sd = 0.5)
TransformedAlphaIntercept <- exp(AlphaIntercept)

FigName <- paste("Results/monoLambda_envAlpha/priors.pdf", sep = "")
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,3))
     plot(density(LambdaIntercept), xlab = "Lambda intercept", ylab = "", main = "")
     plot(density(LambdaSlope), xlab = "Lambda slope", ylab = "", main = "")
     plot(density(TransformedLambdaIntercept), xlab = "exp(Lambda intercept)", ylab = "", main = "")

     plot(density(AlphaIntercept), xlab = "Generic alpha intercept", ylab = "", main = "")
     plot(density(AlphaSlope), xlab = "Generic alpha slope", ylab = "", main = "")
     plot(density(TransformedAlphaIntercept), xlab = "exp(Generic alpha intercept)", ylab = "", main = "")
dev.off()

which(Inclusion_ij == 1)
# 5, 7, 10, 13 (will be 12 here)
which(Inclusion_eij == 1)
# 1, 8, 9, 10

