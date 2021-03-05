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
Focal <- 13
PrelimStanPath <- "StanCode/Prelim_optLambda_envAlpha.stan"
FinalStanPath <- "StanCode/Final_optLambda_envAlpha.stan"

# Load in the appropriate data
FullSim <- read.csv("Simulations/simulation_4.csv")
TrueVals <- read.csv("Simulations/parameters_4.csv")
TrueAlphaMeans <- TrueVals$alpha.13
TrueAlphaSlopes <- TrueVals$alpha.env.gen + TrueVals$alpha.env.spec
TrueLambdaOpt <- TrueVals$z.env[Focal]
TrueLambdaMax <- TrueVals$lambda.max[Focal]
TrueLambdaWidth <- TrueVals$sigma.env[Focal]

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
ChainInitials <- list(lambda_max = 1, lambda_opt = 0, lambda_width = 0.5, alpha_generic_tilde = c(0,0), 
                      alpha_hat_ij_tilde = rep(0, S), local_shrinkage_ij = rep(5, S), c2_tilde = 1.25, 
                      tau_tilde = 15, alpha_hat_eij_tilde = rep(0, S), local_shrinkage_eij = rep(5, S),
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
                  chains = 3, init = InitVals)
PrelimPosteriors <- extract(PrelimFit)
FitFileName <- paste("StanFits/optLambda_envAlpha/", FilePrefix, "PrelimFit.rdata", sep = "")
save(PrelimFit, PrelimPosteriors, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
quartz()
pairs(PrelimFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
traceplot(PrelimFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
traceplot(PrelimFit, pars = "alpha_hat_ij")
traceplot(PrelimFit, pars = "alpha_hat_eij")
acf(PrelimPosteriors$lambda_opt)
acf(PrelimPosteriors$lambda_max)
acf(PrelimPosteriors$lambda_width)
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
alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
alpha_intra <- matrix(NA, nrow = PostLength, ncol = N_ppc)
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        lambda_ei[,i] <- PrelimPosteriors$lambda_max * 
                exp(-1 * ( (PrelimPosteriors$lambda_opt - env_ppc[i]) / (2 * PrelimPosteriors$lambda_width) )^2)
        for(s in 1:S){
                alpha_eij[,i,s] <- exp(PrelimPosteriors$alpha_generic[,1] + PrelimPosteriors$alpha_hat_ij[,s] +
                                 (PrelimPosteriors$alpha_generic[,2] + PrelimPosteriors$alpha_hat_eij[,s])*env_ppc[i])
        }
        alpha_intra[,i] <- exp(PrelimPosteriors$alpha_intra[,1] + PrelimPosteriors$alpha_intra[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
        for(j in 1:N_ppc){
                SigmaTerm <- sum(alpha_eij[i,j,] * SpMatrix_ppc[j,]) + alpha_intra[i,j] * Nt_ppc[j]
                Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
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

# Calculate the deviation in the lambda estimates
PrelimLambdaDevs <- matrix(data = NA, nrow = PostLength, ncol = 3)
PrelimLambdaDevs[,1] <- PrelimPosteriors$lambda_opt - TrueLambdaOpt
PrelimLambdaDevs[,2] <- PrelimPosteriors$lambda_max - TrueLambdaMax
PrelimLambdaDevs[,3] <- PrelimPosteriors$lambda_width - TrueLambdaWidth

# Now calculate the alpha estimates
PrelimAlphaEsts <- array(NA, dim = c(3,2,S+1)) 
for(s in 1:S){
        Intercept <- PrelimPosteriors$alpha_generic[,1] + PrelimPosteriors$alpha_hat_ij[,s]
        Slope <-  PrelimPosteriors$alpha_generic[,2] + PrelimPosteriors$alpha_hat_eij[,s]
        InterceptDev <- Intercept - TrueAlphaMeans[OtherSpecies[s]]
        SlopeDev <- Slope - TrueAlphaSlopes[OtherSpecies[s]]
        
        PrelimAlphaEsts[1,1,s] <- mean(InterceptDev)
        PrelimAlphaEsts[2:3,1,s] <- HDInterval::hdi(InterceptDev)
        PrelimAlphaEsts[1,2,s] <- mean(SlopeDev)
        PrelimAlphaEsts[2:3,2,s] <- HDInterval::hdi(SlopeDev)
}
# Now the intraspefic values
Intercept <- PrelimPosteriors$alpha_intra[,1]
Slope <-  PrelimPosteriors$alpha_intra[,2]
InterceptDev <- Intercept - TrueAlphaMeans[Focal]
SlopeDev <- Slope - TrueAlphaSlopes[Focal]

PrelimAlphaEsts[1,1,S+1] <- mean(InterceptDev)
PrelimAlphaEsts[2:3,1,S+1] <- HDInterval::hdi(InterceptDev)
PrelimAlphaEsts[1,2,S+1] <- mean(SlopeDev)
PrelimAlphaEsts[2:3,2,S+1] <- HDInterval::hdi(SlopeDev)

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

# Reset initial conditions to allow faster fitting
ChainInitials <- list(lambda_max = mean(PrelimPosteriors$lambda_max), lambda_opt = mean(PrelimPosteriors$lambda_opt), 
                      lambda_width = mean(PrelimPosteriors$lambda_width), alpha_generic_tilde = colMeans(PrelimPosteriors$alpha_generic_tilde), 
                      alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
                      c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
                      alpha_hat_eij_tilde = colMeans(PrelimPosteriors$alpha_hat_eij_tilde), local_shrinkage_eij = colMeans(PrelimPosteriors$local_shrinkage_eij),
                      alpha_intra_tilde = colMeans(PrelimPosteriors$alpha_intra_tilde))
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals)
FinalPosteriors <- extract(FinalFit)
FitFileName <- paste("StanFits/optLambda_envAlpha/", FilePrefix, "FinalFit.rdata", sep = "")
save(FinalFit, FinalPosteriors, Inclusion_ij, Inclusion_eij, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
quartz()
pairs(FinalFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = c("lambda_opt", "lambda_max", "lambda_width", "alpha_generic", "alpha_intra"))
which(Inclusion_ij == 1)
traceplot(FinalFit, pars = "alpha_hat_ij")
which(Inclusion_eij == 1)
traceplot(FinalFit, pars = "alpha_hat_eij")
acf(FinalPosteriors$lambda_opt)
acf(FinalPosteriors$lambda_max)
acf(FinalPosteriors$lambda_width)
par(mfrow = c(2,2))
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
alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
alpha_intra <- matrix(NA, nrow = PostLength, ncol = N_ppc)
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        lambda_ei[,i] <- FinalPosteriors$lambda_max * 
                exp(-1 * ( (FinalPosteriors$lambda_opt - env_ppc[i]) / (2 * FinalPosteriors$lambda_width) )^2)
        for(s in 1:S){
                if(Inclusion_ij[s] == 1){
                        if(Inclusion_eij[s] == 1){
                                alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        FinalPosteriors$alpha_hat_ij[,s] +
                                        (FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s])*env_ppc[i])
                        }else{
                                alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        FinalPosteriors$alpha_hat_ij[,s] +
                                        FinalPosteriors$alpha_generic[,2]*env_ppc[i])
                        }
                }else{
                        if(Inclusion_eij[s] == 1){
                                alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        (FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s])*env_ppc[i])
                        }else{
                                alpha_eij[,i,s] <- exp(FinalPosteriors$alpha_generic[,1] + 
                                        FinalPosteriors$alpha_generic[,2]*env_ppc[i])
                        }
                }
        }
        alpha_intra[,i] <- exp(FinalPosteriors$alpha_intra[,1] + FinalPosteriors$alpha_intra[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_eij[i,j,] * SpMatrix_ppc[j,]) + alpha_intra[i,j] * Nt_ppc[j]
          Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
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

# Calculate the deviation in the lambda estimates
FinalLambdaDevs <- matrix(data = NA, nrow = PostLength, ncol = 3)
FinalLambdaDevs[,1] <- FinalPosteriors$lambda_opt - TrueLambdaOpt
FinalLambdaDevs[,2] <- FinalPosteriors$lambda_max - TrueLambdaMax
FinalLambdaDevs[,3] <- FinalPosteriors$lambda_width - TrueLambdaWidth

# Now calculate the alpha estimates
FinalAlphaEsts <- array(NA, dim = c(3,2,S+1)) 
for(s in 1:S){
        if(Inclusion_ij[s] == 1){
                Intercept <- FinalPosteriors$alpha_generic[,1] + FinalPosteriors$alpha_hat_ij[,s]
        }else{
                Intercept <- FinalPosteriors$alpha_generic[,1]
        }
        if(Inclusion_eij == 1){
                Slope <-  FinalPosteriors$alpha_generic[,2] + FinalPosteriors$alpha_hat_eij[,s]
        }else{
                Slope <-  FinalPosteriors$alpha_generic[,2]
        }
        InterceptDev <- Intercept - TrueAlphaMeans[OtherSpecies[s]]
        SlopeDev <- Slope - TrueAlphaSlopes[OtherSpecies[s]]
        
        FinalAlphaEsts[1,1,s] <- mean(InterceptDev)
        FinalAlphaEsts[2:3,1,s] <- HDInterval::hdi(InterceptDev)
        FinalAlphaEsts[1,2,s] <- mean(SlopeDev)
        FinalAlphaEsts[2:3,2,s] <- HDInterval::hdi(SlopeDev)
}
# Now the intraspefic values
Intercept <- FinalPosteriors$alpha_intra[,1]
Slope <-  FinalPosteriors$alpha_intra[,2]
InterceptDev <- Intercept - TrueAlphaMeans[Focal]
SlopeDev <- Slope - TrueAlphaSlopes[Focal]

FinalAlphaEsts[1,1,S+1] <- mean(InterceptDev)
FinalAlphaEsts[2:3,1,S+1] <- HDInterval::hdi(InterceptDev)
FinalAlphaEsts[1,2,S+1] <- mean(SlopeDev)
FinalAlphaEsts[2:3,2,S+1] <- HDInterval::hdi(SlopeDev)

# Plot the results from the ppc
FigName <- paste("Results/optLambda_envAlpha/", FilePrefix, "ppc.pdf", sep = "")
Pred_yRange <- c(-6,1.5)
Dev_yRange <- c(-2,3)
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
LambdaCol <- "forestgreen"
OptimumRange <- c(-0.5, 0.5)
MaxRange <- c(-5, 25)
WidthRange <- c(-0.2, 0.2)
FigName <- paste("Results/optLambda_envAlpha/", FilePrefix, "lambdas.pdf", sep = "")
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
        par(mfrow = c(2,3), mar = c(5,4,2,2) + 0.1, oma = c(0, 2, 0, 0))
        # Upper left: Prelim lambda optimum
        plot(density(PrelimLambdaDevs[,1]), main = "", xlab = "", 
             col = LambdaCol, xlim = OptimumRange, las = 1)
        abline(v = 0, lty = 2)
        abline(v = hdi(PrelimLambdaDevs[,1]), lty = 3, col = LambdaCol)
        abline(v = mean(PrelimLambdaDevs[,1]), lty = 1, col = LambdaCol)
        mtext("Preliminary model fit", side = 2, line = 4.5)
        # Upper middle: Prelim lambda maximum
        plot(density(PrelimLambdaDevs[,2]), main = "", xlab = "", 
             col = LambdaCol, xlim = MaxRange, las = 1)
        abline(v = 0, lty = 2)
        abline(v = hdi(PrelimLambdaDevs[,2]), lty = 3, col = LambdaCol)
        abline(v = mean(PrelimLambdaDevs[,2]), lty = 1, col = LambdaCol)
        # Upper left: Prelim lambda width
        plot(density(PrelimLambdaDevs[,3]), main = "", xlab = "", 
             col = LambdaCol, xlim = WidthRange, las = 1)
        abline(v = 0, lty = 2)
        abline(v = hdi(PrelimLambdaDevs[,3]), lty = 3, col = LambdaCol)
        abline(v = mean(PrelimLambdaDevs[,3]), lty = 1, col = LambdaCol)
        # Lower left: Final lambda optimum
        plot(density(FinalLambdaDevs[,1]), main = "", xlab = "Optimum deviation", 
             col = LambdaCol, xlim = OptimumRange, las = 1)
        abline(v = 0, lty = 2)
        abline(v = hdi(FinalLambdaDevs[,1]), lty = 3, col = LambdaCol)
        abline(v = mean(FinalLambdaDevs[,1]), lty = 1, col = LambdaCol)
        mtext("Final model fit", side = 2, line = 4.5)
        # Lower middle: Final lambda maximum
        plot(density(FinalLambdaDevs[,2]), main = "", xlab = "Maximum deviation", 
             col = LambdaCol, xlim = MaxRange, las = 1)
        abline(v = 0, lty = 2)
        abline(v = hdi(FinalLambdaDevs[,2]), lty = 3, col = LambdaCol)
        abline(v = mean(FinalLambdaDevs[,2]), lty = 1, col = LambdaCol)
        # Lower left: Final lambda width
        plot(density(FinalLambdaDevs[,3]), main = "", xlab = "Width deviation", 
             col = LambdaCol, xlim = WidthRange, las = 1)
        abline(v = 0, lty = 2)
        abline(v = hdi(FinalLambdaDevs[,3]), lty = 3, col = LambdaCol)
        abline(v = mean(FinalLambdaDevs[,3]), lty = 1, col = LambdaCol)
dev.off()


# Now create the alpha graph
InterceptSeq <- 1:(S+1) - 0.15
SlopeSeq <- 1:(S+1) + 0.15
pchSeq <- c(rep(1,S), 16)
xRange <- c(0.5, S+1.5)
yRange <- c(-2, 5.5)
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
estCols <- Dark2Cols[1:2]
FigName <- paste("Results/optLambda_envAlpha/", FilePrefix, "alphas.pdf", sep = "")
pdf(file = FigName, width = 10, height = 3, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(1,2), oma = c(0,2,0,0))
     # Preliminary fit results
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange,
          xlab = "", ylab = "", las = 1, xaxt = "n")
     axis(1, at = 1:S, labels = FALSE, tcl = -0.25)
     axis(1, at = seq(3, 15, by = 3))
     points(x = InterceptSeq, y = PrelimAlphaEsts[1,1,], col = estCols[1], pch = pchSeq)
     points(x = SlopeSeq, y = PrelimAlphaEsts[1,2,], col = estCols[2], pch = pchSeq)
     segments(x0 = InterceptSeq, y0 = PrelimAlphaEsts[2,1,], x1 = InterceptSeq,
              y1 = PrelimAlphaEsts[3,1,], col = estCols[1])
     segments(x0 = SlopeSeq, y0 = PrelimAlphaEsts[2,2,], x1 = SlopeSeq,
              y1 = PrelimAlphaEsts[3,2,], col = estCols[2])
     abline(h = 0, lty = 2)
     mtext("Preliminary model fit", side = 3, line = 1)
     mtext("Alpha deviations", side = 2, line = -1, outer = TRUE, adj = 0.6)
     # Final fit results
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange,
          xlab = "", ylab = "", las = 1, xaxt = "n")
     axis(1, at = 1:S, labels = FALSE, tcl = -0.25)
     axis(1, at = seq(3, 15, by = 3))
     points(x = InterceptSeq, y = FinalAlphaEsts[1,1,], col = estCols[1], pch = pchSeq)
     points(x = SlopeSeq, y = FinalAlphaEsts[1,2,], col = estCols[2], pch = pchSeq)
     segments(x0 = InterceptSeq, y0 = FinalAlphaEsts[2,1,], x1 = InterceptSeq,
              y1 = FinalAlphaEsts[3,1,], col = estCols[1])
     segments(x0 = SlopeSeq, y0 = FinalAlphaEsts[2,2,], x1 = SlopeSeq,
              y1 = FinalAlphaEsts[3,2,], col = estCols[2])
     abline(h = 0, lty = 2)
     legend(x = "bottomleft", legend = c("Intercept", "Slope"), pch = 1, col = estCols,
            bty = "n")
     mtext("Final model fit", side = 3, line = 1)

     mtext("Species", side = 1, line = -2, outer = TRUE)     
dev.off()

# Finally, make a plot of the fixed parameter priors
LambdaOpt <- rnorm(n = 10000, mean = 0, sd = 1)
LambdaMax <- rnorm(n = 10000, mean = 0, sd = 7.5)
LambdaWidth <- rnorm(n = 10000, mean = 0, sd = 1)

AlphaIntercept <- rnorm(n = 10000, mean = -2, sd = 0.75)
AlphaSlope <- rnorm(n = 10000, mean = 0, sd = 0.5)
TransformedAlphaIntercept <- exp(AlphaIntercept)

FigName <- paste("Results/optLambda_envAlpha/priors.pdf", sep = "")
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,3))
     plot(density(LambdaOpt), xlab = "Lambda optimum", ylab = "", main = "")
     plot(density(LambdaMax), xlab = "Lambda maximum", ylab = "", main = "", xlim = c(0, 30))
     plot(density(LambdaWidth), xlab = "Lambda width", ylab = "", main = "", xlim = c(0, 4))
     
     plot(density(AlphaIntercept), xlab = "Generic alpha intercept", ylab = "", main = "")
     plot(density(AlphaSlope), xlab = "Generic alpha slope", ylab = "", main = "")
     plot(density(TransformedAlphaIntercept), xlab = "exp(Generic alpha intercept)", ylab = "", main = "")
dev.off()

which(Inclusion_ij == 1)
# 3, 7, 8, 14
which(Inclusion_eij == 1)
# 4, 10, 12, 13

