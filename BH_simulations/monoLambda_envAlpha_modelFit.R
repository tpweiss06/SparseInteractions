# This script will fit the model with a monotonic lambda and constant alphas
#    (with respect to the environment) to the corresponding data, using a series
#    of dataset sizes (10, 20, 50, 100, 200), and then examine the accuracy of
#    parameter estimates and posterior predictive checks on 300 out of sample
#    data points.



##### Confirm with Chhaya about how alpha_eij is encoded
##### Redo ppc calculations with alpha_eij
##### Or a plot with each alpha_eij as a solid line, and intervals for 
#    generic alpha nad sp.specific alpha like Cath's plots (yeah, that one)
# Also make the lambda plot match Cath's for the paper.
#    Then, those can be a single plot (2 panel; left is lambda, right is alpha)
#    Combine that with a 2 panel ppc plot (left is prelim, right is final)
# We don't need prelim(-) anymore





setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Set the current sample size and associated prefix for all graph and result
#    file names
N <- 10
max_N <- 200
FilePrefix <- paste("N", N, "_", sep = "")
AlphaCredInt <- 0.95

# Now assign the focal species and the file paths for the stan models
Focal <- 3
PrelimStanPath <- "StanCode/Prelim_monoLambda_envAlpha.stan"
FinalStanPath <- "StanCode/Final_monoLambda_envAlpha.stan"

# Load in the appropriate data
FullSim <- read.csv("Simulations/simulation_3.csv")
TrueVals <- read.csv("Simulations/parameters_3.csv")
#TrueAlphas <- exp(TrueVals$alpha.3)

# Load necessary libraries
library(rstan)
library(HDInterval)
library(RColorBrewer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# assign some universal values to be used across model fits and graphs
S <- 15
tau0 <- 1
slab_scale <- log(2)
slab_df <- 25
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
ppcCols <- Dark2Cols[1:2]
estCols <- Dark2Cols[3:4]

# Set initial values to avoid initial problems with the random number generator
ChainInitials <- list(lambdas_tilde = c(1,0), alpha_generic_tilde = 0, alpha_hat_ij_tilde = rep(0, S),
                      local_shrinkage_ij = rep(5, S), c2_tilde = 1.25, tau_tilde = 15)
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
     SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) & (time == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

# Create the data vectors to be passed to rstan for subsequent model fits
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "tau0", "slab_scale", "slab_df",
                   "N_ppc", "Nt_ppc", "SpMatrix_ppc", "env_ppc")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Inclusion_ij",
                  "N_ppc", "Nt_ppc", "SpMatrix_ppc", "env_ppc")

# Set the local values to pass to rstan
CurData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0))
Nt <- CurData$pop
env <- CurData$run.env
SpMatrix <- matrix(data = NA, nrow = N, ncol = S)
for(s in 1:S){
     SpMatrix[,s] <- subset(FullSim, (species == s) & (run <= N) & (time == 0))$pop
}
Ntp1 <- subset(FullSim, (species == Focal) & (run <= N) & (time == 1))$pop

# Now run the perliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, init = InitVals, control = list(max_treedepth = 15, adapt_delta = 0.99))
PrelimPosteriors <- extract(PrelimFit)
FitFileName <- paste("StanFits/monoLambda_constAlpha/", FilePrefix, "PrelimFit.rdata", sep = "")
save(PrelimFit, PrelimPosteriors, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
quartz()
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
traceplot(PrelimFit, pars = "lambdas")
traceplot(PrelimFit, pars = "alpha_generic")
traceplot(PrelimFit, pars = "alpha_hat_ij")
traceplot(PrelimFit, pars = c("local_shrinkage_ij", "c2", "tau"))
acf(PrelimPosteriors$alpha_generic)
########### Posterior Predictive Check
PostLength <- length(PrelimPosteriors$alpha_generic)
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
        alpha_ij[,s] <- exp(PrelimPosteriors$alpha_generic + PrelimPosteriors$alpha_hat_ij[,s])
}
# calculate the posterior distributions of lambda_ei
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        lambda_ei[,i] <- exp(PrelimPosteriors$lambdas[,1] + PrelimPosteriors$lambdas[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
PrelimGrowth_rho <- rep(NA, PostLength)
for(i in 1:PostLength){
        for(j in 1:N_ppc){
                SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,])
                Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
                Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
                Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
        }
        PrelimGrowth_rho[i] <- cor(Growth_ppc, Growth_pred[i,])
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

# Calculate the accuracy of parameter estimates from the preliminary fits
PrelimLambdaEsts <- matrix(data = NA, nrow = 3, ncol = 2)
PrelimLambdaEsts[1,1] <- mean(PrelimPosteriors$lambdas[,1])
PrelimLambdaEsts[1,2] <- mean(PrelimPosteriors$lambdas[,2])
PrelimLambdaEsts[2:3,1] <- HDInterval::hdi(PrelimPosteriors$lambdas[,1])
PrelimLambdaEsts[2:3,2] <- HDInterval::hdi(PrelimPosteriors$lambdas[,2])

# Now calculate the alpha estimates
PrelimAlphaEsts <- matrix(data = NA, nrow = 3, ncol = S)
for(s in 1:S){
     AlphaDev <- alpha_ij[,s] - as.numeric(TrueAlphas[s])
     PrelimAlphaEsts[1,s] <- median(AlphaDev)
     PrelimAlphaEsts[2:3,s] <- HDInterval::hdi(AlphaDev, credMass = AlphaCredInt)
}

# Determine the parameters that should be included and run the final model
plot(PrelimFit, pars = "alpha_hat_ij")

Inclusion_ij <- rep(0, S)
IntLevel <- 0.5
for(s in 1:S){
     Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
     if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
          Inclusion_ij[s] <- 1
     }
}
Inclusion_ij

########### Preliminary- Posterior Predictive Check 
# This will use the preliminary model fit for the ppc, but only include alpha_hat_ij values
#    indicated in the inclusion vector
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
     if(Inclusion_ij[s] == 1){
          alpha_ij[,s] <- exp(PrelimPosteriors$alpha_generic + PrelimPosteriors$alpha_hat_ij[,s])
     }else{
          alpha_ij[,s] <- exp(PrelimPosteriors$alpha_generic)
     }
     
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
PrelimMinusGrowth_rho <- rep(NA, PostLength)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,])
          Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
          Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
     }
     PrelimMinusGrowth_rho[i] <- cor(Growth_ppc, Growth_pred[i,])
}


# Calculate PrelimMinus fit results for the ppc
PrelimMinusPredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
PrelimMinusCIwidths <- rep(NA, N_ppc)
PrelimMinusDevVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
     PrelimMinusPredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
     PrelimMinusPredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
     PrelimMinusCIwidths[i] <- PrelimMinusPredVals[3,i] - PrelimMinusPredVals[2,i]
     PrelimMinusDevVals[1,i] <- mean(Growth_dev[,i], na.rm = TRUE)
     PrelimMinusDevVals[2:3,i] <- HDInterval::hdi(Growth_dev[,i])
}
PrelimMinusCIwidth <- mean(PrelimMinusCIwidths)

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals, control = list(max_treedepth = 15))
FinalPosteriors <- extract(FinalFit)
FitFileName <- paste("StanFits/monoLambda_constAlpha/", FilePrefix, "FinalFit.rdata", sep = "")
save(FinalFit, FinalPosteriors, Inclusion_ij, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
quartz()
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = "lambdas")
traceplot(FinalFit, pars = "alpha_generic")
which(Inclusion_ij == 1)
traceplot(FinalFit, pars = "alpha_hat_ij")

# Double check the autocorrelation in a few potentially suspect traceplots
acf(FinalPosteriors$alpha_hat_ij[,4])
acf(FinalPosteriors$alpha_hat_ij[,15])
acf(FinalPosteriors$alpha_generic)

########### Posterior Predictive Check
PostLength <- length(FinalPosteriors$alpha_generic)
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
     if(Inclusion_ij[s] == 1){
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_generic + FinalPosteriors$alpha_hat_ij[,s])
     }else{
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_generic)
     }
     
}
# calculate the posterior distributions of lambda_ei
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
     lambda_ei[,i] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
Growth_dev <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
FinalGrowth_rho <- rep(NA, PostLength)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,])
          Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
          Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          Growth_dev[i,j] <- Growth_pred[i,j] - Growth_ppc[j]
     }
     FinalGrowth_rho[i] <- cor(Growth_ppc, Growth_pred[i,])
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

# Calculate the accuracy of parameter estimates from the preliminary fits
FinalLambdaEsts <- matrix(data = NA, nrow = 3, ncol = 2)
FinalLambdaEsts[1,1] <- mean(FinalPosteriors$lambdas[,1])
FinalLambdaEsts[1,2] <- mean(FinalPosteriors$lambdas[,2])
FinalLambdaEsts[2:3,1] <- HDInterval::hdi(FinalPosteriors$lambdas[,1])
FinalLambdaEsts[2:3,2] <- HDInterval::hdi(FinalPosteriors$lambdas[,2])

# Now calculate the alpha estimates
FinalAlphaEsts <- matrix(data = NA, nrow = 3, ncol = S)
for(s in 1:S){
     AlphaDev <- alpha_ij[,s] - as.numeric(TrueAlphas[s])
     FinalAlphaEsts[1,s] <- median(AlphaDev)
     FinalAlphaEsts[2:3,s] <- HDInterval::hdi(AlphaDev, credMass = AlphaCredInt)
}

# Plot the results from the ppc
FigName <- paste("Results/monoLambda_constAlpha/", FilePrefix, "ppc.pdf", sep = "")
Pred_yRange <- range(PrelimPredVals, PrelimMinusPredVals, FinalPredVals) #PrelimMinusPredVals
Dev_yRange <- range(PrelimDevVals, PrelimMinusDevVals, FinalDevVals) #PrelimMinusDevVals
xRange <- range(Growth_ppc)
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,3))
     # Upper left: Prelim ppc estimates
     plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
          ylab = "Predicted growth", las = 1)
     points(x = Growth_ppc, y = PrelimPredVals[1,], pch = 1, col = ppcCols[1])
     segments(x0 = Growth_ppc, y0 = PrelimPredVals[2,], x1 = Growth_ppc, y1 = PrelimPredVals[3,], col = ppcCols[1])
     abline(a = 0, b = 1, lty = 2)
     RhoMean <- round(mean(PrelimGrowth_rho), digits = 2)
     RhoCI <- round(HDInterval::hdi(PrelimGrowth_rho), digits = 2)
     Message1 <- paste("rho = ", RhoMean, " (", RhoCI[1], " - ", RhoCI[2], ")", sep = "")
     Message2 <- paste("CI width: ", round(PrelimCIwidth, digits = 2), sep = "")
     mtext(Message1, side = 3, adj = 0.1, line = -2, col = ppcCols[1])
     mtext(Message2, side = 3, adj = 0.1, line = -3.2, col = ppcCols[1])
     mtext("Preliminary model fit", side = 3, line = 0.5)
     # Upper middle: Prelim- ppc estimates
     plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = PrelimMinusPredVals[1,], pch = 1, col = ppcCols[1])
     segments(x0 = Growth_ppc, y0 = PrelimMinusPredVals[2,], x1 = Growth_ppc, y1 = PrelimMinusPredVals[3,], col = ppcCols[1])
     abline(a = 0, b = 1, lty = 2)
     RhoMean <- round(mean(PrelimMinusGrowth_rho), digits = 2)
     RhoCI <- round(HDInterval::hdi(PrelimMinusGrowth_rho), digits = 2)
     Message1 <- paste("rho = ", RhoMean, " (", RhoCI[1], " - ", RhoCI[2], ")", sep = "")
     Message2 <- paste("CI width: ", round(PrelimMinusCIwidth, digits = 2), sep = "")
     mtext(Message1, side = 3, adj = 0.1, line = -2, col = ppcCols[1])
     mtext(Message2, side = 3, adj = 0.1, line = -3.2, col = ppcCols[1])
     mtext("Preliminary(-) model fit", side = 3, line = 0.5)
     # Upper right: Final ppc estimates
     plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = FinalPredVals[1,], pch = 1, col = ppcCols[1])
     segments(x0 = Growth_ppc, y0 = FinalPredVals[2,], x1 = Growth_ppc, y1 = FinalPredVals[3,], col = ppcCols[1])
     abline(a = 0, b = 1, lty = 2)
     RhoMean <- round(mean(FinalGrowth_rho), digits = 2)
     RhoCI <- round(HDInterval::hdi(FinalGrowth_rho), digits = 2)
     Message1 <- paste("rho = ", RhoMean, " (", RhoCI[1], " - ", RhoCI[2], ")", sep = "")
     Message2 <- paste("CI width: ", round(FinalCIwidth, digits = 2), sep = "")
     mtext(Message1, side = 3, adj = 0.1, line = -2, col = ppcCols[1])
     mtext(Message2, side = 3, adj = 0.1, line = -3.2, col = ppcCols[1])
     mtext("Final model fit", side = 3, line = 0.5)
     # Lower left: Prelim ppc deviations
     plot(x = NA, y = NA, xlim = xRange, ylim = Dev_yRange, xlab = "",
          ylab = "Deviation", las = 1)
     points(x = Growth_ppc, y = PrelimDevVals[1,], pch = 1, col = ppcCols[1])
     segments(x0 = Growth_ppc, y0 = PrelimDevVals[2,], x1 = Growth_ppc, y1 = PrelimDevVals[3,], col = ppcCols[1])
     abline(h = 0, lty = 2)
     # Lower middle: Prelim- ppc deviations
     plot(x = NA, y = NA, xlim = xRange, ylim = Dev_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = PrelimMinusDevVals[1,], pch = 1, col = ppcCols[1])
     segments(x0 = Growth_ppc, y0 = PrelimMinusDevVals[2,], x1 = Growth_ppc, y1 = PrelimMinusDevVals[3,], col = ppcCols[1])
     abline(h = 0, lty = 2)
     # Lower right: Final ppc deviations
     plot(x = NA, y = NA, xlim = xRange, ylim = Dev_yRange, xlab = "",
          ylab = "", las = 1)
     points(x = Growth_ppc, y = FinalDevVals[1,], pch = 1, col = ppcCols[1])
     segments(x0 = Growth_ppc, y0 = FinalDevVals[2,], x1 = Growth_ppc, y1 = FinalDevVals[3,], col = ppcCols[1])
     abline(h = 0, lty = 2)
     
     mtext("True growth", side = 1, line = -2, outer = TRUE)
dev.off()

# Plot the accuracy of lambda estimates
FigName <- paste("Results/monoLambda_constAlpha/", FilePrefix, "lambdas.pdf", sep = "")
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), mar = c(5,4,2,2) + 0.1)
     # Upper left: Prelim lambda intercept
     plot(density(PrelimPosteriors$lambdas[,1]), main = "", xlab = "Intercept", col = estCols[1])
     abline(v = PrelimLambdaEsts[,1], lty = c(1,2,2), col = estCols[1])
     abline(v = TrueVals$lambda.mean[Focal], col = estCols[2])
     mtext("Preliminary model fit", side = 3, line = 1)
     # Upper right: Final lambda intercept
     plot(density(FinalPosteriors$lambdas[,1]), main = "", xlab = "Intercept", col = estCols[1], ylab = "")
     abline(v = FinalLambdaEsts[,1], lty = c(1,2,2), col = estCols[1])
     abline(v = TrueVals$lambda.mean[Focal], col = estCols[2])
     mtext("Final model fit", side = 3, line = 1)
     # Lower left: Prelim lambda slope
     plot(density(PrelimPosteriors$lambdas[,2]), main = "", xlab = "Slope", col = estCols[1])
     abline(v = PrelimLambdaEsts[,2], lty = c(1,2,2), col = estCols[1])
     abline(v = TrueVals$lambda.env[Focal], col = estCols[2])
     # Lower right: Final lambda slope
     plot(density(FinalPosteriors$lambdas[,2]), main = "", xlab = "Slope", col = estCols[1], ylab = "")
     abline(v = FinalLambdaEsts[,2], lty = c(1,2,2), col = estCols[1])
     abline(v = TrueVals$lambda.env[Focal], col = estCols[2])
dev.off()

# Now plot the alpha estimates versus the true values
xSeq <- 1:S
FigName <- paste("Results/monoLambda_constAlpha/", FilePrefix, "alphas.pdf", sep = "")
pdf(file = FigName, width = 10, height = 3, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(1,2), oma = c(0,2,0,0))
     # Preliminary fit results
     plot(x = NA, y = NA, xlim = c(1,S), ylim = range(PrelimAlphaEsts),
          xlab = "", ylab = "", las = 1)
     points(x = xSeq, y = PrelimAlphaEsts[1,], col = estCols[1])
     segments(x0 = xSeq, y0 = PrelimAlphaEsts[2,], x1 = xSeq, y1 = PrelimAlphaEsts[3,], col = estCols[1])
     abline(h = 0, lty = 2)
     mtext("Preliminary model fit", side = 3, line = 1)
     mtext("Alpha deviations", side = 2, line = 0.5, outer = TRUE, adj = 0.6)
     # Final fit results
     plot(x = NA, y = NA, xlim = c(1,S), ylim = range(FinalAlphaEsts),
          xlab = "", ylab = "", las = 1)
     points(x = xSeq, y = FinalAlphaEsts[1,], col = estCols[1])
     segments(x0 = xSeq, y0 = FinalAlphaEsts[2,], x1 = xSeq, y1 = FinalAlphaEsts[3,], col = estCols[1])
     abline(h = 0, lty = 2)
     mtext("Final model fit", side = 3, line = 1)

     mtext("Species", side = 1, line = -2, outer = TRUE)     
dev.off()

FileName <- paste("StanFits/monoLambda_constAlpha/", FilePrefix, "rhos.rdata", sep = "")
save(PrelimGrowth_rho, PrelimMinusGrowth_rho, FinalGrowth_rho, file = FileName)






