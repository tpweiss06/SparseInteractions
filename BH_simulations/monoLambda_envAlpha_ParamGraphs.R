# This script will create new parameter graphs for the monoLambda_envAlpha model

setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Set the current sample size and associated prefix for all graph names
N <- 100
FilePrefix <- paste("N", N, "_", sep = "")

# Now assign the focal species and load the data
Focal <- 11
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

# Set data for the current graphs
S <- 15
CurData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0))
Nt <- CurData$pop
env <- CurData$run.env
SpMatrix <- matrix(data = NA, nrow = N, ncol = S)
for(s in 1:S){
     SpMatrix[,s] <- subset(FullSim, (species == s) & (run <= N) & (time == 0))$pop
}
Ntp1 <- subset(FullSim, (species == Focal) & (run <= N) & (time == 1))$pop

# Now load the preliminary model fit
FitFileName <- paste("StanFits/monoLambda_envAlpha/", FilePrefix, "PrelimFit.rdata", sep = "")
load(FitFileName)

# Determine the parameters that should be included and run the final model
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
FitFileName <- paste("StanFits/monoLambda_envAlpha/", FilePrefix, "FinalFit.rdata", sep = "")
load(FitFileName)

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

InterceptRange <- c(-0.75, 0.75)
SlopeRange <- c(-0.5, 0.5)
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
PrelimInterceptDevs <- matrix(data = NA, nrow = 3, ncol = S)
PrelimSlopeDevs <- matrix(data = NA, nrow = 3, ncol = S)
FinalInterceptDevs <- matrix(data = NA, nrow = 3, ncol = S)
FinalSlopeDevs <- matrix(data = NA, nrow = 3, ncol = S)
for(s in 1:S){
     PrelimIntercepts <- PrelimPosteriors$alphas[,1] + PrelimPosteriors$alpha_hat_ij[,s]
     PrelimSlopes <- PrelimPosteriors$alphas[,2] + PrelimPosteriors$alpha_hat_ij[,s]
     if(Inclusion_ij[s] == 1){
          FinalIntercepts <- FinalPosteriors$alphas[,1] + FinalPosteriors$alpha_hat_ij[,s]
     }else{
          FinalIntercepts <- FinalPosteriors$alphas[,1]
     }
     if(Inclusion_eij[s] == 1){
          FinalSlopes <- FinalPosteriors$alphas[,2] + FinalPosteriors$alpha_hat_eij[,s]
     }else{
          FinalSlopes <- FinalPosteriors$alphas[,2]
     }
     PrelimInterceptDev <- PrelimIntercepts - TrueAlphaMeans[s]
     PrelimSlopeDev <- PrelimSlopes - TrueAlphaSlopes[s]
     FinalInterceptDev <- FinalIntercepts - TrueAlphaMeans[s]
     FinalSlopeDev <- FinalSlopes - TrueAlphaSlopes[s]
     
     PrelimInterceptDevs[1,s] <- mean(PrelimInterceptDev)
     PrelimInterceptDevs[2:3,s] <- hdi(PrelimInterceptDev)
     PrelimSlopeDevs[1,s] <- mean(PrelimSlopeDev)
     PrelimSlopeDevs[2:3,s] <- hdi(PrelimSlopeDev)
     FinalInterceptDevs[1,s] <- mean(FinalInterceptDev)
     FinalInterceptDevs[2:3,s] <- hdi(FinalInterceptDev)
     FinalSlopeDevs[1,s] <- mean(FinalSlopeDev)
     FinalSlopeDevs[2:3,s] <- hdi(FinalSlopeDev)
}

InterceptSeq <- 1:S - 0.15
SlopeSeq <- 1:S + 0.15
xRange <- c(0.5, S+0.5)
yRange <- range(c(range(PrelimSlopeDevs), range(PrelimInterceptDevs),
                  range(FinalSlopeDevs), range(FinalInterceptDevs)))
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
     points(x = InterceptSeq, y = PrelimInterceptDevs[1,], col = estCols[1])
     points(x = SlopeSeq, y = PrelimSlopeDevs[1,], col = estCols[2])
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
     points(x = InterceptSeq, y = FinalInterceptDevs[1,], col = estCols[1])
     points(x = SlopeSeq, y = FinalSlopeDevs[1,], col = estCols[2])
     segments(x0 = InterceptSeq, y0 = FinalInterceptDevs[2,], x1 = InterceptSeq,
              y1 = FinalInterceptDevs[3,], col = estCols[1])
     segments(x0 = SlopeSeq, y0 = FinalSlopeDevs[2,], x1 = SlopeSeq,
              y1 = FinalSlopeDevs[3,], col = estCols[2])
     abline(h = 0, lty = 2)
     legend(x = "topleft", legend = c("Intercept", "Slope"), pch = 1, col = estCols,
            bty = "n")
     mtext("Final model fit", side = 3, line = 1)

     mtext("Species", side = 1, line = -2, outer = TRUE)     
dev.off()






