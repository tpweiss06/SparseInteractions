# This script will make some summary, result figures for the SparseInteraction
#    model runs.

setwd("~/Desktop/Wyoming/SparseInteractions/")

# First make a figure with 4 panels (Arca and Waitzia rows, phosphorous and shade columns) 
#     with lambda_i on each panel (two lines; 1 for each reserve)
# Then make four more figures. 1 for Arca for each environment and 1 for Waitzia
#    under each environment. Each will have the number of panels corresponding to 
#    the number of species/reserve deviations from the generic alpha.

# Make vectors for the environmental data
FullData <- read.csv("water_full_env.csv")
ObsPhos <- as.vector(scale(FullData$Colwell.P))
ObsShade <- as.vector(scale(FullData$Canopy))
EnvLength <- 1000
PlotPhos <- seq(min(ObsPhos, na.rm = TRUE), max(ObsPhos, na.rm = TRUE), length.out = EnvLength)
PlotShade <- seq(min(ObsShade, na.rm = TRUE), max(ObsShade, na.rm = TRUE), length.out = EnvLength)

ReserveNames <- c("Bendering", "Perenjori")

############################ Lambda figure
LambdaPlotVals <- array(NA, dim = c(4, 2, EnvLength, 3)) 
# load in the final model fit for each species, calculate the posterior for lambda_i,
#    and save the means and 95% credible interval values
# Waitzia Phos
load("Waitzia/Phosphorous/Model Fits/Waitzia_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotPhos[i]
     LambdaPost[2,,i] <- Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotPhos[i]
     LambdaPlotVals[1,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[1,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[1,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[1,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# Waitzia Shade
load("Waitzia/Shade/Model Fits/Waitzia_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotShade[i]
     LambdaPost[2,,i] <- Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotShade[i]
     LambdaPlotVals[2,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[2,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[2,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[2,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# ARCA Phos
load("ARCA/Phosphorous/Model Fits/ARCA_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotPhos[i]
     LambdaPost[2,,i] <- Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotPhos[i]
     LambdaPlotVals[3,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[3,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[3,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[3,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# ARCA Shade
load("ARCA/Shade/Model Fits/ARCA_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotShade[i]
     LambdaPost[2,,i] <- Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotShade[i]
     LambdaPlotVals[4,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[4,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[4,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[4,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# Now remove the superfluous objects from the environment
rm(env, Fecundity, FinalFit, Inclusion_eij, Inclusion_ij, N, Posteriors, reserve,
   S, slab_df, slab_scale, SpMatrix, SpNames, tau0)

library(RColorBrewer)
PairedCols <- brewer.pal(n = 12, name = "Paired")
ReserveCols <- PairedCols[c(6,10)]
xRange <- matrix(NA, nrow = 4, ncol = 2)
xRange[1,] <- range(PlotPhos)
xRange[2,] <- range(PlotShade)
xRange[3,] <- range(PlotPhos)
xRange[4,] <- range(PlotShade)
yRange <- range(LambdaPlotVals)
pdf(file = "Result Figures/Lambdas.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), oma = c(4,4,2,2), mar = c(1,1,1,1))
     for(i in 1:4){
          plot(NA, NA, xlim = xRange[i,], ylim = yRange, main = "", yaxt = "n",
               xaxt = "n", xlab = "", ylab = "")
          if(i == 1 | i == 3){
               lines(x = PlotPhos, y = LambdaPlotVals[i,1,,1], col = ReserveCols[1])
               lines(x = PlotPhos, y = LambdaPlotVals[i,1,,2], col = ReserveCols[1], lty = 2)
               lines(x = PlotPhos, y = LambdaPlotVals[i,1,,3], col = ReserveCols[1], lty = 2)
               
               lines(x = PlotPhos, y = LambdaPlotVals[i,2,,1], col = ReserveCols[2])
               lines(x = PlotPhos, y = LambdaPlotVals[i,2,,2], col = ReserveCols[2], lty = 2)
               lines(x = PlotPhos, y = LambdaPlotVals[i,2,,3], col = ReserveCols[2], lty = 2)
          }else{
               lines(x = PlotShade, y = LambdaPlotVals[i,1,,1], col = ReserveCols[1])
               lines(x = PlotShade, y = LambdaPlotVals[i,1,,2], col = ReserveCols[1], lty = 2)
               lines(x = PlotShade, y = LambdaPlotVals[i,1,,3], col = ReserveCols[1], lty = 2)
               
               lines(x = PlotShade, y = LambdaPlotVals[i,2,,1], col = ReserveCols[2])
               lines(x = PlotShade, y = LambdaPlotVals[i,2,,2], col = ReserveCols[2], lty = 2)
               lines(x = PlotShade, y = LambdaPlotVals[i,2,,3], col = ReserveCols[2], lty = 2)
          }
          if(i == 1){
               mtext("Waitzia", side = 2, line = 3.5)
               axis(1, labels = FALSE)
               axis(2)
               #legend("bottomleft", legend = c("Bendering", "Perenjori"), lty = c(1,1), 
               #       col = ReserveCols, bty = "n")
          } else if(i == 2){
               axis(1, labels = FALSE)
               axis(2, labels = FALSE)
          } else if(i == 3){
               mtext("ARCA", side = 2, line = 3.5)
               mtext("Standardized Phosphorous", side = 1, line = 3)
               axis(1)
               axis(2)
          }else{
               mtext("Standardized canopy cover", side = 1, line = 3)
               axis(1)
               axis(2, labels = FALSE)
          }
          legend("bottomleft", legend = c("Bendering", "Perenjori"), lty = c(1,1), 
                 col = ReserveCols, bty = "n")
     }
     mtext(expression(paste("ln(", lambda[i], ")", sep = "")), side = 2, outer = TRUE, line = 1.5)
dev.off()

############################ Alpha_eij figures
# First Waitzia and shade
load("Waitzia/Shade/Model Fits/Waitzia_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
# There are 2 instances of species/reserve specific deviations here
SpecificAlpha <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
GenericAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
ij_vals <- which(Inclusion_ij == 1, arr.ind = TRUE)

SpecificPlot <- array(NA, dim = c(2, 3, EnvLength))
GenericPlot <- matrix(NA, nrow = 3, ncol = EnvLength)

for(i in 1:EnvLength){
     GenericAlpha[,i] <- Posteriors$alphas[,1] + Posteriors$alphas[,2] * PlotShade[i]
     GenericPlot[1,i] <- mean(GenericAlpha[,i])
     GenericPlot[2:3,i] <- HDInterval::hdi(GenericAlpha[,i])
     for(s in 1:2){
          SpecificAlpha[s,,i] <- GenericAlpha[,i] + Posteriors$alpha_hat_ij[,ij_vals[s,1],ij_vals[s,2]]
          SpecificPlot[s,1,i] <- mean(SpecificAlpha[s,,i])
          SpecificPlot[s,2:3,i] <- HDInterval::hdi(SpecificAlpha[s,,i])
     }
}

PairedCols <- brewer.pal(n = 12, name = "Paired")
AlphaCols <- PairedCols[c(2, 12)]
xRange <- range(PlotShade)
yRange <- range(c(range(GenericPlot), range(SpecificPlot)))

pdf(file = "Result Figures/Alpha_Waitzia_Shade.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(1,2), oma = c(4,4,2,2), mar = c(1,1,1,1))
     for(i in 1:2){
          plot(NA, NA, xlim = xRange, ylim = yRange, main = "", yaxt = "n",
               xaxt = "n", xlab = "", ylab = "")
          xCoords <- c(PlotShade, PlotShade[EnvLength:1])
          yCoords <- c(GenericPlot[2,], GenericPlot[3,EnvLength:1])
          polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
          lines(x = PlotShade, y = GenericPlot[1,], lty = 1, lwd = 2)
          
          lines(x = PlotShade, y = SpecificPlot[i,1,], lwd = 1.5, col = AlphaCols[i])
          lines(x = PlotShade, y = SpecificPlot[i,2,], lwd = 1, col = AlphaCols[i], lty = 2)
          lines(x = PlotShade, y = SpecificPlot[i,3,], lwd = 1, col = AlphaCols[i], lty = 2)
          
          if(i == 1){
               axis(1)
               axis(2)
               mtext(expression(paste("ln(", alpha, ")", sep = "")), side = 2, line = 2.5)
          } else if(i == 2){
               axis(1)
               axis(2, labels = FALSE)
          } 
          mtext(SpNames[ij_vals[i,2]], side = 1, line = -3.5, col = AlphaCols[i])
          mtext(ReserveNames[ij_vals[i,1]], side = 1, line = -2.25, col = AlphaCols[i])
     }
     mtext("Standardized canopy cover", side = 1, outer = TRUE, line = 2)
dev.off()

# Now Waitzia and phosphorous
load("Waitzia/Phosphorous/Model Fits/Waitzia_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
# There are 4 instances of species/reserve specific deviations here
SpecificAlpha <- array(NA, dim = c(4, length(Posteriors$lambdas[,1,1]), EnvLength))
GenericAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
ij_vals <- which(Inclusion_ij == 1, arr.ind = TRUE)

SpecificPlot <- array(NA, dim = c(4, 3, EnvLength))
GenericPlot <- matrix(NA, nrow = 3, ncol = EnvLength)

for(i in 1:EnvLength){
     GenericAlpha[,i] <- Posteriors$alphas[,1] + Posteriors$alphas[,2] * PlotPhos[i]
     GenericPlot[1,i] <- mean(GenericAlpha[,i])
     GenericPlot[2:3,i] <- HDInterval::hdi(GenericAlpha[,i])
     for(s in 1:4){
          SpecificAlpha[s,,i] <- GenericAlpha[,i] + Posteriors$alpha_hat_ij[,ij_vals[s,1],ij_vals[s,2]]
          SpecificPlot[s,1,i] <- mean(SpecificAlpha[s,,i])
          SpecificPlot[s,2:3,i] <- HDInterval::hdi(SpecificAlpha[s,,i])
     }
}

PairedCols <- brewer.pal(n = 12, name = "Paired")
AlphaCols <- PairedCols[c(2, 6, 10, 12)]
xRange <- range(PlotPhos)
yRange <- range(c(range(GenericPlot), range(SpecificPlot)))

pdf(file = "Result Figures/Alpha_Waitzia_Phos.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), oma = c(4,4,2,2), mar = c(1,1,1,1))
     for(i in 1:4){
          plot(NA, NA, xlim = xRange, ylim = yRange, main = "", yaxt = "n",
               xaxt = "n", xlab = "", ylab = "")
          xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
          yCoords <- c(GenericPlot[2,], GenericPlot[3,EnvLength:1])
          polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
          lines(x = PlotPhos, y = GenericAlpha[1,], lty = 1, lwd = 2)
          
          lines(x = PlotPhos, y = SpecificPlot[i,1,], lwd = 1.5, col = AlphaCols[i])
          lines(x = PlotPhos, y = SpecificPlot[i,2,], lwd = 1, col = AlphaCols[i], lty = 2)
          lines(x = PlotPhos, y = SpecificPlot[i,3,], lwd = 1, col = AlphaCols[i], lty = 2)
     
          if(i == 1){
               axis(1, labels = FALSE)
               axis(2)
          } else if(i == 2){
               axis(1, labels = FALSE)
               axis(2, labels = FALSE)
          } else if(i == 3){
               axis(1)
               axis(2)
          }else{
               axis(1)
               axis(2, labels = FALSE)
          }
          mtext(SpNames[ij_vals[i,2]], side = 1, line = -3.5, col = AlphaCols[i])
          mtext(ReserveNames[ij_vals[i,1]], side = 1, line = -2.25, col = AlphaCols[i])
     }
     mtext(expression(paste("ln(", alpha, ")", sep = "")), side = 2, outer = TRUE, line = 1.5)
     mtext("Standardized phosphorous", side = 1, outer = TRUE, line = 2)
dev.off()

# Now ARCA and shade
load("ARCA/Shade/Model Fits/ARCA_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
# There are 0 instances of species/reserve specific deviations here
GenericAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
GenericPlot <- matrix(NA, nrow = 3, ncol = EnvLength)

for(i in 1:EnvLength){
     GenericAlpha[,i] <- Posteriors$alphas[,1] + Posteriors$alphas[,2] * PlotShade[i]
     GenericPlot[1,i] <- mean(GenericAlpha[,i])
     GenericPlot[2:3,i] <- HDInterval::hdi(GenericAlpha[,i])
}

xRange <- range(PlotShade)
yRange <- range(GenericPlot)

pdf(file = "Result Figures/Alpha_ARCA_Shade.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "Standardized canopy cover",
          ylab = expression(paste("ln(", alpha, ")", sep = "")))
     xCoords <- c(PlotShade, PlotShade[EnvLength:1])
     yCoords <- c(GenericPlot[2,], GenericPlot[3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotShade, y = GenericPlot[1,], lty = 1, lwd = 2)
dev.off()

# Now ARCA and phosphorous
load("ARCA/Phosphorous/Model Fits/ARCA_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
# There are 0 instances of species/reserve specific deviations here
GenericAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
GenericPlot <- matrix(NA, nrow = 3, ncol = EnvLength)
for(i in 1:EnvLength){
     GenericAlpha[,i] <- Posteriors$alphas[,1] + Posteriors$alphas[,2] * PlotPhos[i]
     GenericPlot[1,i] <- mean(GenericAlpha[,i])
     GenericPlot[2:3,i] <- HDInterval::hdi(GenericAlpha[,i])
}

xRange <- range(PlotPhos)
yRange <- range(GenericPlot)
pdf(file = "Result Figures/Alpha_ARCA_Phos.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "Standardized phosphorous",
          ylab = expression(paste("ln(", alpha, ")", sep = "")))
     xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
     yCoords <- c(GenericPlot[2,], GenericPlot[3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotPhos, y = GenericPlot[1,], lty = 1, lwd = 2)
dev.off()


################################################################################
####### Now make plots on the natural scale rather than the log scale ##########
################################################################################
library(RColorBrewer)
############################ Lambda figure
LambdaPlotVals <- array(NA, dim = c(4, 2, EnvLength, 3)) 
# load in the final model fit for each species, calculate the posterior for lambda_i,
#    and save the means and 95% credible interval values
# Waitzia Phos
load("Waitzia/Phosphorous/Model Fits/Waitzia_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- exp(Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotPhos[i])
     LambdaPost[2,,i] <- exp(Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotPhos[i])
     LambdaPlotVals[1,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[1,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[1,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[1,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# Waitzia Shade
load("Waitzia/Shade/Model Fits/Waitzia_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- exp(Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotShade[i])
     LambdaPost[2,,i] <- exp(Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotShade[i])
     LambdaPlotVals[2,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[2,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[2,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[2,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# ARCA Phos
load("ARCA/Phosphorous/Model Fits/ARCA_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- exp(Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotPhos[i])
     LambdaPost[2,,i] <- exp(Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotPhos[i])
     LambdaPlotVals[3,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[3,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[3,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[3,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# ARCA Shade
load("ARCA/Shade/Model Fits/ARCA_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
LambdaPost <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
for(i in 1:EnvLength){
     LambdaPost[1,,i] <- exp(Posteriors$lambdas[,1,1] + Posteriors$lambdas[,1,2] * PlotShade[i])
     LambdaPost[2,,i] <- exp(Posteriors$lambdas[,2,1] + Posteriors$lambdas[,2,2] * PlotShade[i])
     LambdaPlotVals[4,1,i,1] <- mean(LambdaPost[1,,i])
     LambdaPlotVals[4,1,i,2:3] <- HDInterval::hdi(LambdaPost[1,,i])
     LambdaPlotVals[4,2,i,1] <- mean(LambdaPost[2,,i])
     LambdaPlotVals[4,2,i,2:3] <- HDInterval::hdi(LambdaPost[2,,i])
}
# Now remove the superfluous objects from the environment
rm(env, Fecundity, FinalFit, Inclusion_eij, Inclusion_ij, N, Posteriors, reserve,
   S, slab_df, slab_scale, SpMatrix, SpNames, tau0)

library(RColorBrewer)
PairedCols <- brewer.pal(n = 12, name = "Paired")
ReserveCols <- PairedCols[c(6,10)]
xRange <- matrix(NA, nrow = 4, ncol = 2)
xRange[1,] <- range(PlotPhos)
xRange[2,] <- range(PlotShade)
xRange[3,] <- range(PlotPhos)
xRange[4,] <- range(PlotShade)
yRange <- matrix(NA, nrow = 4, ncol = 2)
yRange[1,] <- c(0, 25)
yRange[2,] <- c(0, 25)
yRange[3,] <- c(0, 20)
yRange[4,] <- c(0, 20)
LegendLoc <- c("topleft", "topleft", "topleft", "topright")

pdf(file = "Result Figures/Lambdas_NaturalScale.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), oma = c(4,4,2,2), mar = c(2,2,2,2))
     for(i in 1:4){
          plot(NA, NA, xlim = xRange[i,], ylim = yRange[i,], main = "", xlab = "", ylab = "")
          axis(2, at = seq(0, 25, by = 1), labels = FALSE, tcl = -0.25)
          if(i == 1 | i == 3){
               lines(x = PlotPhos, y = LambdaPlotVals[i,1,,1], col = ReserveCols[1])
               lines(x = PlotPhos, y = LambdaPlotVals[i,1,,2], col = ReserveCols[1], lty = 2)
               lines(x = PlotPhos, y = LambdaPlotVals[i,1,,3], col = ReserveCols[1], lty = 2)
          
               lines(x = PlotPhos, y = LambdaPlotVals[i,2,,1], col = ReserveCols[2])
               lines(x = PlotPhos, y = LambdaPlotVals[i,2,,2], col = ReserveCols[2], lty = 2)
               lines(x = PlotPhos, y = LambdaPlotVals[i,2,,3], col = ReserveCols[2], lty = 2)
               
               axis(1, at = seq(-2, 3, by = 0.2), labels = FALSE, tcl = -0.25)
               mtext(expression(lambda["e,i"]), side = 2, line = 2.5)
          }else{
               lines(x = PlotShade, y = LambdaPlotVals[i,1,,1], col = ReserveCols[1])
               lines(x = PlotShade, y = LambdaPlotVals[i,1,,2], col = ReserveCols[1], lty = 2)
               lines(x = PlotShade, y = LambdaPlotVals[i,1,,3], col = ReserveCols[1], lty = 2)
          
               lines(x = PlotShade, y = LambdaPlotVals[i,2,,1], col = ReserveCols[2])
               lines(x = PlotShade, y = LambdaPlotVals[i,2,,2], col = ReserveCols[2], lty = 2)
               lines(x = PlotShade, y = LambdaPlotVals[i,2,,3], col = ReserveCols[2], lty = 2)
               axis(1, at = seq(-3, 3, by = 0.2), labels = FALSE, tcl = -0.25)
          }
          if(i == 1){
               mtext(expression(italic("W. acuminata")), side = 2, line = 4)
          } else if(i == 3){
               mtext(expression(italic("A. calendula")), side = 2, line = 4)
               mtext("Standardized Phosphorous", side = 1, line = 3)
          }else if (i == 4){
               mtext("Standardized canopy cover", side = 1, line = 3)
          }
          legend(LegendLoc[i], legend = c("Bendering", "Perenjori"), lty = c(1,1), 
                 col = ReserveCols, bty = "n")
     }
dev.off()

# Now Waitzia and phosphorous on the natural scale
load("Waitzia/Phosphorous/Model Fits/Waitzia_Phos_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
# There are 4 instances of species/reserve specific deviations here
SpecificAlpha <- array(NA, dim = c(4, length(Posteriors$lambdas[,1,1]), EnvLength))
GenericAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
ij_vals <- which(Inclusion_ij == 1, arr.ind = TRUE)
SpecificPlot <- array(NA, dim = c(4, 3, EnvLength))
GenericPlot <- matrix(NA, nrow = 3, ncol = EnvLength)
for(i in 1:EnvLength){
     GenericAlpha[,i] <- Posteriors$alphas[,1] + Posteriors$alphas[,2] * PlotPhos[i]
     for(s in 1:4){
          SpecificAlpha[s,,i] <- exp(GenericAlpha[,i] + Posteriors$alpha_hat_ij[,ij_vals[s,1],ij_vals[s,2]])
          SpecificPlot[s,1,i] <- mean(SpecificAlpha[s,,i])
          SpecificPlot[s,2:3,i] <- HDInterval::hdi(SpecificAlpha[s,,i])
     }
     GenericPlot[1,i] <- mean(exp(GenericAlpha[,i]))
     GenericPlot[2:3,i] <- HDInterval::hdi(exp(GenericAlpha[,i]))
}

PairedCols <- brewer.pal(n = 12, name = "Paired")
AlphaCols <- PairedCols[c(2, 6, 10, 12)]
xRange <- range(PlotPhos)
yAxisSeqs <- vector(mode = "list", length = 4)
yAxisSeqs[[1]] <- seq(0, 0.4, by = 0.02)
yAxisSeqs[[2]] <- seq(0, 0.04, by = 0.002)
yAxisSeqs[[3]] <- seq(0, 0.04, by = 0.002)
yAxisSeqs[[4]] <- seq(0, 0.2, by = 0.01)

PlotSpNames <- c(expression(italic("Hypochaeris glabra")), expression(italic("Pentaschistis airoides")),
                 expression(italic("Schoenus nanus")), expression(italic("Waitzia acuminata")))

pdf(file = "Result Figures/Alpha_Waitzia_Phos_NaturalScale.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), oma = c(4,4,2,2), mar = c(2,2,2,2))
     for(i in 1:4){
          yRange <- range(c(range(GenericPlot), range(SpecificPlot[i,,])))
          plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "", ylab = "")
          xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
          yCoords <- c(GenericPlot[2,], GenericPlot[3,EnvLength:1])
          polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
          lines(x = PlotPhos, y = GenericPlot[1,], lty = 1, lwd = 2)
          lines(x = PlotPhos, y = SpecificPlot[i,1,], lwd = 1.5, col = AlphaCols[i])
          lines(x = PlotPhos, y = SpecificPlot[i,2,], lwd = 1, col = AlphaCols[i], lty = 2)
          lines(x = PlotPhos, y = SpecificPlot[i,3,], lwd = 1, col = AlphaCols[i], lty = 2)
          axis(1, at = seq(-2, 3, by = 0.2), labels = FALSE, tcl = -0.25)
          axis(2, at = yAxisSeqs[[i]], labels = FALSE, tcl = -0.25)
          mtext(PlotSpNames[i], side = 3, line = -1.25, col = AlphaCols[i])
          mtext(ReserveNames[ij_vals[i,1]], side = 3, line = -2.5, col = AlphaCols[i])
     }
     mtext(expression(alpha["e,i,j"]), side = 2, outer = TRUE, line = 1.5, cex = 1.5)
     mtext("Standardized phosphorous", side = 1, outer = TRUE, line = 2)
dev.off()

# Now Waitzia and shade
load("Waitzia/Shade/Model Fits/Waitzia_Shade_FinalFit.rdata")
Posteriors <- rstan::extract(FinalFit)
# There are 2 instances of species/reserve specific deviations here
SpecificAlpha <- array(NA, dim = c(2, length(Posteriors$lambdas[,1,1]), EnvLength))
GenericAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
ij_vals <- which(Inclusion_ij == 1, arr.ind = TRUE)
SpecificPlot <- array(NA, dim = c(2, 3, EnvLength))
GenericPlot <- matrix(NA, nrow = 3, ncol = EnvLength)
for(i in 1:EnvLength){
     GenericAlpha[,i] <- Posteriors$alphas[,1] + Posteriors$alphas[,2] * PlotShade[i]
     GenericPlot[1,i] <- mean(exp(GenericAlpha[,i]))
     GenericPlot[2:3,i] <- HDInterval::hdi(exp(GenericAlpha[,i]))
     for(s in 1:2){
          SpecificAlpha[s,,i] <- exp(GenericAlpha[,i] + Posteriors$alpha_hat_ij[,ij_vals[s,1],ij_vals[s,2]])
          SpecificPlot[s,1,i] <- mean(SpecificAlpha[s,,i])
          SpecificPlot[s,2:3,i] <- HDInterval::hdi(SpecificAlpha[s,,i])
     }
}

PairedCols <- brewer.pal(n = 12, name = "Paired")
AlphaCols <- PairedCols[c(2, 12)]
xRange <- range(PlotShade)
yAxisSeqs <- vector(mode = "list", length = 2)
yAxisSeqs[[1]] <- seq(0, 1, by = 0.025)
yAxisSeqs[[2]] <- seq(0, 0.5, by = 0.02)

PlotSpNames <- c(expression(italic("Hypochaeris glabra")), expression(italic("Waitzia acuminata")))

pdf(file = "Result Figures/Alpha_Waitzia_Shade_NaturalScale.pdf", width = 8, height = 4, onefile = FALSE, paper = "special")
     par(mfrow = c(1,2), oma = c(4,4,2,2), mar = c(2,2,2,2))
     for(i in 1:2){
          yRange <- range(c(range(GenericPlot), range(SpecificPlot[i,,])))
          plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "", ylab = "")
          xCoords <- c(PlotShade, PlotShade[EnvLength:1])
          yCoords <- c(GenericPlot[2,], GenericPlot[3,EnvLength:1])
          polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
          lines(x = PlotShade, y = GenericPlot[1,], lty = 1, lwd = 2)
          lines(x = PlotShade, y = SpecificPlot[i,1,], lwd = 1.5, col = AlphaCols[i])
          lines(x = PlotShade, y = SpecificPlot[i,2,], lwd = 1, col = AlphaCols[i], lty = 2)
          lines(x = PlotShade, y = SpecificPlot[i,3,], lwd = 1, col = AlphaCols[i], lty = 2)
          axis(1, at = seq(-3, 3, by = 0.2), labels = FALSE, tcl = -0.25)
          axis(2, at = yAxisSeqs[[i]], labels = FALSE, tcl = -0.25)
          
          mtext(PlotSpNames[i], side = 3, line = -1.25, col = AlphaCols[i])
          mtext(ReserveNames[ij_vals[i,1]], side = 3, line = -2.5, col = AlphaCols[i])
     }
     mtext("Standardized canopy cover", side = 1, outer = TRUE, line = 2)
     mtext(expression(alpha["e,i,j"]), side = 2, outer = TRUE, line = 1.5, cex = 1.5)
dev.off()
     
# Now ARCA for both phosphorous and shade
load("ARCA/Phosphorous/Model Fits/ARCA_Phos_FinalFit.rdata")
PhosPosteriors <- rstan::extract(FinalFit)
load("ARCA/Shade/Model Fits/ARCA_Shade_FinalFit.rdata")
ShadePosteriors <- rstan::extract(FinalFit)
PhosAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
PhosPlot <- matrix(NA, nrow = 4, ncol = EnvLength)
ShadeAlpha <- matrix(NA, nrow = length(Posteriors$lambdas[,1,1]), ncol = EnvLength)
ShadePlot <- matrix(NA, nrow = 4, ncol = EnvLength)

for(i in 1:EnvLength){
     PhosAlpha[,i] <- exp(PhosPosteriors$alphas[,1] + PhosPosteriors$alphas[,2] * PlotPhos[i])
     PhosPlot[1,i] <- mean(PhosAlpha[,i])
     PhosPlot[2:3,i] <- HDInterval::hdi(PhosAlpha[,i])
     PhosPlot[4,i] <- median(PhosAlpha[,i])
     
     ShadeAlpha[,i] <- exp(ShadePosteriors$alphas[,1] + ShadePosteriors$alphas[,2] * PlotShade[i])
     ShadePlot[1,i] <- mean(ShadeAlpha[,i])
     ShadePlot[2:3,i] <- HDInterval::hdi(ShadeAlpha[,i])
     ShadePlot[4,i] <- median(ShadeAlpha[,i])
}

pdf(file = "Result Figures/Alpha_ARCA_NaturalScale.pdf", width = 8, height = 4, onefile = FALSE, paper = "special")
     par(mfrow = c(1,2))
     # First, phosphorous
     plot(NA, NA, xlim = range(PlotPhos), ylim = c(0,30), xlab = "Standardized phosphorous",
          ylab = expression(alpha["e,i,j"]))
     xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
     yCoords <- c(PhosPlot[2,], PhosPlot[3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotPhos, y = PhosPlot[1,], lty = 1, lwd = 2)
     lines(x = PlotPhos, y = PhosPlot[4,], lty = 2, lwd = 2)
     # Now make objects to allow plotting an inset later
     p <- c(.10, .75, .20, .7)
     fig.new.1 <- c(grconvertX(p[1:2], from="npc", to="ndc"),
                    grconvertY(p[3:4], from="npc", to="ndc"))
     # Now make the shade figure
     plot(NA, NA, xlim = range(PlotShade), ylim = c(0, 0.06), xlab = "Standardized canopy cover",
          ylab = expression(alpha["e,i,j"]))
     xCoords <- c(PlotShade, PlotShade[EnvLength:1])
     yCoords <- c(ShadePlot[2,], ShadePlot[3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotShade, y = ShadePlot[1,], lty = 1, lwd = 2)
     lines(x = PlotShade, y = ShadePlot[4,], lty = 2, lwd = 2)
     
     # Now plot the inset on the phosphorous figure
     op <- par(fig=fig.new.1, cex=.5, new=TRUE, mar=rep(0, 4))
     plot(NA, NA, xlim = range(PlotPhos), ylim = c(0, 0.06), xlab = "", ylab = "")
     xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
     yCoords <- c(PhosPlot[2,], PhosPlot[3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotPhos, y = PhosPlot[1,], lty = 1, lwd = 2)
     lines(x = PlotPhos, y = PhosPlot[4,], lty = 2, lwd = 2)
dev.off()


# Native = T
# Exotic = H