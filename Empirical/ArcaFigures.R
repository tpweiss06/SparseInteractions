# This script will make the empirical figures for the York Gum Woodland data.
#    There will be two figures (1 per species), each with a plot of lambda and
#    alphas across both environmental gradients

setwd("~/Desktop/Wyoming/SparseInteractions/Empirical/")

# Set the figure name
FigName <- "Arca.pdf"

# Load in the empirical data and make objects for graphing
FullData <- read.csv("water_full_env.csv")
ObsPhos <- as.vector(scale(FullData$Colwell.P))
ObsShade <- as.vector(scale(FullData$Canopy))
EnvLength <- 1000
PlotPhos <- seq(min(ObsPhos, na.rm = TRUE), max(ObsPhos, na.rm = TRUE), length.out = EnvLength)
PlotShade <- seq(min(ObsShade, na.rm = TRUE), max(ObsShade, na.rm = TRUE), length.out = EnvLength)

ReserveNames <- c("Bendering", "Perenjori")

# Load in the model fits for the current species
# Phosphorous
load("StanFits/ARCA_Phos_FinalFit.rdata")
Post <- rstan::extract(FinalFit)
Phos <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
             reserve = reserve, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij,
             Inclusion_eij = Inclusion_eij, Intra = Intra)
rm(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
   Inclusion_eij, tau0, slab_scale, slab_df, Intra, Post)
# Canopy cover
load("StanFits/ARCA_Shade_FinalFit.rdata")
Post <- rstan::extract(FinalFit)
Shade <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
              reserve = reserve, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij,
              Inclusion_eij = Inclusion_eij, Intra = Intra)
rm(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
   Inclusion_eij, tau0, slab_scale, slab_df, Intra, Post)

# First, calculate the lambda values over both environmental variables
LambdaPlotVals <- array(NA, dim = c(2, 2, EnvLength, 3)) # Environmental co-variate, reserve, environmental sequence, metric (mean, lwr, upr)
for(i in 1:EnvLength){
     for(j in 1:2){
          PhosPost <- exp(Phos$Post$lambdas[,j,1] + Phos$Post$lambdas[,j,2] * PlotPhos[i])
          ShadePost <- exp(Shade$Post$lambdas[,j,1] + Shade$Post$lambdas[,j,2] * PlotShade[i])
          LambdaPlotVals[1,j,i,1] <- mean(PhosPost)
          LambdaPlotVals[1,j,i,2:3] <- HDInterval::hdi(PhosPost)
          LambdaPlotVals[2,j,i,1] <- mean(ShadePost)
          LambdaPlotVals[2,j,i,2:3] <- HDInterval::hdi(ShadePost)
     }
}
rm(PhosPost, ShadePost)

# Now, calculate the alpha values to plot for both environmental variables
GenericPlot <- array(NA, dim = c(2,3,EnvLength))
IntraPlot <- array(NA, dim = c(2,3,EnvLength))

PhosLegend <- c("Generic", "Intraspecific")
ShadeLegend <- c("Generic", "Intraspecific")

for(i in 1:EnvLength){
     # First phosphorous
     GenericPost <- exp(Phos$Post$alpha_generic[,1] + Phos$Post$alpha_generic[,2] * PlotPhos[i])
     IntraPost <- exp(Phos$Post$alpha_intra[,1] + Phos$Post$alpha_intra[,2] * PlotPhos[i])
     GenericPlot[1,1,i] <- mean(GenericPost)
     GenericPlot[1,2:3,i] <- HDInterval::hdi(GenericPost)
     IntraPlot[1,1,i] <- mean(IntraPost)
     IntraPlot[1,2:3,i] <- HDInterval::hdi(IntraPost)
     
     # Now shade
     GenericPost <- exp(Shade$Post$alpha_generic[,1] + Shade$Post$alpha_generic[,2] * PlotShade[i])
     IntraPost <- exp(Shade$Post$alpha_intra[,1] + Shade$Post$alpha_intra[,2] * PlotShade[i])
     GenericPlot[2,1,i] <- mean(GenericPost)
     GenericPlot[2,2:3,i] <- HDInterval::hdi(GenericPost)
     IntraPlot[2,1,i] <- mean(IntraPost)
     IntraPlot[2,2:3,i] <- HDInterval::hdi(IntraPost)
}

AlphaRange <- range(c(range(GenericPlot), range(IntraPlot)))
#AlphaRange <- c(0, 0.06)
LambdaRange <- range(LambdaPlotVals)
PhosRange <- range(ObsPhos, na.rm = TRUE)
ShadeRange <- range(ObsShade, na.rm = TRUE)
library(RColorBrewer)
DarkCols <- brewer.pal(n = 8, name = "Dark2")
LambdaCols <- DarkCols[3:4]
IntraCol <- DarkCols[5]
AlphaCols <- DarkCols[6:7]
LambdaLab <- expression(lambda["e,i"])
AlphaLab <- expression(alpha["e,i,j"])

pdf(file = FigName, width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), oma = c(4,4,2,2), mar = c(2,2,2,2))
     # First, plot the lambda results for phosphorous
     plot(NA, NA, xlim = PhosRange, ylim = LambdaRange, main = "", xlab = "", ylab = LambdaLab,
          las = 1, xpd = NA, cex.lab = 1.5)
     axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
     axis(side = 2, at = seq(0, 20, by = 1), labels = FALSE, tcl = -0.25)
     for(i in 1:2){
          lines(x = PlotPhos, y = LambdaPlotVals[1,i,,1], col = LambdaCols[i], lwd = 1.5)
          lines(x = PlotPhos, y = LambdaPlotVals[1,i,,2], col = LambdaCols[i], lty = 2)
          lines(x = PlotPhos, y = LambdaPlotVals[1,i,,3], col = LambdaCols[i], lty = 2)
     }
     legend("topright", legend = c("Bendering", "Perenjori"), lty = 1, col = LambdaCols,
            bty = "n")
     text(x = 0.95*PhosRange[1], y = 0.95*LambdaRange[2], 
          labels = expression(bold("a")), cex = 1.5)

     # Lambda results for canopy cover
     plot(NA, NA, xlim = ShadeRange, ylim = LambdaRange, main = "", xlab = "", ylab = "",
          las = 1, cex.lab = 1.5)
     axis(side = 1, at = seq(-3, 3, by = 0.25), labels = FALSE, tcl = -0.25)
     axis(side = 2, at = seq(0, 20, by = 1), labels = FALSE, tcl = -0.25)
     for(i in 1:2){
          lines(x = PlotShade, y = LambdaPlotVals[2,i,,1], col = LambdaCols[i], lwd = 1.5)
          lines(x = PlotShade, y = LambdaPlotVals[2,i,,2], col = LambdaCols[i], lty = 2)
          lines(x = PlotShade, y = LambdaPlotVals[2,i,,3], col = LambdaCols[i], lty = 2)
     }
     legend("topright", legend = c("Bendering", "Perenjori"), lty = 1, col = LambdaCols,
            bty = "n")
     text(x = 0.95*ShadeRange[1], y = 0.95*LambdaRange[2], 
          labels = expression(bold("b")), cex = 1.5)

     # Alpha results for phosphorous
     plot(NA, NA, xlim = PhosRange, ylim = AlphaRange, main = "", xlab = "Standardized phosphorous", 
          ylab = AlphaLab, las = 1, xpd = NA, cex.lab = 1.5)
     axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
     axis(side = 2, at = seq(0, 0.1, by = 0.005), labels = FALSE, tcl = -0.25)
     # Add the generic polygon and line
     xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
     yCoords <- c(GenericPlot[1,2,], GenericPlot[1,3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotPhos, y = GenericPlot[1,1,], lty = 1, lwd = 1.5)
     # Add the lines for the intraspecific term
     lines(x = PlotPhos, y = IntraPlot[1,1,], lwd = 1.5, col = IntraCol)
     lines(x = PlotPhos, y = IntraPlot[1,2,], lty = 2, col = IntraCol)
     lines(x = PlotPhos, y = IntraPlot[1,3,], lty = 2, col = IntraCol)
     legend("top", legend = PhosLegend, lty = 1, lwd = 1.5, col = c("black", IntraCol, AlphaCols),
            bty = "n")
     text(x = 0.95*PhosRange[1], y = 0.95*AlphaRange[2], 
          labels = expression(bold("c")), cex = 1.5)

     # Alpha results for canopy cover
     plot(NA, NA, xlim = ShadeRange, ylim = AlphaRange, main = "", xlab = "Standardized canopy cover", 
          ylab = "", las = 1, xpd = NA, cex.lab = 1.5)
     axis(side = 1, at = seq(-3, 3, by = 0.25), labels = FALSE, tcl = -0.25)
     axis(side = 2, at = seq(0, 0.1, by = 0.005), labels = FALSE, tcl = -0.25)
     # Add the generic polygon and line
     xCoords <- c(PlotShade, PlotShade[EnvLength:1])
     yCoords <- c(GenericPlot[2,2,], GenericPlot[2,3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotShade, y = GenericPlot[2,1,], lty = 1, lwd = 1.5)
     # Add the lines for the intraspecific term
     lines(x = PlotShade, y = IntraPlot[2,1,], lwd = 1.5, col = IntraCol)
     lines(x = PlotShade, y = IntraPlot[2,2,], lty = 2, col = IntraCol)
     lines(x = PlotShade, y = IntraPlot[2,3,], lty = 2, col = IntraCol)
     legend("topright", legend = ShadeLegend, lty = 1, lwd = 1.5, col = c("black", IntraCol, AlphaCols),
            bty = "n")
     text(x = 0.95*ShadeRange[1], y = 0.95*AlphaRange[2], 
          labels = expression(bold("d")), cex = 1.5)
dev.off()