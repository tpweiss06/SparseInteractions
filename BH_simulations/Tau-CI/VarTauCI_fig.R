# This script will make a figure to demonstrate the effect of using different values
#    of fixed tau and different CI thresholds on parameter estimation. It will be
#    the same layout as the left panels from Figure 1, but with one column corresponding
#    to different values of tau and the other different values of CI thresholds. There
#    will be three rows.

rm(list = ls())
library(RColorBrewer)
library(here)

# Need to decide how I want to structure this. One list? Two lists of 3?
#    One list of 6 and use mfcol = c(3,2) to fill by columns?

# Create lists for the results from different sizes of datasets
Alphas <- vector(mode = "list", length = 6)
Lambdas <- vector(mode = "list", length = 6)
InterceptInclusion <- vector(mode = "list", length = 6)
SlopeInclusion <- vector(mode = "list", length = 6)
GenericDeviations <- vector(mode = "list", length = 6)

# First populate the lists with values from the varying tau fits
load(here("BH_simulations/Tau-CI/TauFits_GraphStuff.rdata"))
load(here("BH_simulations/Tau-CI/TauFits.rdata"))
for(i in 1:3){
        Alphas[[i]] <- AlphaEsts[i,,,]
        Lambdas[[i]] <- LambdaEsts[i,,]
        InterceptInclusion[[i]] <- AllInclusion_ij[,i]
        SlopeInclusion[[i]] <- AllInclusion_eij[,i]
        
        # Calculate the mean and hdi for the posterior of the deviations in the
        #      generic alpha terms
        InterceptDev <- FinalTauFits[[i]]$alpha_generic[,1] - TrueGenericIntercept[i]
        SlopeDev <- FinalTauFits[[i]]$alpha_generic[,2] - TrueGenericSlope[i]
        GenericDeviations[[i]] <- matrix(data = NA, nrow = 2, ncol = 3)
        GenericDeviations[[i]][1,1] <- mean(InterceptDev)
        GenericDeviations[[i]][2,1] <- mean(SlopeDev)
        GenericDeviations[[i]][1,2:3] <- HDInterval::hdi(InterceptDev)
        GenericDeviations[[i]][2,2:3] <- HDInterval::hdi(SlopeDev)
}
# Now populate the lists with values from the varying CI thresholds
load(here("BH_simulations/Tau-CI/CIFits_GraphStuff.rdata"))
load(here("BH_simulations/Tau-CI/CIFits.rdata"))
for(i in 1:3){
        Alphas[[i+3]] <- AlphaEsts[i,,,]
        Lambdas[[i+3]] <- LambdaEsts[i,,]
        InterceptInclusion[[i+3]] <- AllInclusion_ij[,i]
        SlopeInclusion[[i+3]] <- AllInclusion_eij[,i]
        
        # Calculate the mean and hdi for the posterior of the deviations in the
        #      generic alpha terms
        InterceptDev <- FinalCIFits[[i]]$alpha_generic[,1] - TrueGenericIntercept[i]
        SlopeDev <- FinalCIFits[[i]]$alpha_generic[,2] - TrueGenericSlope[i]
        GenericDeviations[[i+3]] <- matrix(data = NA, nrow = 2, ncol = 3)
        GenericDeviations[[i+3]][1,1] <- mean(InterceptDev)
        GenericDeviations[[i+3]][2,1] <- mean(SlopeDev)
        GenericDeviations[[i+3]][1,2:3] <- HDInterval::hdi(InterceptDev)
        GenericDeviations[[i+3]][2,2:3] <- HDInterval::hdi(SlopeDev)
}

# Now make the figure: a 6 panel figure with each panel showing alternating horizontal bands to delineate
#       (1) lambda parameters, (2) intraspecific alphas, (3) model estimated species-specific
#       alpha terms, and (4) generic alpha terms vs. the weighted averages they are
#       estimating. Points will be colored as intercept or slope. Lambda points will be
#       filled squares. Intraspecific alpha points will be filled triangles. Species specific
#       interspecific alpha terms estimated by the model will be filled circles. The generic
#       interspecific alpha terms estimated by the model will be open circles.
#       The first column will show the results for varying tau values (bottom row = low)
#       The second column will show the results for varying CI thresholds (bottom = low)

AllNonGeneric <- union(which(InterceptInclusion[[3]] == 1), which(SlopeInclusion[[3]] == 1))
NumNonGeneric <- length(AllNonGeneric)

Focal <- 8
S <- 15
param_xRange <- c(-7, 7)
param_yRange <- c(0, 6)

# Set some pch and color values for the graph
InterceptPoints <- 15
SlopePoints <- 17
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
InterceptCol <- Dark2Cols[1]
SlopeCol <- Dark2Cols[2]

# Divide up the plot space into appropriate vertical subdivisions
InterceptOffset <- 0.55
SlopeOffset <- 0.45
Generic_y <- c(InterceptOffset, SlopeOffset)
NonGeneric_y <- matrix(data = NA, nrow = NumNonGeneric, ncol = 2)
NonGeneric_y[,1] <- c(1.35, 1.95, 2.55, 3.15, 3.75)
NonGeneric_y[,2] <- c(1.25, 1.85, 2.45, 3.05, 3.65)
Intra_y <- c(4 + InterceptOffset, 4 + SlopeOffset)
Lambda_y <- c(5 + InterceptOffset, 5 + SlopeOffset)

# Define the points for y axis labels
GenericLab_y <- 0.5
GenericLab <- expression(paste("Generic ", alpha["e,i,j"], sep = ""))
IntraLab_y <- 4.5
IntraLab <- expression(paste(alpha["e,i,i"], sep = ""))
LambdaLab_y <- 5.5
LambdaLab <- expression(paste(lambda["e,i"], sep = ""))
NonGenericLab_y1 <- 2.95
NonGenericLab_y2 <- 2.45
NonGenericLab_y3 <- 1.95
NonGenericLab1 <- "Species"
NonGenericLab2 <- "specific"
NonGenericLab3 <- expression(paste(alpha["e,i,j"], sep = ""))

# Make the letters for the subpanels
Letters <- c(expression(bold("a")), expression(bold("b")), expression(bold("c")),
                    expression(bold("d")), expression(bold("e")), expression(bold("f")))

FigName <- here("BH_simulations/Tau-CI/Tau-CI_SuppResults.pdf")
pdf(file = FigName, width = 7, height = 8, onefile = FALSE, paper = "special")
     par(mfcol = c(3,2), mar = c(2, 2.5, 2, 2.5), oma = c(2, 0, 2, 0))
     for(i in 1:6){
          # Plot the parameters
          plot(NA, NA, xlim = param_xRange, ylim = param_yRange, xlab = "",
               ylab = "", axes = FALSE)
          axis(side = 1, at = seq(-6, 6, by = 2), tcl = -0.5)
          axis(side = 1, at = -7:7, tcl = -0.5, labels = FALSE)
          axis(side = 1, at = seq(-8, 8, by = 0.5), tcl = -0.25, labels = FALSE)
          # Add boxes for intraspecifc and generic alpha components
          rect(xleft = -7.5, ybottom = 5, xright = 7.5, ytop = 6, col = "grey", border = NA, xpd = NA)
          
          rect(xleft = -7.5, ybottom = 1, xright = 7.5, ytop = 1.6, col = "grey", border = NA, xpd = NA)
          rect(xleft = -7.5, ybottom = 1.6, xright = 7.5, ytop = 2.2, col = "lightgrey", border = NA, xpd = NA)
          rect(xleft = -7.5, ybottom = 2.2, xright = 7.5, ytop = 2.8, col = "grey", border = NA, xpd = NA)
          rect(xleft = -7.5, ybottom = 2.8, xright = 7.5, ytop = 3.4, col = "lightgrey", border = NA, xpd = NA)
          rect(xleft = -7.5, ybottom = 3.4, xright = 7.5, ytop = 4, col = "grey", border = NA, xpd = NA)
          
          # Add points for the generic alpha terms
          points(x = GenericDeviations[[i]][1,1], Generic_y[1], pch = InterceptPoints, col = InterceptCol)
          points(x = GenericDeviations[[i]][2,1], Generic_y[2], pch = SlopePoints, col = SlopeCol)
          segments(x0 = GenericDeviations[[i]][1,2], y0 = Generic_y[1], x1 = GenericDeviations[[i]][1,3],
                   y1 = Generic_y[1], col = InterceptCol)
          segments(x0 = GenericDeviations[[i]][2,2], y0 = Generic_y[2], x1 = GenericDeviations[[i]][2,3],
                   y1 = Generic_y[2], col = SlopeCol)
          # Add points for the non-generic alpha terms
          for(j in 1:NumNonGeneric){
               if(InterceptInclusion[[i]][AllNonGeneric[j]] == 1){
                    points(x = Alphas[[i]][1,1,AllNonGeneric[j]], y = NonGeneric_y[j,1], pch = InterceptPoints,
                           col = InterceptCol)
                    segments(x0 = Alphas[[i]][2,1,AllNonGeneric[j]], y0 = NonGeneric_y[j,1],
                             x1 = Alphas[[i]][3,1,AllNonGeneric[j]], y1 = NonGeneric_y[j,1], col = InterceptCol)
               }
               if(SlopeInclusion[[i]][AllNonGeneric[j]] == 1){
                    points(x = Alphas[[i]][1,2,AllNonGeneric[j]], y = NonGeneric_y[j,2], pch = SlopePoints,
                           col = SlopeCol)
                    segments(x0 = Alphas[[i]][2,2,AllNonGeneric[j]], y0 = NonGeneric_y[j,2],
                             x1 = Alphas[[i]][3,2,AllNonGeneric[j]], y1 = NonGeneric_y[j,2], col = SlopeCol)
               }
          }
          # Add points for the intraspecific alpha terms
          points(x = Alphas[[i]][1,1,Focal], Intra_y[1], pch = InterceptPoints, col = InterceptCol)
          points(x = Alphas[[i]][1,2,Focal], Intra_y[2], pch = SlopePoints, col = SlopeCol)
          segments(x0 = Alphas[[i]][2,1,Focal], x1 = Alphas[[i]][3,1,Focal], y0 = Intra_y[1], y1 = Intra_y[1], col = InterceptCol)
          segments(x0 = Alphas[[i]][2,2,Focal], x1 = Alphas[[i]][3,2,Focal], y0 = Intra_y[2], y1 = Intra_y[2], col = SlopeCol)
          # Lambda
          points(x = Lambdas[[i]][1,1], Lambda_y[1], pch = InterceptPoints, col = InterceptCol)
          points(x = Lambdas[[i]][1,2], Lambda_y[2], pch = SlopePoints, col = SlopeCol)
          segments(x0 = Lambdas[[i]][2,1], x1 = Lambdas[[i]][3,1], y0 = Lambda_y[1], y1 = Lambda_y[1], col = InterceptCol)
          segments(x0 = Lambdas[[i]][2,2], x1 = Lambdas[[i]][3,2], y0 = Lambda_y[2], y1 = Lambda_y[2], col = SlopeCol)
          # Add the labels
          text(x = -5.5, y = NonGenericLab_y1, labels = NonGenericLab1, xpd = NA, cex = 1.5)
          text(x = -5.5, y = NonGenericLab_y2, labels = NonGenericLab2, xpd = NA, cex = 1.5)
          text(x = -5.5, y = NonGenericLab_y3, labels = NonGenericLab3, xpd = NA, cex = 1.5)
          text(x = 5, y = GenericLab_y, labels = GenericLab, xpd = NA, cex = 1.5)
          text(x = 5, y = IntraLab_y, labels = IntraLab, xpd = NA, cex = 1.5)
          text(x = 5, y = LambdaLab_y, labels = LambdaLab, xpd = NA, cex = 1.5)
          segments(x0 = 0, y0 = 0, x1 = 0, y1 = 6, col = "black", lty = 2)
          if(i == 3 | i == 6){
               mtext("Parameter deviations", side = 1, line = 2.5)
          }
          if(i == 1){
               legend(x = 4.5, y = 7.5, xpd = NA, bty = "n", horiz = TRUE, legend = c("Intercept"),
                      pch = c(InterceptPoints), col = c(InterceptCol), cex = 1.5)
               legend(x = 9.75, y = 7.5, xpd = NA, bty = "n", horiz = TRUE, legend = c("Slope"),
                      pch = c(SlopePoints), col = c(SlopeCol), cex = 1.5)
          }
          # Add the letters to the figure
          text(x = 0.95*param_xRange[1], y = 0.95*param_yRange[2], 
               labels = Letters[i], cex = 1.5)
     }
dev.off()

taus
CIs

