# This script will make the manuscript figures we agreed on for the sparse
#       interactions paper

library(RColorBrewer)
setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Create lists for the results from different sizes of datasets
ppcPreds <- vector(mode = "list", length = 3)
Alphas <- vector(mode = "list", length = 3)
Lambdas <- vector(mode = "list", length = 3)
InterceptInclusion <- vector(mode = "list", length = 3)
SlopeInclusion <- vector(mode = "list", length = 3)
GenericDeviations <- vector(mode = "list", length = 3)

# Now population the objects
FilePrefixes <- c("N10_", "N50_", "N200_")
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

# Now make the simulation figure: a 6 panel figure divided between the posterior predictive
#       graphs and graphs of the parameter deviations. 3 rows and 2 columns. In the
#       parameter estimate graph, there will be alternating horizontal bands to delineate
#       (1) lambda parameters, (2) intraspecific alphas, (3) model estimated species-specific
#       alpha terms, and (4) generic alpha terms vs. the weighted averages they are
#       estimating. Points will be colored as intercept or slope. Lambda points will be
#       filled squares. Intraspecific alpha points will be filled triangles. Species specific
#       interspecific alpha terms estimated by the model will be filled circles. The generic
#       interspecific alpha terms estimated by the model will be open circles.
Focal <- 10
S <- 15
ppc_xRange <- range(Growth_ppc)
ppc_yRange <- c(-3.5, 1)
param_xRange <- c(-7.5, 7.5)
param_yRange <- c(0, 5.65)

# Determine the maximum number of species with non-generic alpha terms
AllNonGeneric <- union(which(InterceptInclusion[[3]] == 1), which(SlopeInclusion[[3]] == 1))
NumNonGeneric <- length(AllNonGeneric)

# Set some pch and color values for the graph
InterceptPoints <- 15
SlopePoints <- 17
# LambdaPoints <- 15
# IntraPoints <- 17
# GenericPoints <- 1
# SpecificPoints <- 16
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
InterceptCol <- Dark2Cols[1]
SlopeCol <- Dark2Cols[2]
ppcCol <- Dark2Cols[3]

# Divide up the plot space into appropriate vertical subdivisions
InterceptOffset <- 0.55
SlopeOffset <- 0.45
Generic_y <- c(InterceptOffset, SlopeOffset)
NonGeneric_y <- matrix(data = NA, nrow = NumNonGeneric, ncol = 2)
for(i in 2:(NumNonGeneric+1)){
     NonGeneric_y[i-1,1] <- 0.5*i + InterceptOffset
     NonGeneric_y[i-1,2] <- 0.5*i + SlopeOffset
}
Intra_y <- c(4 + InterceptOffset, 4 + SlopeOffset)
Lambda_y <- c(5 + InterceptOffset, 5 + SlopeOffset)

# Lambda_y <- c(15 + InterceptOffset, 15 + SlopeOffset)
# Intra_y <- c(14 + InterceptOffset, 14 + SlopeOffset)
# NumBoth <- length(BothNonGeneric)
# Both_y <- matrix(data = NA, nrow = NumBoth, ncol = 2)
# for(i in 1:NumBoth){
#         Both_y[i,1] <- 14-i + InterceptOffset
#         Both_y[i,2] <- 14-i + SlopeOffset
# }
# NumIntercept <- length(OnlyIntercept)
# Intercept_y <- matrix(data = NA, nrow = NumIntercept, ncol = 2)
# for(i in 1:NumIntercept){
#         Intercept_y[i,1] <-  14 - NumBoth - i + InterceptOffset
#         Intercept_y[i,2] <-  14 - NumBoth - i + SlopeOffset
# }
# NumSlope <- length(OnlySlope)
# Slope_y <- matrix(data = NA, nrow = NumSlope, ncol = 2)
# for(i in 1:NumSlope){
#         Slope_y[i,1] <-  14 - NumBoth - NumIntercept - i + InterceptOffset
#         Slope_y[i,2] <-  14 - NumBoth - NumIntercept - i + SlopeOffset
# }
# NumGeneric <- length(BothGeneric)
# Generic_y <- matrix(data = NA, nrow = NumGeneric, ncol = 2)
# for(i in 1:NumGeneric){
#         Generic_y[i,1] <-  14 - NumBoth - NumIntercept - NumSlope - i + InterceptOffset
#         Generic_y[i,2] <-  14 - NumBoth - NumIntercept - NumSlope - i + SlopeOffset
# }

# Define the points for y axis labels
GenericLab_y <- 0.5
GenericLab <- expression(paste("Generic ", alpha["e,i,j"], " terms", sep = ""))
IntraLab_y <- 4.5
IntraLab <- expression(paste(alpha["e,i,i"], " terms", sep = ""))
LambdaLab_y <- 5.5
LambdaLab <- expression(paste(lambda["e,i"], " terms", sep = ""))
NonGenericLab_y1 <- 2.75
NonGenericLab_y2 <- 2.25
NonGenericLab1 <- "Species specific"
NonGenericLab2 <- expression(paste(alpha["e,i,j"], " terms", sep = ""))
     
# BothGenericLab <- (28 - 2*NumBoth - 2*NumIntercept - 2*NumSlope - NumGeneric) / 2
# OnlySlopeLab <- (28 - 2*NumBoth - 2*NumIntercept - NumSlope) / 2
# OnlyInterceptLab <- (28 - 2*NumBoth - NumIntercept) / 2
# BothNonGenericLab <- (28 - NumBoth) / 2
# IntraLab <- 14.5
# LambdaLab <- 15.5


FigName <- "Results/monoLambda_envAlpha/SimResults.pdf"
pdf(file = FigName, width = 7, height = 8, onefile = FALSE, paper = "special")
     par(mfrow = c(3,2), mar = c(2, 2.5, 2, 2.5), oma = c(2, 0, 2, 0))
     for(i in 1:3){
          # Plot the parameters
          plot(NA, NA, xlim = param_xRange, ylim = param_yRange, xlab = "",
               ylab = "", axes = FALSE)
          axis(side = 1, at = -7:7, tcl = -0.5)
          axis(side = 1, at = seq(-8, 8, by = 0.5), tcl = -0.25, labels = FALSE)
          # Add boxes for intraspecifc and generic alpha components
          rect(xleft = -7.5, ybottom = 5, xright = 7.5, ytop = 6, col = "grey", border = NA, xpd = NA)
          rect(xleft = -7.5, ybottom = 1, xright = 7.5, ytop = 4, col = "grey", border = NA, xpd = NA)
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
          # Add points for the ntraspecific alpha terms
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
          text(x = 4, y = NonGenericLab_y1, labels = NonGenericLab1, xpd = NA)
          text(x = 4, y = NonGenericLab_y2, labels = NonGenericLab2, xpd = NA)
          text(x = 4, y = GenericLab_y, labels = GenericLab, xpd = NA)
          text(x = 4, y = IntraLab_y, labels = IntraLab, xpd = NA)
          text(x = 4, y = LambdaLab_y, labels = LambdaLab, xpd = NA)
          abline(v = 0, lty = 2)
          if(i == 3){
               mtext("Parameter deviations", side = 1, line = 2.5)
          }
          if(i == 1){
               legend(x = -5, y = 7, xpd = NA, bty = "n", horiz = TRUE, legend = c("Intercept", "", "Slope"),
                      pch = c(InterceptPoints, NA, SlopePoints), col = c(InterceptCol, NA, SlopeCol))
          }
     
          # Now plot the posterior predictive values
          plot(x = NA, y = NA, xlim = ppc_xRange, ylim = ppc_yRange, xlab = "",
               ylab = "", las = 1)
          mtext("Predicted growth", side = 2, line = 2.5)
          axis(1, at = seq(-4, 1, by = 0.25), tcl = -0.25, labels = FALSE)
          axis(2, at = seq(-4, 1.5, by = 0.25), tcl = -0.25, labels = FALSE)
          points(x = Growth_ppc, y = ppcPreds[[i]][1,], pch = 1, col = ppcCol)
          segments(x0 = Growth_ppc, y0 = ppcPreds[[i]][2,], 
                   x1 = Growth_ppc, y1 = ppcPreds[[i]][3,], col = ppcCol)
          abline(a = 0, b = 1, lty = 2)
          if(i == 3){
               mtext("True growth", side = 1, line = 2.5)
          }
     }
dev.off()


# FigName <- "Results/monoLambda_envAlpha/SimResults.pdf"
# pdf(file = FigName, width = 7, height = 8, onefile = FALSE, paper = "special")
#         par(mfrow = c(3,2), mar = c(2, 2.5, 2, 2.5), oma = c(2, 0, 0, 0))
#         for(i in 1:3){
#                 # First define the point types for interspecific interactions
#                 InterceptPoints <- rep(1, S)
#                 SlopePoints <- rep(1, S)
#                 InterceptPoints[which(InterceptInclusion[[i]] == 1)] <- 16
#                 SlopePoints[which(SlopeInclusion[[i]] == 1)] <- 16
#                 
#                 # Plot the parameters
#                 plot(NA, NA, xlim = param_xRange, ylim = param_yRange, xlab = "",
#                      ylab = "", axes = FALSE)
#                 axis(side = 1, at = -3:3, tcl = -0.5)
#                 axis(side = 1, at = seq(-3, 3, by = 0.25), tcl = -0.25, labels = FALSE)
#                 # Add boxes for intraspecifc, only intercept, and both generic
#                 rect(xleft = -3, ybottom = 14, xright = 3, ytop = 15, col = "grey",
#                      border = NA, xpd = NA)
#                 rect(xleft = -3, ybottom = 14 - NumBoth - NumIntercept, xright = 3, 
#                      ytop = 14 - NumBoth, col = "grey", border = NA, xpd = NA)
#                 rect(xleft = -3, ybottom = 14 - NumBoth - NumIntercept - NumSlope - NumGeneric,
#                      xright = 3, ytop = 14 - NumBoth - NumIntercept - NumSlope, 
#                      col = "grey", border = NA, xpd = NA)
#                 # Both generic
#                 points(x = Alphas[[i]][1,1,BothGeneric], Generic_y[,1], 
#                        pch = InterceptPoints[BothGeneric], col = InterceptCol)
#                 points(x = Alphas[[i]][1,2,BothGeneric], Generic_y[,2],
#                        pch = SlopePoints[BothGeneric], col = SlopeCol)
#                 segments(x0 = Alphas[[i]][2,1,BothGeneric], x1 = Alphas[[i]][3,1,BothGeneric],
#                          y0 = Generic_y[,1], y1 = Generic_y[,1], col = InterceptCol)
#                 segments(x0 = Alphas[[i]][2,2,BothGeneric], x1 = Alphas[[i]][3,2,BothGeneric],
#                          y0 = Generic_y[,2], y1 = Generic_y[,2], col = SlopeCol)
#                 # Only slope
#                 points(x = Alphas[[i]][1,1,OnlySlope], Slope_y[,1], 
#                        pch = InterceptPoints[OnlySlope], col = InterceptCol)
#                 points(x = Alphas[[i]][1,2,OnlySlope], Slope_y[,2], 
#                        pch = SlopePoints[OnlySlope], col = SlopeCol)
#                 segments(x0 = Alphas[[i]][2,1,OnlySlope], x1 = Alphas[[i]][3,1,OnlySlope],
#                          y0 = Slope_y[,1], y1 = Slope_y[,1], col = InterceptCol)
#                 segments(x0 = Alphas[[i]][2,2,OnlySlope], x1 = Alphas[[i]][3,2,OnlySlope],
#                          y0 = Slope_y[,2], y1 = Slope_y[,2], col = SlopeCol)
#                 # Only intercept
#                 points(x = Alphas[[i]][1,1,OnlyIntercept], Intercept_y[,1], 
#                        pch = InterceptPoints[OnlyIntercept], col = InterceptCol)
#                 points(x = Alphas[[i]][1,2,OnlyIntercept], Intercept_y[,2], 
#                        pch = SlopePoints[OnlyIntercept], col = SlopeCol)
#                 segments(x0 = Alphas[[i]][2,1,OnlyIntercept], x1 = Alphas[[i]][3,1,OnlyIntercept],
#                          y0 = Intercept_y[,1], y1 = Intercept_y[,1], col = InterceptCol)
#                 segments(x0 = Alphas[[i]][2,2,OnlyIntercept], x1 = Alphas[[i]][3,2,OnlyIntercept],
#                          y0 = Intercept_y[,2], y1 = Intercept_y[,2], col = SlopeCol)
#                 # Both non-generic
#                 points(x = Alphas[[i]][1,1,BothNonGeneric], Both_y[,1], 
#                        pch = InterceptPoints[BothNonGeneric], col = InterceptCol)
#                 points(x = Alphas[[i]][1,2,BothNonGeneric], Both_y[,2], 
#                        pch = SlopePoints[BothNonGeneric], col = SlopeCol)
#                 segments(x0 = Alphas[[i]][2,1,BothNonGeneric], x1 = Alphas[[i]][3,1,BothNonGeneric],
#                          y0 = Both_y[,1], y1 = Both_y[,1], col = InterceptCol)
#                 segments(x0 = Alphas[[i]][2,2,BothNonGeneric], x1 = Alphas[[i]][3,2,BothNonGeneric],
#                          y0 = Both_y[,2], y1 = Both_y[,2], col = SlopeCol)
#                 # Intraspecific
#                 points(x = Alphas[[i]][1,1,Focal], Intra_y[1], 
#                        pch = IntraPoints, col = InterceptCol)
#                 points(x = Alphas[[i]][1,2,Focal], Intra_y[2], 
#                        pch = IntraPoints, col = SlopeCol)
#                 segments(x0 = Alphas[[i]][2,1,Focal], x1 = Alphas[[i]][3,1,Focal],
#                          y0 = Intra_y[1], y1 = Intra_y[1], col = InterceptCol)
#                 segments(x0 = Alphas[[i]][2,2,Focal], x1 = Alphas[[i]][3,2,Focal],
#                          y0 = Intra_y[2], y1 = Intra_y[2], col = SlopeCol)
#                 # Lambda
#                 points(x = Lambdas[[i]][1,1], Lambda_y[1], 
#                        pch = LambdaPoints, col = InterceptCol)
#                 points(x = Lambdas[[i]][1,2], Lambda_y[2], 
#                        pch = LambdaPoints, col = SlopeCol)
#                 segments(x0 = Lambdas[[i]][2,1], x1 = Lambdas[[i]][3,1],
#                          y0 = Lambda_y[1], y1 = Lambda_y[1], col = InterceptCol)
#                 segments(x0 = Lambdas[[i]][2,2], x1 = Lambdas[[i]][3,2],
#                          y0 = Lambda_y[2], y1 = Lambda_y[2], col = SlopeCol)
#                 # Add the labels
#                 text(x = -2, y = BothGenericLab, labels = "Fully generic", xpd = NA)
#                 text(x = -2, y = OnlySlopeLab, labels = "Non-generic \n slope", xpd = NA)
#                 text(x = -2, y = OnlyInterceptLab, labels = "Non-generic \n intercept", xpd = NA)
#                 text(x = -2, y = BothNonGenericLab, labels = "Fully \n non-generic", xpd = NA)
#                 text(x = -2, y = IntraLab, labels = "Intraspecific", xpd = NA)
#                 text(x = -2, y = LambdaLab, labels = "Lambda", xpd = NA)
#                 abline(v = 0, lty = 2)
#                 if(i == 3){
#                         mtext("Parameter deviations", side = 1, line = 2.5)
#                 }
#                 
#                 # Now plot the posterior predictive values
#                 plot(x = NA, y = NA, xlim = ppc_xRange, ylim = ppc_yRange, xlab = "",
#                      ylab = "", las = 1)
#                 mtext("Predicted growth", side = 2, line = 2.5)
#                 axis(1, at = seq(-2, 1, by = 0.25), tcl = -0.25, labels = FALSE)
#                 axis(2, at = seq(-3, 1.5, by = 0.25), tcl = -0.25, labels = FALSE)
#                 points(x = Growth_ppc, y = ppcPreds[[i]][1,], pch = 1, col = ppcCol)
#                 segments(x0 = Growth_ppc, y0 = ppcPreds[[i]][2,], 
#                          x1 = Growth_ppc, y1 = ppcPreds[[i]][3,], col = ppcCol)
#                 abline(a = 0, b = 1, lty = 2)
#                 if(i == 3){
#                         mtext("True growth", side = 1, line = 2.5)
#                 }
#         }
# dev.off()
# 
# 
# 
# 
