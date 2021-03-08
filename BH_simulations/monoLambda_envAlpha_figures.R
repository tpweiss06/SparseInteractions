# This script will make some possible manuscript figures for the scenario with 
#    a monotonic lambda and alphas varying with the environment.

setwd("~/Desktop/Wyoming/SparseInteractions/BH_simulations/")

# Create a list to hold the results for each size of dataset
AllResults <- vector(mode = "list", length = 5)

load("StanFits/monoLambda_envAlpha/N10_GraphStuff.rdata")
AllResults[[1]] <- GraphStuff
rm(GraphStuff)
load("StanFits/monoLambda_envAlpha/N20_GraphStuff.rdata")
AllResults[[2]] <- GraphStuff
rm(GraphStuff)
load("StanFits/monoLambda_envAlpha/N50_GraphStuff.rdata")
AllResults[[3]] <- GraphStuff
rm(GraphStuff)
load("StanFits/monoLambda_envAlpha/N100_GraphStuff.rdata")
AllResults[[4]] <- GraphStuff
rm(GraphStuff)
load("StanFits/monoLambda_envAlpha/N200_GraphStuff.rdata")
AllResults[[5]] <- GraphStuff
rm(GraphStuff)

# First, try a ppc graph
Pred_yRange <- c(-3.5,2.75)
Dev_yRange <- c(-1.2,1.2)
xRange <- range(AllResults[[1]]$ppc_x)
ppcCol <- "purple"
Sizes <- c(10, 20, 50, 100, 200)
pdf(file = "Results/monoLambda_envAlpha/ppc_all.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,5))
     # ppc estimates
     for(i in 1:5){
          plot(x = NA, y = NA, xlim = xRange, ylim = Pred_yRange, xlab = "",
               ylab = "", las = 1)
          points(x = AllResults[[i]]$ppc_x, y = AllResults[[i]]$ppc_pred[1,], pch = 1, col = ppcCol)
          segments(x0 = AllResults[[i]]$ppc_x, y0 = AllResults[[i]]$ppc_pred[2,], 
                   x1 = AllResults[[i]]$ppc_x, y1 = AllResults[[i]]$ppc_pred[3,], col = ppcCol)
          abline(a = 0, b = 1, lty = 2)
          mtext(Sizes[i], side = 3, line = 1)
          if(i == 1){
               mtext("Predicted growth", side = 2, line = 2)
          }
     }
     # ppc deviations
     for(i in 1:5){
          plot(x = NA, y = NA, xlim = xRange, ylim = Dev_yRange, xlab = "",
               ylab = "", las = 1)
          points(x = AllResults[[i]]$ppc_x, y = AllResults[[i]]$ppc_dev[1,], pch = 1, col = ppcCol)
          segments(x0 = AllResults[[i]]$ppc_x, y0 = AllResults[[i]]$ppc_dev[2,], 
                   x1 = AllResults[[i]]$ppc_x, y1 = AllResults[[i]]$ppc_dev[3,], col = ppcCol)
          abline(h = 0, lty = 2)
          if(i == 1){
               mtext("Deviation", side = 2, line = 2)
          }
     }
     mtext("True growth", side = 1, line = -2, outer = TRUE)
dev.off()

# Now make a graph for the parameter estimates (alpha and lambda together)
AlphaRange <- c(-3,3)
yRange <- c(0.5, S+1.5)
InterceptSeq <- 1:S - 0.15
SlopeSeq <- 1:S + 0.15

LambdaRange <- c(-3, 3)
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
estCols <- Dark2Cols[1:2]
Sizes <- c(10, 20, 50, 100, 200)
S <- 14
pdf(file = "Results/monoLambda_envAlpha/params_all.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,5))
     # alpha estimates
     for(i in 1:5){
          plot(x = NA, y = NA, xlim = AlphaRange, ylim = yRange, xlab = "", ylab = "", las = 1, yaxt = "n")
          points(y = InterceptSeq, x = AllResults[[i]]$AlphaIntercept[1,OrderedIntercepts$species], col = estCols[1])
          points(y = SlopeSeq, x = AllResults[[i]]$AlphaSlope[1,OrderedIntercepts$species], col = estCols[2])
          segments(y0 = InterceptSeq, x0 = AllResults[[i]]$AlphaIntercept[2,OrderedIntercepts$species], y1 = InterceptSeq,
                   x1 = AllResults[[i]]$AlphaIntercept[3,OrderedIntercepts$species], col = estCols[1])
          segments(y0 = SlopeSeq, x0 = AllResults[[i]]$AlphaSlope[2,OrderedIntercepts$species], y1 = SlopeSeq,
                   x1 = AllResults[[i]]$AlphaSlope[3,OrderedIntercepts$species], col = estCols[2])
          # Add the intra values
          points(y = S+1 - 0.15, x = AllResults[[i]]$AlphaIntercept[1,S+1], col = estCols[1], pch = 16)
          points(y = S+1 + 0.15, x = AllResults[[i]]$AlphaSlope[1,S+1], col = estCols[2], pch = 16)
          segments(y0 = S+1 - 0.15, x0 = AllResults[[i]]$AlphaIntercept[2,S+1], y1 = S+1 - 0.15,
                   x1 = AllResults[[i]]$AlphaIntercept[3,S+1], col = estCols[1])
          segments(y0 = S+1 + 0.15, x0 = AllResults[[i]]$AlphaSlope[2,S+1], y1 = S+1 + 0.15,
                   x1 = AllResults[[i]]$AlphaSlope[3,S+1], col = estCols[2])
          abline(v = 0, lty = 2)
          mtext(Sizes[i], side = 3, line = 1)
          if(i == 1){
               mtext("Species ordered by \n competitive effect", side = 2, line = 2) 
          }
     }
     # lambda estimates
     for(i in 1:5){
          plot(density(AllResults[[i]]$LambdaIntercept), main = "", xlab = "", 
               col = estCols[1], xlim = LambdaRange, las = 1, ylab = "")
          lines(density(AllResults[[i]]$LambdaSlope), col = estCols[2])
          abline(v = hdi(AllResults[[i]]$LambdaIntercept), lty = 3, col = estCols[1])
          abline(v = mean(AllResults[[i]]$LambdaIntercept), lty = 1, col = estCols[1])
          abline(v = hdi(AllResults[[i]]$LambdaSlope), lty = 3, col = estCols[2])
          abline(v = mean(AllResults[[i]]$LambdaSlope), lty = 1, col = estCols[2])
          abline(v = 0, lty = 2)
          if(i == 1){
               mtext("Density", side = 2, line = 2)
          }
     }
     mtext("Parameter deviations", side = 1, line = -2, outer = TRUE)
dev.off()

# Now make a graph for the parameter estimates (alphas vs. true values)
InterceptRange
SlopeRange
pdf(file = "Results/monoLambda_envAlpha/alphas_all.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(5,4,2,2) + 0.1, mfrow = c(2,5))
     # alpha intercepts
     for(i in 1:5){
          plot(x = NA, y = NA, xlim = range(TrueAlphaMeans), ylim = InterceptRange, 
               xlab = "", ylab = "", las = 1)
          points(x = TrueAlphaMeans[OrderedIntercepts$species], 
                 y = AllResults[[i]]$AlphaIntercept[1,OrderedIntercepts$species],
                 col = estCols[1])
          segments(x0 = TrueAlphaMeans[OrderedIntercepts$species], 
                   y0 = AllResults[[i]]$AlphaIntercept[2,OrderedIntercepts$species],
                   y1 = AllResults[[i]]$AlphaIntercept[3,OrderedIntercepts$species], col = estCols[1])
          # Add the intraspecific values
          points(x = TrueAlphaMeans[Focal], y = AllResults[[i]]$AlphaIntercept[1,S+1], 
                 pch = 16, col = estCols[1])
          segments(x0 = TrueAlphaMeans[Focal], y0 = AllResults[[i]]$AlphaIntercept[2,S+1],
                   y1 = AllResults[[i]]$AlphaIntercept[3,S+1], col = estCols[1])
          abline(a = 0, b = 1, lty = 2)
          mtext(Sizes[i], side = 3, line = 1)
          if(i == 1){
               mtext("Estimated alpha intercept", side = 2, line = 2) 
          }
     }
     # alpha slopes
     for(i in 1:5){
          plot(x = NA, y = NA, xlim = range(TrueAlphaSlopes), ylim = SlopeRange, 
               xlab = "", ylab = "", las = 1)
          points(x = TrueAlphaSlopes[OrderedIntercepts$species], 
                 y = AllResults[[i]]$AlphaSlope[1,OrderedIntercepts$species],
                 col = estCols[1])
          segments(x0 = TrueAlphaSlopes[OrderedIntercepts$species], 
                   y0 = AllResults[[i]]$AlphaSlope[2,OrderedIntercepts$species],
                   y1 = AllResults[[i]]$AlphaSlope[3,OrderedIntercepts$species], col = estCols[1])
          # Add the intraspecific values
          points(x = TrueAlphaSlopes[Focal], y = AllResults[[i]]$AlphaSlope[1,S+1], 
                 pch = 16, col = estCols[1])
          segments(x0 = TrueAlphaSlopes[Focal], y0 = AllResults[[i]]$AlphaSlope[2,S+1],
                   y1 = AllResults[[i]]$AlphaSlope[3,S+1], col = estCols[1])
          abline(a = 0, b = 1, lty = 2)
          mtext(Sizes[i], side = 3, line = 1)
          if(i == 1){
               mtext("Estimated alpha slope", side = 2, line = 2) 
          }
     }
     mtext("True value", side = 1, line = -2, outer = TRUE)
dev.off()

