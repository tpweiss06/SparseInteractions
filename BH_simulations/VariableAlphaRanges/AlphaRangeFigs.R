# This script will create the figure showing model performance across 150 independent
#   simulations for both a narrower and wider range of alpha values compared to
#   the results shown in the main simulation.

rm(list = ls())
library(RColorBrewer)
library(here)

# First load in the low range alphas and combine them
Round1 <- read.csv(here("BH_simulations/VariableAlphaRanges/LowAlphaTestResults.csv"))
Round2 <- read.csv(here("BH_simulations/VariableAlphaRanges/LowAlphaTestResults2.csv"))
LowRangeAlpha <- rbind(Round1, Round2)
# Randomly select 150 from here
PossibleRows <- which(!is.na(LowRangeAlpha$LambdaIntDev))
RowsToUse <- sample(PossibleRows, size = 150, replace = FALSE)
LowAlphas <- LowRangeAlpha[RowsToUse,]

# Then I'll load in the full results, subset to 50, and randomly select 150
MidRangeAlpha <- read.csv(here("BH_simulations/Main/StanFits/MainSimFits.csv"))
MidRangeAlpha <- subset(MidRangeAlpha, Size == 50)
RowsToUse <- sample(1:nrow(MidRangeAlpha), size = 150, replace = FALSE)
MidAlphas <- MidRangeAlpha[RowsToUse,]

# Finally, I'll load in the high range alphas and randomly select 150
HighRangeAlpha <- read.csv(here("BH_simulations/VariableAlphaRanges/HighAlphaResults.csv"))
PossibleRows <- which(!is.na(HighRangeAlpha$LambdaIntDev))
RowsToUse <- sample(PossibleRows, size = 150, replace = FALSE)
HighAlphas <- HighRangeAlpha[RowsToUse,]

# Combine the data into a single data frame
LowAlphas$Range <- rep("Low", 150)
MidAlphas$Range <- rep("Medium", 150)
HighAlphas$Range <- rep("High", 150)

results <- rbind(LowAlphas, MidAlphas, HighAlphas)

# Create useful objects for graphing
ranges <- unique(results$Range)

xSeqs <- vector(mode = "list", length = 3)
xSeqs[[1]] <- matrix(NA, nrow = 3, ncol = 8)
xSeqs[[1]][1,] <- c(0.6, 0.7, 0.85, 0.95, 1.1, 1.2, 1.35, 1.45)
xSeqs[[1]][2,] <- xSeqs[[1]][1,] + 1
xSeqs[[1]][3,] <- xSeqs[[1]][1,] + 2

xSeqs[[2]] <- matrix(NA, nrow = 3, ncol = 3)
xSeqs[[2]][1,] <- c(0.75, 1, 1.25)
xSeqs[[2]][2,] <- xSeqs[[2]][1,] + 1
xSeqs[[2]][3,] <- xSeqs[[2]][1,] + 2

xSeqs[[3]] <- c(1,2,3)

xRange <- range(xSeqs[[1]])

# Set colors
Dark2Cols <- brewer.pal(n = 8, name = "Dark2")
InterceptCol <- Dark2Cols[1]
SlopeCol <- Dark2Cols[2]
RMSECol <- Dark2Cols[3]

# Make vectors for the colors for the first plot
ParamCols <- rep(c(InterceptCol, SlopeCol), 4)

# Define the points for parameter labels in the first plot
LambdaCoords <- c(0.65, 1.65, 2.65)
IntraCoords <- c(0.9, 1.9, 2.9)
GenericCoords <- c(1.15, 2.15, 3.15)
NonGenericCoords1 <- c(1.35, 2.35, 3.35)
NonGenericCoords2 <- c(1.45, 2.45, 3.45)

# Define the labels for the parameters
LabCoords <- rep(5.65, 3)
LambdaLab <- expression(paste(lambda["e,i"], sep = ""))
GenericLab <- "Generic"
IntraLab <- "Intraspecific"
NonGenericLab1 <- "Species"
NonGenericLab2 <- "specific"

# Set the yRanges for the different graphs
yRanges <- vector(mode = "list", length = 3)
yRanges[[1]] <- c(-8,8)
yRanges[[2]] <- c(0, 6)
yRanges[[3]] <- c(0, 0.8)

# Make an x axis label
xLab <- expression(paste(alpha, " range", sep = ""))

# Make the letters for the subpanels
Letters <- c(expression(bold("a")), expression(bold("b")), expression(bold("c")))

FigName <- here("BH_simulations/VariableAlphaRanges/AlphaRangeResults.pdf")
pdf(file = FigName, width = 14, height = 4, onefile = FALSE, paper = "special")
     layout(mat = matrix(data = c(1,1,1,1,2,2,3,3), nrow = 1, ncol = 8))
     # First make the plot for the parameter deviations
     plot(x = NA, y = NA, xlim = xRange, ylim = yRanges[[1]], main = "", xlab = xLab,
          ylab = "Parameter deviation", las = 1, axes = FALSE, cex.lab = 1.5)
     box()
     axis(side = 1, at = 1:3, tcl = -0.5, labels = ranges, cex.axis = 1.5)
     axis(side = 2, at = seq(-8,8, by = 2), tcl = -0.5, las = 1, cex.axis = 1.5)
     axis(side = 2, at = seq(-7,7, by = 2), tcl = -0.5, labels = FALSE)
     axis(side = 2, at = seq(-9,9, by = 0.25), tcl = -0.25, labels = FALSE)
     # Add the boxes
     for(p in 1:8){
          boxplot(results[,2+p] ~ results$Range, add = TRUE, at = xSeqs[[1]][,p], boxwex = 0.075,
                  axes = FALSE, col = ParamCols[p])
     }
     # Add the lines for clarity
     abline(v = 1.525, col = "black")
     abline(v = 2.525, col = "black")
     abline(h = 0, lty = 2)
     abline(v = c(0.775, 1.775, 2.775), col = "grey30", lty = 3)
     abline(v = c(1.025, 2.025, 3.025), col = "grey30", lty = 3)
     abline(v = c(1.275, 2.275, 3.275), col = "grey30", lty = 3)
     
     # Add the parameter labels 
     text(x = LambdaCoords, y = LabCoords, labels = LambdaLab, cex = 1.5)
     text(x = IntraCoords, y = LabCoords, labels = IntraLab, srt = 90, cex = 1.5)
     text(x = GenericCoords, y = LabCoords, labels = GenericLab, srt = 90, cex = 1.5)
     text(x = NonGenericCoords1, y = LabCoords, labels = NonGenericLab1, srt = 90, cex = 1.5)
     text(x = NonGenericCoords2, y = LabCoords, labels = NonGenericLab2, srt = 90, cex = 1.5)
     
     legend(x = "top", legend = c("Intercept", "Slope"), horiz = TRUE, xpd = NA,
            bty = "n", pch = 15, col = c(InterceptCol, SlopeCol), inset = -0.125, cex = 1.5)
     # Add the letters to the figure
     text(x = 0.5, y = 9.35, labels = Letters[1], cex = 1.5, xpd = NA)
     
    # Now add the plot for the number of non-generic species
     plot(x = NA, y = NA, xlim = xRange, ylim = yRanges[[2]], main = "", xlab = xLab,
          ylab = "Number of species", las = 1, axes = FALSE, cex.lab = 1.5)
     box()
     axis(side = 1, at = 1:3, tcl = -0.5, labels = ranges, cex.axis = 1.5)
     axis(side = 2, at = 0:9, tcl = -0.5, las = 1, cex.axis = 1.5)
     boxplot(results$NumNonGenericInt ~ results$Range, add = TRUE, at = xSeqs[[2]][,1], boxwex = 0.125,
             axes = FALSE, col = InterceptCol)
     boxplot(results$NumNonGenericSlope ~ results$Range, add = TRUE, at = xSeqs[[2]][,2], boxwex = 0.125,
             axes = FALSE, col = SlopeCol)
     boxplot(results$NumNonGenericSpecies ~ results$Range, add = TRUE, at = xSeqs[[2]][,3], boxwex = 0.125,
             axes = FALSE)
     legend(x = "topleft", legend = c("Species-specific intercept", "Species-specific slope", "Either"),
            bty = "n", pch = 15, col = c(InterceptCol, SlopeCol, "grey"), cex = 1.5)
     # Add the letters to the figure
     text(x = 0.5, y = 9.75, labels = Letters[2], cex = 1.5, xpd = NA)
     
     # Finally add the plot for RMSE
     boxplot(results$RMSE ~ results$Range, col = RMSECol, xlab = xLab, ylab = "RMSE", main = "",
             ylim = yRanges[[3]], cex.lab = 1.5, cex.axis = 1.5, las = 1)
     axis(side = 2, at = seq(0, 1.4, by = 0.05), tcl = -0.25, labels = FALSE)
     # Add the letters to the figure
     text(x = 0.5, y = 1.4, labels = Letters[3], cex = 1.5, xpd = NA)
dev.off()

