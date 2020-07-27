# Use this script to evaluate the current fit of the model from the server
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("ModelFits/Current_BH_FH_Phos_fit.rdata")

# Evaluate the model fit for convergence, mixing, autocorrelation, etc.
PrelimFit
PrelimPosteriors <- rstan::extract(PrelimFit)
# Diagnostic plots
NonShrinkageParams <- c("lambdas[1]", "lambdas[2]", "alphas[1]", "alphas[2]", "alpha_sp_intercept[1]",
                        "alpha_sp_intercept[2]", "alpha_sp_intercept[3]", "alpha_sp_intercept[4]")
pairs(PrelimFit, pars = NonShrinkageParams) 
traceplot(PrelimFit, pars = NonShrinkageParams[1:2])
traceplot(PrelimFit, pars = NonShrinkageParams[3:4])
traceplot(PrelimFit, pars = NonShrinkageParams[5:8])

# autocorrelation of the MCMC samples
quartz()
par(mfrow = c(2,2))
for(i in 1:2){
        acf(PrelimPosteriors$lambdas[,i])
        acf(PrelimPosteriors$alphas[,i])
}

for(i in 1:4){
        acf(PrelimPosteriors$alpha_sp_intercept[,i])
}

acf(PrelimPosteriors$c2_tilde)
acf(PrelimPosteriors$tau_tilde)


# Visually examine the model estimates for key parameters
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "salmon", pars = c("alpha_sp_phos"))
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "salmon", pars = "lambdas")
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "salmon", pars = "alphas")
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "salmon", pars = "alpha_sp_intercept")



# Using these values, Carvalho et al. (2010) and Datta and Ghosh (2013) showed that a simple
#    rule of thumb applies where if 1 - kappa_i > 0.5, then we reject the null hypothesis that
#    alpha_sp_i is 0
# Side note: I don't love this approach due to the frequentist philosophy and language driving it,
#    but, I think it should work as a rule of thumb for including "important" interactions in our
#    downstream analyses.
# Plot the paramters to be included by using the IQR of the shrinkage coefficients
S <- 45
QuantVals <- matrix(data = NA, nrow = S, ncol = 2)
for(s in 1:S){
     QuantVals[s,] <- HDInterval::hdi(1-ShrinkCoef[,s], credMass = 0.9)
}
ShrinkageInclusionCols <- ifelse(QuantVals[,1] > 0.5, "green", "red")
quartz()
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     pars = c("alpha_sp_phos"), fill_color = ShrinkageInclusionCols)

# Now take a look at an alternative inclusion criteria based on the posterior distributions directly
PosteriorInclusionCols <- rep(NA, S)
for(s in 1:S){
     IntVals <- HDInterval::hdi(PrelimPosteriors$alpha_sp_phos[,s], credMass = 0.75)
     PosteriorInclusionCols[s] <- ifelse(IntVals[1] > 0 | IntVals[2] < 0, "green", "red")
}
quartz()
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = PosteriorInclusionCols, pars = c("alpha_sp_phos"))

# If we combine these approaches, here are the parameters we would include
InclusionCols <- rep(NA, S)
for(s in 1:S){
     if(ShrinkageInclusionCols[s] == "green" | PosteriorInclusionCols[s] == "green"){
          InclusionCols[s] <- "green"
     }else{
          InclusionCols[s] <- "red"
     }
}
quartz()
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = InclusionCols, pars = c("alpha_sp"))

################### Now run a beverton-holt model with no shrinkage, but only
#       including the previously identified species to obtain a final fit
SpData <- read.csv("water_full_env.csv")
SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
SpData <- na.omit(SpData) 
SpDataFocal <- subset(SpData, Focal.sp.x == "W")

# From here we need to calculate and create the Intra vector, the
#   SpMatrix of abundances for all heterospecifics, and other necessary
#   objects like N, S, and Fecundity
shade <- as.vector(scale(SpDataFocal$Canopy))
phos <- as.vector(scale(SpDataFocal$Colwell.P))
Intra <- as.integer(SpDataFocal$Waitzia.acuminata)
N <- as.integer(dim(SpDataFocal)[1])
Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
Species <- names(SpDataFocal[10:69]) 
Species <- setdiff(Species, "Waitzia.acuminata")
TempSpMatrix <- subset(SpDataFocal, select = Species)
# Now discount any columns with 0 abundance
SpTotals <- colSums(TempSpMatrix)
SpToKeep <- SpTotals > 0
SpMatrix <- matrix(NA, nrow = nrow(TempSpMatrix), ncol = sum(SpToKeep)+1)
SpMatrix[,1] <- Intra
s <- 2
for(i in 1:ncol(TempSpMatrix)){
        if(SpToKeep[i]){
                SpMatrix[,s] <- TempSpMatrix[,i]
                s <- s + 1
        }
}
S <- ncol(SpMatrix)

load("ModelFits/FH_NoEnv_FinalFit.rdata")
InclusionIntercept <- Inclusion
InclusionSlope <- ifelse(PosteriorInclusionCols == "green", 1, 0)
S_intercept <- sum(InclusionIntercept)
S_slope <- sum(InclusionSlope)
Posteriors <- extract(FinalFit)
intercept_mean <- rep(NA, S_intercept)
intercept_sd <- rep(NA, S_intercept)
for(s in 1:S_intercept){
        intercept_mean[s] <- mean(Posteriors$alpha_sp[,s])
        intercept_sd[s] <- 5 * sd(Posteriors$alpha_sp[,s])
}
DataVec <- c("N", "S", "Fecundity", "SpMatrix", "phos", "InclusionIntercept",
             "InclusionSlope", "S_intercept", "S_slope", "intercept_mean",
             "intercept_sd")
FinalFit <- stan(file = "BH_PhosFinal.stan", data = DataVec, iter = 6000,
                 chains = 3, control = list(adapt_delta = 0.9))
# Evaluate the model fit for convergence, mixing, autocorrelation, etc.
FinalFit
FinalPosteriors <- rstan::extract(FinalFit)
# Diagnostic plots
pairs(FinalFit, pars = c("lambdas", "alphas")) 
traceplot(FinalFit, pars = "lambdas")
traceplot(FinalFit, pars = "alphas")
traceplot(FinalFit, pars = "alpha_sp_intercept")
traceplot(FinalFit, pars = "alpha_sp_slope")

# autocorrelation of the MCMC samples
quartz()
par(mfrow = c(2,2))
for(i in 1:2){
        acf(FinalPosteriors$lambdas[,i])
        acf(FinalPosteriors$alphas[,i])
}

for(i in 1:4){
        acf(FinalPosteriors$alpha_sp_intercept[,i])
}

acf(FinalPosteriors$alpha_sp_slope[,1])

# Visually examine the model estimates for key parameters
quartz(width = 10)
par(mfrow = c(1,1))
plot(x = NA, y = NA, xlim = c(-4, 4), ylim = c(0.5, 4.5), main = "", xlab = "species' specific intercept deviations",
     ylab = "Species")
for(i in 1:4){
     PrelimInts <- HDInterval::hdi(PrelimPosteriors$alpha_sp_intercept[,i], credMass = 0.95)
     FinalInts <- HDInterval::hdi(FinalPosteriors$alpha_sp_intercept[,i], credMass = 0.95)
     segments(x0 = PrelimInts[1], y0 = i-0.1, x1 = PrelimInts[2], y1 = i-0.1, col = "darkred")
     segments(x0 = FinalInts[1], y0 = i+0.1, x1 = FinalInts[2], y1 = i+0.1, col = "darkgreen")
}
abline(v = 0, lty = 2)


quartz(width = 10)
par(mfrow = c(2,2))
hist(PrelimPosteriors$lambdas[,1], main = "", xlab = "lambdas[1]", col = "darkred")
hist(FinalPosteriors$lambdas[,1], add = TRUE, col = "darkgreen")

hist(PrelimPosteriors$lambdas[,2], main = "", xlab = "lambdas[2]", col = "darkred")
hist(FinalPosteriors$lambdas[,2], add = TRUE, col = "darkgreen")

hist(PrelimPosteriors$alphas[,1], main = "", xlab = "alphas[1]", col = "darkred")
hist(FinalPosteriors$alphas[,1], add = TRUE, col = "darkgreen")

hist(PrelimPosteriors$alphas[,2], main = "", xlab = "alphas[2]", col = "darkred")
hist(FinalPosteriors$alphas[,2], add = TRUE, col = "darkgreen")

quartz()
par(mfrow = c(1,1))
hist(PrelimPosteriors$alpha_sp_phos[,which(InclusionSlope == 1)], main = "", xlab = "alpha_sp_slope",
     col = "darkred")
hist(FinalPosteriors$alpha_sp_slope[,1], add = TRUE, col = "darkgreen")

HDInterval::hdi(PrelimPosteriors$alpha_sp_phos[,18])
HDInterval::hdi(FinalPosteriors$alpha_sp_slope[,1])

## Save the final fit object along with the inclusion vectors and a vector with
#    species' names to map between the inclusion vectors and the data
SpeciesNames <- c("Intra", Species[SpToKeep])
save(FinalFit, InclusionIntercept, InclusionSlope, SpeciesNames, file = "ModelFits/Final_Waitzia_acuminata_phos.rdata")

#################################### Finally, make a figure of how alpha varies 
#                                       with phosphorous according to this model
PostLength <- nrow(FinalPosteriors$lambdas)
PhosLength <- 1000
PhosSeq <- seq(min(phos), max(phos), length.out = PhosLength)
SpSpecificIntercept <- which(InclusionIntercept == 1)
SpSpecificSlope <- which(InclusionSlope == 1)
TotalSpSpecific <- length(union(SpSpecificIntercept, SpSpecificSlope))
# Calculate the full posterior for the ln(alpha) values
GenericAlpha <- matrix(NA, nrow = PhosLength, ncol = PostLength)
SpSpecificAlphas <- array(NA, dim = c(TotalSpSpecific, PhosLength, PostLength))
for(i in 1:PhosLength){
     GenericAlpha[i,] <- FinalPosteriors$alphas[,1] + FinalPosteriors$alphas[,2]*PhosSeq[i]
     for(j in 1:4){
          if(j != 2){
               SpSpecificAlphas[j,i,] <- GenericAlpha[i,] + FinalPosteriors$alpha_sp_intercept[,j]
          } else{
               SpSpecificAlphas[j,i,] <- GenericAlpha[i,] + FinalPosteriors$alpha_sp_intercept[,j] +
                    FinalPosteriors$alpha_sp_slope[,1] * PhosSeq[i]
          }
     }
}

# Now calculate the mean and 95% CI for the posteriors
GenericInts <- matrix(NA, nrow = PhosLength, ncol = 3)
SpSpecificInts <- array(NA, dim = c(4, PhosLength, 3))
for(i in 1:PhosLength){
     GenericInts[i,1] <- mean(GenericAlpha[i,])
     GenericInts[i,2:3] <- HDInterval::hdi(GenericAlpha[i,])
     for(j in 1:4){
          SpSpecificInts[j,i,1] <- mean(SpSpecificAlphas[j,i,])
          SpSpecificInts[j,i,2:3] <- HDInterval::hdi(SpSpecificAlphas[j,i,])
     }
}

library(RColorBrewer)
ColSeq <- brewer.pal(n = 4, "Dark2")
xRange <- range(PhosSeq)
yRange <- range(c(range(GenericInts), range(SpSpecificInts)))
pdf(file = "PhosphorousAlphas.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), oma = c(4,4,1,1), mar = c(0.5, 0.5, 0.5, 0.5))
     for(j in 1:4){
          plot(NA, NA, xlim = xRange, ylim = yRange, main = "", ylab = "", xlab = "",
               xaxt = "n", yaxt = "n")
          xCoords <- c(PhosSeq, PhosSeq[PhosLength:1])
          yCoords <- c(GenericInts[,2], GenericInts[PhosLength:1,3])
          polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
          lines(x = PhosSeq, y = GenericInts[,1], lwd = 2)
          lines(x = PhosSeq, y = SpSpecificInts[j,,1], col = ColSeq[j], lwd = 2)
          lines(x = PhosSeq, y = SpSpecificInts[j,,2], col = ColSeq[j], lwd = 0.75, lty = 2)
          lines(x = PhosSeq, y = SpSpecificInts[j,,3], col = ColSeq[j], lwd = 0.75, lty = 2)
          if(j == 1 | j == 3){
               axis(2)
          } else{
               axis(2, labels = FALSE)
          }
          if(j == 3 | j == 4){
               axis(1)
          } else{
               axis(1, labels = FALSE)
          }
     }
     mtext("Standardized phosphorous", side = 1, outer = TRUE, line = 2)
     mtext(expression(paste("ln(", alpha, ")", sep = "")), side = 2, outer = TRUE, line = 2)
dev.off()


