# Use this script to evaluate the current fit of the model from the server



# Evaluate the model fit for convergence, mixing, autocorrelation, etc.
PrelimFit
PrelimPosteriors <- rstan::extract(PrelimFit)
# Diagnostic plots
pairs(PrelimFit, pars = c("log_lambda", "alpha_hat", "tau_tilde", "c2_tilde")) 
traceplot(PrelimFit, pars = c("log_lambda", "alpha_hat", "tau_tilde", "c2_tilde"))
traceplot(PrelimFit, pars = "alpha_sp_tilde")
traceplot(PrelimFit, pars = "local_shrinkage")

# autocorrelation of the MCMC samples
acf(PrelimPosteriors$log_lambda)
acf(PrelimPosteriors$tau_tilde)
acf(PrelimPosteriors$alpha_hat)
acf(PrelimPosteriors$c2_tilde)
for(s in 1:S){
     acf(posteriors$alpha_sp_tilde[,s])
     acf(posteriors$local_shrinkage[,s])
}

# Visually examine the model estimates for key parameters
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "salmon", pars = c("alpha_sp"))
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "salmon", pars = c("log_lambda", "alpha_hat", "tau"))



# Using these values, Carvalho et al. (2010) and Datta and Ghosh (2013) showed that a simple
#    rule of thumb applies where if 1 - kappa_i > 0.5, then we reject the null hypothesis that
#    alpha_sp_i is 0
# Side note: I don't love this approach due to the frequentist philosophy and language driving it,
#    but, I think it should work as a rule of thumb for including "important" interactions in our
#    downstream analyses.
# Plot the paramters to be included by using the IQR of the shrinkage coefficients
QuantVals <- matrix(data = NA, nrow = 2, ncol = S)
for(s in 1:S){
     QuantVals[,s] <- HDInterval::hdi(1-ShrinkageCoef[,s], credMass = 0.9)
}
ShrinkageInclusionCols <- ifelse(QuantVals[1,] > 0.5, "green", "red")
quartz()
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = ShrinkageInclusionCols, pars = c("alpha_sp"))

# Now take a look at an alternative inclusion criteria based on the posterior distributions directly
PosteriorInclusionCols <- rep(NA, S)
for(s in 1:S){
     IntVals <- HDInterval::hdi(PrelimPosteriors$alpha_sp[,s], credMass = 0.8)
     PosteriorInclusionCols[s] <- ifelse(IntVals[1] > 0 | IntVals[2] < 0, "green", "red")
}
quartz()
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = PosteriorInclusionCols, pars = c("alpha_sp"))

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

# Because I'm curious, let's plot these estimates against the observed abundance
#    of each species
Abunds <- colSums(SpMatrix)
MeanEsts <- rep(NA, S)
IntVals <- matrix(NA, nrow = 2, ncol = S)
for(s in 1:S){
     MeanEsts[s] <- mean(PrelimPosteriors$alpha_sp[,s])
     IntVals[,s] <- HDInterval::hdi(PrelimPosteriors$alpha_sp[,s], credMass = 0.95)
}
quartz()
plot(x = Abunds, y = MeanEsts, pch = 20, col = InclusionCols, ylim = c(-4, 4),
     xlab = "Abundance", ylab = "Deviation from generic competitor")
abline(v = 15, lty = 2)


################### Try running a Beverton-Holt model with no shrinkage on only
#    the previously identified parameters. See how it effects the model fit and 
#    if it seems like an appropriate next step...
Inclusion <- ifelse(InclusionCols == "green", 1, 0)
S_total <- S
S_eff <- sum(Inclusion)
DataVec <- c("N", "S_total", "S_eff", "Inclusion", "Fecundity", "SpMatrix")
FinalFit <- stan(file = "BevertonHolt_NoEnv.stan", data = DataVec, iter = 6000,
                 chains = 3, control = list(adapt_delta = 0.9))
# Evaluate the model fit for convergence, mixing, autocorrelation, etc.
FinalFit
FinalPosteriors <- rstan::extract(FinalFit)
# Diagnostic plots
pairs(FinalFit, pars = c("log_lambda", "alpha_hat")) 
traceplot(FinalFit, pars = c("log_lambda", "alpha_hat"))
traceplot(FinalFit, pars = "alpha_sp")

# autocorrelation of the MCMC samples
acf(FinalPosteriors$log_lambda)
acf(FinalPosteriors$alpha_hat)
for(s in 1:S_eff){
     acf(FinalPosteriors$alpha_sp[,s])
}

# Visually examine the model estimates for key parameters
quartz(width = 10)
par(mfrow = c(1,1))
ToPlot <- which(Inclusion == 1)
plot(x = NA, y = NA, xlim = c(-4, 4), ylim = c(0.5, 4.5), main = "", xlab = "species' specific deviations",
     ylab = "Species")
for(i in 1:4){
     PrelimInts <- HDInterval::hdi(PrelimPosteriors$alpha_sp[,ToPlot[i]], credMass = 0.95)
     FinalInts <- HDInterval::hdi(FinalPosteriors$alpha_sp[,i], credMass = 0.95)
     segments(x0 = PrelimInts[1], y0 = i-0.1, x1 = PrelimInts[2], y1 = i-0.1, col = "darkred")
     segments(x0 = FinalInts[1], y0 = i+0.1, x1 = FinalInts[2], y1 = i+0.1, col = "darkgreen")
}
abline(v = 0, lty = 2)


quartz(width = 10)
par(mfrow = c(1,2))
hist(PrelimPosteriors$log_lambda, main = "", xlab = "log(lambda)", col = "darkred")
hist(FinalPosteriors$log_lambda, add = TRUE, col = "darkgreen")

hist(PrelimPosteriors$alpha_hat, main = "", xlab = expression(hat(alpha)), col = "darkred")
hist(FinalPosteriors$alpha_hat, add = TRUE, col = "darkgreen")

HDInterval::hdi(PrelimPosteriors$log_lambda)
HDInterval::hdi(FinalPosteriors$log_lambda)
# HDI for log_lambda: 2.27-2.44 (Prelim)
#                     2.30-2.46 (Final)
HDInterval::hdi(PrelimPosteriors$alpha_hat)
HDInterval::hdi(FinalPosteriors$alpha_hat)
# HDI for alpha_hat: -5.26 - -4.00 (Prelim)
#                    -4.78 - -3.92

