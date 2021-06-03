library(tidyverse)
library(patchwork)
library(HDInterval)
setwd("~/Documents/Work/Current Papers/SparseInteractions/BH_simulations/")


# plot theme
theme_cw <- function () { 
     theme_bw(base_size=12) %+replace% 
          theme(
               panel.background = element_blank(), 
               plot.background = element_blank(), 
               axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
               panel.grid.minor = element_blank(), 
               panel.grid.major = element_blank(),
               legend.background = element_blank(), 
               legend.key = element_blank(),
               strip.background = element_blank(), 
               # axis.text=element_text(size=12),
               #     strip.text=element_text(size=12),
               complete = TRUE
          )
}

# Create lists for the results from different sizes of datasets
sample.size <- c(10, 20, 50, 80, 100, 200) # need the final fits still for some of these
file.prefixes <- paste('N', sample.size, '_', sep = "")
n.samples <- length(file.prefixes)
alpha.type <- c('envAlpha', 'constAlpha')
Focal <- 1
S <- 15

alpha.fit <- tibble(sample.size = rep(sample.size, times = length(alpha.type)),
                    alpha.type = rep(alpha.type, each = n.samples),
                    n.ij = -1, n.eij = -1, n.total = -1,
                    dev.mean = -1, dev.sd = -1, dev.se = -1)

## env alpha fits
FullSim <- read.csv("Simulations/simulation_perturb2.csv")
ppc_data <- subset(FullSim, (species == Focal) & (run > 200) & (time == 0) & (thinned == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
     SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) & (time == 0) & (thinned == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

for(sn in 1:n.samples){
     ## number of samples included
     file.prelim <- paste("StanFits/monoLambda_envAlpha/",
                          file.prefixes[sn],"PrelimFit_b.rdata", sep = "")
     load(file.prelim)
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
          if(s == Focal){
               Inclusion_ij[s] <- 0
               Inclusion_eij[s] <- 0
          }
     }
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'envAlpha', 'n.ij'] <- sum(Inclusion_ij)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'envAlpha', 'n.eij'] <- sum(Inclusion_eij)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'envAlpha', 'n.total'] <- sum(Inclusion_ij + Inclusion_eij)
     
     ## ppc deviations
     file.final <- paste("StanFits/monoLambda_envAlpha/",
                      file.prefixes[sn],"FinalFit_b.rdata", sep = "")
     load(file.final)
     
     PostLength <- length(Posteriors$alpha_generic[,1])
     # calculate the posterior distributions of the interaction coefficients and lambdas
     alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
     lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
     for(i in 1:N_ppc){
          lambda_ei[,i] <- exp(Posteriors$lambdas[,1] + Posteriors$lambdas[,2]*env_ppc[i])
          for(s in 1:S){
               if(s == Focal){
                    alpha_eij[,i,s] <- exp(Posteriors$alpha_intra[,1] + Posteriors$alpha_intra[,2] * env_ppc[i])
               }else{
                    alpha_eij[,i,s] <- exp(Posteriors$alpha_generic[1] + Inclusion_ij[s] * Posteriors$alpha_hat_ij[,s] +
                                                (Posteriors$alpha_generic[,2] + Inclusion_eij[s] * Posteriors$alpha_hat_eij[,s]) * env_ppc[i])
               }
          }
     }
     
     # use the above quantities to calculate the posterior prediction intervals for the new data
     Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
     GrowthRMSE <- numeric(length = PostLength)
     for(i in 1:PostLength){
          for(j in 1:N_ppc){
               SigmaTerm <- sum(alpha_eij[i,j,] * SpMatrix_ppc[j,])
               Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
               Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          }
       deviation_sq <- (Growth_pred[i,] - Growth_ppc)^2
       GrowthRMSE[i] <- sqrt(sum(deviation_sq) / N_ppc)
     }
     
     # Calculate final fit results for the ppc
     PredVals <- colMeans(Growth_pred)
     pred.dev <- Growth_ppc - PredVals
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'envAlpha', 'dev.mean'] <- mean(abs(pred.dev))
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'envAlpha', 'dev.sd'] <- sd(abs(pred.dev))
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'envAlpha', 'dev.se'] <- sd(abs(pred.dev))/sqrt(length(pred.dev))
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'envAlpha', 'rmse'] <- mean(GrowthRMSE)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'envAlpha', 'rmse.low'] <- HDInterval::hdi(GrowthRMSE)[1]   
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'envAlpha', 'rmse.high'] <- HDInterval::hdi(GrowthRMSE)[2] 
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'envAlpha', 'tau'] <- mean(PrelimPosteriors$tau)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'envAlpha', 'tau.low'] <- HDInterval::hdi(PrelimPosteriors$tau)[1]   
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'envAlpha', 'tau.high'] <- HDInterval::hdi(PrelimPosteriors$tau)[2] 
}


## const alpha fits

FullSim <- read.csv("Simulations/simulation_perturb2_const.csv")
ppc_data <- subset(FullSim, (species == Focal) & (run > 200) & (time == 0) & (thinned == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
     SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) & (time == 0) & (thinned == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)
## number of samples included
for(sn in 1:n.samples){
     file.prelim <- paste("StanFits/monoLambda_constAlpha/",
                          file.prefixes[sn],"PrelimFit_b.rdata", sep = "")
     load(file.prelim)
     Inclusion_ij <- rep(0, S)
     IntLevel <- 0.5
     for(s in 1:S){
          Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
          if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
               Inclusion_ij[s] <- 1
          }
          if(s == Focal){
               Inclusion_ij[s] <- 0
          }
     }
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'constAlpha', 'n.ij'] <- sum(Inclusion_ij)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'constAlpha', 'n.total'] <- sum(Inclusion_ij)
     
     ## ppc deviations
     file.final <- paste("StanFits/monoLambda_constAlpha/",
                         file.prefixes[sn],"FinalFit_b.rdata", sep = "")
     load(file.final)
     
     PostLength <- length(FinalPosteriors$alpha_generic)
     # calculate the posterior distributions of the interaction coefficients and lambdas
     alpha_ij <- array(NA, dim = c(PostLength, N_ppc, S))
     lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
     for(i in 1:N_ppc){
          lambda_ei[,i] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[i])
          for(s in 1:S){
               if(s == Focal){
                    alpha_ij[,i,s] <- exp(FinalPosteriors$alpha_intra)
               }else{
                    alpha_ij[,i,s] <- exp(FinalPosteriors$alpha_generic + Inclusion_ij[s] * FinalPosteriors$alpha_hat_ij[,s])
               }
          }
     }
     
     # use the above quantities to calculate the posterior prediction intervals for the new data
     Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
     GrowthRMSE <- numeric(length = PostLength)
     
     for(i in 1:PostLength){
          for(j in 1:N_ppc){
               SigmaTerm <- sum(alpha_ij[i,j,] * SpMatrix_ppc[j,])
               Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
               Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
          }
       deviation_sq <- (Growth_pred[i,] - Growth_ppc)^2
       GrowthRMSE[i] <- sqrt(sum(deviation_sq) / N_ppc)
     }
     
     # Calculate final fit results for the ppc
     PredVals <- colMeans(Growth_pred)
     pred.dev <- Growth_ppc - PredVals
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'constAlpha', 'dev.mean'] <- mean(abs(pred.dev))
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'constAlpha', 'dev.sd'] <- sd(abs(pred.dev))
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                    alpha.fit$alpha.type == 'constAlpha', 'dev.se'] <- sd(abs(pred.dev))/sqrt(length(pred.dev))
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'constAlpha', 'rmse'] <- mean(GrowthRMSE)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'constAlpha', 'rmse.low'] <- HDInterval::hdi(GrowthRMSE)[1]   
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'constAlpha', 'rmse.high'] <- HDInterval::hdi(GrowthRMSE)[2]   
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'constAlpha', 'tau'] <- mean(PrelimPosteriors$tau)
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'constAlpha', 'tau.low'] <- HDInterval::hdi(PrelimPosteriors$tau)[1]   
     alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                 alpha.fit$alpha.type == 'constAlpha', 'tau.high'] <- HDInterval::hdi(PrelimPosteriors$tau)[2] 
}

## figure
library(inauguration)
InterceptCol <- inauguration("inauguration_2021")[3]
SlopeCol <- inauguration("inauguration_2021")[4]
#ppcCol <- Dark2Cols[3]

alpha.fit[alpha.fit == -1] <- NA
alpha.fit$sample.size.2 <- factor(alpha.fit$sample.size)
alpha.fit.long <- alpha.fit %>% 
  pivot_longer(cols = starts_with('n'), names_to = 'specific', values_to = 'number') %>%
  filter(!(alpha.type == 'constAlpha' & specific == 'n.total') & !(specific == 'n.eij'))

a.terms <- ggplot(alpha.fit.long, 
       aes(x = sample.size, y = number, 
           color = interaction(alpha.type, specific),
           linetype = interaction(alpha.type, specific))) + 
      geom_line() +
    #geom_point(position = position_dodge(5), aes(color = alpha.type, shape = specific)) + 
     theme_cw() + 
  theme(legend.position = c(0.75, 0.2)) + 
     scale_color_manual(values = c(InterceptCol, SlopeCol, SlopeCol), name = "",
                        labels = c(expression(Simple~context:~alpha[ij]~pairs), 
                                   expression(Complex~context:~alpha[ij]~pairs),
                                   expression(Complex~context:~all~terms)),
                        guide = guide_legend(label.hjust = 0)) +
     scale_linetype_manual(values = c('solid','solid', 'dashed'), name = "",
                           labels = c(expression(Simple~context:~alpha[ij]~pairs), 
                                      expression(Complex~context:~alpha[ij]~pairs),
                                      expression(Complex~context:~all~terms)),
                           guide = guide_legend(label.hjust = 0)) +
     ylab('Number of non-generic terms') +
     xlab('Input data sample size')
#ggsave(filename = 'Results/Box/alpha_terms_updated.pdf', width = 8, height = 5, units = 'in')

a.post <- ggplot(alpha.fit, aes(x = sample.size, y = rmse, 
                      ymin = rmse.low, ymax = rmse.high,
                      color = alpha.type)) + 
     geom_point(position = position_dodge(5)) + 
     geom_errorbar(position = position_dodge(5), width = 0) +
     theme_cw() + 
     scale_color_manual(values = c(InterceptCol, SlopeCol), name = '',
                        labels = c('Simple context','Complex context'),
                        guide = guide_legend(label.hjust = 0)) +
     ylab('Root mean squared error') +
     xlab('Input data sample size') +
     theme(legend.position = c(0.8, 0.9))
#ggsave(filename = 'Results/Box/alpha_post.pdf', width = 6, height = 5, units = 'in')

a.post/a.terms
# ggsave(filename = 'Results/Box/alpha_combined_4.pdf', width = 5, height = 8, units = 'in')

patchwork <- (p.mono + p.opt) / (a.post + a.terms)

patchwork + 
  plot_annotation(tag_levels = 'a')  & 
  theme(plot.tag = element_text(face = "bold"))
ggsave(filename = 'Results/Box/box_figure.pdf', width = 10, height = 8, units = 'in')



a.tau <- ggplot(alpha.fit, aes(x = sample.size, y = tau, 
                                ymin = tau.low, ymax = tau.high,
                                color = alpha.type)) + 
  geom_point(position = position_dodge(5)) + 
  geom_errorbar(position = position_dodge(5), width = 0) +
  theme_cw() + 
  scale_color_manual(values = c(InterceptCol, SlopeCol), name = '',
                     labels = c('Simple context','Complex context'),
                     guide = guide_legend(label.hjust = 0)) +
  ylab('Tau') +
  xlab('Input data sample size') +
  theme(legend.position = c(0.8, 0.8))
#ggsave(filename = 'Results/Box/alpha_tau_mean.pdf', width = 6, height = 5, units = 'in')
