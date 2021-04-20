
library(tidyverse)
# monotonic lambda data
setwd("~/Documents/Work/Current Papers/SparseInteractions/BH_simulations/")
load("StanFits/monoLambda_constAlpha/N50_FinalFit.rdata")
PostLength <- length(FinalPosteriors$alpha_generic)
Focal <- 1

# individual draws from the posterior
n = 50
env.seq <- seq(-2, 2, length.out = n)
lambda.ei <- tibble(env.seq = rep(env.seq, times = PostLength), 
                    post = rep(1:PostLength, each = n), 
                    lambda.int = rep(FinalPosteriors$lambdas[,1], each = n),
                    lambda.slope = rep(FinalPosteriors$lambdas[,2], each = n))
lambda.ei$lambda.mono <- lambda.ei %>% with(exp(lambda.int + lambda.slope*env.seq))
lambda.ei$post <- as.factor(lambda.ei$post)

# sample 1000 from this to make graphing faster
sample.use <- sample(1:PostLength, size = 1000)
lambda.ei.small <- filter(lambda.ei, post %in% sample.use)
lambda.ei.small$ind <- 'ind'

# average from the posterior (could leave this out)
post.mean <- mean(FinalPosteriors$lambdas[,1])
post.env <- mean(FinalPosteriors$lambdas[,2])
lambda.post <- tibble(env.seq, post = as.factor(0), ind = 'post.mean', 
                      lambda.int = rep(post.mean, each = n),
                      lambda.slope = rep(post.env, each = n))
lambda.post$lambda.mono <- lambda.post %>% with(exp(lambda.int + lambda.slope*env.seq))


# true value
TrueVals <- read.csv("Simulations/parameters_perturb2.csv")
lambda.mean <- with(TrueVals, lambda.mean[species == Focal])
lambda.env <- with(TrueVals, lambda.env[species == Focal])
lambda.true <- tibble(env.seq, post = as.factor(0), ind = 'true', 
                      lambda.int = rep(lambda.mean, each = n),
                      lambda.slope = rep(lambda.env, each = n))
lambda.true$lambda.mono <- lambda.true %>% with(exp(lambda.int + lambda.slope*env.seq))

# lambda.ei.all <- left_join(lambda.ei.small, lambda.true)

ggplot(lambda.ei.small, aes(x = env.seq, y = lambda.mono)) + 
     geom_line(aes(group = post), alpha = 0.02, color = 'grey50') + 
     geom_line(data = lambda.true, color = 'red') +
     geom_line(data = lambda.post, color = 'black', linetype = 'dashed') + 
     theme_classic() + 
     ylab(expression(lambda[ei])) +
     xlab('Environment')
#ggsave(filename = 'Results/Box/monoLambda_constAlpha_N50.pdf', width = 5, height = 5, units = 'in')

## Same process with optimum lambda
load("StanFits/optLambda_constAlpha/N200_FinalFit.rdata")
PostLength <- length(FinalPosteriors$alpha_generic)
Focal <- 8

# individual draws from the posterior
n = 50
env.seq <- seq(-2, 2, length.out = n)
lambda.ei <- tibble(env.seq = rep(env.seq, times = PostLength), 
                    post = rep(1:PostLength, each = n), 
                    lambda.max = rep(FinalPosteriors$lambda_max, each = n),
                    z.env = rep(FinalPosteriors$lambda_opt, each = n),
                    sigma.env = rep(FinalPosteriors$lambda_width, each = n))
lambda.ei$lambda.opt <- lambda.ei %>% 
     with(lambda.max * exp(-((z.env - env.seq)/(2*sigma.env))^2))
lambda.ei$post <- as.factor(lambda.ei$post)

# sample 1000 from this to make graphing faster
sample.use <- sample(1:PostLength, size = 1000)
lambda.ei.small <- filter(lambda.ei, post %in% sample.use)
lambda.ei.small$ind <- 'ind'

# average from the posterior (could leave this out)
post.max <- mean(FinalPosteriors$lambda_max)
post.z.env <- mean(FinalPosteriors$lambda_opt)
post.sigma.env <- mean(FinalPosteriors$lambda_width)
lambda.post <- tibble(env.seq, post = as.factor(0), ind = 'post.mean', 
                      lambda.max = rep(post.max, each = n),
                      z.env = rep(post.z.env, each = n),
                      sigma.env = rep(post.sigma.env, each = n))
lambda.post$lambda.opt <- lambda.post %>% 
     with(lambda.max * exp(-((z.env - env.seq)/(2*sigma.env))^2))


# true value
TrueVals <- read.csv("Simulations/parameters_perturb_opt.csv")
lambda.max <- with(TrueVals, lambda.max[species == Focal])
z.env <- with(TrueVals, z.env[species == Focal])
sigma.env <- with(TrueVals, sigma.env[species == Focal])
lambda.true <- tibble(env.seq, post = as.factor(0), ind = 'true', 
                      lambda.max = rep(lambda.max , each = n),
                      z.env = rep(z.env , each = n),
                      sigma.env = rep(sigma.env, each = n))
lambda.true$lambda.opt <- lambda.true %>% 
     with(lambda.max * exp(-((z.env - env.seq)/(2*sigma.env))^2))


ggplot(lambda.ei.small, aes(x = env.seq, y = lambda.opt)) + 
     geom_line(aes(group = post), alpha = 0.02, color = 'grey50') + 
     geom_line(data = lambda.true, color = 'red') +
     geom_line(data = lambda.post, color = 'black', linetype = 'dashed') + 
     theme_classic() + 
     ylab(expression(lambda[ei])) +
     xlab('Environment')
#ggsave(filename = 'Results/Box/optLambda_constAlpha_N200.pdf', width = 5, height = 5, units = 'in')
