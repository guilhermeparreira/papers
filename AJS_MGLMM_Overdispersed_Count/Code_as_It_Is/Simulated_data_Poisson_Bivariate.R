# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html

# Packages --------------------------------------------------------------------
rm(list = ls())
library(mcglm)
library(tictoc)
library(MCMCglmm) # asreml idea
library(TMB)
library(brms)
library(mvtnorm)
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/unify_packages.R")
# Simulated Dataset --------------------------------------------------------------------
beta1 <- c(2, 0.8)
beta2 <- c(0.5, -1)
n <- 1500
# n <- 10000
set.seed(1234)
x1 <- rnorm(n)
X <- model.matrix(~x1)

## Random effects covariance matrix
s1 <- log(0.3) #.3 é s2;   .54 é o s
s2 <- log(0.15) #.15 é s2; .38 é o s
# Transformação Z de Fisher
rho <- 0.5*log((1 + 0.5)/ (1-0.5)) # .5 é o verdadeira valor de rho
Sigma <- matrix(NA, 2, 2)
Sigma[1,1] <- exp(s1)
Sigma[2,2] <- exp(s2)
Sigma[1,2] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
Sigma[2,1] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
set.seed(1234)
U <- rmvnorm(n, mean = c(0,0), sigma = Sigma)
mu1 <- exp(X%*%beta1 + U[,1])
mu2 <- exp(X%*%beta2 + U[,2])

## Response variables
set.seed(1234)
Y1 <- rpois(n, lambda = mu1)
set.seed(1234)
Y2 <- rpois(n, lambda = mu2)

# True <- c(2,0.8,0.5,-1,rho, .54, .38)
True <- c(2,0.5,0.8,-1, .54, .38, rho)
ds <- data.frame(Y1 = Y1,
                 Y2 = Y2,
                 X = x1,
                 id = 1:n)

# Descriptive
# First Estimation: Poisson data
# Brms -----------------------------------------------------------------------
# tic() #4.48min
brms_poisson <- brm(mvbind(Y1, Y2) ~ X + (1|p|id),
            family = poisson(link = "log"),
            data = ds, 
            chains = 4, cores = 4,
            iter = 6000)
# toc()
# Plota as posterioris - Diagnóstico
# x11()
# plot(fit1, N = 7)
# plot(fit1, pars = c("shape*"))
# pairs(fit1)

# Summary to comparison
(ab <- summary(brms_poisson))

# brms_summary(ab)

# MCMCglmm -------------------------------------------------------------------
mcmc_poisson <- MCMCglmm::MCMCglmm(cbind(Y1, Y2) ~  trait:(X) + trait
                                 -1, # Surpress an overall intercept
                      # random=~ us(trait):id,
                      family = rep("poisson",2),
                      data= ds,
                      rcov=~ us(trait):units)

# Diagnóstico
# x11()
# plot(MMLM1.fit2)
# coda::autocorr(MMLM1.fit2$VCV)

# Summary
(ab1 <- summary(mcmc_poisson))


# TMB ------------------------------------------------------------------------
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "02_poisson_multivariate" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

data <- list(Y = matrix(cbind(Y1, Y2), ncol = 2), 
             X = X)
# Sigma
# cov2cor(Sigma)
## Parameters
parameters <- list(beta = matrix(c(log(mean(Y1)), 0,
                                   log(mean(Y2)), 0),
                                 ncol = 2),
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = .549,
                   sigma = c(sqrt(.3),
                             sqrt(.15)))


# openmp(n=3)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 # method = "BFGS",
                 hessian = T,
                 silent = F
                 ,random = c("U"))
# Rinterface(paste0(model, ".cpp"))
# obj$fn()
# opt <- do.call("optim", obj)
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(abs.tol = 1e-4))
opt <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

rep <- sdreport(obj)
rep_fe <- summary(rep, "fixed")
rep_cor <- summary(rep, "report")

# Joining all together
# (theta <- opt$par[7:9])
# (sigma <- opt$par[10:12])   #s
# cov2cor(getSigma(theta, sigma)) #Valores das variâncias batem)
TMB_summary(rep_fe, rep_cor, 2)

# fixef <- rep_fe
# cor <- rep_cor



# Joining all results -------------------------------------------------------

cbind(True = round(True,2),
      all_summaries(ab, 
              mcmc_poisson, 
              rep_fe, rep_cor, 2)
      )



# MCMCglmm study ------------------------------------------------------------

# Materiais:
# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf
# http://www.wildanimalmodels.org/tiki-download_wiki_attachment.php?attId=24
# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf
# https://github.com/JonBrommer/Multivariate-Mixed-Models-in-R/wiki/MCMCglmm-examples
# MCMCglmm automatically melts the data for us (and assigns the name trait the same way we did manually above), that's why we use trait
# units which index rows of the response variable  (automatically)
# trait which index columns of the response variable (letting it know we want to fit a multivariate mixed model)
# prioris (not using right now)
prior.model.1 <- list(R = list(V=diag(2)/2, nu=0.04),         # Residual Variance
                      G = list(G1=list(V=diag(.84)/2, nu=0.04)) # Random Effect
                      # B would be the fixed effect
)
# The priors for the variance structures (R and G) are lists with the:
# expected (co)variances (V) and degree of belief parameter (nu) for the inverse-Wishart, 
# and also the mean vector (alpha.mu) and covariance matrix (alpha.V) for the redundant working parameters. 
# The defaults are nu=0, V=1, alpha.mu=0, and alpha.V=0. When alpha.V is non-zero, parameter expanded algorithms are used.

# BRMS study -------------------------------------------------------------------

# Test studies
methods(class = "brmsfit")
stancode(brms_poisson)
str(brms_poisson)
brms_poisson$data
brms_poisson$fit
VarCorr(brms_poisson)
p_cor <- posterior_samples(brms_poisson, "^cor")
sd(p_cor$cor_id__Nadm_Intercept__Ndoc_Intercept)
p_sd <- posterior_samples(brms_poisson, "^sd")

sd(p_sd$sd_id__Nadm_Intercept)
sd(p_sd$sd_id__Ndoc_Intercept)

# TMB study ---------------------------------------------------------------------
# For 3 responses
data <- list(Y = matrix(c(ahs$Nadm,ahs$Nhosp, ahs$Ndoc), ncol = 3),
             X = model.matrix(~actdays, ahs))
n <- nrow(ahs)
parameters <- list(beta = matrix(c(log(mean(ahs$Nadm)), 0,
                                   log(mean(ahs$Nhosp)), 0,
                                   log(mean(ahs$Ndoc)), 0),
                                 ncol = 3),
                   U = matrix(0, ncol = 3, nrow = n),
                   rho = c(with(ahs, cor(Nadm, Nhosp)),
                           with(ahs, cor(Nadm, Ndoc)),
                           with(ahs, cor(Nhosp, Ndoc))), 
                   sigma = c(sd(ahs$Nadm),
                             sd(ahs$Nhosp), 
                             sd(ahs$Ndoc)))
