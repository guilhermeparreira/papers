## Compoisson distribution -----------------------------------------------------
## Authors: Multivariate Statistical Modelling Group ---------------------------
## Date: June, 18 2020 ---------------------------------------------------------
devtools::install_github("thomas-fung/mpcmp")
rm(list=ls())
require(MASS)
require(mvtnorm)
require(TMB)
require(compoisson)
## Data simulation -------------------------------------------------------------
set.seed(1234)
beta1 <- c(.2, 0.8)
beta2 <- c(0.5, -1)
n <- 500
x1 <- rnorm(n)
X <- model.matrix(~x1)
## Random effects covariance matrix
s1 <- log(0.5)
s2 <- log(0.3)
true_rho <- 0
rho <- 0.5*log((1 + true_rho)/ (1-true_rho))
Sigma <- matrix(NA, 2, 2)
Sigma[1,1] <- exp(s1)
Sigma[2,2] <- exp(s2)
Sigma[1,2] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
Sigma[2,1] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
set.seed(1234)
U <- rmvnorm(n, mean = c(0,0), sigma = Sigma)
cor(U)
mu1 <- exp(X%*%beta1 + U[,1])
mu2 <- exp(X%*%beta2 + U[,2])
## Response variables
nu1 <- .7
nu2 <- .7
l_nu1 <- log(nu1)
l_nu2 <- log(nu2)
set.seed(1234)
(Y1 <- mpcmp::rcomp(n, mu = mu1, nu = nu1))
set.seed(1234)
(Y2 <- mpcmp::rcomp(n, mu = mu2, nu = nu2))
# data_export_R <- list(Y1 = Y1,
#                       Y2 = Y2,
#                       X = X)
# saveRDS(data_export_R, file = "com_poisson_data_R.rds")
## Data set
data <- list(Y1 = Y1,
             Y2 = Y2,
             X = X)
## Parameters
parameters <- list(beta1 = beta1,
                   beta2 = beta2,
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = true_rho,
                   sigma = c(sqrt(exp(s1)), sqrt(exp(s2))),
                   nu = c(l_nu1, l_nu2))

## Compilation
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "04_compoisson_bi_full"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
openmp(n=2)
obj <- MakeADFun(data, parameters, DLL=model, random = "U")
system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr))
## Obtaining parameters in a logical way
ss <- sdreport(obj)
rep_fixed <- summary(ss, "fixed")
rep_report <- summary(ss, "report")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/unify_packages.R")
round(TMB_summary(rep_fixed, rep_report, 2), 4)
