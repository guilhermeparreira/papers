## Negative Binomial distribution ----------------------------------------------
## Authors: Multivariate Statistical Modelling Group ---------------------------
## Date: June, 04 2020 ---------------------------------------------------------
rm(list=ls())
require(MASS)
require(mvtnorm)
require(TMB)
## Data simulation ------------------------------------------------------------
beta1 <- c(2, 0.8)
beta2 <- c(0.5, -1)
n <- 200
# n <- 10000
set.seed(1234)
x1 <- rnorm(n)
X <- model.matrix(~x1)

## Random effects covariance matrix
s1 <- log(0.5) #.5 é s2
s2 <- log(0.3) #.3 é s2
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
phi0 <- log(1)
phi1 <- log(1)
set.seed(1234)
Y1 <- rnbinom(n, size = exp(phi0), mu = mu1)
set.seed(1234)
Y2 <- rnbinom(n, size = exp(phi1), mu = mu2)

## Data set
data <- list(Y1 = Y1, Y2 = Y2, X = X, phi = c(phi0, phi1))
## Parameters
True <- c(2,0.8,0.5,-1,rho,log(0.5),log(0.3), phi0,phi1)
parameters <- list(beta1 = c(2, 0.8), beta2 = c(0.5, -1),
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = 0, sigma = c(0.1,0.1))

## Compilation
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "06_negbin_bivariate_phi_fixed"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

obj <- MakeADFun(data, parameters, DLL=model, random = "U")
obj$fn()     #1073.514  valor original
# obj$gr()


system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)) #44seg
system.time(fit_NM <- optim(par = obj$par, fn = obj$fn, gr = obj$gr)) #2.8min
system.time(fit_BFGS <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")) #35seg
system.time(fit_CG <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "CG")) # 2.083seg

dput(round(fit_NLB$par,3))
# Otimizado
pars_phi_fixed <- c(beta1 = 2.057, beta1 = 0.814, beta2 = 0.478, beta2 = -1.174, 
                    rho = 0.503, sigma = 0.705, sigma = 0.707)
nll_phi_fixed <- 1064.1


## Maximized log-likelihood values
cbind("Nelder-Mead" = fit_NM$value, "BFGS" = fit_BFGS$value,
      "CG" = fit_CG$value, "NLB" = fit_NLB$objective)

## Estimates
cbind("True" =  True, "Nelder-Mead" = fit_NM$par, "BFGS" = fit_BFGS$par,
      "CG" = fit_CG$par, "NLB" = fit_NLB$par)

inv.transf.fisher <- function(z){ # Inverso da Transf Fisher Z
  r <- (exp(2*z)-1)/(exp(2*z)+1)
  return(r)
}

## Correlation value
ss <- sdreport(obj)
summary(ss, "report")
