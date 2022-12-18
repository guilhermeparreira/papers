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
n <- 1500
# n <- 10000
set.seed(1234)
x1 <- rnorm(n)
X <- model.matrix(~x1)

## Random effects covariance matrix
s1 <- log(1) #1 é s2
s2 <- log(1) #1 é s2
# Transformação Z de Fisher
rho <- 0.5*log((1 + .5)/ (1-.5)) # .5 é o verdadeira valor de rho
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
l_phi0 <- log(1)
l_phi1 <- log(1)
set.seed(1234)
Y1 <- rnbinom(n, size = exp(l_phi0), mu = mu1)
set.seed(1234)
Y2 <- rnbinom(n, size = exp(l_phi1), mu = mu2)

## Data set
data <- list(Y1 = Y1, Y2 = Y2, X = X,
             sigma = c(sqrt(exp(s1)),sqrt(exp(s2))))
## Parameters
# True <- c(2,0.8,0.5,-1,rho,log(0.5),log(0.3),log(phi0),log(phi1))
parameters <- list(beta1 = beta1, beta2 = beta2,
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = 0, phi = c(l_phi0, l_phi1))

## Compilation
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "06_negbin_bivariate_var_fixed"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

obj <- MakeADFun(data, parameters, DLL=model, random = "U")
obj$fn()     #1073.514  valor original
# obj$gr()

#NLB se perdeu com rho = .75; com rho = .5 estimou bem;
#NM se perdeu com rho = .75;


system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)) #44seg
# system.time(fit_NM <- optim(par = obj$par, fn = obj$fn, gr = obj$gr)) #2.8min
# system.time(fit_BFGS <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")) #35seg
# system.time(fit_CG <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "CG")) # 2.083seg
# 
# # Otimizado
# pars_var_fixed <- c(beta1 = 2.23, beta1 = 0.754, beta2 = 0.762, beta2 = -1.06, 
#                     rho = 0.87, phi = -0.295, phi = -0.484)
# nll_var_fixed <- 1148.459
# 
# ## Maximized log-likelihood values
# cbind("Nelder-Mead" = fit_NM$value, "BFGS" = fit_BFGS$value,
#       "CG" = fit_CG$value, "NLB" = fit_NLB$objective)
# 
# ## Estimates
# cbind("True" =  True, "Nelder-Mead" = fit_NM$par, "BFGS" = fit_BFGS$par,
#       "CG" = fit_CG$par, "NLB" = fit_NLB$par)
# 
# inv.transf.fisher <- function(z){ # Inverso da Transf Fisher Z
#   r <- (exp(2*z)-1)/(exp(2*z)+1)
#   return(r)
# }
TMB_summary <- function(fixef, cor, nresp=2){
  # TMB, UNSTRUCTURRED AND VSCALE PARAMETRIZATION
  # Beta (Fixed Effects)
  beta <- fixef[grep("beta", rownames(fixef), invert = F),]
  ncovs <- nrow(beta)/nresp
  rownames(beta) <- rep(c(rep("beta0"), rep("beta1...",ncovs-1)), nresp)
  beta <- beta[order(rownames(beta)), ]
  # Sigma (Standard error of random effect)
  sigma <- fixef[grep("sigma", rownames(fixef), invert = F),]
  # Phi Dispersion parameter (2nd parameter of a distribution)
  sigma <- fixef[grep("phi", rownames(fixef), invert = F),]
  # Correlation of random effects
  cor <- cor[as.vector(lower.tri(diag(nresp))), ]
  out <- rbind(beta, sigma, cor)
  return(out)
}

# Parameters
rep <- sdreport(obj)
rep_fe <- summary(rep, "fixed")
rep_cor <- summary(rep, "report")
TMB_summary(rep_fe, rep_cor)
# Correlation parameters
# Profile na mão
par(mfrow=c(2,2), mai = rep(.3,4))
# Rho
with(p5[1:70,], plot(value~rho, type = "l"))
p5 <- tmbprofile(obj,5)
plot(p5, main = "rho", parm.range = c(-1,10))
# Beta1
p1 <- tmbprofile(obj,1)
plot(p1, main = "beta1")
# Phi
p6 <- tmbprofile(obj,6)
plot(p6, main = "phi1")
# PHi
p7 <- tmbprofile(obj,7)
plot(p7, main = "phi2")



# Profile
l <- length(fit_NLB$par)
par(mfrow=c(4,2), mai = rep(.2,4))
for (i in 1:l){
  paux <- tmbprofile(obj,i)
  print(i)
  print(plot(paux))
}
p5
plot(p1)
assign()
# p1 <- tmbprofile(obj, i)
# paste0("p",i) <- i
# assign(paste0("p",i), i)

# Em paralelo
plan(multiprocess)
future_map(1:l,
        function(i) 
        paux <- tmbprofile(obj,i)
        print(i)
        print(plot(paux)))
# Profile
l <- length(fit_NLB$par)
par(mfrow=c(4,2), mai = rep(.1,4))
for (i in 1:l){
  paux <- tmbprofile(obj,i)
  print(i)
  print(plot(paux))
}
