rm(list=ls())
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/unify_packages.R")
require(mvtnorm)
require(TMB)
require(mpcmp)
# Dado gerado pela binomial negativa
data <- readRDS("data_cor0.rds")
set.seed(2321)
vecs <- sample(1:10000, 4000)
with(data, cor(Y1, Y2))
data$Y1 <- data$Y1[vecs]
data$Y2 <- data$Y2[vecs]
data$X <- data$X[vecs, ]
n <- length(data$Y1)
true_rho <- 0
nu0 <- .7
nu1 <- .7
parameters <- list(beta1 = c(2, 0.8), beta2 = c(0.5, -1),
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = true_rho, sigma = c(0.1,0.1),
                   nu = c(log(nu0), log(nu1)))
True <- c(beta1 = c(2, 0.8), beta2 = c(0.5, -1),
          sigma = c(0.5,0.3),
          rho = true_rho, 
          nu = c(NA, NA))
# Estimação pelo COM-Poisson
model <- "04_compoisson_bi_full"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

obj <- MakeADFun(data, parameters, DLL=model, random = "U")
obj$fn()     #1073.514  valor original (1098,039 agora, último)
# 17h40min
system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr))
(46152/60)/60
#13 hours
cbind(True = True,
      Com_Poisson = fit_NLB$par)
## Obtaining parameters in a logical way
ss <- sdreport(obj)
rep_fixed <- summary(ss, "fixed")
rep_report <- summary(ss, "report")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/unify_packages.R")
saveRDS(cbind(True = True,
              TMB = TMB_summary(rep_fixed, rep_report)),
        "NB_biv_4000.rds")
