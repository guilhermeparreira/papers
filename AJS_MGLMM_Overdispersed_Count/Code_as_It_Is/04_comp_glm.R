library(TMB)
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/unify_packages.R")
# Simulating Data GLM ------------------------------------------------------------------------
seed <- 1123
beta1 <- c(0, .15)
n <- 2000
set.seed(seed)
x1 <- rnorm(n)
X <- model.matrix(~x1)
mu1 <- exp(X%*%beta1)
## Response variables
nu1 <- 1
l_nu1 <- log(nu1)
True <- c(beta1, nu1)
set.seed(seed)
Y1 <- mpcmp::rcomp(n, mu = mu1, nu = nu1)
# Input 
data <- list(Y1 = Y1,
             X = X)
var(Y1);mean(Y1)
parameters <- list(beta1 = beta1,
                   nu = c(l_nu1))
model <- "04_comp_glm" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
openmp(n=2)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F)
obj$fn()
# opt <- do.call("optim", obj)
(opt <- nlminb(obj$par, obj$fn, obj$gr,
               control = list(abs.tol = 1e-4, rel.tol = 1e-4)))
rep <- sdreport(obj)
fixef <- summary(rep, "fixed")
cor <- summary(rep, "report")
# Joining all together
(pars_used <- TMB_summary(fixef, cor, 1))
