library(TMB)
setwd("~/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels")
source("/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/unify_packages.R")
# Simulating Data GLM ------------------------------------------------------------------------
seed <- 789
beta1 <- c(-3.01)
n <- 50000*2
set.seed(seed)
X <- matrix(rep(1,n))
mu1 <- exp(X%*%beta1)
## Response variables
phi1 <- .5
l_phi1 <- log(phi1)
True <- c(beta1, phi1)
set.seed(seed)
Y1 <- rnbinom(n, size = exp(l_phi1), mu = mu1)
# Input 
data <- list(Y1 = Y1,
             X = X)
var(Y1);mean(Y1)
parameters <- list(beta1 = beta1,
                   phi = c(l_phi1))
model <- "03_nb_glm" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
# openmp(n=3)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T)
obj$fn()
# opt <- do.call("optim", obj)
(opt <- nlminb(obj$par, obj$fn, obj$gr,
               control = list(eval.max = 1e9, iter.max = 1e9,
                              abs.tol = 1e-4, rel.tol = 1e-4)))
rep <- sdreport(obj)
fixef <- summary(rep, "fixed")
cor <- summary(rep, "report")
fixef;cor
var(Y1)
mean(Y1)
# var
mean(Y1)+(mean(Y1)^2)/cor[,1]
# Joining all together
# (pars_used <- TMB_summary(fixef, cor, 1))
