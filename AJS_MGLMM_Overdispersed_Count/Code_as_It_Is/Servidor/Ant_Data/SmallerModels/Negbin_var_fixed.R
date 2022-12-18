######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
# load("Negbin_var_fixed.RData")
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
model <- "negbin_multivariate_var_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Obtaining data #######################################################
########################################################################
# Data
library(mvabund)
data("antTraits")
y <- antTraits$abund # Selecting response variables
X <- antTraits$env # Selecting covariates
nr <- ncol(y)
# First attempt: Fix at 1
data <- list(X = as.matrix(cbind(1, X)),
             Y = as.matrix(y),
             sigma = rep(1, nr))
n_params <- 6
# Parameters
U <- matrix(0, nrow = nrow(y), ncol = ncol(y))
parameters <- list(beta = matrix(0, nrow = n_params, ncol = nr),
                   U = U,
                   rho = rep(0, (nr*(nr-1))/2),
                   phi = rep(log(1), nr))
gc()
threads <- 1
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                    control = list(eval.max = 1e7, iter.max = 1e7,
                                                   abs.tol = 1e-04, rel.tol = 1e-04,
                                                   x.tol = 1e-04, xf.tol = 1e-04)))
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
save.image("Negbin_var_fixed.RData")
load("Negbin_var_fixed.RData")
# # Second Round
# save.image("Negbin_var_fixed.RData")
threads <- 12
openmp(n=threads)
parameters <- list(beta = matrix(0, nrow = n_params, ncol = nr),
                   U = U,
                   rho = estimates_negbin$rho,
                   phi = estimates_negbin$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()

# # Beta too high (I will correct!)
# # parameters$beta[which.max(parameters$beta)] <- 1.5
# parameters$beta[parameters$beta>2] <- 1.5
# apply(parameters$beta, 1, max)
# range(parameters$rho)
# range(parameters$sigma)
# range(parameters$phi)
# 
# 
FreeADFun(obj)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e7, iter.max = 1e7,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# gc()
# a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
#                                     method = "BFGS",
#                                     control = list(maxit = 1e8,
#                                                    abstol = 1e-04, reltol = 1e-04)))
# estimates_negbin <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr,
#                          model = model, method = "bfgs", threads = threads)
# estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, 
#                        model = model, method = "bfgs", threads = threads)
# 
# save.image("Negbin_var_fixed.RData")
# load("Negbin_var_fixed.RData")
# 
# gc()
