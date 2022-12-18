######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
model <- "negbin_multivariate_phi_fixed" #1st try
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
data <- list(X = as.matrix(cbind(1, X)),
             Y = as.matrix(y))
data$phi <- rep(log(1), nr)
n_params <- 6
# Parameters
U <- matrix(0, nrow = nrow(y), ncol = ncol(y))
parameters <- list(beta = matrix(0, nrow = ncol(data$X), ncol = nr),
                   U = U,
                   rho = rep(0, (nr*(nr-1))/2),
                   sigma = sapply(y, sd))
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
                                    control = list(eval.max = 1e9, iter.max = 1e9, 
                                                   abs.tol = 1e-04, rel.tol = 1e-04,
                                                   x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
save.image("Negbin_phi_fixed.RData")
###############################
#### 2º PORT
###############################
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma)
gc()
threads <- 12
openmp(n=threads)
obj2 <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_NLB2 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                             model = model, method = "nlimnb", threads = threads)
save.image("Negbin_phi_fixed.RData")
# threads <- 12
# openmp(n=threads)
# parameters <- list(beta = estimates_negbin$beta,
#                    U = U,
#                    rho = estimates_negbin$rho,
#                    sigma = estimates_negbin$sigma)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
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
# save.image("Negbin_phi_fixed.RData")
# gc()
