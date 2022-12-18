######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
# load("Negbin_var_comum.RData")
model <- "compoisson_multi_var_comum" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# ########################################################################
# # Obtaining data #######################################################
# ########################################################################
# # Data
# # Parameters
# parameters <- list(beta = estimates_negbin$beta,
#                    U = U,
#                    rho = estimates_negbin$rho,
#                    sigma = mean(sapply(y, sd)),
#                    nu = rep(log(1.5), nr))
# gc()
# threads <- 10
# openmp(n=threads)
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T,
#                  random = c("U"))
# gc()
# a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
#                                    control = list(eval.max = 1e7, iter.max = 1e7,
#                                                   abs.tol = 1e-04, rel.tol = 1e-04,
#                                                   x.tol = 1e-04, xf.tol = 1e-04)))
# # Second Round
# estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
#                        model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_var_comum.RData")
load("compoisson_var_comum.RData")
# threads <- 1
# openmp(n=threads)
# parameters <- list(beta = estimates_compoisson$beta,
#                    U = U,
#                    rho = estimates_compoisson$rho,
#                    sigma = estimates_compoisson$sigma,
#                    nu = estimates_compoisson$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
# gc()
# a2 <- system.time(fit_CMP2 <- optim(obj$par, obj$fn, obj$gr,
#                                     method = "BFGS",
#                                     control = list(maxit = 1e8,
#                                                    abstol = 1e-04, reltol = 1e-04)))
# estimates_compoisson <- sv(fit = fit_CMP2, time = a2, n_betas = n_params, nr = nr, 
#                        model = model, method = "bfgs", threads = threads)
# fit_CMP2[2]
# save.image("compoisson_var_comum.RData")
openmp(n=1)
# estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP2, time = a2, nr = nr, 
#                              model = model, method = "bfgs", threads = threads)
# save.image("compoisson_var_comum.RData")
# gc()

### PIOROU DE VEZ!!


# threads <- 1
# openmp(n=threads)
# parameters <- list(beta = estimates_compoisson$beta,
#                    U = U,
#                    rho = estimates_compoisson$rho,
#                    sigma = estimates_compoisson$sigma,
#                    nu = estimates_compoisson$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
# gc()
# a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
#                                    control = list(eval.max = 1e9, iter.max = 1e9,
#                                                   abs.tol = 1e-04, rel.tol = 1e-04,
#                                                   x.tol = 1e-04, xf.tol = 1e-04)))
# # Second Round
# estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
#                        model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_var_comum.RData")