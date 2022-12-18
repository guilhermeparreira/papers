# Package --------------------------------------------------------------------
rm(list = ls())
model <- "compoisson_multi" #1st try
###########
# PC SERVIDOR
###########
# library(Rcpp, lib.loc = "/home/est/bonat/nobackup/")
# library(TMB, lib.loc = "/home/est/bonat/nobackup/")
# library(RcppEigen, lib.loc = "/home/est/bonat/nobackup/")
# source(paste0("unify_packages.R"))
# gc()
# compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
# dyn.load(dynlib(paste0(model)))

###########
# PC MEU
###########
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "AHS_DATA/Full/"))
# load("negbin.RData")
source(paste0(mainpath, "unify_packages.R"))
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))

#####################
# compoisson Multivariate ------------------------------------------------------------------------
#####################
#####################
# SAMPLE ------------------------------------------------------------------------
#####################
################ PORT
# threads <- 8
# TMB::openmp(n=threads)
# parameters <- list(beta = estimates_negbin$beta,
#                    U = U_am,
#                    rho = estimates_negbin$rho,
#                    sigma = estimates_negbin$sigma,
#                    nu = log(rep(1.5, nr)))
# obj <- MakeADFun(data = data_am,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T,
#                  random = "U")
# b <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
#                                            control = list(eval.max = 1e7, iter.max = 1e7,
#                                                           abs.tol = 1e-04, rel.tol = 1e-04)))
# estimates_compoisson_am <- sv(fit = fit_CMP, time = b, n_betas = n_params,
# nr = nr, model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_sample.RData")
# save.image("compoisson_sample.RData", version = 2)
################ BFGS
load("compoisson_sample.RData")
threads <- 6
TMB::openmp(n=threads)
config(tape.parallel = FALSE, DLL = model)
parameters <- list(beta = estimates_compoisson_am$beta,
                   U = U_am,
                   rho = estimates_compoisson_am$rho,
                   sigma = estimates_compoisson_am$sigma,
                   nu = estimates_compoisson_am$disp)
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
# FreeADFun(obj)
obj$fn()
b <- system.time(fit_CMP <- optim(obj$par, obj$fn, obj$gr,
                                  method = "BFGS",
                                  control = list(maxit = 1e9, 
                                                 abstol = 1e-04, reltol = 1e-04)))
estimates_compoisson_am <- sv(fit = fit_CMP, time = b, n_betas = n_params,
                              nr = nr, model = model, method = "nlimnb", threads = threads)
# Check
save.image("compoisson_sample.RData", version = 2)
#####################
# FULL ------------------------------------------------------------------------
#####################
################ PORT
# load("compoisson_sample.RData")
# threads <- 4
openmp(n=threads)
parameters <- list(beta = estimates_compoisson_am$beta,
                   U = matrix(0, nrow = 1e3, ncol = 5),
                   rho = estimates_compoisson_am$rho,
                   sigma = estimates_compoisson_am$sigma,
                   nu = estimates_compoisson_am$disp)
dyn.load(dynlib(paste0(mainpath, model)))
# dyn.load(dynlib(paste0(model)))
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
estimates_compoisson_s <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr, 
                             model = model, method = "nlimnb", threads = threads)
save.image("compoisson_1.RData", version = 2)
FreeADFun(obj)
# load("compoisson.RData")
################ BFGS
# threads <- 4
openmp(n=threads)
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U,
                   rho = estimates_compoisson_s$rho,
                   sigma = estimates_compoisson_s$sigma,
                   nu = estimates_compoisson_s$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
dyn.load(dynlib(paste0(model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_CMP2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e9,
                                                   abstol = 1e-04, reltol = 1e-04)))
save.image("compoisson_2.RData", version = 2)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP2, time = a2, nr = nr, 
                          model = model, method = "BFGS", threads = threads)
gc()
save.image("compoisson_2.RData", version = 2)
FreeADFun(obj)