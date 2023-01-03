# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels/"))
source(paste0(mainpath, "unify_packages.R"))
load("negbin_rho_fixed.RData")
gc()
#####################
# Sample port------------------------------------------------------------------------
#####################
model <- "compoisson_multi_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
#### Sample
threads <- 12
TMB::openmp(n=threads)
parameters <- list(beta = estimates_negbin0$beta,
                   U = U_am,
                   sigma = estimates_negbin0$sigma,
                   nu = log(rep(1, nr)))
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
b <- system.time(fit_compoisson_s <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
# Rapidão mesmo.
save.image("comp_sample_rho_fixed.RData")
load("comp_sample_rho_fixed.RData")
estimates_comp_am <- sv(fit = fit_compoisson_s, time = b, n_betas = n_params,
nr = nr, model = model, method = "nlimnb", threads = threads)
dyn.load(dynlib(paste0(mainpath, model)))
#####################
# Full - PORT ------------------------------------------------------------------------
#####################
load("comp_sample_rho_fixed.RData")
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_comp_am$beta,
                   U = U,
                   sigma = estimates_comp_am$sigma,
                   nu = estimates_comp_am$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_compoisson <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
estimates_comp_s <- sv(fit = fit_compoisson, time = a, n_betas = n_params, nr = nr, model = model, method = "nlimnb", threads = threads)
# From std
estimates_comp_s_2 <- sv_sd(obj = obj, fit = fit_compoisson, time = a, nr = nr, model = model, method = "nlimnb", threads = threads)
save.image("comp_rho_fixed.RData")
load("comp_rho_fixed.RData")
####################
# Full - BFGS ------------------------------------------------------------------------
####################
threads <- 6
openmp(n=threads)
parameters <- list(beta = estimates_comp_s$beta,
                   U = U,
                   sigma = estimates_comp_s$sigma,
                   nu = estimates_comp_s$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_compoisson2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
gc()
estimates_comp <- sv_sd(obj = obj, fit = fit_compoisson2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
save.image("comp_rho_fixed.RData")
load("comp_rho_fixed.RData")
####################
# Full - PORT ------------------------------------------------------------------------
####################
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_comp$beta,
                   U = U,
                   sigma = estimates_comp$sigma,
                   nu = estimates_comp$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_compoisson03 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
# estimates_comp é o escolhido
estimates_comp03 <- sv(fit = fit_compoisson03, time = a, n_betas = n_params, nr = nr,
                     model = model, method = "nlimnb", threads = threads)
# From std
estimates_comp03_sd <- sv_sd(obj = obj, fit = fit_compoisson03, time = a, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
save.image("comp_rho_fixed.RData")
load("comp_rho_fixed.RData")