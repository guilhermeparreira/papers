# Package --------------------------------------------------------------------
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load(paste0(mainpath, "NHANES_DATA/Full/double_sample.RData"))
library(TMB)
model <- "doublepoisson_multi_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
#####################
# Sample port------------------------------------------------------------------------
#####################
#### Sample
threads <- 10
TMB::openmp(n=threads)
parameters <- list(beta = estimates_double_am_2$beta,
                   U = U_am,
                   sigma = estimates_double_am_2$sigma,
                   nu = estimates_double_am_2$disp)
data_am$rho <- rep(0, nr)
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
b <- system.time(fit_double_s <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
# Rapidão mesmo.
save.image("doublepoisson_sample_rho_fixed.RData")
load("doublepoisson_sample_rho_fixed.RData")
estimates_double_am <- sv(fit = fit_double_s, time = b, n_betas = n_params,
                        nr = nr, model = model, method = "nlimnb", threads = threads)
dyn.load(dynlib(paste0(mainpath, model)))
#####################
# Full - PORT ------------------------------------------------------------------------
#####################
load("doublepoisson_sample_rho_fixed.RData")
parameters <- list(beta = estimates_double_am$beta,
                   U = U,
                   sigma = estimates_double_am$sigma,
                   nu = estimates_double_am$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# First Round
data$rho <- rep(0, 3)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_double <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
estimates_double_s <- sv(fit = fit_double, time = a, n_betas = n_params, nr = nr, model = model, method = "nlimnb", threads = threads)
# From std
estimates_double_s_2 <- sv_sd(obj = obj, fit = fit_double, time = a, nr = nr, model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_rho_fixed.RData")
load("doublepoisson_rho_fixed.RData")
####################
# Full - BFGS ------------------------------------------------------------------------
####################
parameters <- list(beta = estimates_double_s$beta,
                   U = U,
                   sigma = estimates_double_s$sigma,
                   nu = estimates_double_s$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_double2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
gc()
estimates_double_s_1 <- sv(fit = fit_double2, time = a, n_betas = n_params, nr = nr, model = model, method = "bfgs", threads = threads)

estimates_double_sd_2 <- sv_sd(obj = obj, fit = fit_double2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
save.image("doublepoisson_rho_fixed.RData")
load("doublepoisson_rho_fixed.RData")
####################
# Full - PORT ------------------------------------------------------------------------
####################
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_double_s_1$beta,
                   U = U,
                   sigma = estimates_double_s_1$sigma,
                   nu = estimates_double_s_1$disp)
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
a <- system.time(fit_double03 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
# estimates_double é o escolhido
estimates_double03 <- sv(fit = fit_double03, time = a, n_betas = n_params, nr = nr,
                     model = model, method = "nlimnb", threads = threads)
# From std
estimates_double03_sd <- sv_sd(obj = obj, fit = fit_double03, time = a, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_rho_fixed.RData")
load("doublepoisson_rho_fixed.RData")