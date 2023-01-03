######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load(paste0(mainpath, "NHANES_DATA/Full/double_sample.RData"))
library(TMB)
model <- "doublepoisson_multi_var_comum" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Full data ############################################################
########################################################################
############ PORT
parameters <- list(beta = estimates_double_am_2$beta,
                   U = U,
                   rho = estimates_double_am_2$rho,
                   sigma = .1,
                   nu = estimates_double_am_2$disp)
threads <- 12
openmp(n=threads)
TMB::config(tape.parallel = F, DLL = model)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_double_v1 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e7, iter.max = 1e7,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# fit_double_v1
# Second Round
estimates_double_1 <- sv(fit = fit_double_v1, time = a, n_betas = n_params, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
estimates_double_1_sd <- sv_sd(obj = obj, fit = fit_double_v1, time = a, nr = nr,
                             model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_var_comum.RData")
load("doublepoisson_var_comum.RData")
gc()
########################################################################
# PORT again ###########################################################
########################################################################
parameters <- list(beta = estimates_double_1$beta,
                   U = U,
                   rho = estimates_double_1$rho,
                   sigma = estimates_double_1$sigma,
                   nu = estimates_double_1$disp)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_double_v2 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_double_2 <- sv(fit = fit_double_v2, time = a, n_betas = n_params, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
estimates_double_2_sd <- sv_sd(obj = obj, fit = fit_double_v2, time = a, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_var_comum.RData")
########################################################################
# OPTIM ###########################################################
########################################################################
parameters <- list(beta = estimates_double_1$beta,
                   U = U,
                   rho = estimates_double_1$rho,
                   sigma = estimates_double_1$sigma,
                   nu = estimates_double_1$disp)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a2 <- system.time(fit_double_v3 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e7,
                                                   abstol = 1e-04, reltol = 1e-04)))
fit_double_v3$objective
estimates_double_3 <- sv(fit = fit_double_v3, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "optim", threads = threads)
estimates_double_3_sd <- sv_sd(obj = obj, fit = fit_double_v3, time = a, nr = nr,
                               model = model, method = "optim", threads = threads)
save.image("doublepoisson_var_comum.RData")