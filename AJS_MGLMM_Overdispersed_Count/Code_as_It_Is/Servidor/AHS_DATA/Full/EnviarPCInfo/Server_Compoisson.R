# Package --------------------------------------------------------------------
rm(list = ls())
model <- "compoisson_multi" #1st try
###########
# PC SERVIDOR
###########
# library(Rcpp, lib.loc = "/home/est/bonat/nobackup/")
# library(TMB, lib.loc = "/home/est/bonat/nobackup/")
# library(RcppEigen, lib.loc = "/home/est/bonat/nobackup/")
library(TMB)
source(paste0("unify_packages.R"))
# gc()
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(model)))

###########
# PC MEU
###########
# library(TMB)
# mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
# setwd(paste0(mainpath, "AHS_DATA/Full/"))
# load("negbin.RData")
# source(paste0(mainpath, "unify_packages.R"))
# compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
# dyn.load(dynlib(paste0(mainpath, model)))
#####################
# compoisson Multivariate ------------------------------------------------------------------------
#####################
#####################
# FULL PORT ------------------------------------------------------------------------
#####################
# load("compoisson_sample.RData")
# threads <- 4
load("compoisson_sample.RData")
threads <- 5
TMB::openmp(n=threads)
config(tape.parallel = FALSE, DLL = model)
parameters <- list(beta = estimates_compoisson_am$beta,
                   U = U,
                   rho = estimates_compoisson_am$rho,
                   sigma = estimates_compoisson_am$sigma,
                   nu = estimates_compoisson_am$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
dyn.load(dynlib(paste0(model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP_1 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
fit_CMP_1
estimates_compoisson_s <- sv(fit = fit_CMP_1, time = a, n_betas = n_params, nr = nr, 
                             model = model, method = "nlimnb", threads = threads)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP_1, time = a, nr = nr, 
                                 model = model, method = "nlminb", threads = threads)
save.image("compoisson_1.RData", version = 2)
FreeADFun(obj)
###########################
################ BFGS
###########################
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
fit_CMP_2
save.image("compoisson_2.RData", version = 2)
estimates_compoisson_bfgs_s <- sv(fit = fit_CMP2, time = a2, n_betas = n_params, nr = nr, 
                                      model = model, method = "BFGS", threads = threads)
estimates_compoisson_bfgs_sd <- sv_sd(obj = obj, fit = fit_CMP2, time = a2, nr = nr, 
                          model = model, method = "BFGS", threads = threads)
gc()
save.image("compoisson_2.RData", version = 2)
FreeADFun(obj)
###########################
################ PORT (PELO bfgs)
###########################
parameters <- list(beta = estimates_compoisson_bfgs_s$beta,
                   U = U,
                   rho = estimates_compoisson_bfgs_s$rho,
                   sigma = estimates_compoisson_bfgs_s$sigma,
                   nu = estimates_compoisson_bfgs_s$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
dyn.load(dynlib(paste0(model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP_3 <- nlminb(obj$par, obj$fn, obj$gr,
                                     control = list(eval.max = 1e9, iter.max = 1e9,
                                                    abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
fit_CMP_3
estimates_compoisson_port_2 <- sv(fit = fit_CMP_3, time = a, n_betas = n_params, nr = nr, 
                             model = model, method = "nlimnb", threads = threads)
estimates_compoisson_port_2_sd <- sv_sd(obj = obj, fit = fit_CMP_3, time = a, nr = nr, 
                                 model = model, method = "nlminb", threads = threads)
save.image("compoisson_3.RData", version = 2)
FreeADFun(obj)
###########################
################ PORT (PELO PORT)
###########################
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U,
                   rho = estimates_compoisson_s$rho,
                   sigma = estimates_compoisson_s$sigma,
                   nu = estimates_compoisson_s$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
dyn.load(dynlib(paste0(model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP_4 <- nlminb(obj$par, obj$fn, obj$gr,
                                     control = list(eval.max = 1e9, iter.max = 1e9,
                                                    abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
fit_CMP_4
estimates_compoisson_port_3 <- sv(fit = fit_CMP_4, time = a, n_betas = n_params, nr = nr, 
                                  model = model, method = "nlimnb", threads = threads)
estimates_compoisson_port_3_sd <- sv_sd(obj = obj, fit = fit_CMP_4, time = a, nr = nr, 
                                        model = model, method = "nlminb", threads = threads)
save.image("compoisson_4.RData", version = 2)
FreeADFun(obj)
