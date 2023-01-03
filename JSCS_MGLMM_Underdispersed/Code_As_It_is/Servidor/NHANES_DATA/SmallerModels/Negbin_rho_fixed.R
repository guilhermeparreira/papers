# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels/"))
source(paste0(mainpath, "unify_packages.R"))
load("Poisson_rho_fixed.RData")
gc()
######################
# Sample ------------------------------------------------------------------------
######################
model <- "negbin_multivariate_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
#### Sample
threads <- 6
TMB::openmp(n=threads)
parameters <- list(beta = estimates_poisson$beta,
                   U = U_am,
                   sigma = estimates_poisson$sigma,
                   phi = log(rep(1, nr)))
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = "U")
b <- system.time(fit_NLB_neg_bin <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
# Reduz o uso de memória!!
config(tape.parallel = 0, DLL = model)
#### MUITO LENTO AQUI!!!
# From std
estimates_negbin_am_2 <- sv_sd(obj = obj, fit = fit_NLB_neg_bin, time = b, nr = nr, model = model, method = "nlimnb", threads = threads)
# Rapidão mesmo.
estimates_negbin_am <- sv(fit = fit_NLB_neg_bin, time = b, n_betas = n_params,
nr = nr, model = model, method = "nlimnb", threads = threads)
# Check
save.image("negbin_sample_rho_fixed.RData")
load("negbin_sample_rho_fixed.RData")
######################
# Full - PORT ------------------------------------------------------------------------
######################
load("negbin_sample.RData")
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_negbin_am$beta,
                   U = U,
                   sigma = estimates_negbin_am$sigma,
                   phi = estimates_negbin_am$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
TMB::FreeADFun(obj)
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
estimates_negbin0 <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
# Não rdou
estimates_negbin0_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                          model = model, method = "nlimnb", threads = threads)
save.image("negbin_rho_fixed.RData")
# load("negbin.RData")
# estimates_negbin_am
# estimates_negbin_s
# estimates_negbin
######################
# Full - BFGS ------------------------------------------------------------------------
######################
load("negbin_rho_fixed.RData")
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_negbin0$beta,
                   U = U,
                   sigma = estimates_negbin0$sigma,
                   phi = estimates_negbin0$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
# Não rodou
estimates_negbin1_sd_2 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
gc()
save.image("negbin_rho_fixed.RData")
# ls()
######################
# PORT - AGAIN (WORSE THAN THE FIRST) ------------------------------------------------------------------------
######################
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_negbin0$beta,
                   U = U,
                   sigma = estimates_negbin0$sigma,
                   phi = estimates_negbin0$disp)
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
a <- system.time(fit_NLB3 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
# Modelo final, já é o suficiente.
estimates_negbin <- sv(fit = fit_NLB3, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
# Não rodou
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB3, time = a, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
save.image("negbin_rho_fixed.RData")
