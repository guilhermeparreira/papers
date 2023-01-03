######################
# DoublePoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load(paste0(mainpath, "NHANES_DATA/Full/double_sample.RData"))
library(TMB)
model <- "doublepoisson_multi_nu_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Full data ############################################################
########################################################################
rm(fit_NLB);rm(fit_NLB2);rm(fit_double_v1);rm(fit_double_v2)
parameters <- list(beta = estimates_double_am_2$beta,
                   U = U,
                   rho = estimates_double_am_2$rho,
                   sigma = estimates_double_am_2$sigma)
data$nu <- rep(log(.4), nr)
threads <- 10
openmp(n=threads)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
obj$fn()
a <- system.time(fit_double_v1_nu_fixed <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e7, iter.max = 1e7, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Works
# rep <- sdreport(obj)
# summary(rep, 'fixed')

# Second Round
estimates_double_nu_fixed <- sv(fit = fit_double_v1_nu_fixed, time = a, n_betas = n_params, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_nu_fixed.RData")
load("doublepoisson_nu_fixed.RData")
parameters <- list(beta = estimates_double_nu_fixed$beta,
                   U = U,
                   rho = estimates_double_nu_fixed$rho,
                   sigma = estimates_double_nu_fixed$sigma)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_double_v2_nu_fixed <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_double_nu_fixed_v2 <- sv(fit = fit_double_v2_nu_fixed, time = a2, n_betas = n_params, nr = nr, 
                       model = model, method = "bfgs", threads = threads)
save.image("doublepoisson_nu_fixed.RData")
estimates_double_nu_fixed_v2 <- sv_sd(obj = obj, fit = fit_double_v2_nu_fixed, time = a2, nr = nr, 
                          model = model, method = "bfgs", threads = threads)
save.image("doublepoisson_nu_fixed.RData")
gc()
########################################################################
# PORT - SECOND ############################################################
########################################################################
parameters <- list(beta = estimates_double_nu_fixed_v2$beta,
                   U = U,
                   rho = estimates_double_nu_fixed_v2$rho,
                   sigma = estimates_double_nu_fixed_v2$sigma)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_double_v3_nu_fixed <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_double_nu_fixed_v3_sv <- sv(fit = fit_double_v3_nu_fixed, time = a, n_betas = n_params, nr = nr, 
                           model = model, method = "nlimnb", threads = threads)
# mesmo com 1e-04 não deu certo o sdreport() e melhor otimizador
estimates_double_nu_fixed_v3 <- sv_sd(obj = obj, fit = fit_double_v3_nu_fixed, time = a2, nr = nr, 
                              model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_nu_fixed.RData")