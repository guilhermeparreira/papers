######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load("Negbin_var_comum_sample.RData")
library(TMB)
model <- "negbin_multivariate_var_comum" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Full data ############################################################
########################################################################
rm(fit_NLB);rm(fit_NLB2)
parameters <- list(beta = estimates_negbin_s$beta,
                   U = U,
                   rho = estimates_negbin_s$rho,
                   sigma = estimates_negbin_s$sigma,
                   phi = estimates_negbin_s$disp)
threads <- 4
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
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
save.image("Negbin_var_comum.RData")
load("Negbin_var_comum.RData")
threads <- 4
openmp(n=threads)
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma,
                   phi = estimates_negbin$disp)
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
estimates_negbin <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
estimates_negbin <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, 
                       model = model, method = "bfgs", threads = threads)
save.image("Negbin_var_comum.RData")
gc()
########################################################################
# PORT AGAIN ############################################################
########################################################################
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma,
                   phi = estimates_negbin$disp)
threads <- 4
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
# TENTOU 1E-02 EM VÃO
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
save.image("Negbin_var_comum.RData")
