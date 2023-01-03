######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load("compoisson_nu_fixed_sample_15.RData")
library(TMB)
model <- "compoisson_multi_nu_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Full data ############################################################
########################################################################
rm(fit_NLB);rm(fit_NLB2);rm(res_cov)
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U,
                   rho = estimates_compoisson_s$rho,
                   sigma = estimates_compoisson_s$sigma)
data$nu <- c(log(15),log(15),log(15))
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
estimates_compoisson <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
save.image("compoisson_nu_fixed_15.RData")
load("compoisson_nu_fixed_15.RData")
threads <- 1
openmp(n=threads)
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma)
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
                                    control = list(maxit = 1e9,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_compoisson <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr, 
                       model = model, method = "bfgs", threads = threads)
# TERMINEI AQUI!!
save.image("compoisson_nu_fixed_15.RData")
load("compoisson_nu_fixed_15.RData")
openmp(n=12)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, 
                          model = model, method = "bfgs", threads = threads)
save.image("compoisson_nu_fixed_15.RData")
gc()
########################################################################
# PORT - SECOND ############################################################
########################################################################
threads <- 3
openmp(n=threads)
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_compoisson <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                           model = model, method = "nlimnb", threads = threads)
# mesmo com 1e-04 não deu certo o sdreport() e melhor otimizador
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, 
                              model = model, method = "nlimnb", threads = threads)
save.image("compoisson_nu_fixed_15.RData")