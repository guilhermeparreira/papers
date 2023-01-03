######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load("Negbin_phi_fixed_sample.RData")
library(TMB)
model <- "negbin_multivariate_phi_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Full data ############################################################
########################################################################
rm(fit_NLB);rm(fit_NLB2);rm(res_cov)
parameters <- list(beta = estimates_negbin_s$beta,
                   U = U,
                   rho = estimates_negbin_s$rho,
                   sigma = estimates_negbin_s$sigma)
# parameters <- list(beta = matrix(0, nrow = 5, ncol = 3),
#                    U = U,
#                    rho = rep(0,3),
#                    sigma = rep(1,3))
data$phi <- rep(log(1), nr)
threads <- 6
openmp(n=threads)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
#Com eval e iter 1e9, objective = 4025,195
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                    control = list(eval.max = 1e9, 
                                                   iter.max = 1e9, 
                                                   abs.tol = 1e-04, rel.tol = 1e-04,
                                                   x.tol = 1e-04, xf.tol = 1e-04)))
# fit_NLB
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
save.image("Negbin_phi_fixed.RData")
####### BFGS OPTIM, JÁ ERA! DEU CONVERGE IGUAL A 0, COM A MESMA COISA QUE O PORT
# estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                       # model = model, method = "nlimnb", threads = threads)
# save.image("Negbin_phi_fixed_port.RData")
# load("Negbin_phi_fixed.RData")
threads <- 6
openmp(n=threads)
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma)
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
# estimates_negbin <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, 
#                        model = model, method = "bfgs", threads = threads)
# save.image("Negbin_phi_fixed.RData")
# gc()
######################################
### PORT - SECOND TIME
######################################
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma)
threads <- 6
openmp(n=threads)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, 
                                                  iter.max = 1e9, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
save.image("Negbin_phi_fixed.RData")
######################################
### PORT - THIRD TIME
######################################
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma)
threads <- 6
openmp(n=threads)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, 
                                                  iter.max = 1e9, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
# Error in optimHess - gradient in optim evaluated to length 1 not 21 
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
save.image("Negbin_phi_fixed.RData")
#######################################
### OPTIM - BFGS
#######################################
# threads <- 6
# openmp(n=threads)
# parameters <- list(beta = estimates_negbin$beta,
#                    U = U,
#                    rho = estimates_negbin$rho,
#                    sigma = estimates_negbin$sigma)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
# gc()
# a2 <- system.time(fit_NLB <- optim(obj$par, obj$fn, obj$gr,
#                                    method = "BFGS",
#                                    control = list(maxit = 1e9,
#                                                   abstol = 1e-04, reltol = 1e-04)))
# estimates_negbin <- sv(fit = fit_NLB, time = a2, n_betas = n_params, nr = nr,
#                        model = model, method = "bfgs", threads = threads)
# estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a2, nr = nr,
#                              model = model, method = "BFGS", threads = threads)
