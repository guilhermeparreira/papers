######################
# DoublePoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load(paste0(mainpath, "NHANES_DATA/Full/double_sample.RData"))
library(TMB)
model <- "doublepoisson_multi_var_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
### Sample
parameters <- list(beta = estimates_double_am_2$beta,
                   U = U,
                   rho = estimates_double_am_2$rho,
                   nu = estimates_double_am_2$disp)
data$sigma <- rep(1, nr)
gc()
threads <- 12
openmp(n=threads)
##################
## First Attempt
##################
TMB::config(tape.parallel = F, DLL = model)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_double_v1 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_doublepoisson_1 <- sv(fit = fit_double_v1, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
estimates_doublepoisson_sd_1 <- sv_sd(obj = obj, fit = fit_double_v1, time = a, nr = nr,
                             model = model, method = "nlimnb", threads = threads)
save.image("doublepoisson_var_fixed.RData")
load("doublepoisson_var_fixed.RData")

parameters <- list(beta = estimates_doublepoisson_1$beta,
                   U = U,
                   rho = estimates_doublepoisson_1$rho,
                   nu = estimates_doublepoisson_1$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_double_v2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_doublepoisson_2 <- sv(fit = fit_double_v2, time = a2, n_betas = n_params, nr = nr,
                         model = model, method = "bfgs", threads = threads)
estimates_doublepoisson_sd_2 <- sv_sd(obj = obj, fit = fit_double_v2, time = a2, nr = nr,
                             model = model, method = "bfgs", threads = threads)
# gc()
save.image("doublepoisson_var_fixed.RData")

#####################################
## PORT AGAIN
#####################################
parameters <- list(beta = estimates_doublepoisson_2$beta,
                   U = U,
                   rho = estimates_doublepoisson_2$rho,
                   nu = estimates_doublepoisson_2$disp)
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_double_v3 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_doublepoisson_3 <- sv(fit = fit_double_v3, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
estimates_doublepoisson_sd_3 <- sv_sd(obj = obj, fit = fit_double_v3, time = a, nr = nr,
                                 model = model, method = "nlminb", threads = threads)
save.image("compoisson_var_fixed.RData")