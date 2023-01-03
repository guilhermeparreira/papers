######################
# compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load("compoisson_var_fixed_sample.RData")
source(paste0(mainpath, "unify_packages.R")) # Update functions
library(TMB)
model <- "compoisson_multi_var_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
### Sample
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U,
                   rho = estimates_compoisson_s$rho,
                   nu = estimates_compoisson_s$disp)
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
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
fit_CMP
estimates_compoisson_2 <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
estimates_compoisson_sd_2 <- sv_sd(obj = obj, fit = fit_CMP, time = a, nr = nr,
                             model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_fixed.RData")



load("compoisson_var_fixed.RData")
threads <- 12
openmp(n=threads)
parameters <- list(beta = estimates_compoisson_2$beta,
                   U = U,
                   rho = estimates_compoisson_2$rho,
                   nu = estimates_compoisson_2$disp)
dyn.load(dynlib(paste0(mainpath, model)))
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
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
fit_CMP2
estimates_compoisson <- sv(fit = fit_CMP2, time = a2, n_betas = n_params, nr = nr,
                         model = model, method = "bfgs", threads = threads)
# estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP2, time = a2, nr = nr, 
#                              model = model, method = "bfgs", threads = threads)
# gc()
# save.image("compoisson_var_fixed.RData")

#####################################
## PORT AGAIN
#####################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   nu = estimates_compoisson$disp)
gc()
threads <- 12
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP, time = a, nr = nr,
                                 model = model, method = "nlminb", threads = threads)
save.image("compoisson_var_fixed2.RData")
##### Comparando modelos
load("compoisson_var_fixed.RData")
estimates_compoisson$beta
load("compoisson_var_fixed2.RData")
estimates_compoisson$beta
load("compoisson_var_fixed.RData")
estimates_compoisson$sigma
load("compoisson_var_fixed2.RData")
estimates_compoisson$sigma
load("compoisson_var_fixed.RData")
estimates_compoisson$rho
load("compoisson_var_fixed2.RData")
estimates_compoisson$rho
load("compoisson_var_fixed.RData")
estimates_compoisson$disp
load("compoisson_var_fixed2.RData")
estimates_compoisson$disp
load("compoisson_var_fixed.RData")
estimates_compoisson$summary
load("compoisson_var_fixed2.RData")
estimates_compoisson$summary
