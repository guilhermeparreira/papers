######################
# compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load(paste0(mainpath, "nhanes_mcglm.RData"))
source(paste0(mainpath, "unify_packages.R")) # Update functions
library(TMB)
library(mcglm)
model <- "compoisson_multi_var_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
#### Sample
mcglm_s <- betas_mcglm(fit_refe)
parameters <- list(beta = mcglm_s$beta,
                   U = U_am,
                   rho = rho_start,
                   nu = rep(log(1), nr))
data_am$sigma <- rep(1, nr)
gc()
threads <- 6
openmp(n=threads)
##################
## First Attempt
##################
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e7, iter.max = 1e7,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_compoisson_s <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_fixed_sample.RData")
threads <- 2
openmp(n=threads)
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U_am,
                   rho = estimates_compoisson_s$rho,
                   nu = estimates_compoisson_s$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data_am,
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
estimates_compoisson <- sv(fit = fit_CMP2, time = a2, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP2, time = a2, nr = nr, 
                             model = model, method = "bfgs", threads = threads)
gc()
save.image("compoisson_var_fixed_sample.RData")

