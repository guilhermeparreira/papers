######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load("Negbin_var_comum.RData")
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
model <- "compoisson_multi_var_comum" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
#### Sample
parameters <- list(beta = estimates_negbin$beta,
                   U = U_am,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma,
                   nu = rep(log(1.5), nr))
threads <- 6
openmp(n=threads)
##################
## First Attempt - 918,4854
##################
obj <- MakeADFun(data = data_am,
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
# Second Round - 918,4823
estimates_compoisson_s <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_comum_sample.RData")
# load("compoisson_var_comum_sample.RData")
threads <- 6
openmp(n=threads)
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U_am,
                   rho = estimates_compoisson_s$rho,
                   sigma = estimates_compoisson_s$sigma,
                   nu = estimates_compoisson_s$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_compoisson_s <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
gc()
save.image("compoisson_var_comum_sample.RData")
