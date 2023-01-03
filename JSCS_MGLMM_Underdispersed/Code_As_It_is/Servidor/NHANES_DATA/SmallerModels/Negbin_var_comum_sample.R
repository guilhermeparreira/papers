######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load(paste0(mainpath, "nhanes_mcglm.RData"))
source(paste0(mainpath, "unify_packages.R")) # Update functions
library(TMB)
library(mcglm)
model <- "negbin_multivariate_var_comum" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
#### Sample
mcglm_s <- betas_mcglm(fit_refe)
parameters <- list(beta = mcglm_s$beta,
                   U = U_am,
                   rho = rho_start,
                   sigma = mean(mcglm_s$sigma), #Var comum to all data
                   phi = rep(log(1), nr))
gc()
threads <- 8
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
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e7, iter.max = 1e7, 
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_negbin_s <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_negbin_s <- sv_sd(obj = obj, fit = fit_NLB, time = a, 
                            nr = nr, model = model, method = "nlimnb", threads = threads)
save.image("Negbin_var_comum_sample.RData")
threads <- 8
openmp(n=threads)
parameters <- list(beta = estimates_negbin_s$beta,
                   U = U_am,
                   rho = estimates_negbin_s$rho,
                   sigma = estimates_negbin_s$sigma,
                   phi = estimates_negbin_s$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data_am,
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
estimates_negbin_s <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
gc()
save.image("Negbin_var_comum_sample.RData")

