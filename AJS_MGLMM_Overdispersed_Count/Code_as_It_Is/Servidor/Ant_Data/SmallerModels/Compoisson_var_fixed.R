######################
# compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
# load("Negbin_var_fixed.RData")
# load("compoisson_var_fixed.RData")
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
model <- "compoisson_multi_var_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Obtaining data #######################################################
########################################################################
# Parameters
# parameters <- list(beta = estimates_negbin$beta,
#                    U = U,
#                    rho = estimates_negbin$rho,
#                    nu = rep(log(1), nr))
# gc()
# threads <- 6
# openmp(n=threads)
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T,
#                  random = c("U"))
# gc()
# a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
#                                     control = list(eval.max = 1e7, iter.max = 1e7,
#                                                    abs.tol = 1e-04, rel.tol = 1e-04,
#                                                    x.tol = 1e-04, xf.tol = 1e-04)))
# estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
# model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_var_fixed.RData")
# # Second Round
load("compoisson_var_fixed.RData")
openmp(n=1)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP, time = a, nr = nr, 
                       model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_fixed.RData")
gc()

##################
# Do PORT para o PORT (PIOROU)
##################
# Parameters
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   nu = estimates_compoisson$disp)
gc()
# TO RUN NOW!!!!!
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
                                    control = list(eval.max = 1e7, iter.max = 1e7,
                                                   abs.tol = 1e-04, rel.tol = 1e-04,
                                                   x.tol = 1e-04, xf.tol = 1e-04)))
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
model = model, method = "nlimnb", threads = threads)
