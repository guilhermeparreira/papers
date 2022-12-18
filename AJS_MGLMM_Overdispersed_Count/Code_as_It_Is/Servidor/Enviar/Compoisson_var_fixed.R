######################
# compoisson Multivariate ------------------------------------------------------------------------
######################
# TMB::compile("/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/compoisson_multi_var_fixed.cpp")
rm(list=ls())
model <- "compoisson_multi_var_fixed" #1st try
########### Meu PC
# mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
# setwd(paste0(mainpath, "AHS_DATA/SmallerModels"))
# source(paste0(mainpath, "unify_packages.R"))
# library(TMB)
# compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
# dyn.load(dynlib(paste0(mainpath, model)))
########### SERVER
library(Rcpp, lib.loc = "/home/est/bonat/nobackup/")
library(TMB, lib.loc = "/home/est/bonat/nobackup/")
library(RcppEigen, lib.loc = "/home/est/bonat/nobackup/")
source("unify_packages.R")
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
mainpath <- ""
########### BOTH/COMUM
load("compoisson_var_fixed_sample.RData")
gc()
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U,
                   rho = estimates_compoisson_s$rho,
                   nu = estimates_compoisson_s$disp)

parameters <- list(beta = matrix(0, nrow = n_params, ncol = nr),
                   U = U,
                   rho = rep(0, length(estimates_compoisson_s$rho)),
                   nu = rep(0, nr))

data$sigma <- rep(1, nr)
gc()
threads <- 1
openmp(n=threads)
TMB::FreeADFun(obj)
##################
## First Attempt
##################
# Rinterface(paste0(mainpath, model, ".cpp"))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
obj$fn()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_compoisson <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlimnb", threads = threads)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr,
                                 model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_fixed_1.RData", version = 2)
# load("compoisson_var_fixed.RData")

###############################################
######### BFGS
###############################################

threads <- 1
openmp(n=threads)
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   nu = estimates_compoisson$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj2 <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a2 <- system.time(fit_NLB2 <- optim(obj2$par, obj2$fn, obj2$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e9,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_compoisson2 <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
# estimates_compoisson2_sd <- sv_sd(obj = obj2, fit = fit_NLB2, time = a2, nr = nr, 
#                        model = model, method = "bfgs", threads = threads) 
gc()
save.image("compoisson_var_fixed_2.RData", version = 2)
