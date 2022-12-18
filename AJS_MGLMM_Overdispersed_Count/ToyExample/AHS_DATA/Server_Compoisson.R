###
# DISCLAIMER: THIS CODE TAKES FOREVER TO RUN
###
# Package --------------------------------------------------------------------
rm(list = ls())
model <- "compoisson_multi" #1st try
###########
# PC MEU
###########
library(TMB)
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
setwd(paste0(mainpath, "AHS_DATA/"))
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
load("negbin.RData")
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
model <- 'compoisson_multi'
#####################
# Compoisson Multivariate ------------------------------------------------------------------------
#####################

#####################
# SAMPLE ------------------------------------------------------------------------
#####################
################ PORT
threads <- 8
TMB::openmp(n=threads)
parameters <- list(beta = estimates_negbin_s$beta,
                   U = U_am,
                   rho = estimates_negbin_s$rho,
                   sigma = estimates_negbin_s$sigma,
                   nu = log(rep(1.5, nr)))
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
b <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
# From standard functions of TMB
rep <- sdreport(obj)
summary(rep, select = 'fixed')
summary(rep, select = 'random')
# My customized function to obtain estimates in the way I need to use as input to the next model
estimates_compoisson_am <- sv(fit = fit_CMP, time = b, n_betas = n_params, nr = nr, model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_sample.RData", version = 2)

#####################
# FULL ------------------------------------------------------------------------
#####################
################ PORT
load("compoisson_sample.RData")
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
threads <- 4
openmp(n=threads)
parameters <- list(beta = estimates_compoisson_am$beta,
                   U = matrix(0, nrow = 5190, ncol = 5),
                   rho = estimates_compoisson_am$rho,
                   sigma = estimates_compoisson_am$sigma,
                   nu = estimates_compoisson_am$disp)
dyn.load(dynlib(paste0(mainpath, model)))
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
estimates_compoisson_s <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr, 
                             model = model, method = "nlimnb", threads = threads)