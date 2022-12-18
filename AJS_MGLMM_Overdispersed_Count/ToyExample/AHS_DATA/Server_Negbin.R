# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
setwd(paste0(mainpath, "AHS_DATA/"))
load("Poisson.RData")
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
gc()
######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
model <- "negbin_multivariate" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
#### Sample
threads <- 6
TMB::openmp(n=threads)
parameters <- list(beta = estimates_poisson$beta,
                   U = U_am,
                   rho = estimates_poisson$rho,
                   sigma = estimates_poisson$sigma,
                   phi = log(rep(1, nr)))
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = "U")
b <- system.time(fit_NLB_neg_bin <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
# From std
estimates_negbin_am_2 <- sv_sd(obj = obj, fit = fit_NLB_neg_bin, time = b, nr = nr, model = model, method = "nlimnb", threads = threads)
# Rapidão mesmo.
# estimates_negbin_am <- sv(fit = fit_NLB_neg_bin, time = b, n_betas = n_params,
                          # nr = nr, model =model, method = "nlimnb", threads = threads)
# Check
# save.image("negbin_sample.RData")
#### 5190 full
load("negbin_sample.RData")
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
threads <- 1 #Where to use 1 or maximum number of threads of your PC is hard to define. Try and error
openmp(n=threads)
parameters <- list(beta = estimates_negbin_am$beta,
                   U = U,
                   rho = estimates_negbin_am$rho,
                   sigma = estimates_negbin_am$sigma,
                   phi = estimates_negbin_am$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e8, iter.max = 1e8,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
estimates_negbin_s <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                             model = model, method = "nlimnb", threads = threads)
load("negbin.RData")
obj$fn()