# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
# MEU PC
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels/"))
load("Negbin_rho_fixed.RData")
source(paste0(mainpath, "unify_packages.R"))
model <- "compoisson_multi_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# SERVIDOR
source(paste0("unify_packages.R"))
model <- "compoisson_multi" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(model)))
gc()

######################
# Negbin Multivariate ------------------------------------------------------------------------
######################
# gc()
# nrow(y)*nr #1230 observations in total
# #1148 parameters to estimate (I am not taking into account the prediction of Random effects)
# prod(dim(start_beta))+length(rep(0, (nr*(nr-1))/2))+length(start_sigma)+length(rep(1,nr))
parameters <- list(beta = estimates_negbin_1$beta,
                   U = matrix(0, nrow(y), nr),
                   sigma = estimates_negbin_1$sigma,
                   nu = rep(log(1), nr))

gc()
threads <- 12
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
obj$fn()
################################
# 1ª Tentativa
################################
# gc()
a <- system.time(fit_comp0 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
ss <- sv(fit = fit_comp0, time = a, n_betas = 6,
         nr = nr, model = model, method = "PORT", threads = threads)
save.image("Compoisson_full_rho_fixed.RData")
# gc()
################################
# 2ª Tentativa
################################
parameters <- list(beta = ss$beta,
                   U = matrix(0, nrow(y), nr),
                   sigma = as.numeric(ss$sigma),
                   nu = as.numeric(ss$disp))
gc()
# dyn.load(dynlib(model))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))

b <- system.time(fit_comp1 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
ss <- sv(fit = fit_comp1, time = a, n_betas = 6,
         nr = nr, model = model, method = "PORT", threads = threads)

save.image("Compoisson_full_rho_fixed.RData")
################################
# 3ª Tentativa - Optim
################################
# library(TMB)
threads <- 12
openmp(n = threads)
parameters <- list(beta = ss$beta,
                   U = matrix(0, nrow(y), nr),
                   sigma = ss$sigma,
                   nu = ss$disp)
gc()
dyn.load(dynlib(paste0(mainpath, model)))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
## Optim
a <- system.time(fit_comp2 <- optim(obj$par, obj$fn, obj$gr,
                                  method = "BFGS",
                                  control = list(maxit = 1e9,
                                                 abstol = 1e-04, reltol = 1e-04)))
# save.image("Compoisson_full.RData")
ss <- sv(fit = fit_comp2, time = b, n_betas = 6,
         nr = nr, model = model, method = "BFGS", threads = threads)
# load("Compoisson_full.RData")
save.image("Compoisson_full_rho_fixed.RData", version = 3)
threads <- 12
openmp(n = threads)
ss_sd <- sv_sd(obj = obj, fit = fit_comp2, time = a, nr = nr,
               model = model, method = "BFGS", threads = threads)
save.image("Compoisson_full.RData", version = 2)
load("Compoisson_full.RData")
################################
# TUNNED - PORT (FROM BFGS) - DIDN'T WORK
################################
parameters <- list(beta = ss$beta,
                   U = matrix(0, nrow(y), nr),
                   sigma = ss$sigma,
                   nu = ss$disp)
gc()
dyn.load(dynlib(paste0(mainpath, model)))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
b <- system.time(fit_comp1 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
ss <- sv(fit = fit_comp1, time = b, n_betas = 6,
         nr = nr, model = model, method = "PORT", threads = threads)
ss_sd <- sv_sd(obj = obj, fit = fit_comp, time = a, nr = nr,
               model = model, method = "PORT", threads = threads)
# fit_comp1[2]
save.image("Compoisson_rho_fixed.RData")
