# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
source(paste0(mainpath, "unify_packages.R"))
load("Poisson_rho_fixed.RData")
gc()
# Já que as estimativas do beta deram ruim, vou rodar com valor inicial 0.
parameters <- list(U = U,
                   # beta = matrix(0, nrow = ncol(data$X), ncol = nr),
                   beta = estimates_poisson_s$beta, # Didn't work
                   sigma = estimates_poisson_s$sigma,
                   # sigma = sapply(y, sd),
                   phi = rep(log(1), nr))
# Complexity
# nrow(y)*nr #1230 observations in total
# prod(dim(estimates_poisson_s$beta))+length(rep(0, (nr*(nr-1))/2))+length(estimates_poisson_s$sigma)+length(rep(1,nr))

######################
# Full - PORT ------------------------------------------------------------------------
######################
model <- "negbin_multivariate_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
threads <- 12
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
# gc()
# obj$fn()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
estimates_negbin_1 <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlminb", threads = threads)
save.image("Negbin_rho_fixed.RData")
######################
# Full - 2º PORT ------------------------------------------------------------------------
######################
parameters <- list(U = U,
                   beta = estimates_negbin_1$beta, # Didn't work
                   sigma = estimates_negbin_1$sigma,
                   phi = estimates_negbin_1$disp)
threads <- 12
openmp(n=threads)
obj2 <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
# gc()
# obj$fn()
a <- system.time(fit_NLB <- nlminb(obj2$par, obj2$fn, obj2$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
estimates_negbin_2 <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlminb", threads = threads)
estimates_negbin_2_sd <- sv_sd(obj = obj2, fit = fit_NLB, time = a, nr = nr,
                         model = model, method = "nlminb", threads = threads)
save.image("Negbin_rho_fixed.RData")
load("Negbin_rho_fixed.RData")

######################
# Full - BFGS ------------------------------------------------------------------------
######################
parameters <- list(beta = estimates_negbin_2$beta,
                   U = U,
                   sigma = estimates_negbin_2$sigma,
                   phi = estimates_negbin_2$disp)
gc()
threads <- 12
openmp(n=threads)
dyn.load(dynlib(paste0(mainpath, model)))
# Rinterface(paste0(mainpath, model, ".cpp"))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))

a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
# Não rodou
estimates_negbin_s2 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr,
                              model = model, method = "bfgs", threads = threads)
estimates_negbin_1 <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr,
                         model = model, method = "BFGS", threads = threads)
save.image("Negbin_rho_fixed.RData")
# #### Second time
parameters <- list(beta = estimates_negbin_1$beta,
                   U = U,
                   sigma = estimates_negbin_1$sigma,
                   phi = estimates_negbin_1$disp)
# gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
# 
a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
# ###### Third time
# # 1753
estimates_negbin_1 <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr,
                         model = model, method = "BFGS", threads = threads)
estimates_negbin_s2 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr,
                             model = model, method = "bfgs", threads = threads)
save.image("Negbin_rho_fixed.RData")
