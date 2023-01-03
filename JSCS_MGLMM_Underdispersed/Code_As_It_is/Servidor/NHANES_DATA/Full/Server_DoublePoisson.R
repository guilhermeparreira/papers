# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/Full/"))
# load("Poisson.RData")
source(paste0(mainpath, "unify_packages.R"))
gc()
######################
# Sample ------------------------------------------------------------------------
######################
model <- "doublepoisson_multi" #1st try
# system(paste0('rm ', mainpath, model, '.so'))
# system(paste0('rm ', mainpath, model, '.o'))
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
### Sample
threads <- 10
TMB::config(tape.parallel = FALSE, DLL = model)
TMB::openmp(n=threads, DLL = model)

# data_am$Y <- data_am$Y+20 # In order to correct numerical instability
# Usar o modelo fixo para cada resposta para pegar os chutes iniciais.
sigma <- c(sigma1 = -1.0662, sigma2 = -0.9495, sigma3 = -0.1748)
mi <- structure(c(0.381, -0.0943, -0.0017, 0.4191, -0.0078, -0.1777, 
                  -0.4062, -0.0637, 0.013, -0.0054, -0.8105, 0.0472, -0.0769, 0.469, 
                  -0.0016), dim = c(5L, 3L), dimnames = list(c("Beta", "Beta", 
                                                               "Beta", "Beta", "Beta"), c("beta", "beta2", "beta3")))
sigma <- c(0,0,0)
mi <- matrix(0, nrow = 5, ncol = 3)
parameters <- list(beta = mi,
                   U = U_am,
                   rho = rep(0,3),
                   sigma = rep(.1,3),
                   nu = sigma)

obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")

obj$fn()
# b <- system.time(fit_double_v1 <- optim(obj$par, obj$fn, obj$gr,
#                                           method = 'BFGS'))

b <- system.time(fit_double_v1 <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
ad <- sdreport(obj)
summary(ad, select = 'fixed')
# From std
estimates_double_am_1_sd <- sv_sd(obj = obj, fit = fit_double_v1, time = b, nr = nr, model = model, method = "nlimnb", threads = threads)
# Rapidão mesmo.

estimates_double_am_1 <- sv(fit = fit_double_v1, time = b, n_betas = n_params,
                          nr = nr, model = model, method = "nlimnb", threads = threads)
# Check
save.image("double_sample.RData")
# load("double_sample.RData")

### Sample again
parameters <- list(beta = estimates_double_am_1$beta,
                   U = U_am,
                   rho = estimates_double_am_1$rho,
                   sigma = estimates_double_am_1$sigma,
                   nu = estimates_double_am_1$disp)

obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
# Valor inicial muito ruim
# Usar o modelo fixo para cada resposta para pegar os chutes iniciais.
# Roda só uma vez e avalia. Passa rho = 0 e variãncia da NM pequeno
obj$fn()

b <- system.time(fit_double_v2 <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))

# b <- system.time(fit_double_v2 <- optim(obj$par, obj$fn, obj$gr,
                                          # method = 'BFGS'))

# From std
estimates_double_am_2_sd <- sv_sd(obj = obj, fit = fit_double_v2, time = b, nr = nr, model = model, method = "nlimnb", threads = threads)
# Rapidão mesmo.
estimates_double_am_2 <- sv(fit = fit_double_v2, time = b, n_betas = n_params,
                          nr = nr, model = model, method = "nlimnb", threads = threads)
estimates_double_am_2_sd
save.image("double_sample.RData")
# Sample again III - Running with optim gave the same result

######################
# Full - PORT ------------------------------------------------------------------------
######################

parameters <- list(beta = estimates_double_am_2$beta,
                   U = U,
                   rho = estimates_double_am_2$rho,
                   sigma = estimates_double_am_2$sigma,
                   nu = estimates_double_am_2$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()

# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
obj$fn()
gc()
a <- system.time(fit_double_v3 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
# ad <- sdreport(obj)
# summary(ad, select = 'fixed')

gc()
estimates_double_1 <- sv(fit = fit_double_v3, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_double_1_sd <- sv_sd(obj = obj, fit = fit_double_v3, time = a, nr = nr, 
                          model = model, method = "nlimnb", threads = threads)
save.image("double.RData")
######################
# Full - BFGS ------------------------------------------------------------------------
######################
parameters <- list(beta = estimates_double_1$beta,
                   U = U,
                   rho = estimates_double_1$rho,
                   sigma = estimates_double_1$sigma,
                   nu = estimates_double_1$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
obj$fn()
a2 <- system.time(fit_double_v4 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_double_2 <- sv(fit = fit_double_v4, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
estimates_negbin_2_sd <- sv_sd(obj = obj, fit = fit_double_v4, time = a2, nr = nr, 
                               model = model, method = "bfgs", threads = threads)

gc()
save.image("double.RData")
# ls()
######################
# PORT - AGAIN (WORSE THAN THE FIRST) ------------------------------------------------------------------------
######################
######################
# Full - BFGS ------------------------------------------------------------------------
######################
parameters <- list(beta = estimates_double_2$beta,
                   U = U,
                   rho = estimates_double_2$rho,
                   sigma = estimates_double_2$sigma,
                   nu = estimates_double_2$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
obj$fn()
a <- system.time(fit_double_v5 <- nlminb(obj$par, obj$fn, obj$gr,
                                         control = list(eval.max = 1e9, iter.max = 1e9,
                                                        abs.tol = 1e-04, rel.tol = 1e-04)))
estimates_double_3 <- sv(fit = fit_double_v5, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_negbin_3_sd <- sv_sd(obj = obj, fit = fit_double_v5, time = a2, nr = nr, 
                               model = model, method = "nlimnb", threads = threads)

gc()
save.image("double.RData")
load("double.RData")


######################
# Full - PORT ------------------------------------------------------------------------
######################

parameters <- list(beta = matrix(0, nrow = 5, ncol = 3),
                   U = U,
                   rho = rep(0, 3),
                   sigma = rep(0.1, 3),
                   nu = rep(log(.5), 3))
dyn.load(dynlib(paste0(mainpath, model)))
gc()

# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
obj$fn()
gc()
a <- system.time(fit_double_v6 <- nlminb(obj$par, obj$fn, obj$gr,
                                         control = list(eval.max = 1e9, iter.max = 1e9,
                                                        abs.tol = 1e-04, rel.tol = 1e-04)))
# ad <- sdreport(obj)
# summary(ad, select = 'fixed')

gc()
estimates_double_4 <- sv(fit = fit_double_v6, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_double_4_sd <- sv_sd(obj = obj, fit = fit_double_v6, time = a, nr = nr, 
                               model = model, method = "nlimnb", threads = threads)
save.image("double.RData")

######################
# Full - OPTIM ------------------------------------------------------------------------
######################

parameters <- list(beta = estimates_double_4$beta,
                   U = U,
                   rho = estimates_double_4$rho,
                   sigma = estimates_double_4$sigma,
                   nu = estimates_double_4$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
obj$fn()
a <- system.time(fit_double_v6 <- nlminb(obj$par, obj$fn, obj$gr,
                                         control = list(eval.max = 1e9, iter.max = 1e9,
                                                        abs.tol = 1e-04, rel.tol = 1e-04)))
estimates_double_5 <- sv(fit = fit_double_v6, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "bfgs", threads = threads)
estimates_negbin_5_sd <- sv_sd(obj = obj, fit = fit_double_v6, time = a2, nr = nr, 
                               model = model, method = "bfgs", threads = threads)
save.image("double.RData")

gc()

