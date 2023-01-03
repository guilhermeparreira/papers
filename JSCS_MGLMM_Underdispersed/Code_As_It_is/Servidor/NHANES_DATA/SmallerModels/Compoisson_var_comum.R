######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/SmallerModels"))
load("compoisson_var_comum_sample.RData")
library(TMB)
model <- "compoisson_multi_var_comum" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
########################################################################
# Full data ############################################################
########################################################################
############ PORT
parameters <- list(beta = estimates_compoisson_s$beta,
                   U = U,
                   rho = estimates_compoisson_s$rho,
                   sigma = estimates_compoisson_s$sigma,
                   nu = estimates_compoisson_s$disp)
threads <- 12
openmp(n=threads)
TMB::config(tape.parallel = F, DLL = model)
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
fit_CMP
# Second Round
estimates_compoisson_2 <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
estimates_compoisson_2_sd <- sv_sd(obj = obj, fit = fit_CMP, time = a, nr = nr,
                             model = model, method = "nlimnb", threads = threads)

save.image("compoisson_var_comum.RData")
load("compoisson_var_comum.RData")
# estimates_compoisson <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr,
#                           model = model, method = "BFGS", threads = threads)
# save.image("compoisson_var_comum.RData")
gc()
# load("compoisson_var_comum.RData")
########################################################################
# PORT again ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                       model = model, method = "nlimnb", threads = threads)
# estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP, time = a, nr = nr,
#                            model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_var_comum.RData")
########################################################################
# PORT again II ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson_2$beta,
                   U = U,
                   rho = estimates_compoisson_2$rho,
                   sigma = estimates_compoisson_2$sigma,
                   nu = estimates_compoisson_2$disp)
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
                                   control = list(eval.max = 1e6, iter.max = 1e6,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
fit_CMP$objective
# save.image("compoisson_var_comum.RData")
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
########################################################################
# PORT again III ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
########################################################################
# PORT again IV ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_var_comum.RData")
########################################################################
# PORT again V ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
# threads <- 6
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
# save.image("compoisson_var_comum.RData")
########################################################################
# PORT again VI ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
# threads <- 6
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a <- system.time(fit_CMP <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_comum.RData")
# load("compoisson_var_comum.RData")
########################################################################
# PORT again VII ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_comum.RData")
########################################################################
# PORT again VIII ###########################################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_comum.RData")
########################################################################
# PORT again IX (Final model) ##########################################
########################################################################
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

# Second Round
estimates_compoisson <- sv(fit = fit_CMP, time = a, n_betas = n_params, nr = nr,
                           model = model, method = "nlimnb", threads = threads)
estimates_compoisson_sd <- sv_sd(obj = obj, fit = fit_CMP, time = a, nr = nr,
                                 model = model, method = "nlimnb", threads = threads)
save.image("compoisson_var_comum.RData")
load("compoisson_var_comum.RData")
########################################################################
# BFGS IX ###########################################################
########################################################################
# Acho que não rodou
FreeADFun(obj)
TMB::config(tape.parallel = F)
parameters <- list(beta = estimates_compoisson$beta,
                   U = U,
                   rho = estimates_compoisson$rho,
                   sigma = estimates_compoisson$sigma,
                   nu = estimates_compoisson$disp)
threads <- 1
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
gc()
a2 <- system.time(fit_CMP2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e7,
                                                   abstol = 1e-04, reltol = 1e-04)))
fit_CMP2

