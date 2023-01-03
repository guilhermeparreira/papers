######################
# Full - PORT (PELO VAR COMUM) ------------------------------------------------------------------------
######################
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
locals <- "/NHANES_DATA/SmallerModels"
setwd(paste0(path, "/NHANES_DATA/Full"))
load(paste0(path, locals,"/doublepoisson_var_comum.RData"))
model <- "doublepoisson_multi"
parameters <- list(beta = estimates_double_3$beta,
                   U = U,
                   rho = estimates_double_3$rho,
                   sigma = rep(estimates_double_3$sigma, 3),
                   nu = estimates_double_3$disp)
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
save.image("double_v2.RData")
######################
# Full - PORT AGAIN (PELO VAR COMUM) ------------------------------------------------------------------------
######################
parameters <- list(beta = estimates_double_4$beta,
                   U = U,
                   rho = estimates_double_4$rho,
                   sigma = estimates_double_4$sigma,
                   nu = estimates_double_4$disp)
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
a <- system.time(fit_double_v7 <- nlminb(obj$par, obj$fn, obj$gr,
                                         control = list(eval.max = 1e9, iter.max = 1e9,
                                                        abs.tol = 1e-04, rel.tol = 1e-04)))
# ad <- sdreport(obj)
# summary(ad, select = 'fixed')

gc()
estimates_double_5 <- sv(fit = fit_double_v7, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
estimates_double_5_sd <- sv_sd(obj = obj, fit = fit_double_v7, time = a, nr = nr, 
                               model = model, method = "nlimnb", threads = threads)

save.image("double_v2.RData")
######################
# Full - optim (MESMA COISA) ------------------------------------------------------------------------
######################
# parameters <- list(beta = estimates_double_5$beta,
#                    U = U,
#                    rho = estimates_double_5$rho,
#                    sigma = estimates_double_5$sigma,
#                    nu = estimates_double_5$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# # First Round
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
# obj$fn()
# gc()
# a2 <- system.time(fit_double_v8 <- optim(obj$par, obj$fn, obj$gr,
#                                          method = "BFGS",
#                                          control = list(maxit = 1e8,
#                                                         abstol = 1e-04, reltol = 1e-04)))
# # ad <- sdreport(obj)
# # summary(ad, select = 'fixed')
# 
# gc()
# estimates_double_6 <- sv(fit = fit_double_v8, time = a, n_betas = n_params, nr = nr, 
#                          model = model, method = "nlimnb", threads = threads)
# estimates_double_6_sd <- sv_sd(obj = obj, fit = fit_double_v8, time = a, nr = nr, 
#                                model = model, method = "nlimnb", threads = threads)
# 
# save.image("double_v2.RData")
