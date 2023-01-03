# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
load("Poisson.RData") # I loaded the initial guess to the parameters from Poisson model.
mainpath <- "~/Dropbox/Underdispersed_Count/Code/ExampleWithRealData"
setwd(mainpath)
gc()
######################
# Sample ------------------------------------------------------------------------
######################
model <- "negbin_multivariate" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
#### Sample
threads <- 6
TMB::openmp(n=threads)
parameters <- list(beta = estimates_poisson_s2$beta,
                   U = U_am,
                   rho = estimates_poisson_s2$rho,
                   sigma = estimates_poisson_s2$sigma,
                   phi = log(rep(1, nr)))
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
b <- system.time(fit_NLB_neg_bin <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))
# It calculates the standard error of the parameters
rep <- sdreport(obj)
summary(rep, c("fixed"), p.value = T)
summary(rep, c("report"), p.value = T)
######################
# Full - PORT ------------------------------------------------------------------------
######################
load("negbin_sample.RData")
mainpath <- "~/Dropbox/Underdispersed_Count/Code/"
threads <- 6
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
# It calculates the standard error of the parameters from default TMB function
rep <- sdreport(obj)
summary(rep, c("fixed"), p.value = T)
summary(rep, c("report"), p.value = T)

# It obtain the pontual estimates of the parameters in order to use as input to the next model run
estimates_negbin <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                         model = model, method = "nlimnb", threads = threads)
# It calculates the SE of parameters from my customized function sv_sd. 
estimates_negbin_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                          model = model, method = "nlimnb", threads = threads)
# save.image("negbin.RData")
# load("negbin.RData")

######################
# Full - BFGS ------------------------------------------------------------------------
######################
load("negbin.RData")
mainpath <- "~/Dropbox/Underdispersed_Count/Code/"
threads <- 12
openmp(n=threads)
parameters <- list(beta = estimates_negbin$beta,
                   U = U,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma,
                   phi = estimates_negbin$disp)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
obj <- MakeADFun(data = data,
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
estimates_negbin_sd_2 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
gc()
save.image("negbin.RData")