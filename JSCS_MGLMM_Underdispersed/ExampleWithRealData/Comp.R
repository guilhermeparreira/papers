# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "~/Dropbox/Underdispersed_Count/Code/ExampleWithRealData"
setwd(paste0(mainpath))
load("negbin.RData")
gc()
#####################
# Sample port------------------------------------------------------------------------
#####################
model <- "compoisson_multi" #1st try
mainpath <- "~/Dropbox/Underdispersed_Count/Code/ExampleWithRealData"
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
# COM-Poisson is sooooo time and memory consuming. Using tape.parallel = FALSE means don't run AD (Automatic differentiation) in parallel. If your computer is strong enough you can set it to TRUE
config(tape.parallel = FALSE, DLL = model)
#### Sample
threads <- 6
TMB::openmp(n=threads)
parameters <- list(beta = estimates_negbin$beta,
                   U = U_am,
                   rho = estimates_negbin$rho,
                   sigma = estimates_negbin$sigma,
                   nu = log(rep(1, nr)))
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = "U")
b <- system.time(fit_compoisson_s <- nlminb(obj$par, obj$fn, obj$gr,
                                           control = list(eval.max = 1e7, iter.max = 1e7,
                                                          abs.tol = 1e-04, rel.tol = 1e-04)))

# It calculates the standard error of the parameters
rep <- sdreport(obj)
summary(rep, c("fixed"), p.value = T)
summary(rep, c("report"), p.value = T)


#####################
# Sample bfgs (too much memory consuming)------------------------------------------------------------------------
#####################
load("compoisson_sample.RData")
mainpath <- "~/Dropbox/Underdispersed_Count/Code/"
threads <- 2
TMB::openmp(n=threads)
parameters <- list(beta = estimates_compoisson_am$beta,
                   U = U_am,
                   rho = estimates_compoisson_am$rho,
                   sigma = estimates_compoisson_am$sigma,
                   nu = estimates_compoisson_am$disp)
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = "U")
gc()
b <- system.time(fit_compoisson_s2 <- optim(obj$par, obj$fn, obj$gr,
                                            method = "BFGS",
                                            control = list(maxit = 1e8,
                                            abstol = 1e-04, reltol = 1e-04)))

# It calculates the standard error of the parameters
rep <- sdreport(obj)
summary(rep, c("fixed"), p.value = T)
summary(rep, c("report"), p.value = T)