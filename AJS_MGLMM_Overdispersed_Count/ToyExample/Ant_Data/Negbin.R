# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "~/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
setwd(paste0(mainpath, "Ant_Data"))
load("Poisson.RData")
gc()

# Complexity
# nrow(y)*nr #1230 observations in total
# prod(dim(estimates_poisson_s$beta))+length(rep(0, (nr*(nr-1))/2))+length(estimates_poisson_s$sigma)+length(rep(1,nr))

######################
# Full - PORT ------------------------------------------------------------------------
######################
model <- "negbin_multivariate" #1st try
mainpath <- "~/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# With all "standard parameters, and beta from Poisson didn't work
# With all "standard parameters, and sigma from Poisson worked
# With all "standard parameters, and rho from Poisson worked
# With all "standard parameters, rho and sigma from Poisson worked
# Once beta estimates didn't work, I will use 0
parameters <- list(U = U,
                   beta = matrix(0, nrow = ncol(data$X), ncol = nr),
                   # beta = estimates_poisson_s$beta, # Didn't work
                   # rho = rep(0, nr*(nr-1)/2),
                   rho = estimates_poisson_s$rho,
                   sigma = estimates_poisson_s$sigma,
                   # sigma = sapply(y, sd),
                   phi = rep(log(1), nr))
# gc()
threads <- 1
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = c("U"))
# gc()
# obj$fn()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
estimates_negbin_1 <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr,
                         model = model, method = "nlminb", threads = threads)
save.image("Negbin_full.RData")
######################
# Full - 2º PORT - Example using input just obtained ------------------------------------------------------------------------
######################
load("Negbin_full.RData")
mainpath <- "~/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
dyn.load(dynlib(paste0(mainpath, model)))

parameters <- list(U = U,
                   beta = matrix(0, nrow = ncol(data$X), ncol = nr), # Didn't work
                   rho = estimates_negbin_1$rho,
                   sigma = estimates_negbin_1$sigma,
                   phi = estimates_negbin_1$disp)
threads <- 1
openmp(n=threads)
obj2 <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
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
save.image("Negbin_full.RData")