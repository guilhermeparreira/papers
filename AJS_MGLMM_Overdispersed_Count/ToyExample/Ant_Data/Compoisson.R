# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
# MEU PC
mainpath <- "~/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
setwd(paste0(mainpath, "Ant_Data/"))
model <- "compoisson_multi" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
load("Poisson.RData")
model <- "compoisson_multi" #1st try
###############
# In this code, you can try the first attempt, or go straigth to the 3º Try.
###############
######################
# COM-Poisson Multivariate ------------------------------------------------------------------------
######################
# gc()
# #### 5190 Full
parameters <- list(beta = estimates_poisson_s$beta,
                   U = matrix(0, nrow(y), nr),
                   # rho = start_rho,
                   rho = estimates_poisson_s$rho,
                   sigma = estimates_poisson_s$sigma,
                   nu = rep(log(1), nr))

gc()
threads <- 1
openmp(n=threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = c("U"))
obj$fn()
################################
# 1ª TRY - Didn't converge
################################
# gc()
a <- system.time(fit_comp0 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))

ss <- sv(fit = fit_comp0, time = a, n_betas = 6,
         nr = nr, model = model, method = "PORT", threads = threads)

save.image("Compoisson_full.RData")
# gc()
################################
# 2ª TRY
################################
load("Compoisson_full.RData")
mainpath <- "~/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
parameters <- list(beta = ss$beta,
                   U = matrix(0, nrow(y), nr),
                   rho = as.numeric(ss$rho),
                   sigma = as.numeric(ss$sigma),
                   nu = as.numeric(ss$disp))
gc()
dyn.load(dynlib(model))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = c("U"))

b <- system.time(fit_comp1 <- nlminb(fit_comp0$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04,
                                                  x.tol = 1e-04, xf.tol = 1e-04)))
ss <- sv(fit = fit_comp1, time = a, n_betas = 6,
         nr = nr, model = model, method = "PORT", threads = threads)

save.image("Compoisson_full.RData")
# Ficou na mesma coisa, não andou
# Optim

################################
# 3ª TRY - Optim
################################
# library(TMB)
load("Compoisson_full.RData")
mainpath <- "~/Dropbox/SuperDisperso_AJS/Code/ToyExample/"

threads <- 1
openmp(n = threads)
parameters <- list(beta = ss$beta,
                   U = matrix(0, nrow(y), nr),
                   rho = ss$rho,
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
save.image("Compoisson_full.RData", version = 2)
threads <- 12
openmp(n = threads)
ss_sd <- sv_sd(obj = obj, fit = fit_comp2, time = a, nr = nr,
               model = model, method = "BFGS", threads = threads)
save.image("Compoisson_full.RData", version = 2)
load("Compoisson_full.RData")