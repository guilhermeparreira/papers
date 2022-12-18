# Poisson Multivariate ------------------------------------------------------------------------
######################
rm(list = ls())
# # Packages
library(mvabund) #Ant data
data("antTraits")
library(TMB)
# # External Codes
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
source(paste0(mainpath, "unify_packages.R"))
# # .cpp File
model <- "02_poisson_multivariate_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
# # Data
y <- antTraits$abund # Selecting response variables
X <- antTraits$env # Selecting covariates
nr <- ncol(y)
data <- list(X = as.matrix(cbind(1, X)),
             Y = as.matrix(y),
             rho = rep(0, (nr*(nr-1))/2))
n_params <- 6
# Parameters
U <- matrix(0, nrow = nrow(y), ncol = ncol(y))
parameters <- list(beta = matrix(0, nrow = ncol(data$X), ncol = nr),
                   U = U,
                   sigma = sapply(y, sd))
gc()
threads <- 12
openmp(n=threads)
# Object
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
gc()
# First Round
a <- system.time(fit_poisson <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
estimates_poisson_s <- sv(fit = fit_poisson, time = a, n_betas = n_params, nr = nr,
                          model = model, method = "nlimnb", threads = threads)
# nÃO RODOU
# estimates_poisson_s2 <- sv_sd(obj = obj, fit = fit_poisson, time = a2, nr = nr, model = model, 
                              # method = "nlminb", threads = threads)
# estimates_poisson_s$summary
save.image("Poisson_rho_fixed.RData")
# Second Round
threads <- 12
openmp(n=threads)
parameters <- list(beta = estimates_poisson_s$beta,
                   U = U,
                   sigma = estimates_poisson_s$sigma)
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
                                    control = list(maxit = 1e9,
                                                   abstol = 1e-04, reltol = 1e-04)))
# Second Round
save.image("Poisson_rho_fixed.RData")
# Rodar dps só com ela na máquina (e nada mais)
threads <- 6
openmp(threads)
load("Poisson_rho_fixed.RData")
estimates_poisson_s1 <- sv(fit = fit_NLB2, time = a, n_betas = n_params, nr = nr,
                          model = model, method = "bfgs", threads = threads)
# não rodou
estimates_poisson_s12 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
# estimates_poisson_s12$summary
save.image("Poisson_rho_fixed.RData")
gc()
load("Poisson_rho_fixed.RData")
####################################
######### PORT 2º (demorou demais, e eu matei). Timing stopped at: 1.395e+06 903.6 1.18e+05
####################################
parameters <- list(beta = estimates_poisson_s$beta,
                   U = U,
                   sigma = estimates_poisson_s$sigma)
gc()
threads <- 12
openmp(n=threads)
# Object
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,
                 random = "U")
gc()
# First Round
a <- system.time(fit_poisson2 <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
fit_poisson2
estimates_poisson <- sv(fit = fit_poisson, time = a, n_betas = n_params, nr = nr,
                          model = model, method = "nlimnb", threads = threads)
estimates_poisson_sd <- sv_sd(obj = obj, fit = fit_poisson, time = a, nr = nr,
                        model = model, method = "nlimnb", threads = threads)
save.image("Poisson_rho_fixed.RData")
load("Poisson_rho_fixed.RData")
