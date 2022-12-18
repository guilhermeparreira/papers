# Poisson Multivariate ------------------------------------------------------------------------
######################
rm(list = ls())
# # Packages
library(mvabund) #Ant data
data("antTraits")
library(TMB)
# # External Codes
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/Full"))
source(paste0(mainpath, "unify_packages.R"))
# # .cpp File
model <- "02_poisson_multivariate" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
# # Data
y <- antTraits$abund # Selecting response variables
X <- antTraits$env # Selecting covariates
# nr <- ncol(y)
data <- list(X = as.matrix(cbind(1, X)),
Y = as.matrix(y))
# n_params <- 6
# # Parameters
# U <- matrix(0, nrow = nrow(y), ncol = ncol(y))
# parameters <- list(beta = matrix(0, nrow = ncol(data$X), ncol = nr),
#                    U = U,
#                    rho = rep(0, (nr*(nr-1))/2),
#                    sigma = sapply(y, sd))
# gc()
# threads <- 6
# openmp(n=threads)
# # Object
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = F,
#                  random = "U")
# gc()
# # First Round
# a <- system.time(fit_poisson <- nlminb(obj$par, obj$fn, obj$gr,
#                                    control = list(eval.max = 1e7, iter.max = 1e7,
#                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
# estimates_poisson_s <- sv(fit = fit_poisson, time = a, n_betas = n_params, nr = nr, 
#                           model = model, method = "nlimnb", threads = threads)
# estimates_poisson_s$summary
# save.image("Poisson.RData")
# Second Round
# threads <- 12
# openmp(n=threads)
# parameters <- list(beta = estimates_poisson_s$beta,
#                    U = U,
#                    rho = estimates_poisson_s$rho,
#                    sigma = estimates_poisson_s$sigma)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = F
#                  ,random = c("U"))
# gc()
# a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
#                                     method = "BFGS",
#                                     control = list(maxit = 1e8,
#                                                    abstol = 1e-04, reltol = 1e-04)))
# Second Round
# save.image("Poisson.RData")
# Rodar dps só com ela na máquina (e nada mais)
# threads <- 6
# openmp(threads)
# load("Poisson.RData")
# estimates_poisson_s2 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
# estimates_poisson_s2$summary
# save.image("Poisson.RData")
# gc()
# load("Poisson.RData")
####################################
######### PORT 2º
####################################
parameters <- list(beta = estimates_poisson_s$beta,
                   U = U,
                   rho = estimates_poisson_s$rho,
                   sigma = estimates_poisson_s$sigma)
gc()
threads <- 1
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
                                   control = list(eval.max = 1e7, iter.max = 1e7,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
estimates_poisson <- sv(fit = fit_poisson, time = a, n_betas = n_params, nr = nr,
                          model = model, method = "nlimnb", threads = threads)
estimates_poisson_sd <- sv_sd(obj = obj, fit = fit_poisson, time = a, nr = nr,
                        model = model, method = "nlimnb", threads = threads)
save.image("Poisson.RData")
load("Poisson.RData")
