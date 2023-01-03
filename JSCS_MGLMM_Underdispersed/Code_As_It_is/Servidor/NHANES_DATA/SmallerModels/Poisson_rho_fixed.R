# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
library(mcglm)
library(gmailr)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
complement <- "NHANES_DATA/SmallerModels/"
pathlocal <- paste0(mainpath, complement)
setwd(paste0(pathlocal))
source(paste0(mainpath, "unify_packages.R"))
# # Dataset --------------------------------------------------------------------
nhanes <- read.table(paste0(mainpath, "NHANES_DATA/nhanes.txt"), header = T)
# Sample
set.seed(2390)
nam <- 350
ss <- sample(1:nrow(nhanes), nam, replace = F)
nhanesam <- nhanes[ss, ]
## Full
n <- nrow(nhanes)
# General variables
vars_resp <- c("Nmsp", "Nmosp", "Nspfy")
nr <- length(vars_resp)
units <- 60*60 # Second to hours
threads <- 6
n_params <- 5
## Fixed Effects -------------------------------------------------------------
form_Nmsp <- Nmsp ~ Race + Education + Marital + Age
form_Nmosp <- Nmosp ~ Race + Education + Marital + Age
form_Nspfy <- Nspfy ~ Race + Education + Marital + Age
# All Data sets for TMB --------------------------------------------------------------------
## Poisson
data_am <- list(Y = as.matrix(nhanesam[, vars_resp]),
                X = model.matrix(form_Nmsp, nhanesam),
                rho = rep(0, (nr*(nr-1))/2)) #rho fixed at 0
data <- list(Y = as.matrix(nhanes[, vars_resp]),
             X = model.matrix(form_Nmsp, nhanes),
             rho = rep(0, (nr*(nr-1))/2))
# Start Parameters ------------------------------------------------------------------------
U_am = matrix(0, ncol = nr, nrow = nam)
U = matrix(0, ncol = nr, nrow = n)
######
## Start Model - MCGLM -------------------------------------------------------------
######
Z0 <- mc_id(nhanes)
system.time(fit_refe <- mcglm(linear_pred = c(form_Nmsp, form_Nmosp, form_Nspfy),
                              matrix_pred = list(Z0,Z0,Z0),
                              link = rep("log", 3), variance = rep("tweedie", 3),
                              power_fixed = rep(TRUE, 3), data = nhanes,
                              control_algorithm = list(correct = FALSE, max_iter = 100)))
gc()
######################
# Poisson Multivariate ------------------------------------------------------------------------
######################
model <- "02_poisson_multivariate_rho_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath,model)))
gc()
#### am sample
mcglm_s <- betas_mcglm(fit_refe)
parameters <- list(beta = mcglm_s$beta,
                   U = U_am,
                   sigma = mcglm_s$sigma)
openmp(n=threads)
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = "U")
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e7, iter.max = 1e7,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
# From std (não rodou)
# estimates_poisson_am_2 <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, model = model, method = "nlimnb", threads = threads)
# Rapidão mesmo. (pois não calcula SE)
estimates_poisson_am <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, model = model, method = "nlimnb", threads = threads)
save.image("Poisson_sample_rho_fixed.RData")
#### 5190 full
load("Poisson_sample_rho_fixed.RData")
threads <- 12
openmp(n=threads)
parameters <- list(beta = estimates_poisson_am$beta,
                   U = U,
                   sigma = estimates_poisson_am$sigma)
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
                                   control = list(eval.max = 1e8, iter.max = 1e8,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
# Second Round
estimates_poisson_s <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                          model = model, method = "nlminb", threads = threads)
estimates_poisson_s2 <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                              model = model, method = "nlminb", threads = threads)
save.image("Poisson_rho_fixed.RData")
# load("Poisson_p.RData")
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
                 silent = F
                 ,random = c("U"))
gc()
a2 <- system.time(fit_NLB2 <- optim(obj$par, obj$fn, obj$gr,
                                    method = "BFGS",
                                    control = list(maxit = 1e8,
                                                   abstol = 1e-04, reltol = 1e-04)))
estimates_poisson_s12 <- sv(fit = fit_NLB2, time = a2, n_betas = n_params, nr = nr, model = model, method = "bfgs", threads = threads)
# estimates_poisson_sd_s12 <- sv_sd(obj = obj, fit = fit_NLB2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
gc()

###################################
### PORT AGAIN! (from bfgs), which is almost the same as from port (same llk)
###################################
parameters <- list(beta = estimates_poisson_s12$beta,
                   U = U,
                   sigma = estimates_poisson_s12$sigma)
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# First Round
threads <- 2
openmp(threads)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 1e8, iter.max = 1e8,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
# Second Round
estimates_poisson <- sv(fit = fit_NLB, time = a, n_betas = n_params, nr = nr, 
                          model = model, method = "nlimnb", threads = threads)
# Não compilou
# estimates_poisson_sd <- sv_sd(obj = obj, fit = fit_NLB, time = a, nr = nr, 
                          # model = model, method = "nlimnb", threads = threads)
save.image("Poisson_rho_fixed.RData")
# load("Poisson_rho_fixed.RData")
