# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "~/Dropbox/Underdispersed_Count/Code/ExampleWithRealData/"
setwd(paste0(mainpath))
load("negbin.RData") #Load real data
gc()
#####################
# Example of Initial Estimates from GAMLSS------------------------------------------------------------------------
#####################
library(gamlss)
library(gamlss.dist)
Nmsp <- data$Y[,1]
Nmosp <- data$Y[,2]
Nspfy <- data$Y[,3]

a0 <- gamlss::gamlss(Nmsp ~ data$X -1, family = gamlss.dist::DPO())
a0 <- summary(a0)
a1 <- gamlss::gamlss(Nmosp ~ data$X -1, family = 	gamlss.dist::DPO())
a1 <- summary(a1)
a2 <- gamlss::gamlss(Nspfy ~ data$X -1, family = 	gamlss.dist::DPO())
a2 <- summary(a2)
betas <- cbind(a0[1:5,1],a1[1:5,1],a2[1:5,1])
colnames(betas) <- paste0("beta", 1:3)
initial_estimates <- list()
initial_estimates$beta <- betas
initial_estimates$rho <- rep(0, nr)
initial_estimates$sigma <- c(a0[6,1],a1[6,1],a2[6,1])
initial_estimates$disp <- c(log(.5), log(.5), log(1.2))


#####################
# Initial Estimates from fitting double a poisson for each outcome individually with sample data
#####################

sigma <- c(sigma1 = -1.0662, sigma2 = -0.9495, sigma3 = -0.1748)
mi <- structure(c(0.381, -0.0943, -0.0017, 0.4191, -0.0078, -0.1777, 
                  -0.4062, -0.0637, 0.013, -0.0054, -0.8105, 0.0472, -0.0769, 0.469, 
                  -0.0016), dim = c(5L, 3L), dimnames = list(c("Beta", "Beta", 
                                                               "Beta", "Beta", "Beta"), c("beta", "beta2", "beta3")))
sigma <- c(0,0,0)
mi <- matrix(0, nrow = 5, ncol = 3)
parameters <- list(beta = mi,
                   U = U,
                   rho = rep(0,3),
                   sigma = rep(.1,3),
                   nu = sigma)
initial_estimates <- parameters

#####################
# Running the model with port ------------------------------------------------------------------------
#####################
model <- "doublepoisson_multi" #1st try
mainpath <- "~/Dropbox/Underdispersed_Count/Code/ExampleWithRealData/"
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
# Using tape.parallel = FALSE means don't run AD (Automatic differentiation) in parallel. If your computer is strong enough you can set it to TRUE
config(tape.parallel = FALSE, DLL = model)

parameters <- list(beta = initial_estimates$beta,
                   U = U,
                   rho = initial_estimates$rho,
                   sigma = initial_estimates$sigma,
                   nu = initial_estimates$disp)
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
a <- system.time(fit_double <- nlminb(obj$par, obj$fn, obj$gr,
                                      control = list(eval.max = 1e9, iter.max = 1e9,
                                      abs.tol = 1e-04, rel.tol = 1e-04)))
ad <- sdreport(obj)
summary(ad, select = 'fixed', p.value = T)
summary(ad, select = 'fixed', p.value = T)