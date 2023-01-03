rm(list = ls())
library(TMB)
library(mcglm)
# library(mvtnorm)
setwd("~/Dropbox/Underdispersed_Count/Code")
# Start values from MCGLM --------------------------------------------------------------------
betas_mcglm <- function(fit){
  df <- coef(fit, type = "beta")
  nbetas <- nrow(df[df$Response==1,])
  df$betas <- rep(1:nbetas, times = max(df$Response))
  dfnew <- reshape(df[, -c(2,3)], 
                   idvar = "betas",
                   timevar = "Response",
                   direction = "wide")[, -1]
  names(dfnew) <- paste0("beta", 1:ncol(dfnew))
  dfnew <- as.matrix(dfnew)
  
  res_cov <- as.matrix(residuals(fit))
  sigma <- sqrt(diag(cov(res_cov)))
  return(list(beta = dfnew,
              sigma = sigma))
}
# Dataset --------------------------------------------------------------------
# Sample
set.seed(2390)
nam <- 300                                  # Run first with an sample
ss <- sample(1:nrow(ahs), nam, replace = F) # ahs data from mcglm package
ahsam <- ahs[ss, ]
ahsam$id <- 1:nam
## Full Data
ahs$id <- 1:nrow(ahs)
n <- nrow(ahs)
# General variables
vars_resp <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
nr <- length(vars_resp)
n_params <- 11 # Number of regressors
units <- 60*60 # Second to hours
threads <- 6
## Fixed Effects -------------------------------------------------------------
form_Ndoc <- Ndoc ~ sex + age + income + levyplus + freepoor + freerepa + illness + actdays + chcond
form_Nndoc = Nndoc ~ sex + age + income + levyplus + freepoor + freerepa + illness + actdays + chcond
form_Nmed = Nmed ~ sex + age + income + levyplus + freepoor + freerepa + illness + actdays + chcond
form_Nhosp = Nhosp ~ sex + age + income + levyplus + freepoor + freerepa + illness + actdays + chcond
form_Nadm = Nadm ~ sex + age + income + levyplus + freepoor + freerepa + illness + actdays + chcond
# All Data sets for TMB --------------------------------------------------------------------
## Poisson
data_am <- list(Y = as.matrix(ahsam[, vars_resp]),
                X = model.matrix(form_Ndoc, ahsam))
data <- list(Y = as.matrix(ahs[, vars_resp]),
             X = model.matrix(form_Ndoc, ahs))
# Start Parameters ------------------------------------------------------------------------
U_am = matrix(0, ncol = nr, nrow = nam)
U = matrix(0, ncol = nr, nrow = n)
rho_start = rep(0, (nr*(nr-1))/2)
######
## Start Model - MCGLM -------------------------------------------------------------
######
Z0 <- mc_id(ahs)
system.time(fit_refe <- mcglm(linear_pred = c(form_Ndoc, form_Nndoc, form_Nmed, form_Nhosp, form_Nadm),
                              matrix_pred = list(Z0,Z0,Z0,Z0,Z0),
                              link = rep("log", 5), variance = rep("tweedie", 5),
                              power_fixed = rep(TRUE, 5), data = ahs,
                              control_algorithm = list(correct = FALSE, max_iter = 100)))
# plogLik(fit_refe)
# summary(fit_refe)
# gc()
######################
# Poisson Multivariate ------------------------------------------------------------------------
######################
# It compiles the TMB template file
model <- "02_poisson_multivariate" 
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
#### Sample
mcglm_s <- betas_mcglm(fit_refe) # Function that obtain the parameter estimates from MCGLM
parameters <- list(beta = mcglm_s$beta,
                   U = U_am,
                   rho = rho_start,
                   sigma = mcglm_s$sigma) # Provide the initial values of the parameters
openmp(n=threads) # Run in n threads the code (in parallel)
#### Just creates the objective and differentiated functions.
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T,   # Don't print% http://www.ajs.or.at/
                 random = "U") # Name of random variable in .cpp
gc()
# Optimize
a <- system.time(fit_NLB <- nlminb(start = obj$par, 
                                   objective = obj$fn, 
                                   gradient = obj$gr,
                                   control = list(eval.max = 7000, iter.max = 7000,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
# It calculates the standard error of the parameters
rep <- sdreport(obj)
summary(rep, c("fixed"), p.value = T)
summary(rep, c("report"), p.value = T)

## SLOW due to the high number of parameters in the model -------------------
## Profiling  ---------------------------------------------------------------
par(mfrow=c(2,2), mai = rep(.7,4))
# Beta1 as log(Y)
p1 <- tmbprofile(obj,1, trace = F)
plot(p1, main = "beta1")
# Sigma5 is already in natural scale
p3 <- tmbprofile(obj, 70, trace = F)
plot(p3, main = "Sigma5")
# Rho_4_5_ in the TMB scale
p4 <- tmbprofile(obj, 65, trace = F)
plot(p4, main = "Rho_4_5")