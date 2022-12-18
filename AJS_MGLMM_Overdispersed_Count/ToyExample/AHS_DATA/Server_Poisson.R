# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html
# Package --------------------------------------------------------------------
rm(list = ls())
library(mcglm)
library(TMB)
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
setwd(paste0(mainpath, "AHS_DATA/"))
# source(paste0(mainpath, "unify_packages.R"))
# Dataset --------------------------------------------------------------------
# Sample
set.seed(2390)
nam <- 300
ss <- sample(1:nrow(ahs), nam, replace = F)
ahsam <- ahs[ss, ]
ahsam$id <- 1:nam
## Full
# ahs <- ahs[ss[-1], ] # COMMENT LATER
ahs$id <- 1:nrow(ahs)
n <- nrow(ahs)
# General variables
vars_resp <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
nr <- length(vars_resp)
n_params <- 11
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
#########################################################################
# Start values from MCGLM
#########################################################################
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
######
## Start Model - MCGLM -------------------------------------------------------------
######
Z0 <- mc_id(ahs)
system.time(fit_refe <- mcglm(linear_pred = c(form_Ndoc, form_Nndoc, form_Nmed, form_Nhosp, form_Nadm),
                              matrix_pred = list(Z0,Z0,Z0,Z0,Z0),
                              link = rep("log", 5), variance = rep("tweedie", 5),
                              power_fixed = rep(TRUE, 5), data = ahs,
                              control_algorithm = list(correct = FALSE, max_iter = 100)))
plogLik(fit_refe)
summary(fit_refe)
gc()
######################
# Poisson Multivariate ------------------------------------------------------------------------
######################
model <- "02_poisson_multivariate" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath,model)))
#### am sample
mcglm_s <- betas_mcglm(fit_refe)
# save.image("mcglm_ahs.RData")
parameters <- list(beta = mcglm_s$beta,
                   U = U_am,
                   rho = rho_start,
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
                                   control = list(eval.max = 7000, iter.max = 7000,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
rep <- sdreport(obj)
summary(rep, select = 'fixed')
summary(rep, select = 'random')

# Check
# identical(round(estimates_poisson_am_2$rho, 3),round(estimates_poisson_am$rho, 3))
# identical(round(estimates_poisson_am_2$beta, 3),round(estimates_poisson_am$beta, 3))
# identical(round(estimates_poisson_am_2$sigma, 3),round(estimates_poisson_am$sigma, 3))
# identical(round(estimates_poisson_am_2$disp, 3),round(estimates_poisson_am$disp, 3))
# identical(estimates_poisson_am$summary, estimates_poisson_am_2$summary)
#### 5190 full
load("Poisson_sample.RData")
mainpath <- "/home/guilherme/Dropbox/SuperDisperso_AJS/Code/ToyExample/"
threads <- 12
openmp(n=threads)
parameters <- list(beta = estimates_poisson_am$beta,
                   U = U,
                   rho = estimates_poisson_am$rho,
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
                                   control = list(eval.max = 1e9, iter.max = 1e9,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
gc()
#-----------------------------------------
# Obtain estimates

rep <- sdreport(obj)
summary(rep, select = 'fixed')
summary(rep, select = 'random')
